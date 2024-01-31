#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Identify partitions that contain cells."""

import dataclasses
import json
import os
from collections import defaultdict

import martian
import numpy as np
import pandas as pd
from six import ensure_binary, ensure_str

import cellranger.cell_calling_helpers as helpers
import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
import cellranger.rna.matrix as rna_matrix
import cellranger.utils as cr_utils
import tenkit.safe_json as tk_safe_json
from cellranger.barcodes.utils import load_probe_barcode_map
from cellranger.cell_barcodes import InvalidCellBarcode
from cellranger.cell_calling_helpers import (
    FilterMethod,
    HighOccupancyGemSummary,
    MetricGroups,
    get_filter_method_from_string,
)
from cellranger.csv_utils import combine_csv, write_filtered_barcodes
from cellranger.feature.antibody.analysis import FRACTION_CORRECTED_READS, FRACTION_TOTAL_READS
from cellranger.metrics import BarcodeFilterResults
from cellranger.multi.config import CrMultiGraph
from cellranger.reference_paths import get_reference_genomes
from cellranger.targeted.rtl_multiplexing import get_probe_bc_defn, get_probe_bc_whitelist
from cellranger.targeted.simple_utils import determine_targeting_type_from_csv
from cellranger.targeted.targeted_constants import TARGETING_METHOD_HC

FILTER_BARCODES_MIN_MEM_GB = 2.0

PROBE_BC_SAMPLE_ID = "id"
PROBE_BC_SEQS = "sequence"
PROBE_BC_OFFSET = "offset"
PROBE_BC_LEN = "length"

np.random.seed(0)


# TODO: We can get rid of all the gem-group logic here,
# because we run this only on data from a single gem-group

__MRO__ = """
struct ProbeBCDef(
    string   id,
    string[] sequence,
    int      offset,
    int      length,
)

struct CellCallingParam(
    int      per_gem_well,
    map<int> per_sample,
)

struct CellCalling(
    CellCallingParam recovered_cells,
    CellCallingParam force_cells,
    CellCallingParam emptydrops_minimum_umis,
    json             cell_barcodes,
    string           override_mode,
    string[]         override_library_types,
    bool             disable_ab_aggregate_detection,
    bool             disable_high_occupancy_gem_detection,
)

stage FILTER_BARCODES(
    in  string       sample_id,
    in  h5           matrices_h5,
    in  csv          barcode_correction_csv,
    in  bool         is_antibody_only,
    in  path         reference_path,
    in  int[]        gem_groups,
    in  string       chemistry_description,
    in  CellCalling  config,
    in  csv          target_set,
    in  ChemistryDef chemistry_def,
    in  json         multi_graph,
    in  csv          per_barcode_metrics,
    in  bool         is_spatial,
    out json         summary,
    out csv          filtered_barcodes,
    out csv          aggregate_barcodes,
    out h5           filtered_matrices_h5,
    out path         filtered_matrices_mex,
    out csv          nonambient_calls,
    src py           "stages/counter/filter_barcodes",
) split (
    in  ProbeBCDef   probe_bc_def,
    out json         filtered_metrics_groups,
    out json         filtered_bcs_groups,
) using (
    mem_gb   = 8,
    volatile = strict,
)
"""


def split(args):
    # We need to store one full copy of the matrix.
    mem_gb = 2.3 * cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.matrices_h5)
    mem_gb = max(mem_gb, FILTER_BARCODES_MIN_MEM_GB)

    offset, length = get_probe_bc_defn(args.chemistry_def["barcode"])
    chunks = []
    # Perform cell-calling per sample when input is multiplexed using probe bcs
    if length is not None and args.multi_graph is not None:
        probe_bcs = parse_multi_graph(args.multi_graph, args.chemistry_def["barcode"])
        for sample_id, probe_bc_seqs in probe_bcs.items():
            chunks.append(
                {
                    "probe_bc_def": {
                        "id": sample_id,
                        "sequence": probe_bc_seqs,
                        "offset": offset,
                        "length": length,
                    },
                    "__mem_gb": mem_gb,
                }
            )
    else:
        chunks.append({"probe_bc_def": None, "__mem_gb": mem_gb})
    return {
        "chunks": chunks,
        "join": {
            "__mem_gb": mem_gb,
        },
    }


def main(args, outs):
    filter_barcodes(args, outs)


def join(args, outs, chunk_defs, chunk_outs):
    filtered_metrics_groups = parse_chunked_metrics(chunk_outs)
    filtered_bcs, genome_filtered_bcs = parse_chunked_filtered_bcs(chunk_outs)
    summary = parse_chunked_summary(chunk_outs)

    # high occupancy GEM filtering to mitigate any GEMs that appear overloaded
    # in initial cell calls. Only applies if input has probe bcs and not disabled.
    config = helpers.CellCalling(args.config)
    probe_bc_offset, _ = get_probe_bc_defn(args.chemistry_def["barcode"])
    if probe_bc_offset and not config.disable_high_occupancy_gem_detection:
        (
            filtered_bcs,
            genome_filtered_bcs,
            high_occupancy_gem_metrics,
        ) = helpers.remove_bcs_from_high_occupancy_gems(
            filtered_bcs, genome_filtered_bcs, args.per_barcode_metrics, probe_bc_offset
        )
        high_occupancy_gem_metrics = dataclasses.asdict(high_occupancy_gem_metrics)
    else:
        high_occupancy_gem_metrics = dataclasses.asdict(HighOccupancyGemSummary())

    # Select cell-associated barcodes
    raw_matrix = cr_matrix.CountMatrix.load_h5_file(args.matrices_h5)
    filtered_matrix = raw_matrix.select_barcodes_by_seq(filtered_bcs)

    # If non-spatial targeted ensure all filtered barcodes have nonzero targeted UMI counts
    is_targeted = raw_matrix.feature_ref.has_target_features()
    if (not args.is_spatial) and is_targeted:
        filtered_bcs, genome_filtered_bcs = helpers.remove_cells_with_zero_targeted_counts(
            raw_matrix, filtered_bcs, genome_filtered_bcs, summary
        )
        filtered_matrix = raw_matrix.select_barcodes_by_seq(filtered_bcs)

    # subset the filtered matrix to only targeted genes
    if is_targeted:
        target_features = raw_matrix.feature_ref.get_target_feature_indices()
        filtered_matrix = filtered_matrix.remove_genes_not_on_list(target_features)

    # Write the filtered barcodes file
    write_filtered_barcodes(outs.filtered_barcodes, genome_filtered_bcs)

    matrix_attrs = cr_matrix.make_matrix_attrs_count(
        args.sample_id, args.gem_groups, args.chemistry_description
    )
    filtered_matrix.save_h5_file(
        outs.filtered_matrices_h5,
        extra_attrs=matrix_attrs,
        sw_version=martian.get_pipelines_version(),
    )

    rna_matrix.save_mex(
        filtered_matrix, outs.filtered_matrices_mex, martian.get_pipelines_version()
    )

    parse_filtered_bcs_method(filtered_metrics_groups, summary)

    genomes = get_reference_genomes(args.reference_path)
    summary = helpers.combine_initial_metrics(
        genomes, filtered_metrics_groups, genome_filtered_bcs, summary
    )

    # Add keys in from high occupancy GEM filtering in if it was performed
    summary.update(high_occupancy_gem_metrics)

    # Write metrics json
    with open(ensure_binary(outs.summary), "w") as f:
        tk_safe_json.dump_numpy(summary, f, indent=4, sort_keys=True)

    chunk_aggregate_barcodes = [
        chunk.aggregate_barcodes
        for chunk in chunk_outs
        if chunk.aggregate_barcodes is not None and os.path.exists(chunk.aggregate_barcodes)
    ]
    combine_csv(chunk_aggregate_barcodes, outs.aggregate_barcodes)

    chunk_nonambient_calls = [
        chunk.nonambient_calls
        for chunk in chunk_outs
        if chunk.nonambient_calls is not None and os.path.exists(chunk.nonambient_calls)
    ]
    combine_csv(chunk_nonambient_calls, outs.nonambient_calls)


def filter_barcodes(args, outs):
    """Identify cell-associated partitions."""

    # Apply the cell calling config in args.config
    config = helpers.CellCalling(args.config)
    force_cells = helpers.get_force_cells(
        config.force_cells,
        args.probe_bc_def.get(PROBE_BC_SAMPLE_ID, None) if args.probe_bc_def else None,
    )
    recovered_cells = helpers.get_recovered_cells(
        config.recovered_cells,
        args.probe_bc_def.get(PROBE_BC_SAMPLE_ID, None) if args.probe_bc_def else None,
    )

    try:
        correction_data = pd.read_csv(
            ensure_str(args.barcode_correction_csv), converters={"barcode": ensure_binary}
        )
    except pd.errors.EmptyDataError:
        correction_data = None  # will circumvent aggregate detection below
    if correction_data is not None:
        correction_data["barcode"] = correction_data["barcode"].astype("S")

    raw_matrix = cr_matrix.CountMatrix.load_h5_file(args.matrices_h5)
    has_cmo_data = raw_matrix.feature_ref.has_feature_type(rna_library.MULTIPLEXING_LIBRARY_TYPE)

    # If there is a probe barcode def, restrict to data from this specific sample
    if args.probe_bc_def and args.probe_bc_def[PROBE_BC_SEQS] != ["all"]:
        probe_bc_seq = [seq.encode() for seq in args.probe_bc_def[PROBE_BC_SEQS]]
        offset = args.probe_bc_def[PROBE_BC_OFFSET]
        ending_idx = offset + args.probe_bc_def[PROBE_BC_LEN]
        mask = [bc[offset:ending_idx] in probe_bc_seq for bc in raw_matrix.bcs]
        bcs_idx = [i for i, x in enumerate(mask) if x]
        raw_matrix = raw_matrix.select_barcodes(bcs_idx)

    is_targeted = raw_matrix.feature_ref.has_target_features()
    is_rtl = (
        determine_targeting_type_from_csv(args.target_set) != TARGETING_METHOD_HC
        if is_targeted
        else False
    )

    if config.cell_barcodes is not None:
        method = FilterMethod.MANUAL
    elif force_cells is not None:
        method = FilterMethod.TOP_N_BARCODES
    else:
        # Override from MRO takes precedence
        if config.override_mode is not None:
            method = get_filter_method_from_string(config.override_mode)
        elif args.is_antibody_only:
            method = FilterMethod.ORDMAG
        elif is_targeted and not is_rtl:
            method = FilterMethod.TARGETED
        else:
            method = FilterMethod.ORDMAG_NONAMBIENT

    # Get unique gem groups
    unique_gem_groups = sorted(list(set(args.gem_groups)))

    if config.override_library_types is None:
        if args.is_antibody_only:
            feature_types = [rna_library.ANTIBODY_LIBRARY_TYPE]
        else:
            feature_types = [rna_library.GENE_EXPRESSION_LIBRARY_TYPE]
    else:
        feature_types = config.override_library_types

    # For targeted GEX, retrieve target gene indices for cell calling
    if is_targeted:
        target_features = raw_matrix.feature_ref.get_target_feature_indices()
        assert all(
            x in range(raw_matrix.features_dim) for x in target_features
        ), "Invalid index value found in target_features."
    else:
        target_features = None

    lib_types = raw_matrix.get_library_types()
    # In case of SC_MULTI determine library_types from input config.csv
    if args.multi_graph is not None:
        multi_graph = CrMultiGraph.from_json_file(args.multi_graph)
        lib_types = {
            feature
            for lib in multi_graph.libraries
            for feature in lib.library_features.library_features()
        }

    if (
        method not in [FilterMethod.MANUAL]
        and correction_data is not None
        and any(
            x in correction_data.library_type.unique()
            for x in [rna_library.ANTIBODY_LIBRARY_TYPE, rna_library.ANTIGEN_LIBRARY_TYPE]
        )
    ):

        (matrix, summary, removed_bcs_df,) = helpers.remove_antibody_antigen_aggregates(
            correction_data, raw_matrix, lib_types, config.disable_ab_aggregate_detection
        )
        # report all identified aggregate barcodes, together with their umis,
        # umi corrected reads, fraction of corrected reads, and fraction of total reads
        removed_bcs_df = removed_bcs_df.round(
            {FRACTION_CORRECTED_READS: 3, FRACTION_TOTAL_READS: 3}
        )
        if len(removed_bcs_df) != 0:
            removed_bcs_df.to_csv(outs.aggregate_barcodes)
    else:
        matrix = raw_matrix
        summary = {}

    if config.cell_barcodes is not None:
        assert method == FilterMethod.MANUAL
    elif force_cells is not None:
        assert method == FilterMethod.TOP_N_BARCODES

    total_diversity_key = (
        args.probe_bc_def[PROBE_BC_SAMPLE_ID] + "_total_diversity"
        if args.probe_bc_def
        else "total_diversity"
    )
    summary[total_diversity_key] = matrix.bcs_dim

    genomes = get_reference_genomes(args.reference_path)

    # Get per-gem group cell load
    if recovered_cells is not None:
        gg_recovered_cells = int(float(recovered_cells) / float(len(unique_gem_groups)))
    elif method == FilterMethod.ORDMAG or method == FilterMethod.ORDMAG_NONAMBIENT:
        gg_recovered_cells = None
    else:
        gg_recovered_cells = cr_constants.DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP

    ### Get the initial cell calls
    probe_barcode_sample_id = args.probe_bc_def[PROBE_BC_SAMPLE_ID] if args.probe_bc_def else None
    num_probe_barcodes = len(args.probe_bc_def[PROBE_BC_SEQS]) if args.probe_bc_def else None
    try:
        filtered_metrics_groups, filtered_bcs_groups = helpers.call_initial_cells(
            matrix,
            genomes,
            probe_barcode_sample_id,
            unique_gem_groups,
            method,
            gg_recovered_cells,
            config.cell_barcodes,
            force_cells,
            feature_types,
            args.chemistry_description,
            target_features,
            has_cmo_data,
            num_probe_barcodes=num_probe_barcodes,
        )
    except InvalidCellBarcode as ex:
        # This is thrown deeper in the code, but caught here with a martian.exit call
        msg = (
            "Cell Barcodes did not match list of valid barcodes with observed reads.  "
            "If attempting manual cell calling, please be sure all input "
            "barcodes are expected for that chemistry.\n\nError: {}".format(ex)
        )
        martian.exit(msg)

    ### Do additional cell calling
    if method == FilterMethod.ORDMAG_NONAMBIENT:
        filtered_bcs_groups, nonambient_summary = helpers.call_additional_cells(
            matrix,
            unique_gem_groups,
            genomes,
            filtered_bcs_groups,
            feature_types,
            args.chemistry_description,
            num_probe_barcodes=num_probe_barcodes,
            emptydrops_minimum_umis=helpers.get_emptydrops_minimum_umis(
                config.emptydrops_minimum_umis, probe_barcode_sample_id
            ),
        )
        if nonambient_summary.empty:
            outs.nonambient_calls = None
        else:
            nonambient_summary.to_csv(outs.nonambient_calls)
    else:
        outs.nonambient_calls = None

    def remap_keys_metrics(groups):
        return [{"key": k, "value": v.__dict__} for k, v in groups.items()]

    with open(outs.filtered_metrics_groups, "w") as f:
        tk_safe_json.dump_numpy(
            remap_keys_metrics(filtered_metrics_groups),
            f,
            indent=4,
            sort_keys=True,
        )

    def remap_keys_bcs(groups):
        return [{"key": k, "value": v} for k, v in groups.items()]

    with open(outs.filtered_bcs_groups, "w") as f:
        tk_safe_json.dump_numpy(
            remap_keys_bcs(filtered_bcs_groups),
            f,
            indent=4,
            sort_keys=True,
        )

    with open(outs.summary, "w") as f:
        tk_safe_json.dump_numpy(summary, f, indent=4, sort_keys=True)


def parse_chunked_metrics(chunk_outs):
    filtered_metrics_groups_list = []
    for chunk_out in chunk_outs:
        if chunk_out.filtered_metrics_groups is not None:
            with open(chunk_out.filtered_metrics_groups) as infile:
                filtered_metrics_groups_list.append(json.load(infile))

    filtered_metrics_groups = defaultdict(set)
    for filtered_metrics in filtered_metrics_groups_list:
        for groups in filtered_metrics:
            key_tuple = MetricGroups(
                groups["key"][0], groups["key"][1], groups["key"][2], groups["key"][3]
            )
            result = BarcodeFilterResults(0)
            result.filtered_bcs = groups["value"]["filtered_bcs"]
            result.filtered_bcs_var = groups["value"]["filtered_bcs_var"]
            result.filtered_bcs_cv = groups["value"]["filtered_bcs_cv"]
            result.filtered_bcs_lb = groups["value"]["filtered_bcs_lb"]
            result.filtered_bcs_ub = groups["value"]["filtered_bcs_ub"]
            filtered_metrics_groups[key_tuple] = result
    return filtered_metrics_groups


def parse_chunked_filtered_bcs(chunk_outs):
    filtered_bcs_groups_list = []
    for chunk_out in chunk_outs:
        if chunk_out.filtered_bcs_groups is not None:
            with open(chunk_out.filtered_bcs_groups) as infile:
                filtered_bcs_groups_list.append(json.load(infile))

    # collapse out gem_group and (gem_group, genome)
    genome_filtered_bcs = defaultdict(set)
    filtered_bcs = set()
    for filtered_bcs_groups in filtered_bcs_groups_list:
        for groups in filtered_bcs_groups:
            genome = groups["key"][1]
            bcs = groups["value"]
            bcs = [bc.encode() for bc in bcs]
            genome_filtered_bcs[genome].update(bcs)
            filtered_bcs.update(bcs)

    # Deduplicate and sort filtered barcode sequences
    # Sort by (gem_group, barcode_sequence)
    def barcode_sort_key(x):
        return cr_utils.split_barcode_seq(x)[::-1]

    for genome, bcs in genome_filtered_bcs.items():
        genome_filtered_bcs[genome] = sorted(list(set(bcs)), key=barcode_sort_key)
    filtered_bcs = sorted(list(set(filtered_bcs)), key=barcode_sort_key)

    return filtered_bcs, genome_filtered_bcs


def parse_chunked_summary(chunk_outs):
    summary = defaultdict(int)
    for chunk_out in chunk_outs:
        if chunk_out.summary is not None:
            with open(chunk_out.summary) as infile:
                per_chunk_summary = json.load(infile)
                for k, v in per_chunk_summary.items():
                    summary[k] += v
    return dict(summary)


def parse_filtered_bcs_method(filtered_metrics_groups, summary):
    for key_tuple in filtered_metrics_groups:
        if key_tuple.sample is not None:
            summary.update({key_tuple.sample + "_filter_barcodes_method": key_tuple.method})
        else:
            summary.update({"filter_barcodes_method": key_tuple.method})


def parse_multi_graph(path, barcode_def):
    probe_bc_wl = get_probe_bc_whitelist(barcode_def)
    wl_map = load_probe_barcode_map(probe_bc_wl)
    assert wl_map is not None
    multi_graph = CrMultiGraph.from_json_file(path)
    tags_by_sample = {}
    for sample in multi_graph.samples:
        seen_tags = []
        for fingerprint in sample.fingerprints:
            if fingerprint.tag_name is not None:
                seen_tags.append(wl_map[fingerprint.tag_name])
        assert len(seen_tags) == len(set(seen_tags))
        # To account for inputs with MFRP chemistry where no sample is defined in the
        # multi config.csv and the implicit assumption is that all observed probe bcs belong
        # to the input library == sample
        tags_by_sample[sample.sample_id] = seen_tags if len(seen_tags) > 0 else ["all"]
    return tags_by_sample
