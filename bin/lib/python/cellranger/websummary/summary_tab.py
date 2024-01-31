#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#

"""Here, you could find everything shown in the summary tab."""


from __future__ import annotations

from typing import TYPE_CHECKING

import cellranger.constants as cr_constants
import cellranger.reference as cr_reference
import cellranger.rna.library as rna_library
import cellranger.vdj.constants as vdj_constants
import cellranger.webshim.common as cr_webshim
import cellranger.webshim.constants.shared as shared_constants
import cellranger.websummary.plotly_tools as pltly
import cellranger.websummary.sample_properties as wsp
from cellranger.feature.throughputs import (
    HT_THROUGHPUT,
    INCONSISTENT_THROUGHPUT_METRIC,
    LT_THROUGHPUT,
    MT_THROUGHPUT,
    THROUGHPUT_INFERRED_METRIC,
)
from cellranger.h5_constants import H5_CHEMISTRY_DESC_KEY
from cellranger.targeted.targeted_constants import TARGETING_METHOD_HC, TARGETING_METHOD_TL
from cellranger.webshim.data import FILTERED_BCS_TRANSCRIPTOME_UNION
from cellranger.websummary.metrics import INFO_THRESHOLD, WARNING_THRESHOLD

if TYPE_CHECKING:
    from cellranger.webshim.data import SampleData

ALARMS = shared_constants.ALARMS

_ANTIBODY_reads_lost_to_aggregate_GEMs = "ANTIBODY_reads_lost_to_aggregate_GEMs"
_ANTIGEN_reads_lost_to_aggregate_GEMs = "ANTIGEN_reads_lost_to_aggregate_GEMs"
RTL_GENES_DETECTED = "num_genes_detected_on_target"
RTL_MEDIAN_GENES_BARCODE = "median_genes_per_cell_on_target"
ARC_V1_MULTIOME_DESCRIPTION = "Single Cell Multiome ATAC + Gene Expression v1"

# Feature barcoding internel name <-> display name
FEATURE_BARCODINGS = [
    rna_library.CRISPR_METRIC_PREFIX,
    rna_library.ANTIBODY_METRIC_PREFIX,
    rna_library.ANTIGEN_METRIC_PREFIX,
    rna_library.CUSTOM_METRIC_PREFIX,
]

FB_DISPLAY_NAME = {
    rna_library.CRISPR_METRIC_PREFIX: "CRISPR",
    rna_library.ANTIBODY_METRIC_PREFIX: "Antibody",
    rna_library.ANTIGEN_METRIC_PREFIX: "Antigen",
    rna_library.CUSTOM_METRIC_PREFIX: "Custom Feature",
}

RANK_PLOT_HELP = [
    [
        "Barcode Rank Plot",
        [
            "The plot shows the count of filtered UMIs mapped to each barcode. Barcodes are not determined to be cell-associated strictly based on their UMI count. Instead, they could be determined based on their expression profile, or removed via Protein Aggregate Detection and Filtering and/or High Occupancy GEM Filtering. Therefore, some regions of the graph contain both cell-associated and background-associated barcodes. The color of the graph in these regions is based on the local density of barcodes that are cell-associated."
        ],
    ]
]

ZOOM_IMAGE_HELP = [
    [
        "Tissue Detection and Fiducial Alignment",
        [
            (
                "Shows the tissue image in gray tones with an overlay of the aligned fiducial frame (open blue circles) and the capture area spots (gray circles). "
                "For the latter, the circles filled in red denote selected tissue-associated spots and the remaining open gray circles denote unselected spots. "
                "Hover mouse cursor over the image to magnify the view. "
                "Confirm fiducial frame aligns well with fiducial spots, e.g. the corner shapes match, and confirm selection of tissue-covered spots. "
                "If the result shows poor fiducial alignment or tissue detection, consider sharing the image with support@10xgenomics.com so we can improve the algorithm. "
                "Otherwise, perform manual alignment and spot selection with Loupe Browser."
            )
        ],
    ],
]

REGISTRATION_IMAGE_HELP = [
    "CytAssist Image Alignment",
    [
        (
            "Shows the CytAssist image aligned onto the microscope image. "
            "Click-drag the opacity slider to blend the two images to confirm good alignment, i.e. tissue boundaries and features should remain fixed in place. "
            "For QC purposes, fluorescence microscopy images are inverted to have a light background. "
            "If alignment is poor, rerun with Loupe Browser manual alignment.",
        )
    ],
]


CELL_CALLING_METRIC_KEYS = [
    "filtered_bcs",
    "filtered_bcs_conf_mapped_barcoded_reads_cum_frac",
    "multi_transcriptome_total_raw_reads_per_filtered_bc",
    "filtered_reads_per_filtered_bc",
    "filtered_bcs_median_counts",
    "filtered_bcs_median_unique_genes_detected",
    "filtered_bcs_total_unique_genes_detected",
]

TARGETED_CELL_CALLING_METRIC_KEYS = [
    "filtered_bcs",
    "filtered_bcs_conf_mapped_barcoded_reads_cum_frac",
    "multi_transcriptome_total_raw_reads_per_filtered_bc",
    "filtered_reads_per_filtered_bc",
    "total_targeted_reads_per_filtered_bc",
    "median_genes_per_cell_on_target",
    "num_genes_detected_on_target",
    "median_umis_per_cell_on_target",
]

ANTIBODY_filtered_bcs_transcriptome_union = "ANTIBODY_filtered_bcs_transcriptome_union"
ANTIBODY_CELL_CALLING_METRIC_KEYS = [
    ANTIBODY_filtered_bcs_transcriptome_union,
    "ANTIBODY_multi_transcriptome_total_raw_reads_per_filtered_bc",
]

CELL_CALLING_ALARM_KEYS = ["filtered_bcs_conf_mapped_barcoded_reads_cum_frac"]
TARGETED_CELL_CALLING_ALARM_KEYS = [
    "filtered_bcs_conf_mapped_barcoded_reads_cum_frac",
    FILTERED_BCS_TRANSCRIPTOME_UNION,
    "filtered_reads_per_filtered_bc",
]

# metric keys for sequencing (GEX and feature barcoding)
SEQUENCING_METRIC_KEYS = [
    "total_read_pairs",  # CR Number of reads
    "sequenced_reads_count",  # SR Number of reads.
    "unprocessed_read_pairs",
    "good_bc_frac",
    "good_umi_frac",
    "multi_cdna_pcr_dupe_reads_frac",
    "bc_bases_with_q30_frac",
    "read_bases_with_q30_frac",
    "read2_bases_with_q30_frac",
    "umi_bases_with_q30_frac",
]

TARGETED_SEQUENCING_METRIC_KEYS = [
    "total_read_pairs",  # CR Number of reads
    "sequenced_reads_count",
    "unprocessed_read_pairs",
    "subsampled_frac",
    "good_bc_frac",
    "good_umi_frac",
    "multi_cdna_pcr_dupe_reads_frac_on_target",
    "bc_bases_with_q30_frac",
    "read_bases_with_q30_frac",
    "read2_bases_with_q30_frac",
    "umi_bases_with_q30_frac",
]

SEQUENCING_ALARM_KEYS = [
    "good_bc_frac",
    "good_umi_frac",
    "bc_bases_with_q30_frac",
    "read_bases_with_q30_frac",
    "umi_bases_with_q30_frac",
]

AGGREGATION_METRIC_KEYS = [
    "frac_reads_kept",
    "pre_normalization_raw_reads_per_filtered_bc",
    "pre_normalization_cmb_reads_per_filtered_bc",
    "pre_normalization_targeted_cmb_reads_per_filtered_bc",
]

# metric keys for feature barcoding application
FB_APP_METRIC_KEYS = {
    "CRISPR": [
        "CRISPR_feature_bc_extracted_frac",
        "CRISPR_recognized_feature_bc_frac",
        "CRISPR_frac_feature_reads_usable",
        "CRISPR_feature_reads_usable_per_cell",
        "CRISPR_unrecognized_feature_bc_frac",
        "CRISPR_feature_reads_in_cells",
        "CRISPR_frac_cells_with_protospacer",
        "CRISPR_frac_cells_with_multiple_protospacer",
        "CRISPR_multi_filtered_bcs_median_counts",
    ],
    "ANTIBODY": [
        "ANTIBODY_recognized_feature_bc_frac",
        "ANTIBODY_frac_feature_reads_usable",
        "ANTIBODY_feature_reads_usable_per_cell",
        _ANTIBODY_reads_lost_to_aggregate_GEMs,
        "ANTIBODY_unrecognized_feature_bc_frac",
        "ANTIBODY_feature_reads_in_cells",
        "ANTIBODY_multi_filtered_bcs_median_counts",
    ],
    "ANTIGEN": [
        "ANTIGEN_recognized_feature_bc_frac",
        "ANTIGEN_frac_feature_reads_usable",
        "ANTIGEN_feature_reads_usable_per_cell",
        _ANTIGEN_reads_lost_to_aggregate_GEMs,
        "ANTIGEN_unrecognized_feature_bc_frac",
        "ANTIGEN_feature_reads_in_cells",
        "ANTIGEN_multi_filtered_bcs_median_counts",
    ],
    "Custom": [
        "Custom_recognized_feature_bc_frac",
        "Custom_frac_feature_reads_usable",
        "Custom_feature_reads_usable_per_cell",
        "Custom_unrecognized_feature_bc_frac",
        "Custom_feature_reads_in_cells",
        "Custom_multi_filtered_bcs_median_counts",
    ],
}

# targeted hero metric keys for targeted GEX -- hacky way to override hero metrics
# dict keys are the targeted hero metrics to display, values are the original GEX
# hero metrics, which are also the keys used to display the values in the html template

HERO_METRIC_MAPPING = {
    TARGETING_METHOD_TL: {
        FILTERED_BCS_TRANSCRIPTOME_UNION: FILTERED_BCS_TRANSCRIPTOME_UNION,
        "multi_transcriptome_total_raw_reads_per_filtered_bc": "multi_transcriptome_total_raw_reads_per_filtered_bc",
        "median_genes_per_cell_on_target": "filtered_bcs_median_unique_genes_detected",
    },
    TARGETING_METHOD_HC: {
        FILTERED_BCS_TRANSCRIPTOME_UNION: FILTERED_BCS_TRANSCRIPTOME_UNION,
        "multi_transcriptome_total_raw_reads_per_filtered_bc": "multi_transcriptome_total_raw_reads_per_filtered_bc",
        "multi_transcriptome_targeted_conf_mapped_reads_frac": "filtered_bcs_median_unique_genes_detected",
    },
}


MAPPING_KEYS = [
    "genome_mapped_reads_frac",
    "genome_conf_mapped_reads_frac",
    "intergenic_conf_mapped_reads_frac",
    "intronic_conf_mapped_reads_frac",
    "exonic_conf_mapped_reads_frac",
    "transcriptome_conf_mapped_reads_frac",
    "multi_transcriptome_targeted_conf_mapped_reads_frac",
    "antisense_reads_frac",
]

TEMP_LIG_MAPPING_KEYS = [
    "genome_mapped_reads_frac",
    "genome_conf_mapped_reads_frac",
    "multi_transcriptome_targeted_conf_mapped_reads_frac",
]

MAPPING_ALARMS = [
    "transcriptome_conf_mapped_reads_frac",
    "antisense_reads_frac",
]

TARGETED_MAPPING_ALARMS = [
    "multi_transcriptome_targeted_conf_mapped_reads_frac",
    "targeted_unsupported_panel",
] + MAPPING_ALARMS

TEMP_LIG_MAPPING_ALARMS = [
    "genome_conf_mapped_reads_frac",
    "multi_transcriptome_targeted_conf_mapped_reads_frac",
]


def get_empty_rank_plot():
    """Generates a template for the barcode rank plot.

    The template can be used by components making different types of rank plots.
    config/layout are pre-filled, the consumer is responsible for adding the data.

    Returns:
        dict: data for a plotly plot
    """
    return {
        "config": pltly.PLOT_CONFIG,
        "layout": {
            "title": "Barcode Rank Plot",
            "xaxis": {
                "title": "Barcodes",
                "type": "log",
                "showline": True,
                "zeroline": False,
                "fixedrange": False,
            },
            "yaxis": {
                "title": "UMI counts",
                "type": "log",
                "showline": True,
                "zeroline": False,
                "fixedrange": False,
            },
            "font": pltly.DEFAULT_WEB_FONT,
            "hovermode": "closest",
        },
        "data": [],
    }


def add_data(websummary_json, alarm_list, input_data):
    """Adds data to global dictionary."""
    if input_data is None:
        return

    if ALARMS in input_data:
        alarm_list.extend(input_data[ALARMS])
        del input_data[ALARMS]
    websummary_json.update(input_data)
    return


def hero_metrics(metadata, sample_data: SampleData, species_list):
    if sample_data is None or sample_data.summary is None:
        return None

    def _generate_data_with_alarms(data, alarm_keys):
        alarms = metadata.gen_metric_list(sample_data.summary, alarm_keys, species_list)
        new_alarms = [metric.alarm_dict for metric in alarms if metric.alarm_dict]
        if new_alarms:
            data[ALARMS] = new_alarms
        return data

    # For FB-only web summaries, only use cell counts and total antibody reads
    if sample_data.is_antibody_only():
        data = {}
        for key in ANTIBODY_CELL_CALLING_METRIC_KEYS:
            metrics = metadata.gen_metric_list(sample_data.summary, [key], species_list)
            for metric in metrics:
                # remove the prefix, so the key name matches with the input json expects
                new_key = metric.key.replace("ANTIBODY_", "")
                data[new_key] = metric.gen_metric_dict()

        alarm_keys = [ANTIBODY_filtered_bcs_transcriptome_union]
        data_with_alarms = _generate_data_with_alarms(data, alarm_keys)

        return data_with_alarms

    if sample_data.targeting_method is not None:
        data = {}
        for key, new_key in HERO_METRIC_MAPPING[sample_data.targeting_method].items():
            metrics = metadata.gen_metric_list(sample_data.summary, [key], species_list)
            for metric in metrics:
                # translate the key, so the key name matches with the input json expects
                data[new_key] = metric.gen_metric_dict()

        data_with_alarms = _generate_data_with_alarms(data, [])
        return data_with_alarms

    data = {}
    for key in [
        FILTERED_BCS_TRANSCRIPTOME_UNION,
        "multi_transcriptome_total_raw_reads_per_filtered_bc",
    ]:
        metrics = metadata.gen_metric_list(sample_data.summary, [key], species_list)
        for metric in metrics:
            data[metric.key] = metric.gen_metric_dict()

    is_barnyard = len(species_list) > 1
    if not is_barnyard:
        median_unique_genes = "filtered_bcs_median_unique_genes_detected"
        metrics = metadata.gen_metric_list(sample_data.summary, [median_unique_genes], species_list)
        for metric in metrics:
            data[median_unique_genes] = metric.gen_metric_dict()

    alarm_keys = [FILTERED_BCS_TRANSCRIPTOME_UNION]
    data_with_alarms = _generate_data_with_alarms(data, alarm_keys)

    return data_with_alarms


def _get_chemistry_description(sample_data, sample_properties):
    chemistry = sample_data.summary.get("chemistry_description")

    if hasattr(sample_properties, "is_spatial") and sample_properties.is_spatial:
        if sample_data.targeting_method == TARGETING_METHOD_TL and hasattr(
            sample_properties, "barcode_whitelist"
        ):
            whitelist = sample_properties.barcode_whitelist
            if whitelist == "visium-v1" or whitelist == "visium-v2":
                assay = "FFPE v1"
            elif whitelist is not None and whitelist in ["visium-v3", "visium-v4", "visium-v5"]:
                assay = "FFPE v2"
            else:
                assay = "FFPE"
        elif sample_data.targeting_method == TARGETING_METHOD_HC:
            assay = "Targeted"
        else:
            assay = "3'"
        chemistry += f" - {assay}"

    return chemistry


def pipeline_info_table(
    sample_data: SampleData, sample_properties, pipeline, metadata=None, species_list=None
):
    """Generates a table of general pipeline information.

    Args:
        sample_data:
        sample_properties:
        pipeline:
        metadata:
        species_list:

    Returns:
    """
    assert isinstance(sample_properties, wsp.SampleProperties)

    if sample_data is None or sample_data.summary is None:
        return None

    alarms = []

    throughput_inferred = sample_data.summary.get(
        THROUGHPUT_INFERRED_METRIC,
    )

    chemistry = _get_chemistry_description(sample_data, sample_properties)

    chemistry_with_throughput = chemistry
    if throughput_inferred == HT_THROUGHPUT and chemistry[-2:] not in [
        LT_THROUGHPUT,
        MT_THROUGHPUT,
        HT_THROUGHPUT,
    ]:
        chemistry_with_throughput = f"{chemistry} {HT_THROUGHPUT}"

    if (
        chemistry.endswith(HT_THROUGHPUT)
        and throughput_inferred
        and throughput_inferred != HT_THROUGHPUT
    ):
        sample_data.summary[INCONSISTENT_THROUGHPUT_METRIC] = HT_THROUGHPUT

    rows = [
        ["Sample ID", sample_properties.sample_id],
        ["Sample Description", sample_properties.sample_desc],
        ["Chemistry", chemistry_with_throughput],
    ]

    if sample_data.summary.get("spatial_slide_info", None) is not None:
        rows.append(["Slide Serial Number", sample_data.summary["spatial_slide_info"]])

    if pipeline in shared_constants.PIPELINE_COUNT and not sample_properties.is_spatial:
        rows.append(["Include introns", str(sample_properties.include_introns)])

    if isinstance(sample_properties, wsp.ExtendedCountSampleProperties):

        if sample_properties.reference_path:
            rows.append(["Reference Path", sample_properties.reference_path])
        if sample_properties.target_set:
            if sample_data.targeting_method == TARGETING_METHOD_TL:
                rows.append(["Probe Set Name", sample_properties.target_set])
            else:
                rows.append(["Target Panel Name", sample_properties.target_set])
            rows.append(["Number of Genes Targeted", sample_data.summary["num_genes_on_target"]])

        # This was meant to be enabled in 3.1 but due to a bug was not included./
        # if sample_properties.barcode_whitelist:
        #    rows.append(
        #        ['Barcode Whitelist', sample_properties.barcode_whitelist])

    # Find references in the summary
    if (
        isinstance(sample_properties, wsp.AggrCountSampleProperties)
        and not sample_properties.genomes
    ):
        rows.append(
            [
                cr_constants.REFERENCE_TYPE,
                "Not applicable for aggr with feature barcoding-only samples",
            ]
        )
    elif pipeline in [shared_constants.PIPELINE_AGGR, shared_constants.PIPELINE_REANALYZE]:
        genomes = sample_properties.genomes
        if genomes is not None:
            rows.append(
                [cr_constants.REFERENCE_TYPE, cr_reference.get_ref_name_from_genomes(genomes)]
            )
    else:
        reference_metric_prefixes = [
            cr_constants.REFERENCE_METRIC_PREFIX,
            vdj_constants.REFERENCE_METRIC_PREFIX,
        ]
        # Find all references in the summary
        for prefix in reference_metric_prefixes:
            ref_type_key = f"{prefix}{cr_constants.REFERENCE_TYPE_KEY}"
            if ref_type_key in sample_data.summary:
                ref_type = sample_data.summary[ref_type_key]

                ref_version_key = f"{prefix}{cr_constants.REFERENCE_VERSION_KEY}"
                if ref_version_key in sample_data.summary:
                    ref_version = "-%s" % sample_data.summary.get(ref_version_key)
                else:
                    ref_version = ""

                ref_name_key = f"{prefix}{cr_constants.REFERENCE_GENOMES_KEY}"
                if ref_name_key in sample_data.summary:
                    ref_name = sample_data.summary.get(ref_name_key)
                    if isinstance(ref_name, list):
                        ref_name = cr_reference.get_ref_name_from_genomes(ref_name)

                    rows.append([ref_type, f"{ref_name}{ref_version}"])

    # add pipeline version
    rows.append(["Pipeline Version", sample_properties.version])

    # add reorientation mode description to Sample table but not for aggr
    if (
        hasattr(sample_properties, "is_spatial")
        and sample_properties.is_spatial
        and pipeline not in shared_constants.PIPELINE_AGGR
    ):
        if (
            sample_properties.reorientation_mode == "rotation"
            or sample_properties.reorientation_mode == "rotation+mirror"
        ):
            rows.append(["Image Reorientation", "On"])
        else:
            rows.append(["Image Reorientation", "Off"])

    # add filter_probes mode description to the Sample table
    if sample_data.targeting_method == TARGETING_METHOD_TL:
        rows.append(["Filter Probes", "Off" if sample_properties.filter_probes is False else "On"])

    pipeline_info = {
        "header": ["Sample"],
        "rows": rows,
    }
    to_return = {"pipeline_info_table": pipeline_info}
    # We want to alarm if ARC is used in GEX chemistry.
    if metadata is not None and species_list is not None:
        alarms = metadata.gen_metric_list(
            sample_data.summary, [H5_CHEMISTRY_DESC_KEY, "inconsistent_throughput"], species_list
        )
        new_alarms = [metric.alarm_dict for metric in alarms if metric.alarm_dict]
        if new_alarms:
            to_return[ALARMS] = new_alarms
    if _get_chemistry_description(sample_data, sample_properties) == ARC_V1_MULTIOME_DESCRIPTION:
        to_return.setdefault(ALARMS, []).append(
            {
                "formatted_value": None,
                "title": "Unsupported workflow used",
                "message": "Multiome Gene Expression only analysis is not a supported workflow. Results may vary.",
                "level": WARNING_THRESHOLD,
            }
        )

    if (
        isinstance(sample_properties, wsp.CountSampleProperties)
        and sample_properties.include_introns
        and not sample_properties.is_targeted
        and not sample_data.is_antibody_only()
    ):
        if ALARMS not in to_return:
            to_return[ALARMS] = []

        program = "Space Ranger" if sample_properties.is_spatial else "Cell Ranger"

        intron_info_alarm = {
            "formatted_value": None,
            "title": "Intron mode used",
            "message": f"""This data has been analyzed with intronic reads included in the count matrix. This behavior is different from previous {program} versions. If you would not like to count intronic reads, please rerun with the "include-introns" option set to "false". Please contact support@10xgenomics.com for any further questions.""",
            "level": INFO_THRESHOLD,
        }
        to_return[ALARMS].extend([intron_info_alarm])

    # If --aligner was used to force the aligner warn the user
    if (
        isinstance(sample_properties, wsp.CountSampleProperties)
        and sample_properties.aligner
        and sample_properties.is_spatial
    ):
        if ALARMS not in to_return:
            to_return[ALARMS] = []
        aligner_info_alarm = {
            "formatted_value": None,
            "title": "Force Aligner Used",
            "message": f"The --aligner option was used to set the sequencing read aligner to {sample_properties.aligner}. Incorrect usage of this option will lead to unusable data and low mapping metrics. Please contact support@10xgenomics.com for any further questions.",
            "level": INFO_THRESHOLD,
        }
        to_return[ALARMS].extend([aligner_info_alarm])

    # If loupe alignment file contains redundant information: image registration info when only Cytassist image provided
    to_return = _display_loupe_warning(sample_properties, to_return)

    return to_return


def create_table_with_alarms(
    table_key, title, metric_keys, alarm_keys, metadata, sample_data, species_list
):
    """Sequencing info for GEX."""
    if sample_data is None or sample_data.summary is None:
        return None

    data_dict = {}
    # Account for imaging block not having metrics
    if len(metric_keys):
        metrics = metadata.gen_metric_list(sample_data.summary, metric_keys, species_list)
        if metrics:
            # Not all requested metrics will appear, and we only generate help text for those
            # that show up in the table
            observed_keys = {x.parent_metric_info.name for x in metrics}
            filtered_keys = [x for x in metric_keys if x in observed_keys]
            data_dict["help"] = {
                "title": title,
                "data": metadata.gen_metric_helptext(filtered_keys),
            }
            data_dict["table"] = {
                "rows": [[metric.name, metric.value_string] for metric in metrics]
            }
    else:
        data_dict["help"] = {"title": title, "data": []}
        data_dict["table"] = {}
    if not data_dict:
        return None

    result = {table_key: data_dict}

    # Alerts.
    if alarm_keys:
        alarms = metadata.gen_metric_list(sample_data.summary, alarm_keys, species_list)
        # If a metric is from a barnyard and the cumulative version of the metric should be tested
        # we do not check the non-cumulative metrics
        alarms = [
            x
            for x in alarms
            if not (
                x.is_barnyard and x.parent_metric_info.include_cumulative and not x.is_cumulative
            )
        ]
        new_alarms = [metric.alarm_dict for metric in alarms if metric.alarm_dict]

        if new_alarms:
            result[ALARMS] = new_alarms

    return result


def sequencing_table(metadata, sample_data, species_list, is_targeted=False):
    """Sequencing info for GEX."""
    return create_table_with_alarms(
        "sequencing",
        "Sequencing",
        SEQUENCING_METRIC_KEYS if not is_targeted else TARGETED_SEQUENCING_METRIC_KEYS,
        SEQUENCING_ALARM_KEYS,
        metadata,
        sample_data,
        species_list,
    )


def feature_barcode_sequencing_table(metadata, sample_data, species_list, feature_barcode):
    metric_keys = [f"{feature_barcode}_{i}" for i in SEQUENCING_METRIC_KEYS]
    alarm_keys = [f"{feature_barcode}_{i}" for i in SEQUENCING_ALARM_KEYS]

    return create_table_with_alarms(
        f"{feature_barcode.upper()}_sequencing",
        f"{FB_DISPLAY_NAME[feature_barcode]} Sequencing",
        metric_keys,
        alarm_keys,
        metadata,
        sample_data,
        species_list,
    )


def feature_barcode_application_table(metadata, sample_data, species_list, feature_barcode):
    """Feature barcoding application metric."""
    return create_table_with_alarms(
        f"{feature_barcode.upper()}_application",
        f"{FB_DISPLAY_NAME[feature_barcode]} Application",
        FB_APP_METRIC_KEYS[feature_barcode],
        FB_APP_METRIC_KEYS[feature_barcode],
        metadata,
        sample_data,
        species_list,
    )


def mapping_table(metadata, sample_data, species_list):
    """Mapping info table."""
    if sample_data.is_targeted():
        # Post library targeting
        if sample_data.targeting_method == TARGETING_METHOD_HC:
            alarm_keys = TARGETED_MAPPING_ALARMS
        # template ligation targeting
        elif sample_data.targeting_method == TARGETING_METHOD_TL:
            alarm_keys = TEMP_LIG_MAPPING_ALARMS
    else:
        alarm_keys = MAPPING_ALARMS

    return create_table_with_alarms(
        "mapping",
        "Mapping",
        MAPPING_KEYS
        if "targeting_method" not in sample_data.summary
        or sample_data.targeting_method == TARGETING_METHOD_HC
        else TEMP_LIG_MAPPING_KEYS,
        alarm_keys,
        metadata,
        sample_data,
        species_list,
    )


def cell_or_spot_calling_table(
    metadata,
    sample_data,
    sample_properties,
    species_list,
    metric_keys,
    alarm_keys,
):
    """Cell calling data (table and plot)."""
    # TODO: Barnyard not currently in spatial
    is_barnyard = len(species_list) > 1
    if is_barnyard:
        metric_keys.insert(0, FILTERED_BCS_TRANSCRIPTOME_UNION)

    # Replace genes detected metric in "Spots" table if templated ligation assay but only for count
    # not for aggr because the the targeting metric generation has already been taken care of in the count run.
    if sample_data.targeting_method == TARGETING_METHOD_TL and not isinstance(
        sample_properties, wsp.AggrCountSampleProperties
    ):
        RTL_METRIC_TRANSLATION = {
            "filtered_bcs_median_unique_genes_detected": RTL_MEDIAN_GENES_BARCODE,
            "filtered_bcs_total_unique_genes_detected": RTL_GENES_DETECTED,
        }
        metric_keys = [RTL_METRIC_TRANSLATION.get(key, key) for key in metric_keys]
    tbl_name = "Spots" if sample_properties.is_spatial else "Cells"

    table_dict = create_table_with_alarms(
        "cells", tbl_name, metric_keys, alarm_keys, metadata, sample_data, species_list
    )
    if table_dict is None:
        return None

    # add image
    if sample_properties.is_spatial:
        return table_dict
    else:
        # the data we are interested in is in data_dict["cells"]
        data_dict = table_dict["cells"]
        to_return = {}
        # Be sure to bubble up alarms
        if ALARMS in table_dict:
            to_return[ALARMS] = table_dict[ALARMS]
        chart = get_empty_rank_plot()
        knee_plot = cr_webshim.plot_barcode_rank(chart, sample_properties, sample_data)
        if knee_plot:
            data_dict["barcode_knee_plot"] = knee_plot
            data_dict["help"]["data"] = data_dict["help"]["data"] + RANK_PLOT_HELP
        to_return["cells"] = data_dict
        return to_return


def summary_image_table(
    metadata,
    sample_data,
    species_list,
    metric_keys,
    alarm_keys,
    zoom_images,
    regist_images=None,
):
    """Image block for spatial samples."""
    tbl_name = "Image"
    # Keeping this core function just in case we want to throw image based alarms in the future
    table_dict = create_table_with_alarms(
        "image", tbl_name, metric_keys, alarm_keys, metadata, sample_data, species_list
    )
    if table_dict is None:
        return None
    # add tissue detection image
    table_dict["image"]["zoom_images"] = zoom_images
    # add tissue detection help
    table_dict["image"]["help"]["data"] = ZOOM_IMAGE_HELP
    if regist_images:
        # add cytassist tissue registration image
        table_dict["image"]["regist_images"] = regist_images
        # add cytassist tissue registration help
        table_dict["image"]["help"]["data"].append(REGISTRATION_IMAGE_HELP)
    return table_dict


def batch_correction_table(metadata, sample_data, species_list):
    metric_keys = [
        "batch_effect_score_before_correction",
        "batch_effect_score_after_correction",
    ]

    return create_table_with_alarms(
        "batch_correction",
        "Chemistry Batch Correction",
        metric_keys,
        None,
        metadata,
        sample_data,
        species_list,
    )


def aggregation_table(metadata, sample_data, sample_properties):
    """Report normalization metrics in aggr.

    The trick here is to use the batch prefix as the
    a species/genome prefix, and define these metric as species_specific.

    TODO: the above trick doesn't generate good web summary if it's a barnyard aggr sample.
    """
    if not isinstance(sample_properties, wsp.AggrCountSampleProperties):
        return None
    metric_keys = [
        "pre_normalization_total_reads",
        "post_normalization_total_reads",
        "pre_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc",
        "post_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc",
    ] + AGGREGATION_METRIC_KEYS

    alarm_keys = [
        "lowest_frac_reads_kept",
    ]
    batches = sample_properties.agg_batches
    return create_table_with_alarms(
        "aggregation", "Aggregation", metric_keys, alarm_keys, metadata, sample_data, batches
    )


def feature_barcode_aggregation_table(metadata, sample_data, sample_properties, feature_barcode):
    if not isinstance(sample_properties, wsp.AggrCountSampleProperties):
        return None
    metric_keys = [f"{feature_barcode}_{i}" for i in AGGREGATION_METRIC_KEYS]
    batches = sample_properties.agg_batches
    return create_table_with_alarms(
        f"{feature_barcode}_aggregation",
        f"{FB_DISPLAY_NAME[feature_barcode]} Aggregation",
        metric_keys,
        None,
        metadata,
        sample_data,
        batches,
    )


def _display_loupe_warning(sample_properties, to_return):
    """Display warning to websummary.

    if loupe file contains redundant information
    """
    if (
        isinstance(sample_properties, wsp.CountSampleProperties)
        and sample_properties.redundant_loupe_alignment
    ):
        if ALARMS not in to_return:
            to_return[ALARMS] = []
        loupe_info_alarm = {
            "formatted_value": None,
            "title": "Unused Loupe alignment info",
            "message": "Loupe alignment file contains unused tissue registration info because the microscope image used to produce the alignment file was not input to Space Ranger.\
                 If this was intentional, then you can proceed with analyzing your data. Please contact support@10xgenomics.com for any further questions.",
            "level": INFO_THRESHOLD,
        }
        to_return[ALARMS].extend([loupe_info_alarm])
    return to_return
