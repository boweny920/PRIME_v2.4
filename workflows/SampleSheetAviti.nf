def create_metadata_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.libID = row.libID
    meta.refGenome = row.genomeVersion
    meta.species = row.speciesName
    meta.orderType = row.orderType
    // meta.laneID = row.laneId
    meta.nanalysis_path = row.resultPaths
    meta.pi_name = row.requestingDepartment
    meta.molngID = row.prnOrderNo
    meta.requester_name = row.requester
    meta.annotation = row.annotation

    return meta
    }

process samplesheet_make {

    label 'small_mem'

    input:
    path lims_table_csv

    output:
    path "samplesheet.csv"

    script:
    """
    Aviti_samplesheet_Make.py --lims_info_table ${lims_table_csv} --RoboIndex_samplesheet ${params.samplesheet_RoboIndex}
    """
}

workflow SAMPLESHEET {
    take:
    lims_info_csv

    main:
    samplesheet_make(lims_info_csv)
    samplesheet_make.out.splitCsv ( header:true, sep:',' )
        .map { create_metadata_channel(it) }
        .set { data_meta }
    
    emit:
    data_meta
}
