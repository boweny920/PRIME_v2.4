process lims_info_process {
    
    label 'small_mem'

    input:
    val fcid 

    output:
    path "lims_info.csv"

    publishDir "${params.fcpath}", mode: 'copy'

    script:

    //sendMail(to: 'by2747@stowers.org', subject: "${fcid} PrimeAG started", body: "${fcid} PrimeAG WorkDir: ${workflow.workDir}")

    if (!params.molng)
    """
    Aviti_lims_json_Make_V2.py --fcid ${fcid} --samplesheet_Robo ${params.samplesheet_RoboIndex}
    """
    else
    """
    Aviti_lims_json_Make_V2.py --fcid ${fcid} --samplesheet_Robo ${params.samplesheet_RoboIndex} --molng ${params.molng} 
    """
}

process run_base2fastq {

    label 'big_mem'

    input:
    path lims_info_csv

    output:
    val "${params.fcpath}/base2fastq_outs" 

    script:
    """
    bases2fastq ${params.fcpath} ${params.fcpath}/base2fastq_outs 
    """
}

process namechange_and_copy_fastqs {
    
    errorStrategy 'retry'
    maxRetries 5
    label 'lil_mem'

    input:
    val meta
    val demutiplex_dir

    output:
    tuple val("${fastq_path}"), val(meta) //fastq path, meta info

    script:

    if (!params.only_cp_fastq) {
        fastq_path = "${demutiplex_dir}/fastqs"        
    } else {
        fastq_path = "${demutiplex_dir}"       
    }
    
    """
    mkdir -p ${demutiplex_dir}/fastqs
    find ${demutiplex_dir}/Samples -name "*.fastq.gz" -exec mv {} ${demutiplex_dir}/fastqs \\;
    Aviti_nameChange_copyFiles.py --fastq_dir ${fastq_path} --libID ${meta.libID} --nanalysis_dir ${meta.nanalysis_path}
    """
}

process run_cellranger {

    label 'big_mem'

    input:
    tuple val(fastq_dir), val(meta)
    path lims_info_csv

    output:
    tuple val(secondary_path), val(meta.orderType), path("outs")

    publishDir "${output_lib_folder}", mode: 'copy'

    script:
    // If user defined a genome to use, use that genome as reference. Otherwise, use the LIMs fetched reference genome as reference
    
    if (!params.annotation) {
        annotation = meta.annotation
    } else {
        annotation = params.annotation 
    }

    if (!params.genome){
        output_lib_folder = file("${params.outdir}/${meta.pi_name}/${meta.requester_name}/${meta.molngID}.${meta.refGenome}.${annotation}/${meta.libID}/")
        index_genome = "${meta.species}/${meta.refGenome}"
        secondary_path = "${params.outdir}/${meta.pi_name}/${meta.requester_name}/${meta.molngID}.${meta.refGenome}.${annotation}"
    } else {
        output_lib_folder = file("${params.outdir}/${meta.pi_name}/${meta.requester_name}/${meta.molngID}.${params.genome}.${annotation}/${meta.libID}/")
        index_genome = params.genome
        secondary_path = "${params.outdir}/${meta.pi_name}/${meta.requester_name}/${meta.molngID}.${params.genome}.${annotation}"
    }
    output_lib_folder.mkdirs()
    
    if (meta.orderType.contains("10x_RNA_Flex")) //This will be the cellranger FLEX branch, have to make a config csv too
    """
    Aviti_FLEXscRNAseq_configcsv_Make.py --lims_info_csv ${lims_info_csv} --libID ${meta.libID} --fastqDir ${fastq_dir}
    
    cp ${meta.libID}.config.csv ${params.fcpath}

    cellranger multi --id ${meta.libID} --csv ${meta.libID}.config.csv
    
    cp -r ${meta.libID}/outs ./

    """
    else
    """
    cellranger count --id ${meta.libID}-CellrangerOut --nopreflight --fastqs=${fastq_dir} --sample=${meta.libID} --transcriptome=${params.indexDir}/${index_genome}/annotation/${annotation}/10x/cellranger --chemistry=auto --disable-ui 

    cp -r ${meta.libID}-CellrangerOut ${params.fcpath}

    cp -r ${meta.libID}-CellrangerOut/outs ./

    """
}

process combined_report_generate {

    label 'small_mem'

    input:
    tuple val(secondary_path), val(orderType), path(cellranger_outs)

    output:
    val "https://webfs${secondary_path}/analysis_10xscRNASeq.html"

    script:
    if (orderType.contains("10x_RNA_Flex"))
    """
    Aviti_FLEXscRNAseqReport.py --secondary_path ${secondary_path}
    """
    else
    """
    Aviti_10xscRNAseqReport.py --secondary_path ${secondary_path}
    """
}

process namechange_and_copy_fastqs_noRef {

    label 'lil_mem'

    input:
    val meta
    val demutiplex_dir

    output:
    val("No reference genome in index collection for this order. Fastqs copied to: ${meta.nanalysis_path}") //fastq path, meta info

    script:
    if (!params.only_cp_fastq) {
        """
        Aviti_nameChange_copyFiles.py --fastq_dir ${demutiplex_dir}/Samples/${meta.libID} --libID ${meta.libID} --nanalysis_dir ${meta.nanalysis_path}
        """
    } else {
        """
        Aviti_nameChange_copyFiles.py --fastq_dir ${demutiplex_dir}/ --libID ${meta.libID} --nanalysis_dir ${meta.nanalysis_path}
        """
    }
}