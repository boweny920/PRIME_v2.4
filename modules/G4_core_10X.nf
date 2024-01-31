process copy_fastq_10X_G4 {
    label 'small_mem'

    input:
    val meta
    val fastq_dir

    output:
    val "copy complete"

    script:
    file("${meta.nanalysis_path}").mkdirs()

    """
    cp ${fastq_dir}/${meta.libID}*.fastq.gz ${meta.nanalysis_path}/
    """
}

process run_cellranger_G4{

    label 'big_mem'

    input:
    val meta
    path lims_info_csv
    path fastq_dir
    

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
    
    cp -f ${meta.libID}.config.csv ${params.fcpath}

    cellranger multi --id ${meta.libID} --csv ${meta.libID}.config.csv
    
    cp -r ${meta.libID}/outs ./

    """
    else
    """
    cellranger count --id ${meta.libID}-CellrangerOut --nopreflight --fastqs=${fastq_dir} --sample=${meta.libID} --transcriptome=${params.indexDir}/${index_genome}/annotation/${annotation}/10x/cellranger --chemistry=auto --disable-ui 

    cp -r -f ${meta.libID}-CellrangerOut ${params.fcpath}

    cp -r ${meta.libID}-CellrangerOut/outs ./

    """
}