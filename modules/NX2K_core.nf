process nx2k_samplesheet_create {
    label 'small_mem'

    input:
    path lims_info_csv

    output:
    path "SampleSheet.csv"

    publishDir "${params.fcpath}" , mode: 'copy'

    script:

    """
    NextSeq2K_bcl_SampleSheet_V2.py --lims_csv ${lims_info_csv} --run_type ${params.run_type} --novaSeq ${params.novaseq}
    """
}

process bclconvert {

    input:
    path nx2k_samplesheet
    
    output:
    val "${params.fcpath}/PrimeNX2K_fastqs"

    // publishDir "${params.fcpath}" , mode: 'copy' 

    script:
    """
    bcl-convert-3.10.5 --bcl-input-dir ${params.fcpath} --output-dir ${params.fcpath}/PrimeNX2K_fastqs -f  
    """

}

process nx2k_driver_make {

    label 'small_mem'

    input:
    path lims_table_csv
    val fastqs

    publishDir "${params.fcpath}" , mode: 'copy'
    
    output:
    path "PRIME_driver.csv"

    script:
    """
    G4_Driver_Make.py --lims_csv ${lims_table_csv} --fastq_path ${fastqs} --RoboIndex_samplesheet ${params.samplesheet_RoboIndex}
    """
}

process nx2k_cellranger_mkfastq {

    label 'big_mem'

    input:
    path nx2k_samplesheet
    
    output:
    val "${params.fcpath}/PrimeNX2K_fastqs/${params.fcid}"

    script:
    """
    cellranger mkfastq --run=${params.fcpath} --csv=${nx2k_samplesheet} --nopreflight --output-dir ${params.fcpath}/PrimeNX2K_fastqs 
    """
}