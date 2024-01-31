process bowtie2_fastqc_Aviti {
    // No need to worry about slum resources, i have specified them in the config just for this process 
    
    conda "${projectDir}/assets/fqc_bt2_condaEnv.yml"

    input:
    tuple val(fastq_dir), val(driver)

    output:
    tuple val(driver), path("${driver.file_name}.log"), path("${driver.file_name}.bam"), path("*_fastqc.*")

    publishDir "${params.fcpath}/PRIME_output/mapped" , pattern: '*.log',mode: 'copy'
    publishDir "${params.fcpath}/PRIME_output/mapped" , pattern: '*.bam',mode: 'copy'
    publishDir "${params.fcpath}/PRIME_output/qc/fastqc" , pattern: '*_fastqc.*',mode: 'copy'


    script:

    file("${params.fcpath}/PRIME_output/mapped").mkdirs()
    file("${params.fcpath}/PRIME_output/qc/fastqc").mkdirs()
    file("${driver.nanalysis_path}/fastqc").mkdirs()

    """
    bowtie2 --threads 5 -x ${params.indexDir}/${driver.refSpecies}/${driver.refGenome}/bowtie2/${driver.refGenome} -U ${driver.fastq} 2> ${driver.file_name}.log | samtools view -Sbh -o ${driver.file_name}.bam &

    fastqc -t 10 ${driver.fastq} -o ./  --quiet &

    wait

    cp *_fastqc.* ${driver.nanalysis_path}/fastqc 

    """
}

process driver_make_Aviti {

    label 'small_mem'

    input:
    path lims_table_csv
    val fastq_dir
    val place_holder_for_nameChange_toFinish

    publishDir "${params.fcpath}" , mode: 'copy'
    
    output:
    path "PRIME_driver.csv"

    script:

    """
    G4_Driver_Make.py --lims_csv ${lims_table_csv} --fastq_path '${fastq_dir}/fastqs' --RoboIndex_samplesheet ${params.samplesheet_RoboIndex}
    """
    
}
