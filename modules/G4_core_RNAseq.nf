process driver_make {

    label 'small_mem'

    input:
    path lims_table_csv
    val fastq_dir

    publishDir "${params.fcpath}" , mode: 'copy'
    
    output:
    path "PRIME_driver.csv"

    script:
    """
    G4_Driver_Make.py --lims_csv ${lims_table_csv} --fastq_path ${fastq_dir} --RoboIndex_samplesheet ${params.samplesheet_RoboIndex}
    """
}

process copy_fastqs {
    
    label 'small_mem'

    input:
    val driver

    output: 
    val "copy complete"

    script:
    file("${driver.nanalysis_path}").mkdirs()

    """
    cp ${driver.fastq} ${driver.nanalysis_path}/
    """
}

process bowtie2_fastqc {
    // No need to worry about slum resources, i have specified them in the config just for this process 
    
    // conda "${projectDir}/assets/fqc_bt2_condaEnv.yml"
    conda "/home/by2747/miniconda3/envs/fqc_bt2"

    input:
    val driver

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

process bam_stats {

    label 'small_mem'

    input:
    tuple val(driver), path(bt2Log), path(bt2Bam), val(fastqc_outs) //bt2_fqc output will be the input of this process 

    output:
    path "${driver.file_name}_bamstats.txt" 

    publishDir "${params.fcpath}/PRIME_output/mapped" , mode: 'copy'

    script:
    """
    bamtools stats -in ${bt2Bam} > ${driver.file_name}_bamstats.txt
    """
}

process flagstat {

    label 'small_mem'

    input:
    tuple val(driver), path(bt2Log), path(bt2Bam), val(fastqc_outs) //bt2_fqc output will be the input of this process  

    output:
    path "${driver.file_name}_flagstat.tsv" 
    
    publishDir "${params.fcpath}/PRIME_output/mapped" , mode: 'copy'

    script:
    """
    samtools flagstat ${bt2Bam} -@ 4 > ${driver.file_name}_flagstat.tsv
    """
}

process multiQC {

    label 'small_mem'
    // maxForks 10
    errorStrategy 'retry'
    maxRetries 5

    input:
    // tuple val(driver), path(bt2Log), path(bt2Bam), val(fastqc_outs) //bt2_fqc output will be the input of this process
    val all_bt_outputs
    path driver_csv

    output:
    tuple path("multiqc.html"), path("multiqc_data"), val(all_bt_outputs)

    // publishDir "${params.fcpath}/PRIME_output/qc/" , mode: 'copy' 

    script:
    
    """
    multiqc --config ${projectDir}/assets/multiqc_config.yaml\
            --force -o ./ -n multiqc.html ${params.fcpath}/PRIME_output/qc/fastqc ${params.fcpath}/PRIME_output/mapped/*.log\
            --quiet

    copy_multiQC.py -d ${driver_csv} -m multiqc.html 

    cp -f multiqc.html ${params.fcpath}/PRIME_output/qc/
    cp -r -f multiqc_data ${params.fcpath}/PRIME_output/qc/
    """
}

process sample_report_make {

    label 'small_mem'
    // maxForks 10
    errorStrategy 'retry'
    maxRetries 5

    input:
    path driver_csv 
    tuple path(multiQC_html), path(multiqc_data_folder), val(driver)

    output:
    tuple path("Prime_Sample_Report.csv"), val("https://webfs/${driver.nanalysis_path}/multiqc.html")

    publishDir "${params.fcpath}/" , mode: 'copy' 

    script:
    """
    G4_RNAseqSampleReport_Make.py --driver_csv ${driver_csv} --multiqc_bowtie2_txt ${multiqc_data_folder}/multiqc_bowtie2.txt 

    """
}

process multiQC_report_provide {

    label 'small_mem'
    
    input:
    tuple path(sampleReport), val(multiQC_link)

    output:
    val multiQC_link

    script:
    """
    echo ${multiQC_link}
    """
}

process sgdemux_demultiplex {

    input: 
    path lims_info_csv

    output:
    val "${params.fcpath}/PRIME_demux_filtered_fastqs/MOLNG-*"

    script:

    if (params.run_type.contains("RNA-Seq")) {
    """
    mkdir -p ${params.fcpath}/PRIME_demux_filtered_fastqs

    sgdemux --fastqs ${params.fcpath}/unfiltered_fastqs/Undetermined --sample-metadata ${params.fcpath}/samplesheet.csv --output-dir ${params.fcpath}/PRIME_demux_filtered_fastqs
    """
    }
    else if (params.run_type.contains("10X")) {
    """
    sgdemux --fastqs ${params.fcpath}/unfiltered_fastqs/Undetermined --sample-metadata ${params.fcpath}/samplesheet.csv --output-dir ${params.fcpath}/PRIME_demux_filtered_fastqs --output-types TB	
    """    
    }
}