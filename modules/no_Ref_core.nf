process no_Ref_DriverMake {
    
    label 'small_mem'

    input:
    path lims_table_csv
    val fastq_dir

    publishDir "${params.fcpath}/" , mode: 'copy'
    
    output:
    path "PRIME_driver.csv"

    script:

    if (params.machine_type.contains("Aviti")) {

        if (!params.only_cp_fastq) {
        """
        no_Ref_Driver_Make.py --lims_csv ${lims_table_csv} --fastq_path ${fastq_dir} --aviti True
        """
        } else {
            """
            no_Ref_Driver_Make.py --lims_csv ${lims_table_csv} --fastq_path ${fastq_dir}
            """ 
        }

    } else {
        """
        no_Ref_Driver_Make.py --lims_csv ${lims_table_csv} --fastq_path ${fastq_dir}
        """
    }
    
}

process no_Ref_SampleReport_Make {

    label 'small_mem'
    errorStrategy 'retry'
    maxRetries 5
    
    input:
    path driver_csv 

    output:
    path "Prime_Sample_Report.csv"

    publishDir "${params.fcpath}/" , mode: 'copy' 

    script:
    """
    no_Ref_SampleReport_Make.py --driver_csv ${driver_csv}
    """
}

process no_Ref_fastqc {
    // No need to worry about slum resources, i have specified them in the config just for this process 
    
    // conda "${projectDir}/assets/fqc_bt2_condaEnv.yml"
    conda "/home/by2747/miniconda3/envs/fqc_bt2"

    input:
    val driver

    output:
    tuple val(driver), path("*_fastqc.*")

    publishDir "${params.fcpath}/PRIME_output/qc/fastqc" , pattern: '*_fastqc.*',mode: 'copy'


    script:

    file("${params.fcpath}/PRIME_output/qc/fastqc").mkdirs()
    file("${driver.nanalysis_path}/fastqc").mkdirs()

    """

    fastqc -t 10 ${driver.fastq} -o ./  --quiet

    cp *_fastqc.* ${driver.nanalysis_path}/fastqc 

    """
}

process no_Ref_multiQC {

    label 'small_mem'
    // maxForks 10
    errorStrategy 'retry'
    maxRetries 5

    input:
    val fastqc_out1st_driver
    path driver_csv
    val placeholder

    output:
    val "https://webfs/${fastqc_out1st_driver.nanalysis_path}/multiqc.html"

    // publishDir "${params.fcpath}/PRIME_output/qc/" , mode: 'copy' 

    script:
    
    """
    multiqc --config ${projectDir}/assets/multiqc_config.yaml\
            --force -o ./ -n multiqc.html ${params.fcpath}/PRIME_output/qc/fastqc \
            --quiet

    copy_multiQC.py -d ${driver_csv} -m multiqc.html 

    cp -f multiqc.html ${params.fcpath}/PRIME_output/qc/
    cp -r -f multiqc_data ${params.fcpath}/PRIME_output/qc/
    """
}