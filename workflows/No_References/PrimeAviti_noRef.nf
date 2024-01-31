include {run_base2fastq; namechange_and_copy_fastqs_noRef} from "../../modules/Aviti_core_10X.nf"
include {driver_make_Aviti} from "../../modules/Aviti_core_RNAseq.nf"
include {SAMPLESHEET} from "../../workflows/SampleSheetAviti.nf"
include {DRIVERCSV} from "../DriverG4.nf"
include {no_Ref_DriverMake; no_Ref_SampleReport_Make; no_Ref_fastqc; no_Ref_multiQC} from "../../modules/no_Ref_core.nf"

workflow AVITI_noRef {
    take:
    lims_info_csv

    main:
    meta = SAMPLESHEET(lims_info_csv)
   // Determine orgin of the fastq directory
    if (!params.fastq_dir) {
        run_base2fastq(lims_info_csv)
        fastq_dir = run_base2fastq.out
    } else {
        fastq_dir = params.fastq_dir
    }
    //  
    driver_csv = no_Ref_DriverMake(lims_info_csv,fastq_dir)
    no_Ref_SampleReport_Make(driver_csv) 
    
    
    // Add QC process to no_ref orders
    driver = DRIVERCSV(driver_csv)
    no_Ref_fastqc(driver)

    def fastQcOut_1stDriver = no_Ref_fastqc.out
                                    .collect()
                                    .map {it[0]}

    
    
    //moved name changing step to the end becasue the name change will mess up the information for the driver.csv
    namechange_and_copy_fastqs_noRef(meta, fastq_dir) 
    def all_nameChange_outs_placeHolder = namechange_and_copy_fastqs_noRef.out
                                                        .collect()
                                                        .map {it[0]}

    no_Ref_multiQC(fastQcOut_1stDriver, driver_csv, all_nameChange_outs_placeHolder)

    emit:
    no_Ref_multiQC.out

}