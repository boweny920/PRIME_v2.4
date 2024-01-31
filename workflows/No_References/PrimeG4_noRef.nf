include {SAMPLESHEET} from "../../workflows/SampleSheetAviti.nf"
include {DRIVERCSV} from "../DriverG4.nf"
include {copy_fastq_10X_G4} from "../../modules/G4_core_10X.nf"
include {sgdemux_demultiplex; driver_make} from "../../modules/G4_core_RNAseq.nf"
include {no_Ref_DriverMake; no_Ref_SampleReport_Make; no_Ref_fastqc; no_Ref_multiQC} from "../../modules/no_Ref_core.nf"

workflow G4_noRef {
    take:
    lims_info_csv

    main:
    meta = SAMPLESHEET(lims_info_csv)
     // Check orgin of fastqs
    if (!params.fastq_dir) {
        sgdemux_demultiplex(lims_info_csv)
        fastq_dir = sgdemux_demultiplex.out
    } else {
        fastq_dir = params.fastq_dir 
    }
    //
    driver_csv = no_Ref_DriverMake(lims_info_csv, fastq_dir)
    no_Ref_SampleReport_Make(driver_csv)
    copy_fastq_10X_G4(meta, fastq_dir)

    // Add QC process to no_ref orders

    driver = DRIVERCSV(driver_csv)

    no_Ref_fastqc(driver)

    def fastQcOut_1stDriver = no_Ref_fastqc.out
                                    .collect()
                                    .map {it[0]}

    no_Ref_multiQC(fastQcOut_1stDriver, driver_csv, "placeHolder")
    
    emit:
    no_Ref_multiQC.out

}