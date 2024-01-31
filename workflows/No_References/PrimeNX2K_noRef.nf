include {DRIVERCSV} from "../DriverG4.nf"
include {SAMPLESHEET} from "../../workflows/SampleSheetAviti.nf"
include {bowtie2_fastqc; bam_stats; flagstat; copy_fastqs; multiQC; sample_report_make; multiQC_report_provide} from "../../modules/G4_core_RNAseq.nf"
include {nx2k_samplesheet_create; bclconvert; nx2k_driver_make; nx2k_cellranger_mkfastq} from "../../modules/NX2K_core.nf"
include {no_Ref_DriverMake; no_Ref_SampleReport_Make; no_Ref_fastqc; no_Ref_multiQC} from "../../modules/no_Ref_core.nf"

workflow NX2K_noRef {
    take:
    lims_info_csv

    main:
    nx2k_samplesheet_create(lims_info_csv)
    // Check Fastq origin
    if (!params.fastq_dir) {
        bclconvert(nx2k_samplesheet_create.out)
        fastq_dir = bclconvert.out
    } else {
        fastq_dir = params.fastq_dir
    }
    //
    // driver_csv = nx2k_driver_make(lims_info_csv, fastq_dir)
    driver_csv = no_Ref_DriverMake(lims_info_csv, fastq_dir)
    no_Ref_SampleReport_Make(driver_csv)
    driver = DRIVERCSV(driver_csv)
    copy_fastqs(driver)

    // Add QC process to no_ref orders
    no_Ref_fastqc(driver)

    def fastQcOut_1stDriver = no_Ref_fastqc.out
                                    .collect()
                                    .map {it[0]}

    no_Ref_multiQC(fastQcOut_1stDriver, driver_csv, "placeHolder")
    
    emit:
    no_Ref_multiQC.out
}
