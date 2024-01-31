include {run_base2fastq; namechange_and_copy_fastqs; run_cellranger; combined_report_generate} from "../modules/Aviti_core_10X.nf"
include {SAMPLESHEET} from "../workflows/SampleSheetAviti.nf"
include {DRIVERCSV} from "./DriverG4.nf"
include {driver_make; bowtie2_fastqc; bam_stats; flagstat; copy_fastqs; multiQC; sample_report_make; multiQC_report_provide} from "../modules/G4_core_RNAseq.nf"
include {bowtie2_fastqc_Aviti; driver_make_Aviti} from "../modules/Aviti_core_RNAseq.nf"

workflow AVITI_10X {
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
    namechange_and_copy_fastqs(meta, fastq_dir)
    run_cellranger(namechange_and_copy_fastqs.out, lims_info_csv)

    def all_cellranger_outs = run_cellranger.out
                                            .collect()
                                            .map{it[0,1,2]}

    combined_report_generate(all_cellranger_outs)

    emit:
    combined_report_generate.out
}

workflow AVITI_RNASEQ { 
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
    namechange_and_copy_fastqs(meta, fastq_dir)
    def all_nameChange_outs = namechange_and_copy_fastqs.out
                                                        .collect()
                                                        .map {it[0]}
    driver_csv = driver_make_Aviti(lims_info_csv, fastq_dir, all_nameChange_outs) 
    driver = DRIVERCSV(driver_csv)
    bowtie2_fastqc(driver)
    bam_stats(bowtie2_fastqc.out)
    flagstat(bowtie2_fastqc.out)
    
    def all_bowtieOuts = bowtie2_fastqc.out
                                    .collect()
                                    .map {it[0]}

    multiQC(all_bowtieOuts, driver_csv)
    sample_report_make(driver_csv, multiQC.out)
    multiQC_report_provide(sample_report_make.out)

    emit:
    multiQC_report_provide.out

}
