include {DRIVERCSV} from "./DriverG4.nf"
include {SAMPLESHEET} from "../workflows/SampleSheetAviti.nf"
include {driver_make; bowtie2_fastqc; bam_stats; flagstat; copy_fastqs; multiQC; sample_report_make; multiQC_report_provide; sgdemux_demultiplex} from "../modules/G4_core_RNAseq.nf"
include {combined_report_generate} from "../modules/Aviti_core_10X.nf" 
include {copy_fastq_10X_G4; run_cellranger_G4 } from "../modules/G4_core_10X.nf"

workflow G4_RNASEQ {
    take:
    lims_info_csv

    main:
    // Check orgin of fastqs
    if (!params.fastq_dir) {
        sgdemux_demultiplex(lims_info_csv)
        fastq_dir = sgdemux_demultiplex.out
    } else {
        fastq_dir = params.fastq_dir 
    }
    //
    driver_csv = driver_make(lims_info_csv,fastq_dir)
    driver = DRIVERCSV(driver_csv)
    bowtie2_fastqc(driver)
    copy_fastqs(driver)
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

workflow G4_10X { 
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
    copy_fastq_10X_G4(meta, fastq_dir)
    run_cellranger_G4(meta, lims_info_csv, fastq_dir)

    def all_cellranger_outs = run_cellranger_G4.out
                                            .collect()
                                            .map{it[0,1,2]}

    combined_report_generate(all_cellranger_outs)
    

    emit:
    combined_report_generate.out
}