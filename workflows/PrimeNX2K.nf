include {DRIVERCSV} from "./DriverG4.nf"
include {SAMPLESHEET} from "../workflows/SampleSheetAviti.nf"
include {bowtie2_fastqc; bam_stats; flagstat; copy_fastqs; multiQC; sample_report_make; multiQC_report_provide} from "../modules/G4_core_RNAseq.nf"
include {combined_report_generate} from "../modules/Aviti_core_10X.nf" 
include {copy_fastq_10X_G4; run_cellranger_G4} from "../modules/G4_core_10X.nf"
include {nx2k_samplesheet_create; bclconvert; nx2k_driver_make; nx2k_cellranger_mkfastq} from "../modules/NX2K_core.nf"

workflow NX2K_RNASEQ {
    take:
    lims_info_csv

    main:
    // use AAC7NV3HV , /n/ngs/data/NextSeq2K/230919_VH00629_133_AAC7NV3HV as test 
    nx2k_samplesheet_create(lims_info_csv)
    // Check Fastq origin
    if (!params.fastq_dir) {
         bclconvert(nx2k_samplesheet_create.out)
        fastq_dir = bclconvert.out
    } else {
        fastq_dir = params.fastq_dir
    }
    //
    driver_csv = nx2k_driver_make(lims_info_csv, fastq_dir)
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

workflow NX2K_10X {
    take:
    lims_info_csv

    main:
    // use AAAYVMKHV , /n/ngs/data/NextSeq2K/230809_VH00629_123_AAAYVMKHV as test 
    nx2k_samplesheet_create(lims_info_csv)
    meta = SAMPLESHEET(lims_info_csv)
    // Check Fastq orgin
    if (!params.fastq_dir) { 
        fastq_dir = nx2k_cellranger_mkfastq(nx2k_samplesheet_create.out)
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