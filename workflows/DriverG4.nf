include {driver_make} from "../modules/G4_core_RNAseq.nf"

def create_driver_channel(LinkedHashMap row) {
    // create meta map
    def driver = [:]
    driver.fastq = row.Output
    driver.refGenome = row.Reference
    driver.refSpecies = row.Species
    driver.nanalysis_path = row.resultPaths
    driver.molngID = row.Order
    driver.file_name = row.fileNames
    driver.annotation = row.annotation

    return driver
    }

workflow DRIVERCSV {
    take:
    driver_csv

    main:
    driver_csv.splitCsv ( header:true, sep:',' )
        .map { create_driver_channel(it) }
        .set { data_driver }
    
    emit:
    data_driver
}
