#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// Logging functions
def log(message) {
    log_file.withWriterAppend { writer ->
        writer.writeLine("${message}")
    }
}

def tab_log(message) {
    log_file.withWriterAppend { writer ->
        writer.writeLine("\t- ${message}")
    }
}

def date_log(message) {
    def timestamp = new java.text.SimpleDateFormat('yyyy-MM-dd HH:mm:ss').format(new java.util.Date())
    log_file.withWriterAppend { writer ->
        writer.writeLine("[${timestamp}] ${message}")
    }
}
def cmd_args = workflow.commandLine

// SNPRS Main Script
// Params are read in from command line or from nextflow.config and/or conf/profiles.config

timestamp = "${params.timestamp}"

// Base SNPRS Directory
if("${params.out}"==""){
    error "Must provide base directory for output (--out)..."
} else{
    snprs_directory = file(params.out)
    if(!snprs_directory.isDirectory()){
        if(!snprs_directory.getParent().isDirectory()){
            error "Parent directory for output is not a valid directory [${snprs_directory.getParent()}]..."
        } else{
            snprs_directory.mkdirs()
            new_dir = true 
        }
    } else{
        new_dir = false
    }
}

// Log File
log_directory = file("${snprs_directory}/Run_Logs")
if(!log_directory.isDirectory()){
    log_directory.mkdirs()
}

log_file = file("${log_directory}/SNPRS_Log_${timestamp}.txt")

if (log_file.exists()) {
    error "Log file ${log_file} already exists?"
} else {
    log("SNPRS Log File")
    log("${new java.text.SimpleDateFormat('yyyy-MM-dd HH:mm:ss').format(new java.util.Date())}\n")
    log("Command: ${cmd_args}\n")

    if(new_dir){
        tab_log("Created output directory: ${snprs_directory}")
    } else{
        tab_log("Found output directory: ${snprs_directory}")
    }
}

// Pangenome output infromation
if(params.pg_name == ""){
    pg_name = "SNPRS_${params.timestamp}"
    params.pg_name = "${pg_name}"
} else{
    pg_name = "${params.pg_name}"
}

// Pangenome Directory
pangenome_directory = file("${snprs_directory}/SNPRS_Pangenomes")
if(!pangenome_directory.isDirectory()){
    pangenome_directory.mkdirs() 
    tab_log("Created pangenome directory: ${pangenome_directory}...")
} else{
    tab_log("Found pangenome directory: ${pangenome_directory}...")
}

// Pangenome Reads
if("${params.pg_reads}" == ""){
    pg_read_directory = ""
} else{
    pg_read_directory = file("${params.pg_reads}")
}

// Mapping Directory
mapping_directory = file("${snprs_directory}/Mapping")
if(!mapping_directory.isDirectory()){
    mapping_directory.mkdirs() 
    tab_log("Created mapping directory: ${mapping_directory}...")
} else{
    tab_log("Found mapping directory: ${mapping_directory}...")
}

// Parameterize major directories
params.snprs_directory = file(snprs_directory)
params.pangenome_directory = file(pangenome_directory)
params.mapping_directory = file(mapping_directory)
params.log_directory = file(log_directory)
params.log_file = file(log_file)

include {makePangenome} from "./subworkflows/make_pangenome/main.nf"
include {indexGenome} from "./subworkflows/make_pangenome/main.nf"
include {checkGenome} from "./subworkflows/make_pangenome/main.nf"
include {mapReads} from "./subworkflows/mapping/main.nf"
include {fetchBAM} from "./subworkflows/mapping/main.nf"
include {fetchRawParquet} from "./subworkflows/convert_bam/main.nf"
include {bamToParquet} from "./subworkflows/convert_bam/main.nf"
include {callBases} from "./subworkflows/call_bases/main.nf"
include {fetchCalledBases} from "./subworkflows/call_bases/main.nf"

workflow{

    // Assemble pangenome from reads
    if("${params.pg_reads}" != ""){
        pangenome_info = makePangenome(pangenome_directory,pg_name,pg_read_directory).first()
    } 
    // Get FASTA from --fasta (creates fai/ref in needed)
    else if("${params.fasta}" != ""){
        pangenome_info = indexGenome(params.fasta).first()
    }
    // Specify pangenome by name (checks for fasta, fai, and ref in SNPRS_Pangenomes/PG_NAME)
    else if("${pg_name}" != "SNPRS_${params.timestamp}"){
        check_dir = file("${pangenome_directory}/${pg_name}")
        pangenome_info = checkGenome(check_dir).first()
    }
    // No pangenome information provided
    else {
        pangenome_info = Channel.empty()
    }

    // Get BAM files
    new_bam_data = (pangenome_info && params.map_reads) ? mapReads(params.map_reads, pangenome_info, mapping_directory) : Channel.empty()
    existing_bam_data =(params.bam_files) ? fetchBAM(params.bam_files) : Channel.empty()
    bam_data = new_bam_data.concat(existing_bam_data)


    // Get raw parquets
    new_parquet_data = (pangenome_info && bam_data) ? bamToParquet(bam_data,pangenome_info,mapping_directory) : Channel.empty()
    existing_parquet_data = (params.raw_parquet) ? fetchRawParquet(params.raw_parquet) : Channel.empty()
    raw_parquet_data = new_parquet_data.concat(existing_parquet_data)

    // Call bases
    new_called_bases = (pangenome_info && raw_parquet_data) ? callBases(raw_parquet_data,pangenome_info,mapping_directory) : Channel.empty()
    existing_called_bases = (params.called_bases) ? fetchCalledBases(params.called_bases) : Channel.empty()
    called_bases_data = new_called_bases.concat(existing_called_bases)
}   

