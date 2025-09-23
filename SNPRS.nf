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

if(params.out){
    snprs_directory = file(params.out)
    parent_dir = snprs_directory.getParent()
    if(!snprs_directory.isDirectory()){
        if(!parent_dir.isDirectory()){
            error "Parent directory for output is not a valid directory [${parent_dir}]..."
        } else{
            snprs_directory.mkdirs()
            new_dir = true 
        }
    } else{
        new_dir = false
    }
} else{
    snprs_directory = file("SNPRS_${timestamp}")
    snprs_directory.mkdirs()
    new_dir = true 
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
    
    // Cache log file info as params for other processes
    params.log_directory = file(log_directory)
    params.log_file = file(log_file)
    
    log("SNPRS Log File")
    log("${new java.text.SimpleDateFormat('yyyy-MM-dd HH:mm:ss').format(new java.util.Date())}\n")
    log("Command: ${cmd_args}\n")

    if(new_dir){
        tab_log("Created output directory: ${snprs_directory}")
    } else{
        tab_log("Found output directory: ${snprs_directory}")
    }
}

// Major subdirectories

// Pangenome Directory
pangenome_directory = file("${snprs_directory}/SNPRS_Pangenomes")
if(!pangenome_directory.isDirectory()){
    pangenome_directory.mkdirs() 
    tab_log("Created pangenome directory: ${pangenome_directory}...")
} else{
    tab_log("Found pangenome directory: ${pangenome_directory}...")
}

// Mapping Directory
mapping_directory = file("${snprs_directory}/Mapping")
if(!mapping_directory.isDirectory()){
    mapping_directory.mkdirs() 
    tab_log("Created mapping directory: ${mapping_directory}...")
} else{
    tab_log("Found mapping directory: ${mapping_directory}...")
}

// Joined Directory
joined_directory = file("${snprs_directory}/Joined")
if(!joined_directory.isDirectory()){
    joined_directory.mkdirs() 
    tab_log("Created joining directory: ${joined_directory}...")
} else{
    tab_log("Found joining directory: ${joined_directory}...")
}

// Join ID
if(params.validate){
    join_id = "Validation"
} else if(!params.join_id){
    join_id = "SNPRS_${timestamp}"
} else{
    join_id = "${params.join_id}"
}

// Filter ID
if(params.validate){
    filter_id = "Filter_Validation"
} else if(!params.filter_id){
    join_id = "Filter_SNPRS_${timestamp}"
} else{
    filter_id = "${params.filter_id}"
}

// Check for validation mode
if(params.validate){    
    tab_log("Running in validation mode, pangenome reads will be mapped back onto the pangenome and joined into an alignment")
}

// Pangenome output infromation
pg_name = (params.pg_name) ? "${params.pg_name}" : "SNPRS_${timestamp}"

// Set relevant paths based on --pg_name
current_pg_directory = file("${pangenome_directory}/${pg_name}")
validation_directory = file("${current_pg_directory}/Validation")
validation_read_directory = file("${validation_directory}/Reads")

current_joined_directory = (params.validate) ? file("${validation_directory}/Joined") : file("${joined_directory}/${pg_name}")
current_mapping_directory = (params.validate) ? file("${validation_directory}/Mapping") : file("${mapping_directory}/${pg_name}")

include {assembleGenome} from "./subworkflows/prepare_genome/main.nf"
include {prepareGenome} from "./subworkflows/prepare_genome/main.nf"
include {checkSNPRSGenome} from "./subworkflows/prepare_genome/main.nf"

include {mapReads} from "./subworkflows/mapping/main.nf"
include {fetchBAM} from "./subworkflows/mapping/main.nf"

include {bamToParquet} from "./subworkflows/convert_bam/main.nf"
include {fetchRawParquet} from "./subworkflows/convert_bam/main.nf"

include {callBases} from "./subworkflows/call_bases/main.nf"
include {fetchCalledBases} from "./subworkflows/call_bases/main.nf"

include {joinCalledBases} from "./subworkflows/join_parquets/main.nf"
include {fetchJoin} from "./subworkflows/join_parquets/main.nf"

workflow{

    // Assemble pangenome from reads
    if(params.pg_reads){
        pg_read_data = file(params.pg_reads)
        pangenome_info = assembleGenome(pangenome_directory,pg_name,pg_read_data) | first
    } 

    // Get FASTA from --fasta (creates fai/ref in needed)
    else if(params.fasta){
        fasta_file = file(params.fasta)
        pangenome_info = prepareGenome(fasta_file) | first
    }

    // Specify pangenome by name (checks for fasta, fai, and ref in SNPRS_Pangenomes/PG_NAME)
    else if(params.pg_name){
        pangenome_info = checkSNPRSGenome(pangenome_directory,params.pg_name) | first
    }

    // No pangenome information provided
    else {
        pangenome_info = Channel.empty()
    }

    // Map reads, call bases, and join files
    if(!params.joined){

        // Get BAM files    
        existing_bam_data = (params.bam_files && !params.validate) ? fetchBAM(params.bam_files) : Channel.empty()

        if(params.validate){
            new_bam_data = (pangenome_info) ? mapReads(validation_read_directory,pangenome_info,current_mapping_directory) : Channel.empty()
        } else if(params.map_reads){
            mapping_read_data = (params.map_reads) ? file(params.map_reads) : ""
            new_bam_data = (pangenome_info && params.map_reads) ? mapReads(mapping_read_data,pangenome_info,current_mapping_directory) : Channel.empty()
        } else{
            new_bam_data = Channel.empty()
        }

        bam_data = existing_bam_data.concat(new_bam_data)

        // Get raw parquets
        existing_parquet_data = (params.raw_parquet && !params.validate) ? fetchRawParquet(params.raw_parquet) : Channel.empty()
        new_parquet_data = (pangenome_info && bam_data) ? bamToParquet(bam_data,pangenome_info,current_mapping_directory) : Channel.empty()
        raw_parquet_data = existing_parquet_data.concat(new_parquet_data)

        // Get called bases
        existing_called_base_data = (params.called_bases && !params.validate) ? fetchCalledBases(params.called_bases) : Channel.empty()
        new_called_base_data = (pangenome_info && raw_parquet_data) ? callBases(raw_parquet_data,pangenome_info,current_mapping_directory) : Channel.empty()
        called_bases_data = existing_called_base_data.concat(new_called_base_data) | collect  | flatten | collate(2)

        // Join bases
        joined_data = (pangenome_info && called_bases_data) ? joinCalledBases(called_bases_data,pangenome_info,current_joined_directory,join_id) : Channel.empty()
    } 
    
    else{
        // Fetch already joined bases
        join_path = file("${params.joined}")
        joined_data = (pangenome_info && params.joined) ? fetchJoin(join_path) : Channel.empty()
    }
}