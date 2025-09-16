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

workflow{

    // Assemble pangenome from reads
    if("${params.pg_reads}" != ""){
        pangenome_info = makePangenome(pangenome_directory,pg_name,pg_read_directory)
    } 
    // Get FASTA from --fasta (creates fai/ref in needed)
    else if("${params.fasta}" != ""){
        pangenome_info = indexGenome(params.fasta)
    }

    // Specify pangenome by name (checks for fasta, fai, and ref in SNPRS_Pangenomes/PG_NAME)
    else if("${pg_name}" != "SNPRS_${params.timestamp}"){
        check_dir = file("${pangenome_directory}/${pg_name}")
        pangenome_info = checkGenome(check_dir)
    }
    
    // No pangenome information provided
    else {
        pangenome_info = Channel.empty()
    }

    pangenome_info.view()
}