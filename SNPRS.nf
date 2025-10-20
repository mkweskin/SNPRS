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

// Genome Directory
genome_directory = file("${snprs_directory}/Reference_Genome")
if(!genome_directory.isDirectory()){
    genome_directory.mkdirs() 
    tab_log("Created reference genome directory: ${genome_directory}...")
} else{
    tab_log("Found reference genome directory: ${genome_directory}...")
}

// Genome Name
if(params.genome_name){
    genome_name = "${params.genome_name}"
} else if(params.fasta){
    fasta_file = file(params.fasta)
    genome_name = fasta_file.getBaseName().replaceAll(/\.f(ast[an]?)(\.gz)?$/, '')
} else{
    genome_name = "SNPRS_${timestamp}"
}

// Prep Directories
genome_prep_directory = file("${genome_directory}/Prep_${genome_name}")
genome_read_link_directory = file("${genome_directory}/Pangenome_Read_Links")

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

// SNP Directory
snp_directory = file("${snprs_directory}/SNP_Analysis")
if(!snp_directory.isDirectory()){
    snp_directory.mkdirs() 
    tab_log("Created SNP directory: ${snp_directory}...")
} else{
    tab_log("Found SNP directory: ${snp_directory}...")
}

// Fixed Directory
fixed_directory = file("${snp_directory}/Fixed_Sites")
if(!fixed_directory.isDirectory()){
    fixed_directory.mkdirs() 
    tab_log("Created fixed site directory: ${fixed_directory}...")
} else{
    tab_log("Found fixed site directory: ${fixed_directory}...")
}

// Join ID
if(params.validate){
    join_id = "Validation"
} else if(!params.join_id){
    join_id = "Joined_SNPRS_${timestamp}"
} else{
    join_id = "${params.join_id}"
}

// Filter ID
if(params.validate){
    filter_id = "Filter_Validation"
} else if(!params.filter_id){
    filter_id = "Filter_SNPRS_${timestamp}"
} else{
    filter_id = "${params.filter_id}"
}

// Fixed ID
if(!params.fixed_id){
    fixed_id = "Fixed_SNPRS_${timestamp}"
} else{
    fixed_id = "${params.fixed_id}"
}

// Refine ID
if(!params.refine_id){
    refine_id = "Refined_SNPRS_${timestamp}"
} else{
    refine_id = "${params.refine_id}"
}

// Import workflows
include {assembleGenome} from "./subworkflows/prepare_genome/main.nf"
include {useFASTA} from "./subworkflows/prepare_genome/main.nf"
include {checkGenomeDir} from "./subworkflows/prepare_genome/main.nf"

include {mapReads} from "./subworkflows/mapping/main.nf"
include {fetchBAM} from "./subworkflows/mapping/main.nf"

include {bamToParquet} from "./subworkflows/convert_bam/main.nf"
include {fetchRawParquet} from "./subworkflows/convert_bam/main.nf"

include {callBases} from "./subworkflows/call_bases/main.nf"
include {fetchCalledBases} from "./subworkflows/call_bases/main.nf"

include {joinCalledBases} from "./subworkflows/join_parquets/main.nf"
include {fetchJoin} from "./subworkflows/join_parquets/main.nf"

include {filterJoined} from "./subworkflows/filter_joined/main.nf"
include {fetchFiltered} from "./subworkflows/filter_joined/main.nf"

include {getAlignment} from "./subworkflows/alignment_tools/main.nf"
include {generateTree} from "./subworkflows/generate_tree/main.nf"

workflow{

    // Pangenome Info (Genome Name, Directory, FASTA file)

    genome_info = Channel.empty()
    
    // Assemble pangenome from reads
    if(params.pg_reads){
        pg_read_data = file(params.pg_reads)
        genome_info = assembleGenome(genome_directory,genome_name,pg_read_data) | first
    } 

    // Get reference information from an assembly
    else if(params.fasta){
        fasta_file = file(params.fasta)
        genome_info = useFASTA(genome_directory,genome_name,fasta_file) | first
    }

    // Get reference information from a folder
    else{
        genome_dir = (params.genome_dir) ? file(params.genome_dir) : genome_directory
        genome_info = checkGenomeDir(genome_dir) | first
    }

    // Load in joined or filtered datasets if provided
    joined_data = Channel.empty()
    filtered_data = Channel.empty()
    
    if(params.filtered){
        filtered_path = file("${params.filtered}")
        filtered_data = (genome_info) ? fetchFiltered(filtered_path) : Channel.empty()
    } else if(params.joined){
        join_path = file("${params.joined}")
        joined_data = (genome_info) ? fetchJoin(join_path) : Channel.empty()

        // Filter data if requested
        filtered_data = (genome_info && ((params.filter || params.validate))) ? filterJoined(genome_info,joined_data,filter_id) : Channel.empty()
    } 

    // Map, convert, call bases, join, and filter as requested
    
    else{

        // Get BAM files    
        existing_bam_data = (params.bam_files && !params.validate) ? fetchBAM(params.bam_files) : Channel.empty()

        if(params.validate || params.map_reads){
            map_read_dir = (params.validate) ? file("${genome_read_link_directory}") : file("${params.map_reads}")
            new_bam_data = (genome_info) ? mapReads(map_read_dir,genome_info,mapping_directory) : Channel.empty()
        } else{
            new_bam_data = Channel.empty()
        }

        if("${params.runProfile}" == "local"){
            bam_data = existing_bam_data.concat(new_bam_data) | collect | flatten | collate(2)
        } else{
            bam_data = existing_bam_data.concat(new_bam_data) 
        }

        // Get raw parquets
        existing_parquet_data = (params.raw_parquets && !params.validate) ? fetchRawParquet(params.raw_parquets) : Channel.empty()
        new_parquet_data = (genome_info && bam_data) ? bamToParquet(bam_data,genome_info,mapping_directory) : Channel.empty()    

        if("${params.runProfile}" == "local"){
            raw_parquet_data = existing_parquet_data.concat(new_parquet_data) | collect | flatten | collate(2)
        } else{
            raw_parquet_data = existing_parquet_data.concat(new_parquet_data)
        }

        // Get called bases
        existing_called_base_data = (params.called_bases && !params.validate) ? fetchCalledBases(params.called_bases) : Channel.empty()
        new_called_base_data = (genome_info && raw_parquet_data) ? callBases(raw_parquet_data,genome_info,mapping_directory) : Channel.empty()
        called_bases_data = existing_called_base_data.concat(new_called_base_data) | collect  | flatten | collate(2)

        // Join bases
        joined_data = (called_bases_data) ? joinCalledBases(called_bases_data,joined_directory,join_id) : Channel.empty()

        // Filter data if requested
        filtered_data = (genome_info && ((params.filter || params.validate))) ? filterJoined(genome_info,joined_data,filter_id) : Channel.empty()
    }

    // Generate alignment from filtered or joined data if requested
    alignment_file = Channel.empty()
    if (params.alignment || params.validate) {
        if (filtered_data) {
            alignment_file = getAlignment(filtered_data)
        } else if (joined_data) {
            alignment_file = getAlignment(joined_data)
        }
    }

    // Generate tree from filtered data if requested
    tree_file = Channel.empty()
    if (params.tree || params.validate) {
        if (filtered_data) {
            tree_file = filtered_data.combine(alignment_file) | collect | flatten | collate(3) | generateTree
        } else if (joined_data) {
            tree_file = joined_data.combine(alignment_file) | collect | flatten | collate(3) | generateTree
        }
    }
}