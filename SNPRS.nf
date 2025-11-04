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

def validate_dir(path, flag) {
    if (!path.isDirectory()) error "${path} provided by --${flag} does not exist"
}

def validate_file(path, flag) {
    if (!path.exists()) error "${path} provided by --${flag} does not exist"
}

def check_nonos(){

    if(params.genome_dir && params.fasta){
        error "Cannot set --genome_dir and --fasta together"
    } else if(params.genome_dir && params.pg_reads){
        error "Cannot set --genome_dir and --pg_reads together"
    } else if(params.fasta && params.pg_reads){
        error "Cannot set --fasta and --pg_reads together"
    }  else if(params.pg_reads && !params.size){
        error "Cannot set assemble pangenome without --size"
    } else if(params.pangenome && !params.pg_reads){
        error "Cannot set assemble pangenome without --pg_reads"
    } else if((params.joined || params.filtered) && (params.map_reads || params.map_sra)){
        error "Cannot add new mapping data to existing joined/filtered datasets"
    } else if(params.alignment_file && !file(params.alignment_file).exists()){
        error "Alignment file provided by --alignment_file does not exist"
    } else if(params.tree_file && !file(params.tree_file).exists()){
        error "Tree file provided by --tree_file does not exist"
    } else if(params.split_file && !file(params.split_file).exists()){
        error "Split file provided by --split_file does not exist"
    } else if(params.group_file && !file(params.group_file).exists()){
        error "Group file provided by --group_file does not exist"
    } else if(params.classify && !params.snp_dir){
        error "Cannot run in classifier mode (--classify) without providing SNP information (--snp_dir)"
    } else if(params.classify && (!params.map_reads && !params.map_sra && !params.bam_files && !params.raw_parquets && !params.called_bases)){
        error "Cannot run in classifier mode (--classify) without data (--map_reads/--map_sra/--bam_files/--raw_parquets/--called_bases)"
    } else if(params.classify && params.pg_reads){
        error "Cannot run in classifier mode (--classify) and assemble a pangenome (--pg_reads) in a single run"
    }

    return true
}

def cmd_args = workflow.commandLine

// SNPRS Main Script
// Params are read in from command line or from nextflow.config and/or conf/profiles.config

timestamp = "${params.timestamp}"

// Check incompatibilities

if(check_nonos()){

    if(params.out){
        snprs_directory = file(params.out)
        parent_dir = snprs_directory.getParent()

        if(snprs_directory.isDirectory()){
            new_dir = false
        } else{
            validate_dir(parent_dir,"out")
            snprs_directory.mkdirs()
            new_dir = true 
        }

    } else if(params.joined || params.filtered || params.snp_dir ||  params.group_file || params.split_file || params.genome_dir){

        if(params.joined){

            test_dir = file(params.joined)
            validate_dir(test_dir,"joined")
            snprs_directory = test_dir.getParent().getParent()
            new_dir = false


        } else if(params.filtered){

            test_dir = file(params.filtered)
            validate_dir(test_dir,"filtered")
            snprs_directory = test_dir.getParent().getParent()
            new_dir = false

        } else if(params.snp_dir){

            test_dir = file(params.snp_dir)
            validate_dir(test_dir,"snp_dir")
            snprs_directory = test_dir.getParent().getParent()
            new_dir = false

        } else if(params.split_file){

            test_file = file(params.split_file)
            validate_file(test_file,"split_file")
            snprs_directory = test_file.getParent().getParent().getParent()
            new_dir = false

        } else if(params.group_file){
            
            test_file = file(params.group_file)
            validate_file(test_file,"group_file")
            snprs_directory = test_file.getParent().getParent().getParent()
            new_dir = false

        } else if(params.genome_dir){
            
            test_dir = file(params.genome_dir)
            validate_dir(test_dir,"genome_dir")
            snprs_directory = test_dir.getParent()
            new_dir = false

        } else{
            error "Output directory could not be determined"
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
    mapping_directory = file("${snprs_directory}/Mapping")
    joined_directory = file("${snprs_directory}/Joined")
    sra_directory = file("${snprs_directory}/SRA_Reads")
    classified_directory = file("${snprs_directory}/Classification")
    genome_directory = (params.genome_dir) ? file(params.genome_dir) : file("${snprs_directory}/Reference_Genome")

    // Genome Name
    new_genome_name = "SNPRS_${timestamp}"

    if(!params.genome_dir){
        if(params.genome_name){
            new_genome_name = "${params.genome_name}"
        } else if(params.fasta){
            fasta_file = file(params.fasta)
            validate_file(fasta_file,"fasta")
            new_genome_name = fasta_file.getBaseName().replaceAll(/\.f(ast[an]?)(\.gz)?$/, '')
        } 
    } else{
        validate_dir(genome_directory,"genome_dir")
    }

    // Subanalysis IDs
    join_id = (params.join_id) ? "${params.join_id}" : "Joined_${timestamp}"
    filter_id = (params.filter_id) ? "${params.filter_id}" : "Filtered_${timestamp}"

    // SNP Directory
    if (params.snp_dir) {
        snp_dir_path = file(params.snp_dir)
        validate_dir(snp_dir_path,"snp_dir")
        snp_directory = snp_dir_path.getParent()
        snp_id = snp_dir_path.getName()

        check_output_json = file("${snp_directory}/${snp_id}/${snp_id}.json")
        check_output_comparisons = file("${snp_directory}/${snp_id}/${snp_id}_Comparisons.csv")
        check_snp_parquet = file("${snp_directory}/${snp_id}/${snp_id}_SNPs.parquet")

        validate_file(check_output_json,"snp_dir")
        validate_file(check_output_comparisons,"snp_dir")
        validate_file(check_row_numbers,"snp_dir")
        validate_file(check_snp_parquet,"snp_dir")

    } else{
        snp_directory = file("${snprs_directory}/SNP_Analysis")
        snp_id = (params.snp_id) ? "${params.snp_id}" : "SNP_${timestamp}"
    }

    // Parameterize major directories
    params.final_snprs_directory = snprs_directory
    params.final_genome_directory = genome_directory
    params.final_mapping_directory = mapping_directory
    params.final_joined_directory = joined_directory
    params.final_snp_directory = snp_directory
    params.final_sra_directory = sra_directory
    params.final_classified_directory = classified_directory

    // Parameterize IDs
    params.final_join_id = join_id
    params.final_filter_id = filter_id
    params.final_snp_id = snp_id
    params.new_genome_name = new_genome_name
}

// Import subworkflows

include {assembleGenome} from "./subworkflows/prepare_genome/main.nf"
include {useFASTA} from "./subworkflows/prepare_genome/main.nf"
include {checkGenomeDir} from "./subworkflows/prepare_genome/main.nf"

include {fetchBAM} from "./subworkflows/mapping/main.nf"
include {fetchMapReads} from "./subworkflows/mapping/main.nf"
include {fetchSRAReads} from "./subworkflows/mapping/main.nf"
include {mapReads} from "./subworkflows/mapping/main.nf"

include {bamToParquet} from "./subworkflows/convert_bam/main.nf"
include {fetchRawParquet} from "./subworkflows/convert_bam/main.nf"

include {callBases} from "./subworkflows/call_bases/main.nf"
include {fetchCalledBases} from "./subworkflows/call_bases/main.nf"

include {joinCalledBases} from "./subworkflows/join_parquets/main.nf"
include {fetchJoin} from "./subworkflows/join_parquets/main.nf"

include {filterJoined} from "./subworkflows/filter_joined/main.nf"
include {fetchFiltered} from "./subworkflows/filter_joined/main.nf"

include {getAlignment} from "./subworkflows/alignment_tools/main.nf"

include {fetchTree} from "./subworkflows/tree_tools/main.nf"
include {makeSplitTable} from "./subworkflows/tree_tools/main.nf"
include {makeSNPGroups} from "./subworkflows/tree_tools/main.nf"
include {generateTree} from "./subworkflows/tree_tools/main.nf"

include {generateSNPs} from "./subworkflows/snp_tools/main.nf"
include {classifySample} from "./subworkflows/classifier/main.nf"

workflow{

    ///////////////////////////////////////// FETCH GENOME ////////////////////////////////////////////

    genome_info = Channel.empty()

    if(params.pg_reads){
        genome_info = assembleGenome(params.pg_reads) | first
    }  else if(params.fasta){      
        genome_info = useFASTA(params.fasta) | first
    } else{
        genome_info = checkGenomeDir(genome_directory) | first
    }

    ///////////////////////////////////// FETCH RAW PARQUETS /////////////////////////////////////////

    new_parquet_data = Channel.empty()
    raw_parquet_data = Channel.empty()
    called_bases_data = Channel.empty()

    // Get existing files    
    existing_bam_data = (params.bam_files) ? fetchBAM(params.bam_files) : Channel.empty()
    existing_parquet_data = (params.raw_parquets) ? fetchRawParquet(params.raw_parquets) : Channel.empty()
    existing_called_base_data = (params.called_bases) ? fetchCalledBases(params.called_bases) : Channel.empty()

    if((!params.joined && !params.filtered) && (params.bam_files || params.map_reads || params.map_sra)){
        
        // Check for new read data to be mapped
        map_read_data = (params.map_reads) ? fetchMapReads(params.map_reads) : Channel.empty()
        sra_read_data = (params.map_sra) ? fetchSRAReads(params.map_sra): Channel.empty()
        
        // Map reads
        uncollected_bam_data = map_read_data.concat(sra_read_data).combine(genome_info).map{it->tuple(it[0],it[1],it[2],it[5])} | mapReads
        new_bam_data = (params.local) ?  uncollected_bam_data | collect | flatten | collate(2) : uncollected_bam_data
        
        // Convert to parquet
        new_parquet_data = existing_bam_data.concat(new_bam_data).combine(genome_info).map{it->tuple(it[0],it[1],it[4])} | bamToParquet
    }
    
    raw_parquet_data = (params.local) ? existing_parquet_data.concat(new_parquet_data) | collect | flatten | collate(2) : new_parquet_data.concat(existing_parquet_data) 

    //////////////////////////////////////// CALL BASES ///////////////////////////////////////////////

    new_called_base_data = raw_parquet_data | callBases
    called_bases_data = existing_called_base_data.concat(new_called_base_data) | collect | flatten | collate(2)

    //////////////////////////////////////// CLASSIFY /////////////////////////////////////////////////

    if(params.classify && params.snp_dir){
      
        classified_data = genome_info.combine(called_bases_data).map{it->tuple(it[0],it[1],it[3],it[4],snp_directory,snp_id)} | classifySample | collect | flatten | collate(3)
        classified_data | view
        
    } else{

        tree_file = (params.tree_file) ? Channel.fromPath(params.tree_file) : Channel.empty()
        alignment_file = (params.alignment_file) ? Channel.fromPath(params.alignment_file) : Channel.empty()
        split_file = (params.split_file) ? Channel.fromPath(params.split_file) : Channel.empty()
        group_file = (params.group_file) ? Channel.fromPath(params.group_file) : Channel.empty()

    ///////////////////////////////////// JOIN CALLED BASES ///////////////////////////////////////////
    
        joined_data = (params.joined) ? fetchJoin(params.joined) : joinCalledBases(called_bases_data)

    ///////////////////////////////////// FILTER JOINED DATA //////////////////////////////////////////

        filtered_data = (params.filtered) ? fetchFiltered(params.filtered) : filterJoined(joined_data)

    ///////////////////////////////////// GET ALIGNMENT ///////////////////////////////////////////////

        if(!params.filter && !params.snp_dir && !params.split_file && !params.group_file && !params.alignment_file){
            alignment_file = getAlignment(filtered_data) | collect | flatten | collate(1)
        }

    //////////////////////////////////////// GET TREE /////////////////////////////////////////////////

        if(!params.alignment && !params.snp_dir && !params.split_file && !params.group_file && !params.tree_file){
            tree_file = filtered_data.combine(alignment_file) | generateTree | collect | flatten | collate(1)
        }

    //////////////////////////////////////// GET SPLITS ///////////////////////////////////////////////

        if(!params.tree && !params.snp_dir && !params.split_file && !params.group_file){
            split_file = tree_file | makeSplitTable | collect | flatten | collate(1)
        }

    //////////////////////////////////////// GET GROUPS ///////////////////////////////////////////////

        if(!params.tree && !params.split && !params.snp_dir && !params.group_file && params.split_file){
            group_file = tree_file.combine(split_file) | makeSNPGroups | collect | flatten | collate(1)
        }

    //////////////////////////////////////// GET SNPS //////////////////////////////////////////////////

        if(params.group_file && !params.snp_dir){        
            snp_data = filtered_data.combine(tree_file).combine(group_file).map{it->tuple(it[0],it[1],snp_id,snp_directory,it[2],it[3])} | generateSNPs | collect | flatten | collate(2)
        }
    }
}   