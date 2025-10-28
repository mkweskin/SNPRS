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
    }
    
    else{
        return true
    }
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

        } else if(params.filtered){

            test_dir = file(params.filtered)
            validate_dir(test_dir,"filtered")
            snprs_directory = test_dir.getParent().getParent()

        } else if(params.snp_dir){

            test_dir = file(params.snp_dir)
            validate_dir(test_dir,"snp_dir")
            snprs_directory = test_dir.getParent().getParent()

        } else if(params.split_file){

            test_file = file(params.split_file)
            validate_file(test_file,"split_file")
            snprs_directory = test_file.getParent().getParent().getParent()


        } else if(params.group_file){
            
            test_file = file(params.group_file)
            validate_file(test_file,"group_file")
            snprs_directory = test_file.getParent().getParent().getParent()

        } else if(params.genome_dir){
            
            test_dir = file(params.genome_dir)
            validate_dir(test_dir,"genome_dir")
            snprs_directory = test_dir.getParent()

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
    snp_directory = file("${snprs_directory}/SNP_Analysis")
    snp_id = (params.snp_id) ? "${params.snp_id}" : "SNP_${timestamp}"

    if (params.split_file){
        split_path = file(params.split_file)
        validate_file(split_path,"split_file")
        snp_directory = split_path.getParent().getParent()
        snp_id = split_path.getParent().getName()
    } else if (params.group_file){
        group_path = file(params.group_file)
        validate_file(group_path,"group_file")
        snp_directory = group_path.getParent().getParent()
        snp_id = group_path.getParent().getName()
    } else if (params.snp_dir) {
        snp_dir_path = file(params.snp_dir)
        validate_dir(snp_dir_path,"snp_dir")
        snp_directory = snp_dir_path.getParent()
        snp_id = snp_dir_path.getName()
    }

}

// Parameterize major directories
params.final_snprs_directory = snprs_directory
params.final_genome_directory = genome_directory
params.final_mapping_directory = mapping_directory
params.final_joined_directory = joined_directory
params.final_snp_directory = snp_directory
params.final_sra_directory = sra_directory

// Parameterize IDs
params.final_join_id = join_id
params.final_filter_id = filter_id
params.final_snp_id = snp_id
params.new_genome_name = new_genome_name

// Import workflows
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

workflow{

    ///////////////////////////////////// FETCH GENOME ////////////////////////////////////////////////

    genome_info = Channel.empty()

    if(params.pg_reads){
        pg_reads = file(params.pg_reads)
        genome_info = assembleGenome(pg_reads) | first
    }  else if(params.fasta){      
        fasta_file = file(params.fasta)
        genome_info = useFASTA(fasta_file) | first
    } else{
        genome_info = checkGenomeDir(genome_directory) | first
    }





}


/*
    //////////////////////////////////////////////////////////////////////////////////////////////////


    ///////////////////////////////////// FETCH RAW PARQUETS /////////////////////////////////////////

    bam_data = Channel.empty()
    raw_parquet_data = Channel.empty()

    // Get existing BAM/Parquet files    
    existing_bam_data = (params.bam_files) ? fetchBAM(params.bam_files) : Channel.empty()
    existing_parquet_data = (params.raw_parquets) ? fetchRawParquet(params.raw_parquets) : Channel.empty()

    // Check for read data to be mapped
    map_read_data = (params.map_reads) ? fetchMapReads(params.map_reads) : Channel.empty()
    sra_read_data = (params.map_sra) ? fetchSRAReads(params.map_sra,sra_directory): Channel.empty()
    read_data = map_read_data.concat(sra_read_data) | collect | flatten | collate(3)

    // Map reads
    new_bam_data = read_data.combine(genome_info).map{it->(it[0],it[1],it[2],it[5],"${mapping_directory}")} | mapReads : Channel.empty()
    bam_data = ("${params.runProfile}" == "local") ? existing_bam_data.concat(new_bam_data) : existing_bam_data.concat(new_bam_data) | collect | flatten | collate(2)
    
    // Convert to parquet

    new_parquet_data = () bamToParquet(bam_data,genome_info,mapping_directory) : Channel.empty()    


    
    //////////////////////////////////////////////////////////////////////////////////////////////////



        if(params.map_reads || params.map_sra){


            if(genome_info){
                


                all_reads = map_read_data.concat(sra_read_data) | collect | flatten | collate(2)
                
                // Map reads


                // Convert to parquet
                raw_parquet_data = ("${params.runProfile}" == "local") ? existing_parquet_data.concat(new_parquet_data) :  existing_parquet_data.concat(new_parquet_data) | collect | flatten | collate(2)

            }
        }

    }
}



/*


        }




        if(params.validate || params.map_reads){
            map_read_dir = (params.validate) ? file("${genome_read_link_directory}") : file("${params.map_reads}")
            new_bam_data = (genome_info) ? mapReads(map_read_dir,genome_info,mapping_directory) : Channel.empty()
        } else{
            new_bam_data = Channel.empty()
        }



        

        // Check for SRA mapping datas
        sra_bam_data = (params.map_sra) ? mapSRA(params.map_sra,sra_directory,genome_info,mapping_directory) : Channel.empty()        








    }
}


        if(params.no_ref){

            if(params.filtered){
                filtered_path = file("${params.filtered}")
                filtered_data = fetchFiltered(filtered_path)
            }
            
            if(params.joined){
                join_path = file("${params.joined}")
                joined_data = fetchJoin(join_path)
            }
        }

















    }

}





    // Pangenome Info (Genome Name, Directory, FASTA file)




    else{
        

        // Get reference information from an assembly
        else
        // Get reference information from a folder
        else{


        // Load in joined or filtered datasets if provided

        
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



            // Get called bases
            existing_called_base_data = (params.called_bases && !params.validate) ? fetchCalledBases(params.called_bases) : Channel.empty()
            new_called_base_data = (genome_info && raw_parquet_data) ? callBases(raw_parquet_data,genome_info,mapping_directory) : Channel.empty()
            called_bases_data = existing_called_base_data.concat(new_called_base_data) | collect | flatten | collate(2)

            // Join bases
            joined_data = (called_bases_data) ? joinCalledBases(called_bases_data,joined_directory,join_id) : Channel.empty()

            // Filter data if requested
            filtered_data = (genome_info && ((params.filter || params.validate))) ? filterJoined(genome_info,joined_data,filter_id) : Channel.empty()
        }
    }

    // Generate alignment from filtered or joined data if requested
    alignment_file = Channel.empty()
    if(params.alignment_file){
        Channel.fromPath(params.alignment_file).set { alignment_file }
    } else if (params.alignment || params.validate) {
        if (filtered_data) {
            alignment_file = getAlignment(filtered_data) | collect | flatten | collate(1)
        } else if (joined_data) {
            alignment_file = getAlignment(joined_data) | collect | flatten | collate(1)
        }
    }

    // Generate tree from filtered data if requested
    tree_file = Channel.empty()
    if(params.tree_file){
        tree_path = file(params.tree_file) 
        tree_file = fetchTree(tree_path) | collect | flatten | collate(1)
    } else{
        if (params.tree || params.validate) {
            if (filtered_data) {
                tree_file = filtered_data.combine(alignment_file) | collect | flatten | collate(3) | generateTree | collect | flatten | collate(1)
            } else if (joined_data) {
                tree_file = joined_data.combine(alignment_file) | collect | flatten | collate(3) | generateTree | collect | flatten | collate(1)
            }
        }
    }

    // If tree data is available, generate split table if needed
    snp_dir_ch = Channel.of([snp_id, snp_directory]) 
    
    // Intializes SNP directory
    split_file = Channel.empty()
    
    if(!params.group_file){
        split_file = (params.split_file) ? Channel.fromPath(params.split_file) : tree_file.combine(snp_dir_ch) | makeSplitTable | collect | flatten | collate(1)
    }
    
    // Only manual for now to allow for labeling
    split_data_ch = (params.split_file) ? tree_file.combine(split_file).combine(snp_dir_ch) : Channel.empty()
    group_file = (params.group_file) ? Channel.fromPath(params.group_file) : split_data_ch | makeSNPGroups | collect | flatten | collate(1)

    snp_data_ch = (params.group_file) ? snp_dir_ch.combine(tree_file).combine(group_file) : Channel.empty()
    snp_pre_data = (params.filtered) ? filtered_data.combine(snp_data_ch) : joined_data.combine(snp_data_ch) 

    snp_data = snp_pre_data | generateSNPs | collect | flatten | collate(2)
}
*/