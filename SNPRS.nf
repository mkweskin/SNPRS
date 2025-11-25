#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// Logging functions
def log(log_file,message) {
    log_file.withWriterAppend { writer ->
        writer.writeLine("${message}")
    }
}

def tab_log(log_file,message) {
    log_file.withWriterAppend { writer ->
        writer.writeLine("\t- ${message}")
    }
}

def date_log(log_file,message) {
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

    if(!params.out){
        error "Must specify analysis directory via --out"
    } 
    
    // Pangenome assembly checks
    else if(params.no_ref && params.pg_reads){
        "Cannot set --no_ref and --pg_reads together"
    } else if(params.pg_reads && !params.genome_name){
        error "Cannot assemble a pangenome without --genome_name"
    } else if(params.pg_reads && !params.size){
        error "Cannot assemble a pangenome without --size"
    } else if(params.no_ref && params.fasta){
        "Cannot set --no_ref and --fasta together"
    } else if(params.no_ref && params.genome_dir){
        "Cannot set --no_ref and --genome_dir together"
    } else if(params.genome_dir && params.fasta){
        error "Cannot set --genome_dir and --fasta together"
    } else if(params.genome_dir && params.pg_reads){
        error "Cannot set --genome_dir and --pg_reads together"
    } else if(params.fasta && params.pg_reads){
       error "Cannot set --fasta and --pg_reads together"
    } 
    

    // Join/Filter checks
    else if(params.join_csv && params.join_id){
        error("Can't set --join_id and --join_csv together")
    } else if(params.filter_csv && params.filter_id){
        error("Can't set --filter_id and --filter_csv together")
    }
    
    // File checks
    else if(params.pg_reads && !file(params.pg_reads).isDirectory()){
        error "Directory provided by --pg_reads does not exist"
    } else if(params.fasta && !file(params.fasta).exists()){
        error "FASTA provided by --fasta does not exist"
    } else if(params.join_csv && !file(params.join_csv).exists()){
        error "CSV file provided by --join_csv does not exist"
    } else if(params.filter_csv && !file(params.filter_csv).exists()){
        error "CSV file provided by --filter_csv does not exist"
    } else if(params.alignment_file && !file(params.alignment_file).exists()){
        error "Alignment file provided by --alignment_file does not exist"
    } else if(params.tree_file && !file(params.tree_file).exists()){
        error "Tree file provided by --tree_file does not exist"
    } else if(params.manual_counts && !file(params.manual_counts).exists()){
        error "Count file provided by --manual_counts does not exist"
    } else if(params.snp_json && !file(params.snp_json).exists()){
        error("JSON file provided by --snp_json does not exist")
    }

    return true
}

def cmd_args = workflow.commandLine
def snprs_directory 
def log_directory
def log_file 
def genome_directory
def genome_name
def mapping_directory
def joined_directory
def snp_group_directory
def sra_directory

// SNPRS Main Script
// Params are read in from command line or from nextflow.config and/or conf/profiles.config

timestamp = "${params.timestamp}"

// Check incompatibilities

if(check_nonos()){


    snprs_directory = file(params.out)
    
    if(snprs_directory.isDirectory()){
        new_dir = false
    } else{
        parent_dir = snprs_directory.getParent()
        validate_dir(parent_dir,"out")
        snprs_directory.mkdirs()
        new_dir = true 
    }

    // Log File
    log_directory = file("${snprs_directory}/Run_Logs")
    
    if(!log_directory.isDirectory()){
        log_directory.mkdirs()
    }

    log_file = file("${log_directory}/SNPRS_Log_${timestamp}.txt")        
    
    log(log_file,"SNPRS Log File")
    log(log_file,"${new java.text.SimpleDateFormat('yyyy-MM-dd HH:mm:ss').format(new java.util.Date())}\n")
    log(log_file,"Command: ${cmd_args}\n")

    if(new_dir){
        tab_log(log_file,"Created output directory: ${snprs_directory}")
    } else{
        tab_log(log_file,"Found output directory: ${snprs_directory}")
    }

    // Genome Information
    genome_directory = (params.genome_dir) ? file(params.genome_dir) : file("${snprs_directory}/Reference_Genome")

    genome_name = ""
    
    if(params.genome_dir){
        validate_dir(genome_directory,"genome_dir")
    } else if(params.genome_name){
        genome_name = "${params.genome_name}"
    } else if(params.fasta){
        fasta_file = file(params.fasta)
        validate_file(fasta_file,"fasta")
        genome_name = fasta_file.getBaseName().replaceAll(/\.f(ast[an]?)(\.gz)?$/, '')
    } 

    mapping_directory = file("${snprs_directory}/Mapping")
    joined_directory = file("${snprs_directory}/Joined")
    snp_group_directory = file("${snprs_directory}/SNP_Groups")
    sra_directory = file("${snprs_directory}/SRA_Reads")

    params.final_snprs_directory = snprs_directory
    params.final_genome_directory = genome_directory
    params.final_mapping_directory = mapping_directory
    params.final_joined_directory = joined_directory
    params.final_snp_directory = sra_directory
    params.final_sra_directory = sra_directory
    params.final_classified_directory = sra_directory

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

workflow{

    ////////////////////////////////////////// SOLO TASKS /////////////////////////////////////////////

    if(params.join_csv){
        print("Join")        
    }

    else{

        ///////////////////////////////////////// MAIN WORKFLOW ////////////////////////////////////////////


        ///////////////////////////////////////// FETCH GENOME ////////////////////////////////////////////

        genome_info = Channel.empty()

        if(!params.no_ref){
            if(params.pg_reads){
                pg_read_data = file(params.pg_reads)
                genome_info = assembleGenome(pg_read_data,genome_directory,genome_name) | first
            }  else if(params.fasta){     
                fasta_file = file(params.fasta) 
                genome_info = useFASTA(fasta_file,genome_directory,genome_name) | first
            } else{
                genome_info = checkGenomeDir(genome_directory) | first
            }
        }

        ///////////////////////////////////// FETCH RAW PARQUETS /////////////////////////////////////////

        new_parquet_data = Channel.empty()
        
        if(params.bam_files || params.map_reads || params.map_sra){
    
            // Check for new read data to be mapped
            map_read_data = (params.map_reads) ? fetchMapReads(params.map_reads) : Channel.empty()
            sra_read_data = (params.map_sra) ? fetchSRAReads(params.map_sra): Channel.empty()
            
            // Map reads
            uncollected_bam_data = map_read_data.concat(sra_read_data).combine(genome_info).map{it->tuple(it[0],it[1],it[2],it[5])} | mapReads
            new_bam_data = (params.local) ?  uncollected_bam_data | collect | flatten | collate(2) : uncollected_bam_data
            
            // Convert to parquet
            existing_bam_data = (params.bam_files) ? fetchBAM(params.bam_files) : Channel.empty()
            new_parquet_data = existing_bam_data.concat(new_bam_data).combine(genome_info).map{it->tuple(it[0],it[1],it[4])} | bamToParquet
        }
        
        existing_parquet_data = (params.raw_parquets) ? fetchRawParquet(params.raw_parquets) : Channel.empty()
        raw_parquet_data = (params.local) ? existing_parquet_data.concat(new_parquet_data) | collect | flatten | collate(2) : new_parquet_data.concat(existing_parquet_data) 

        //////////////////////////////////////// CALL BASES ///////////////////////////////////////////////

        existing_called_base_data = (params.called_bases) ? fetchCalledBases(params.called_bases) : Channel.empty()
        new_called_base_data = Channel.empty()
        
        if(params.call){
            new_called_base_data = raw_parquet_data | callBases
        }

        called_bases_data = existing_called_base_data.concat(new_called_base_data) | collect | flatten | collate(2)

        ///////////////////////////////////// JOIN CALLED BASES ///////////////////////////////////////////

        if(params.joined){
            joined_data = fetchJoin(params.joined)
        } else if(params.join){
            joined_data = joinCalledBases(called_bases_data)
        }

        ///////////////////////////////////// FILTER JOINED DATA //////////////////////////////////////////

        if(params.filtered){
            filtered_data = fetchFiltered(params.filtered) 
        } else if(params.filter){
            filtered_data = filterJoined(joined_data)
        }

        ///////////////////////////////////// GET ALIGNMENT ///////////////////////////////////////////////

        if(params.alignment_file){
            alignment_file = Channel.fromPath(params.alignment_file) 
        } else if(params.alignment){
            alignment_file = getAlignment(filtered_data) | collect | flatten | collate(1)
        }

        //////////////////////////////////////// GET TREE /////////////////////////////////////////////////

        if(params.tree_file){
            tree_file = Channel.fromPath(params.tree_file) 
        } else if(params.tree){
            tree_file = filtered_data.combine(alignment_file) | generateTree | collect | flatten | collate(1)
        }

        //////////////////////////////////////// GET SPLITS ///////////////////////////////////////////////
        
        if(params.split_file){
            split_file = Channel.fromPath(params.split_file) 
        } else if(params.split){
            split_file = tree_file | makeSplitTable | collect | flatten | collate(1)
        }

        //////////////////////////////////////// GET GROUPS ///////////////////////////////////////////////

        if(params.group_file){
            group_file = Channel.fromPath(params.group_file) 
        } else if(params.split){
            group_file = tree_file.combine(split_file) | makeSNPGroups | collect | flatten | collate(1)
        }

        //////////////////////////////////////// GET SNPS //////////////////////////////////////////////////
        if(params.make_snps){
            snp_data = filtered_data.combine(tree_file).combine(group_file).map{it->tuple(it[0],it[1],snp_id,snp_directory,it[2],it[3])} | generateSNPs | collect | flatten | collate(2)
        }
    }

}   