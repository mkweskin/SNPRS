#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// TO DO:

// SKIP READ COUNTING
// CLEAN UP READ INPUT
// AUTO PE READ REPAIR IF SUBSET FAILS

def validate_dir(path, flag) {
    if (!file(path).isDirectory()) error "${path} provided by --${flag} does not exist"
}

def validate_file(path, flag) {
    if (!file(path).exists()) error "${path} provided by --${flag} does not exist"
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
    } else if(params.joined && (params.join_id || params.join_csv)){
        error("Cannot load joined data via --joined and join other data with --join_id/--join_csv")
    } else if(params.filtered && (params.filter_id || params.filter_csv)){
        error("Cannot load filtered data via --filtered and join other data with --filter_id/--filter_csv")
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
def genome_directory
def genome_name
def mapping_directory
def joined_directory
def snp_group_directory

// SNPRS Main Script
// Params are read in from command line or from nextflow.config and/or conf/profiles.config

timestamp = "${params.timestamp}"

// Check incompatibilities

if(check_nonos()){

    snprs_directory = file(params.out)
    if(!snprs_directory.isDirectory()){
        parent_dir = snprs_directory.getParent()
        validate_dir(parent_dir,"out")
        snprs_directory.mkdirs()
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
    
    // Get strings to full paths
    snprs_directory = snprs_directory.toString()
    genome_directory = genome_directory.toString()
    mapping_directory = file("${snprs_directory}/Mapping").toString()
    joined_directory = file("${snprs_directory}/Joined").toString()
}

// Import subworkflows

include {assembleGenome} from "./subworkflows/prepare_genome/main.nf"
include {useFASTA} from "./subworkflows/prepare_genome/main.nf"
include {checkGenomeDir} from "./subworkflows/prepare_genome/main.nf"

include {fetchBAM} from "./subworkflows/mapping/main.nf"
include {fetchMapReads} from "./subworkflows/mapping/main.nf"
include {mapReads} from "./subworkflows/mapping/main.nf"

include {bamToParquet} from "./subworkflows/convert_bam/main.nf"
include {fetchRawParquet} from "./subworkflows/convert_bam/main.nf"

include {callBases} from "./subworkflows/call_bases/main.nf"
include {fetchCalledBases} from "./subworkflows/call_bases/main.nf"

include {generateScaffold} from "./subworkflows/join_parquets/main.nf"
include {fetchScaffold} from "./subworkflows/join_parquets/main.nf"
include {filterScaffold} from "./subworkflows/join_parquets/main.nf"

workflow{

    // File-ize inputs
    pg_read_data = (params.pg_reads) ? "${file(params.pg_reads)}" : ""
    fasta_file = (params.fasta) ? "${file(params.fasta)}" : ""
    bam_file = (params.bam_files) ? "${file(params.bam_files)}" : ""
    map_file = (params.map_reads) ? "${file(params.map_reads)}" : ""
    raw_parquet_file = (params.raw_parquets) ? "${file(params.raw_parquets)}" : ""
    called_base_file = (params.called_bases) ? "${file(params.called_bases)}" : ""

    ////////////////////////////////////// GENERATE/FETCH GENOME //////////////////////////////////////

    genome_info = Channel.empty()
    if(!params.no_ref){
        if(params.pg_reads){
            genome_info = assembleGenome(pg_read_data,genome_directory,genome_name) | first
        } else if(params.fasta){     
            genome_info = useFASTA(fasta_file,genome_directory,genome_name) | first
        } else{
            genome_info = checkGenomeDir(genome_directory) | first
        }
    }

    ///////////////////////////////// GENERATE/FETCH RAW PARQUETS /////////////////////////////////////

    existing_bam_data = (params.bam_files) ? fetchBAM(bam_file) : Channel.empty()
    existing_parquet_data = (params.raw_parquets) ? fetchRawParquet(raw_parquet_file) : Channel.empty()

    mapped_data = (params.map_reads) ? fetchMapReads(map_file).combine(genome_info).map{it->tuple(it[0],it[1],it[2],it[5])} | mapReads : Channel.empty()
    mapped_data = (params.local) ? mapped_data | collect | flatten | collate(2) : mapped_data
    
    new_parquet_data = mapped_data.concat(existing_bam_data).combine(genome_info).map{it->tuple(it[0],it[1],it[4])} | bamToParquet
    new_parquet_data = (params.local) ? new_parquet_data | collect | flatten | collate(2) : new_parquet_data

    raw_parquet_data = new_parquet_data.concat(existing_parquet_data) 

    /////////////////////////////////// GENERATE/FETCH CALLED BASES /////////////////////////////////////

    existing_called_base_data = (params.called_bases) ? fetchCalledBases(called_base_file) : Channel.empty()
    new_called_base_data = (params.ploidy) ? callBases(raw_parquet_data) : Channel.empty()
    called_bases_data = new_called_base_data.concat(existing_called_base_data) | collect | flatten | collate(2)

    ////////////////////////////////// GENERATE/FETCH JOINED SCAFFOLD ////////////////////////////////////

    scaffold_file = Channel.empty()

    if(params.scaffold){
        scaffold_file = file(params.scaffold)
    } else if(params.join_id){
        scaffold_file = generateScaffold(called_bases_data)
    }

    if(params.filter_id){
        filter_file = filterScaffold(scaffold_file)
    }
}












    ///////////////////////////////////// FILTER JOINED DATA //////////////////////////////////////////

    //if(params.filtered){
    //    filtered_data = fetchFiltered(params.filtered) 
    //} else if(params.filter){
    //    filtered_data = filterJoined(joined_data)
    //}

    ///////////////////////////////////// GET ALIGNMENT ///////////////////////////////////////////////

    //if(params.alignment_file){
    //    alignment_file = Channel.fromPath(params.alignment_file) 
    //} else if(params.alignment){
    //    alignment_file = getAlignment(filtered_data) | collect | flatten | collate(1)
    //}

    //////////////////////////////////////// GET TREE /////////////////////////////////////////////////

    //if(params.tree_file){
    //    tree_file = Channel.fromPath(params.tree_file) 
    //} else if(params.tree){
    //    tree_file = filtered_data.combine(alignment_file) | generateTree | collect | flatten | collate(1)
    //}
    
    //////////////////////////////////////// GET SNPS //////////////////////////////////////////////////




