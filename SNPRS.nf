#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// SNPRS Main Script
// Params are read in from command line or from nextflow.config and/or conf/profiles.config

// Set directory structure
output_directory = file(params.out)
if(!output_directory.isDirectory()){
    if(!output_directory.getParent().isDirectory()){
        error "Parent directory for output is not a valid directory [${output_directory.getParent()}]..."
    } else{
        output_directory.mkdirs() 
    }
}

// Temporary while fixing Ray
if(params.pg_reads != "" && params.ref_path == ""){
    error "Reference path is required for pangenome generation..."
}

// Parameterize major paths
params.snprs_directory = file(output_directory)

include {makePangenome;fetchPangenome} from "./subworkflows/pangenome/main.nf"

workflow{

    // Get pangenome information
    pangenome_info = "${params.pg_reads}" != "" ? makePangenome() : fetchPangenome()
    pangenome_info.subscribe{println(it)}
}