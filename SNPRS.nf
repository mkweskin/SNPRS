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

pangenome_directory = file("${output_directory}/SNPRS_Pangenomes")
mapping_directory = file("${output_directory}/Mapping")

// Parameterize major paths
params.snprs_directory = file(output_directory)
params.pangenome_directory = file(pangenome_directory)
params.mapping_directory = file(mapping_directory)

include {makePangenome;fetchPangenome} from "./subworkflows/pangenome/main.nf"

workflow{

    // Get pangenome information
    if(params.pg_reads != ""){
        pangenome = makePangenome(file(params.pg_reads)) // Build a pangenome if reads are provided
    } else{
        pangenome = fetchPangenome(params.pg_path,params.pg_name) // Fetch an existing pangenome by ID or path
    }
    pangenome.subscribe{println(it)}
}