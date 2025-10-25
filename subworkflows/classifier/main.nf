#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer

workflow classifyCalledBases{

    take:
    snp_data
    called_bases_data

    emit:
    classified_data

    main:

    classified_data = snp_data.join(called_bases_data)
}