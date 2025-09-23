#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer

workflow filterJoined{
    take:
    pangenome_info
    joined_data
    filter_id

    emit:
    filtered_data

    main:

    filtered_data = pangenome_info
    .combine(joined_data)
    .map{it-> tuple(it[0],it[1],it[2],it[3],filter_id)}
    | FILTER_JOINED
    | splitCsv
    | collect | flatten | collate(2)
}

process FILTER_JOINED{

    cpus cpu

    input:
    tuple val(pg_name),val(pg_fasta),val(join_id),val(join_directory),val(filter_id)

    output:
    stdout

    script:

    def filter_script = file("${projectDir}/bin/filter_joined.py")
    def ref_fasta = file("${pg_fasta}")
    def joined_directory = file("${join_directory}")
    def filter_directory = file("${joined_directory}/${filter_id}")

    def site_types = (params.site_types) ? "${params.site_types}" : "btqp"

    // Flags
    def gap_arg = (params.gaps) ? "--gaps" : ""
    def het_arg = (params.het) ? "--het" : ""
    def invalid_arg = (params.invalid) ? "--invalid" : ""
    def nosing_arg = (params.nosing) ? "--nosing" : ""
    def missing_arg = (params.missing != false) ? "--missing ${params.missing}" : ""

    def delete_cmd = (params.overwrite) ? "rm -rf $filter_directory" : ":"
    """
    $delete_cmd &&
    mkdir $filter_directory &&
    python $filter_script --joined $joined_directory --fasta $ref_fasta --out $filter_directory --name $filter_id --alignment --types $site_types $gap_arg $het_arg $invalid_arg $nosing_arg $missing_arg &&
    echo -n "${filter_id},${filter_directory}"
    """
}
