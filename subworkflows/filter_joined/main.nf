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
    | collect | flatten | collate(4)
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
    def output_fasta = file("${filter_directory}/${filter_id}_aln.fasta")
    def output_json = file("${filter_directory}/${filter_id}.json")

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
    python $filter_script --joined $joined_directory --fasta $ref_fasta --out $filter_directory --name $filter_id --types $site_types $gap_arg $het_arg $invalid_arg $nosing_arg $missing_arg &&
    echo -n "${filter_id},${filter_directory},${output_json},${output_fasta}"
    """
}


workflow fetchFiltered{

    take:
    filtered_dir

    emit:
    filtered_info
    
    main:
    filtered_info = FETCH_FILTERED(filtered_dir) | splitCsv()
}

process FETCH_FILTERED{
    
    executor = "local"
    cpus 1

    input:
    val(filtered_dir)

    output:
    stdout

    script:

    def filtered_directory = file(filtered_dir)

    if ( !filtered_directory.isDirectory() ) {
        error "${filtered_directory} does not exist..."
    }

    def filter_id = filtered_directory.name
    def base_file = file("${filtered_directory}/${filter_id}_Bases.parquet")
    def code_file = file("${filtered_directory}/${filter_id}_Codes.parquet")
    def scaffold_file = file("${filtered_directory}/${filter_id}_Scaffold.parquet")
    def site_file = file("${filtered_directory}/${filter_id}_Sites.parquet")
    def output_fasta = file("${filtered_directory}/${filter_id}_aln.fasta")
    def output_json = file("${filtered_directory}/${filter_id}.json")
    
    """
    for f in "${base_file}" "${code_file}" "${scaffold_file}" "${site_file}" "${output_fasta}" "${output_json}" ; do
        if [ ! -s "\$f" ]; then
            echo "ERROR: Missing expected file: \$f" >&2
            exit 1
        fi
    done

    echo -n "${filter_id},${filtered_directory},${output_json},${output_fasta}"  
    """
}
