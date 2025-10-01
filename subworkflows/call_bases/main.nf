#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer
sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu

///// Call bases from raw parquets /////
workflow callBases{

    take:
    raw_parquet_data
    pangenome_info
    mapping_directory
    
    emit:
    base_call_data

    main:

    base_call_data = CALL_BASES(raw_parquet_data,pangenome_info,mapping_directory) | splitCsv()
}

process CALL_BASES{
    
    cpus sample_cpu

    tag "CallBases_${sample_id}"

    input:
    tuple val(sample_id),val(sample_parquet)
    tuple val(pg_name),val(fasta_path)
    val(mapping_directory)

    output:
    stdout

    script:

    def base_call_script = file("${projectDir}/bin/callBases.py")
    def output_directory = file("${mapping_directory}/Base_Calls")
    def output_file = file("${output_directory}/${sample_id}_Called.parquet")

    def delete_cmd = (params.overwrite) ? "rm -rf $output_file" : ":"
    def min_depth = params.min_read as Integer
    def allele_cov = params.min_allele as Integer
    def min_freq = params.min_freq as Float
    def ploidy = params.ploidy as Integer

    """
    mkdir -p ${output_directory} &&
    $delete_cmd &&
    python ${base_call_script} -p ${sample_parquet} -o ${output_file} -min_depth ${min_depth} -min_support ${allele_cov} -min_freq ${min_freq} -max_alleles ${ploidy} &&
    echo -n "${sample_id},${output_file}"
    """
}

///// Fetch Called Bases /////
workflow fetchCalledBases{
    take:
    input_called_bases

    emit:
    called_base_data

    main:
    called_base_data = FETCH_CALLED_BASES(input_called_bases) | splitCsv()
}

process FETCH_CALLED_BASES{
    executor = "local"
    cpus = 1

    input:
    val(input_called_bases)

    output:
    stdout

    script:

    def fetch_called_bases_script = file("${projectDir}/bin/fetchCalledBases.py")
    def full_parquet = file("${input_called_bases}")
    """
    python ${fetch_called_bases_script} -p ${full_parquet}
    """
}