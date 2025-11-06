#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer
sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu
mapping_directory = file(params.final_mapping_directory)


///// Call bases from raw parquets /////
workflow callBases{

    take:
    raw_parquet_data
    
    emit:
    base_call_data

    main:

    pre_base_call_data = CALL_BASES(raw_parquet_data) | splitCsv | collect | flatten | collate(2)
    base_call_data = pre_base_call_data | checkStop | collect | flatten | collate(2)
}

workflow checkStop{
    take:
    pre_base_call_data

    emit:
    base_call_data

    main:

    base_call_data = (params.call) ? Channel.empty() : pre_base_call_data

}

process CALL_BASES{
    
    cpus sample_cpu

    tag "CallBases_${sample_id}"

    input:
    tuple val(sample_id),val(sample_parquet)

    output:
    stdout

    script:

    def base_call_script = file("${projectDir}/bin/callBases.py")
    def output_directory = file("${mapping_directory}/Base_Calls")
    def output_file = file("${output_directory}/${sample_id}_Called.parquet")

    def delete_cmd = (params.overwrite) ? "rm -rf $output_file" 
    : """
if [ -e "$output_file" ] ; then
    echo "❌ Error: ${output_file} file already exists — use --overwrite to replace." >&2
    exit 1
fi"""

    def min_depth = params.min_depth as Integer
    def allele_cov = params.min_allele as Integer
    def min_freq = params.min_freq as Float
    def ploidy
    
    if(!params.ploidy){
        error "Cannot call bases without --ploidy set"
    } else{
        ploidy = params.ploidy as Integer
    } 

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