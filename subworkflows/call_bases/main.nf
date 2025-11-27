#! /usr/bin/env nextflow
nextflow.enable.dsl=2

def cpu = params.cpus as Integer
def sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu
def mapping_directory = file("${params.out}/Mapping")
def called_directory = file("${mapping_directory}/Base_Calls")

def min_depth = params.min_depth as Integer
def allele_cov = params.min_allele as Integer
def min_freq = params.min_freq as Float
def ploidy = (params.ploidy) ? params.ploidy as Integer : 0

///// Call bases from raw parquets /////
workflow callBases{

    take:
    raw_parquet_data
    
    emit:
    base_call_data

    main:

    base_call_data = CALL_BASES(raw_parquet_data) | splitCsv
}


process CALL_BASES{
    
    cpus sample_cpu

    tag "CallBases_${sample_id}"

    input:
    tuple val(sample_id),val(sample_parquet)

    output:
    stdout

    script:

    base_call_script = file("${projectDir}/bin/callBases.py")
    output_file = file("${called_directory}/${sample_id}_Called.parquet")

    delete_cmd = (params.overwrite) ? "rm -rf $output_file" 
    : """
if [ -e "$output_file" ] ; then
    echo "❌ Error: ${output_file} file already exists — use --overwrite to replace." >&2
    exit 1
fi"""

    """
    mkdir -p ${called_directory} &&
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

    fetch_called_bases_script = file("${projectDir}/bin/fetchCalledBases.py")

    """
    python ${fetch_called_bases_script} -p ${input_called_bases}
    """
}