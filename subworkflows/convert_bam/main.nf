#! /usr/bin/env nextflow
nextflow.enable.dsl=2

def cpu = params.cpus as Integer
def sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu
def mapping_directory = file("${params.out}/Mapping")
def raw_parquet_directory = file("${mapping_directory}/Raw_Parquets")

def mapq = params.mapq as Integer
def baseq = params.baseq as Integer
def adj_coef = params.adj_coef as Integer

///// Convert BAM files to Parquet /////
workflow bamToParquet{

    take:
    convert_data
    
    emit:
    raw_parquet_data

    main:
    raw_parquet_data = BAM_TO_PARQUET(convert_data) | splitCsv
}

process BAM_TO_PARQUET{
    
    cpus sample_cpu
    
    tag "BAM2PQ_${sample_id}"

    input:
    tuple val(sample_id),val(sample_bam),val(genome_file)

    output:
    stdout

    script:

    bam_convert_script = file("${projectDir}/bin/bam2parquet.py")

    output_file = file("${raw_parquet_directory}/${sample_id}_Raw.parquet")
    
    delete_cmd = (params.overwrite)
    ? "rm -f $output_file"
    : """
if [ -e "$output_file" ] ; then
    echo "❌ Error: ${output_file} already exists — use --overwrite to replace." >&2
    exit 1
fi"""

    """
    mkdir -p ${raw_parquet_directory} &&
    $delete_cmd &&
    python ${bam_convert_script} --bam ${sample_bam} --fasta ${genome_file} --parquet ${output_file} --mapq ${mapq} --baseq ${baseq} --adj_coef ${adj_coef} &&
    echo -n "${sample_id},${output_file}"
    """
}

///// Fetch Raw Parquets /////
workflow fetchRawParquet{
    take:
    input_raw_parquet

    emit:
    raw_parquet_data

    main:
    raw_parquet_data = FETCH_RAW_PARQUET(input_raw_parquet) | splitCsv()
}

process FETCH_RAW_PARQUET{
    executor = "local"
    cpus = 1

    input:
    val(input_raw_parquet)

    output:
    stdout

    script:

    fetch_raw_parquet_script = file("${projectDir}/bin/fetchRawParquet.py")
    full_parquet = file("${input_raw_parquet}")
    """
    python ${fetch_raw_parquet_script} -p ${full_parquet}
    """
}
