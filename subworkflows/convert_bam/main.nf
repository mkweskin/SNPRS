#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer

///// Convert BAM files to Parquet /////
workflow bamToParquet{

    take:
    bam_data
    pangenome_info
    parquet_directory
    
    emit:
    raw_parquet_data

    main:

    parquet_dir = file(parquet_directory)
    raw_parquet_data = BAM_TO_PARQUET(bam_data,pangenome_info,parquet_dir) | splitCsv()
}

process BAM_TO_PARQUET{
    
    cpus cpu

    input:
    tuple val(sample_id),val(sample_bam)
    tuple val(pg_name),val(fasta_path)
    val(out_dir)

    output:
    stdout

    script:

    def bam_convert_script = file("${projectDir}/bin/bam2parquet.py")
    def output_directory = params.pq_out ? file("${params.pq_out}") : file("${out_dir}/${pg_name}/Raw_Parquet")
    def output_file = file("${output_directory}/${sample_id}_Raw.parquet")
    def full_bam = file("${sample_bam}")
    def full_fasta = file("${fasta_path}")

    def mapq = params.mapq as Integer
    def baseq = params.baseq as Integer
    def adj_coef = params.adj_coef as Integer

    def dup = params.markdup ? "--dup" : ""
    """
    
    if [[ -f "${output_file}" ]]; then
        echo "${output_file} already exists" >&2
        exit 1
    fi

    mkdir -p ${output_directory} &&
    python ${bam_convert_script} -b ${full_bam} -f ${full_fasta} -p ${output_file} --mapq ${mapq} --baseq ${baseq} --adj_coef ${adj_coef} ${dup} &&
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

    def fetch_raw_parquet_script = file("${projectDir}/bin/fetchRawParquet.py")
    def full_parquet = file("${input_raw_parquet}")
    """
    python ${fetch_raw_parquet_script} -p ${full_parquet}
    """
}

// Base workflow/processes if running separate. Must provide --bam_files and --fasta
workflow{
    if (params.bam_files && params.fasta){

        pq_out = params.pq_out ? file("${params.pq_out}") : file("SNPRS_${params.timestamp}")
        
        bam_data = FETCH_BAM(params.bam_files) | splitCsv()
        pangenome_data = CHECK_FASTA(params.fasta) | splitCsv() | first()

        raw_parquet_data = BAM_TO_PARQUET(bam_data,pangenome_data,pq_out)
    }
}

process FETCH_BAM{
    executor = "local"
    cpus = 1

    input:
    val(input_bam_files)

    output:
    stdout

    script:

    def fetchBAMScript = file("${projectDir}/bin/fetchBAMs.py")
    full_bam = file("${input_bam_files}")
    """
    python ${fetchBAMScript} -b ${full_bam}
    """
}

process CHECK_FASTA{
    executor = "local"
    cpus 1

    input:
    val(fasta_path)

    output:
    stdout

    script:
    
    def fasta_file = file("${fasta_path}") 
    def fai_file = file("${fasta_path}.fai") 
    def fasta_dir = fasta_file.getParent()
    def fasta_name = fasta_file.getName()
    def fasta_basename = fasta_file.getBaseName()
    
    """
    if [[ ! -f "${fasta_file}" ]]; then
        echo "${fasta_file} does not exist" >&2
        exit 1
    fi

    if [[ ! -f "${fai_file}" ]]; then
        echo "${fai_file} does not exist" >&2
        exit 1
    fi

    echo -n "${fasta_basename},${fasta_file}"
    """
}
