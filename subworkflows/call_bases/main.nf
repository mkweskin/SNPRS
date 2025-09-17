#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer

///// Call bases from raw parquets /////
workflow callBases{

    take:
    raw_parquet_data
    pangenome_info
    base_call_directory
    
    emit:
    base_call_data

    main:

    base_call_dir = file(base_call_directory)
    base_call_data = CALL_BASES(raw_parquet_data,pangenome_info,base_call_directory) | splitCsv()
}

process CALL_BASES{
    
    cpus cpu

    input:
    tuple val(sample_id),val(sample_parquet)
    tuple val(pg_name),val(fasta_path)
    val(out_dir)

    output:
    stdout

    script:

    def base_call_script = file("${projectDir}/bin/call_bases.py")
    def output_directory = params.cb_out ? file("${params.cb_out}") : file("${out_dir}/${pg_name}/Base_Calls")
    def output_file = file("${output_directory}/${sample_id}_Called.parquet")
    def full_parquet = file("${sample_parquet}")

    def min_depth = params.rd as Integer
    def allele_cov = params.ad as Integer
    def min_freq = params.min_freq as Float
    def ploidy = params.ploidy as Integer

    """
    
    if [[ -f "${output_file}" ]]; then
        echo "${output_file} already exists" >&2
        exit 1
    fi

    mkdir -p ${output_directory} &&
    python ${base_call_script} -p ${full_parquet} -o ${output_file} -min_depth ${min_depth} -min_support ${allele_cov} -min_freq ${min_freq} -max_alleles ${ploidy} &&
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

// Base workflow/processes if running separate. Must provide --raw_parquet and --fasta
workflow{
    if (params.raw_parquet && params.fasta){

        cb_out = params.cb_out ? file("${params.cb_out}") : file("SNPRS_${params.timestamp}")
        
        raw_parquet_data = FETCH_RAW_PARQUET(params.raw_parquet) | splitCsv()
        pangenome_data = CHECK_FASTA(params.fasta) | splitCsv() | first()

        called_base_data = CALL_BASES(raw_parquet_data,pangenome_data,cb_out)
    }
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
