#! /usr/bin/env nextflow
nextflow.enable.dsl=2

def cpu = params.cpus as Integer
def sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu
def chunk_size = params.chunk_size as Integer
def memory_factor = params.mem_factor as Integer

def mapping_directory = file("${params.out}/Mapping")
def called_base_directory = file("${mapping_directory}/Base_Calls")

def join_directory = file("${params.out}/Joined")
def join_id = (params.join_id) ? "${params.join_id}" : ""

def filter_id = (params.filter_id) ? "${params.filter_id}" : ""

workflow generateScaffold{
    
    take:
    called_bases_data

    emit:
    joined_data

    main:

    join_info = called_bases_data | first | INITIATE_JOIN | splitCsv | collect | flatten | collate(2)
        
    called_base_file = called_bases_data.combine(join_info).map { sample_id, called_base_path, joined_id, joined_dir ->
        def out = file("${joined_dir}/${joined_id}_Called_Bases.txt")
        out << (called_base_path + "\n")        
        return out
    }.last()
    
    chunk_parquets = join_info.combine(called_base_file) | INITIATE_SCAFFOLD | splitText | SCORE_SCAFFOLD_CHUNK | collect | flatten | collate(1)

    chunk_parquet_file = chunk_parquets.combine(join_info).map { chunk_parquet, joined_id, joined_dir ->
        def out = file("${joined_dir}/${joined_id}_Chunk_Parquets.txt")
        out << (chunk_parquet + "\n")        
        return out
    }.last()

    joined_data = join_info.combine(chunk_parquet_file) | COMPILE_SCAFFOLD | collect | flatten | collate(1)
}

process INITIATE_JOIN{
    executor = "local"
    cpus 1

    input:
    val(dummy_val)

    output:
    stdout

    script:

    join_dir = file("${join_directory}/${join_id}")

    """
    if [ -d "${join_dir}" ]; then
        echo "ERROR: directory ${join_dir} already exists!" >&2
        exit 1
    fi
    mkdir -p ${join_directory} &&
    mkdir ${join_dir} &&
    echo -n ${join_id},${join_dir}
    """
}

process INITIATE_SCAFFOLD{
    executor = "local"
    cpus 1
        
    input:
    tuple val(joined_id),val(joined_dir),val(called_base_file)

    output:
    stdout

    script:

    initiate_scaffold_script = file("${projectDir}/bin/manual/official/generateChunks.py")

    """
    python ${initiate_scaffold_script} --called ${called_base_file} --name ${joined_id} --output ${joined_dir} --batch ${chunk_size}
    """
}

process SCORE_SCAFFOLD_CHUNK{
    cpus 5
        
    input:
    val(scaffold_chunk)

    output:
    stdout

    script:

    score_scaffold_script = file("${projectDir}/bin/manual/official/processChunk.py")

    """
    python ${score_scaffold_script} --txt ${scaffold_chunk}
    """
}

process COMPILE_SCAFFOLD{
    cpus cpu

    input:
    tuple val(joined_id),val(joined_dir),val(chunk_parquet_file)

    output:
    stdout

    script:

    compile_scaffold_script = file("${projectDir}/bin/manual/official/compileChunks.py")
    output_file = file("${joined_dir}/${join_id}_Scaffold.parquet")

    """
    python ${compile_scaffold_script} --parquets ${chunk_parquet_file} --out ${output_file} --batch ${chunk_size} --join_id ${joined_id}
    echo -n "${output_file}"
    """
}

workflow fetchScaffold{

    take:
    scaffold_file

    emit:
    scaffold_info
    
    main:
    scaffold_info = file(scaffold_file) | collect | flatten | collate(1)
}