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

workflow joinCalledBases{
    
    take:
    called_bases_data

    emit:
    joined_data

    main:
    
    join_info = called_bases_data.first().map{it->tuple(join_directory,join_id)} | PREP_JOIN_DIR | splitCsv | collect | flatten | collate(2)
    called_base_file = called_bases_data.combine(join_info) | SAVE_CALLED_BASE_FILE | collect | map { it[0] }
    scaffold_info = join_info.combine(called_base_file) | CREATE_SCAFFOLD | splitCsv | collect | flatten | collate(2)
    chunk_directories = called_bases_data.combine(join_info).combine(scaffold_info) | SCAFFOLD_SAMPLE | collect | map { it[0] } | VALIDATE_CHUNKS | map {it.split('\n')} | collect | flatten | collate(1)
    processed_chunk_dirs = chunk_directories.combine(scaffold_info) | SCORE_CHUNK | collect | map { it[0] }
    joined_data = join_info.combine(scaffold_info).combine(processed_chunk_dirs) | COMPILE_DATA | splitCsv | collect | flatten | collate(2)
}

process PREP_JOIN_DIR{
    
    executor = "local"
    cpus 1

    input:
    tuple val(joined_dir), val(joined_id)

    output:
    stdout

    script:
    
    join_dir = file("${joined_dir}/${joined_id}")
    temp_dir = file("${join_dir}/Temp_${joined_id}")
    
    delete_cmd = (params.overwrite) ? "rm -rf ${join_dir}"
    : """
if [ -d "${join_dir}" ] ; then
    echo "❌ Error: ${join_dir} already exists! Use --overwrite to replace." >&2
    exit 1
fi"""   

    """
    $delete_cmd &&
    mkdir -p $joined_dir &&
    mkdir $join_dir && 
    mkdir $temp_dir &&
    echo -n "${joined_id},${join_dir}"
    """
}

process SAVE_CALLED_BASE_FILE{
    
    executor = "local"
    cpus 1
    maxForks 0

    input:
    tuple val(sample_id),val(called_base_path),val(joined_id),val(joined_dir)
    
    output:
    stdout

    script:
    
    called_base_file = file("${joined_dir}/Temp_${joined_id}/${joined_id}_Called_Bases.txt")

    """
    echo "$called_base_path" >> $called_base_file
    echo -n $called_base_file 
    """
}

process CREATE_SCAFFOLD{
    cpus cpu
    
    input:
    tuple val(joined_id),val(joined_dir),val(called_base_file)

    output:
    stdout

    script:

    scaffold_script = file("${projectDir}/bin/helper_scripts/create_scaffold.py")

    """
    python ${scaffold_script} --called_bases ${called_base_file} --join_id ${joined_id} --out_dir ${joined_dir} --mem_factor $memory_factor
    """
}

process SCAFFOLD_SAMPLE{

    tag "Scaffold_${sample_id}"

    cpus sample_cpu

    input:
    tuple val(sample_id),val(called_base_path),val(joined_id),val(joined_dir),val(scaffold_parquet),val(chunk_tsv)

    output:
    stdout

    script:

    scaffold_sample_script = file("${projectDir}/bin/helper_scripts/scaffold_sample.py")

    """
    python $scaffold_sample_script --called_bases $called_base_path --join_id $joined_id --scaffold_parquet $scaffold_parquet --out_dir $joined_dir --chunk_tsv $chunk_tsv
    """
}

process VALIDATE_CHUNKS {
    
    executor = "local"
    cpus 1
    maxForks 0

    input:
    val(chunk_tsv)

    output:
    stdout

    script:
    validate_chunk_script = file("${projectDir}/bin/helper_scripts/validate_chunks.py")

    """
    python $validate_chunk_script $chunk_tsv
    """
}

process SCORE_CHUNK{
    cpus cpu

    input:
    tuple val(chunk_dir),val(scaffold_parquet),val(chunk_tsv)

    output:
    stdout

    script:

    score_chunk_script = file("${projectDir}/bin/helper_scripts/score_chunk.py")

    """
    python $score_chunk_script --chunk_dir $chunk_dir --scaffold_parquet $scaffold_parquet --chunk_tsv $chunk_tsv
    """
}

process COMPILE_DATA{

    cpus cpu

    input:
    tuple val(joined_id),val(joined_dir),val(scaffold_file),val(chunk_tsv),val(random_dir)

    output:
    stdout

    script:

    compile_join_script = file("${projectDir}/bin/helper_scripts/compile_join.py")
    summarize_join_script = file("${projectDir}/bin/helper_scripts/summarize_join.py")

    """
    python $compile_join_script --join_id $joined_id --out_dir $joined_dir &&
    echo -n ${join_id},${joined_dir}
    """
}

//////////////////////////////////////////////////////////////////////

workflow joinFromCSV{
    
    take:
    join_csv

    emit:
    joined_data

    main:

    joined_data = join_csv

}

//////////////////////////////////////////////////////////////////////



workflow fetchJoin{

    take:
    join_dir

    emit:
    join_info
    
    main:
    join_info = FETCH_JOIN(join_dir) | splitCsv()
}

/// FIX ///
workflow fetchScaffold{

    take:
    scaffold_file

    emit:
    scaffold_info
    
    main:
    scaffold_info = scaffold_file
}


process FETCH_JOIN{
    
    executor = "local"
    cpus 1

    input:
    val(join_dir)

    output:
    stdout

    script:

    join_path = file(join_dir)

    if(!join_path.isDirectory()){
        error "${join_path} does not exist..."
    }

    """
    cd $join_path

    suffixes=(
      "_Scaffold.parquet"
      "_Bases.parquet"
      "_Codes.parquet"
      "_Sites.parquet"
      "_Missing.tsv"
    )

    prefixes=()

    for suf in "\${suffixes[@]}"; do
        files=( *"\$suf" )
        if [ \${#files[@]} -ne 1 ] || [ ! -f "\${files[0]}" ]; then
            echo "ERROR: Expected exactly one file ending with \$suf in $join_path" >&2
            exit 1
        fi
        fname="\${files[0]}"
        prefix="\${fname%\$suf}"
        prefixes+=( "\$prefix" )
    done

    first="\${prefixes[0]}"
    for p in "\${prefixes[@]}"; do
        if [ "\$p" != "\$first" ]; then
            echo "ERROR: Prefix mismatch in $join_path" >&2
            exit 1
        fi
    done

    echo -n "\$first,$join_path"    
    """
}
