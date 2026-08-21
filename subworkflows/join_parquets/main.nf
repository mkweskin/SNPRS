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

// Filtering
def filter_id = (params.filter_id) ? "${params.filter_id}" : ""

// Distances
def dist_id = (params.dist_id) ? "${params.dist_id}" : ""

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
    scaffold_info = scaffold_file
}

workflow filterScaffold{

    take:
    scaffold_file

    emit:
    filtered_file
    
    main:
    filtered_file = FILTER_SCAFFOLD(scaffold_file) | collect | flatten | collate(1)
}

process FILTER_SCAFFOLD{

    cpus 1

    input:
    val(scaffold_file)

    output:
    stdout

    script:
    
    filter_script = file("${projectDir}/bin/manual/semioffiicial/join_tools/filter_scaffold.py")

    scaf_file = file("${scaffold_file}")
    cov_arg = (params.min_cov) ? " --covered ${params.min_cov}" : ""
    fixed_arg = (params.min_fix) ? " --fixed ${params.min_fix}" : ""
    het_arg = (params.min_het) ? " --het ${params.min_het}" : ""
    ploidy_arg = (params.max_pf) ? " --ploidy ${params.max_pf}" : ""
    clade_arg = (params.min_clade) ? " --min_clade ${params.min_clade}" : ""
    allele_arg = ("${params.alleles}") ? " --alleles ${params.alleles}" : ""
    no_gap_arg = (params.no_gaps) ? " --no_gap": ""
    no_sing_arg = (params.no_sing) ? " --no_sing": ""
    no_het_arg = (params.no_het) ? " --no_het": ""

    """
    python $filter_script --in $scaf_file --filter_id $filter_id $cov_arg $fixed_arg $het_arg $ploidy_arg $clade_arg $allele_arg $no_gap_arg $no_sing_arg $no_het_arg
    """
}

workflow getDistance{

    take:
    called_bases_data
    scaffold_file

    emit:
    final_phylip
    
    main:

    sample_count = called_bases_data.map { it.size() }

    blank_data = CREATE_MATRIX(scaffold_file,sample_count,dist_id) 
    | splitCsv | collect | flatten | collate(4)

    called_base_file = blank_data.combine(called_bases_data)
    .map { output_directory, scaffold, matrix, list, sample_id, called_base_path ->
        def out = file("${list}")
        out << "${called_base_path}\n"
        return out
    }.last()
    
    chunk_file = called_base_file.combine(blank_data) | CHUNK_DATA | splitCsv | collect | flatten | collate(1)

    chunk_ids = chunk_file.combine(called_bases_data).combine(blank_data) | POPULATE_MATRIX | unique | collect | flatten | collate(1)

    chunked_dist = chunk_ids.combine(blank_data).combine(chunk_file) | CHUNK_DIST | collect | flatten | collate(1) | last()

    final_phylip = blank_data.combine(chunked_dist).combine(chunk_file) | CREATE_PHY

    final_tree = RUN_RAPIDNJ(final_phylip)
}

process CREATE_MATRIX{

    cpus 1
    executor "slurm"
    memory  "1G"

    input:
    val(scaffold_file)
    val(sample_count)
    val(dist_id)

    output:
    stdout

    script:
    
    create_matrix_script = file("${projectDir}/bin/manual/semioffiicial/join_tools/createMatrix.py")
    def scaffold_dir = file("${scaffold_file}").getParent()
    def output_dir = file("${scaffold_dir}/${dist_id}")

    """
    mkdir -p $output_dir &&
    python $create_matrix_script -s $scaffold_file -n $dist_id -c $sample_count
    """
}

process CHUNK_DATA{

    executor "local"
    cpus 1

    input:
    tuple val(called_base_file), val(output_directory),val(scaffold_file),val(matrix_file),val(list)

    output:
    stdout

    script:
    
    chunk_sample_script = file("${projectDir}/bin/manual/semioffiicial/join_tools/chunkSamples.py")

    """
    python $chunk_sample_script -c $called_base_file -s $chunk_size -n $dist_id
    """
}

process POPULATE_MATRIX{

    executor "slurm"
    cpus 1
    memory "1G"

    input:
    tuple val(chunk_file),val(sample_id),val(called_base_file),val(output_directory),val(scaffold_file),val(matrix_file),val(list)

    output:
    stdout

    script:
    
    pop_matrix_script = file("${projectDir}/bin/manual/semioffiicial/join_tools/populateMatrix.py")
    """
    python $pop_matrix_script -i $sample_id -b $called_base_file -m $matrix_file -c $chunk_file -s $scaffold_file
    """
}

process CHUNK_DIST{

    executor "slurm"
    cpus cpu

    input:
    tuple val(chunk_id),val(output_directory),val(scaffold_file),val(matrix_file),val(list),val(chunk_file)

    output:
    stdout

    script:
    
    chunk_dist_script = file("${projectDir}/bin/manual/semioffiicial/join_tools/batch_dist_worker.py")
    output_file = file("${output_directory}/Chunk_Dist_${chunk_id}.npy")

    raw_arg = (params.raw) ? " --raw" : ""

    """
    python $chunk_dist_script --chunk_id $chunk_id --chunk_size $chunk_size --chunk_file $chunk_file --out $output_file --processors $cpu --scaffold $scaffold_file --matrix $matrix_file $raw_arg --out_dir $output_directory
    """
}

process CREATE_PHY{

    executor "slurm"
    cpus 1

    input:
    tuple val(output_directory),val(scaffold_file),val(matrix_file),val(list),val(random_dist),val(chunk_file)

    output:
    stdout

    script:
    
    combine_script = file("${projectDir}/bin/manual/semioffiicial/join_tools/combine_phy.py")
    phylip_file = file("${output_directory}/${dist_id}.phylip")
    """
    python $combine_script -i $output_directory -o $phylip_file -c $chunk_file
    """
}

process RUN_RAPIDNJ{

    executor "slurm"
    cpus cpu

    input:
    val(phylip_file)

    output:
    stdout

    script:
    
    newick_file = "${file(phylip_file).getParent()}/${dist_id}.nwk"
    """
    mkdir -p MEM
    rapidnj -i pd -d ./MEM -c 48 -n -x $newick_file $phylip_file
    """
}