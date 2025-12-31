#! /usr/bin/env nextflow
nextflow.enable.dsl=2

def cpu = params.cpus as Integer
def sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu

def mapping_directory = file("${params.out}/Mapping")
def called_base_directory = file("${mapping_directory}/Base_Calls")

def join_directory = file("${params.out}/Joined")
def join_id = (params.join_id) ? "${params.join_id}" : ""

workflow joinCalledBases{
    
    take:
    called_bases_data

    emit:
    joined_data

    main:
    
    join_info = called_bases_data.first().map{it->tuple(join_directory,join_id)} | PREP_JOIN_DIR | splitCsv | collect | flatten | collate(2)
    
    called_base_file = called_bases_data.combine(join_info) | SAVE_CALLED_BASE_FILE | collect | map { it[0] }
    
    scaffold_parquet = CREATE_SCAFFOLD(called_base_file,join_info) | collect | map { it[0] }
    
    random_sample = called_bases_data.combine(join_info).combine(scaffold_parquet) | SCAFFOLD_SAMPLE | collect | map { it[0] }
    
    base_parquet = join_info.combine(scaffold_parquet).combine(random_sample) | COMPILE_BASES | collect | map { it[0] }
    
    joined_data = join_info.combine(scaffold_parquet).combine(base_parquet) | SCORE_SITES | splitCsv | collect | flatten | collate(2)

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
    
    delete_cmd = (params.overwrite) ? "rm -rf ${join_dir}"
    : """
if [ -d "${join_dir}" ] ; then
    echo "❌ Error: ${join_dir} already exists! Use --overwrite to replace." >&2
    exit 1
fi"""   

    """
    $delete_cmd &&
    mkdir -p $joined_dir &&
    mkdir -p "${join_dir}" &&
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
    
    called_base_file = file("${joined_dir}/${joined_id}_Called_Bases.txt")

    """
    echo "$called_base_path" >> $called_base_file
    echo -n $called_base_file 
    """
}

process CREATE_SCAFFOLD{
    cpus cpu

    input:
    val(called_base_file)
    tuple val(joined_id),val(joined_dir)

    output:
    stdout

    script:

    scaffold_script = file("${projectDir}/bin/helper_scripts/create_scaffold.py")
    scaffold_parquet = file("${joined_dir}/${joined_id}_Scaffold.parquet")
    base_call_summary = file("${joined_dir}/${joined_id}_Site_Counts.tsv")

    """
    python ${scaffold_script} --called_bases ${called_base_file} --join_id ${joined_id} --out_dir ${joined_dir} &&
    echo "Sample_ID\tFixed\tHeterozygous\tPloidy_Fail\tUncovered" > ${base_call_summary} &&
    echo -n "${scaffold_parquet}"
    """
}

process SCAFFOLD_SAMPLE{

    tag "Scaffold_${sample_id}"

    cpus sample_cpu

    input:
    tuple val(sample_id),val(called_base_path),val(joined_id),val(joined_dir),val(scaffold_parquet)

    output:
    stdout

    script:

    scaffold_sample_script = file("${projectDir}/bin/helper_scripts/scaffold_sample.py")
    sample_scaffold_file = file("${joined_dir}/Temp_${joined_id}/Scaffolded_${sample_id}.parquet")

    """
    python $scaffold_sample_script --called_bases $called_base_path --join_id $joined_id --out_dir $joined_dir --scaffold $scaffold_parquet &&
    echo -n "${sample_id}"
    """
}

process COMPILE_BASES{

    cpus cpu

    input:
    tuple val(joined_id),val(joined_dir),val(scaffold_file),val(random_sample)

    output:
    stdout

    script:

    compile_base_script = file("${projectDir}/bin/helper_scripts/compile_bases.py")
    base_parquet_file = file("${joined_dir}/${joined_id}_Bases.parquet")

    """
    python $compile_base_script --join_id $joined_id --out_dir $joined_dir &&
    echo -n "${base_parquet_file}"
    """
}

process SCORE_SITES{
    cpus cpu

    input:
    tuple val(joined_id),val(joined_dir),val(scaffold_file),val(base_file)

    output:
    stdout

    script:

    score_site_script = file("${projectDir}/bin/helper_scripts/score_sites.py")
    mem_arg = (params.mem_mode) ? "--mem_mode" : ""


    """
    python $score_site_script --out_dir $joined_dir --join_id $joined_id --scaffold $scaffold_file --base $base_file $mem_arg &&
    echo -n "${joined_id},${joined_dir}"
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
