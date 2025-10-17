#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer

workflow joinCalledBases{
    
    take:
    called_bases_data
    pangenome_info
    joined_dir
    join_id

    emit:
    joined_data

    main:
    
    join_info = PREP_JOIN_DIR(joined_dir,join_id) | splitCsv | collect | flatten | collate(2)
    
    called_base_file = called_bases_data.combine(join_info) | SAVE_CALLED_BASE_FILE | collect | map { it[0] }
    
    scaffold_parquet = CREATE_SCAFFOLD(called_base_file,join_info) | collect | map { it[0] }
    
    prep_scaffold = join_info.combine(scaffold_parquet)
    random_sample = called_bases_data.combine(prep_scaffold) | SCAFFOLD_SAMPLE | collect | map { it[0] }
    
    prep_base = join_info.combine(random_sample)
    base_parquet = CREATE_BASE_PARQUET(prep_base) | collect | map { it[0] }

    parquet_files = scaffold_parquet.combine(base_parquet)
    prep_code = join_info.combine(parquet_files)
    
    joined_data = SCORE_SITES(prep_code) | splitCsv
}

process PREP_JOIN_DIR{
    
    executor = "local"
    cpus 1

    input:
    val(joined_dir)
    val(join_id)

    output:
    stdout

    script:
    def join_dir = file("${joined_dir}/${join_id}")
    def delete_cmd = (params.overwrite) ? "rm -rf $join_dir" : ":"

    if(!params.overwrite && file(join_dir).isDirectory()){
        error "$join_dir already exists, use --overwrite to remove existing directory..."
    }

    """
    $delete_cmd &&
    mkdir -p $join_dir &&
    echo -n $join_id,$join_dir
    """
}

process SAVE_CALLED_BASE_FILE{
    
    executor = "local"
    cpus 1

    input:
    tuple val(sample_id),val(called_base_path),val(join_id),val(join_dir)

    output:
    stdout

    script:
    def called_base_file = file("${join_dir}/${join_id}_Called_Bases.txt")

    """
    echo "$called_base_path" >> $called_base_file
    echo -n $called_base_file 
    """
}

process CREATE_SCAFFOLD{
    cpus cpu

    input:
    val(called_base_file)
    tuple val(join_id),val(join_dir)

    output:
    stdout

    script:

    def scaffold_script = file("${projectDir}/bin/helper_scripts/create_scaffold.py")
    def scaffold_parquet = file("${join_dir}/${join_id}_Scaffold.parquet")
    def base_call_summary = file("${join_dir}/${join_id}_Site_Counts.tsv")

    """
    python $scaffold_script --called_bases $called_base_file --join_id $join_id --out_dir $join_dir &&
    echo "Sample_ID\tFixed_Bases\tFixed_Gaps\tHet_Bases\tHet_Gap\tUncovered\tFiltered" > $base_call_summary &&
    echo -n "$scaffold_parquet"
    """
}

process SCAFFOLD_SAMPLE{
    cpus cpu

    input:
    tuple val(sample_id),val(called_base_path),val(join_id),val(join_dir),val(scaffold_parquet)

    output:
    stdout

    script:

    def scaffold_sample_script = file("${projectDir}/bin/helper_scripts/scaffold_sample.py")
    def sample_scaffold_file = file("${join_dir}/Scaffolded_${sample_id}.parquet")

    def delete_cmd = (params.overwrite) ? "rm -f $sample_scaffold_file" : ":"

    """
    $delete_cmd &&
    python $scaffold_sample_script --called_bases $called_base_path --join_id $join_id --out_dir $join_dir --scaffold $scaffold_parquet &&
    echo -n "${sample_id}"
    """
}

process CREATE_BASE_PARQUET{
    cpus cpu

    input:
    tuple val(join_id),val(join_dir),val(random_sample)

    output:
    stdout

    script:

    def base_parquet_script = file("${projectDir}/bin/helper_scripts/compile_bases.py")
    def base_parquet = file("${join_dir}/${join_id}_Bases.parquet")

    def delete_cmd = (params.overwrite) ? "rm -f $base_parquet" : ":"

    """
    $delete_cmd &&
    python $base_parquet_script --out_dir $join_dir --join_id $join_id &&
    echo -n "${base_parquet}"
    """
}

process SCORE_SITES{
    cpus cpu

    input:
    tuple val(join_id),val(join_dir),val(scaffold_file),val(base_file)

    output:
    stdout

    script:

    def score_site_script = file("${projectDir}/bin/helper_scripts/score_sites.py")

    def site_parquet = file("${join_dir}/${join_id}_Sites.parquet")
    def code_parquet = file("${join_dir}/${join_id}_Codes.parquet")
    def missing_tsv = file("${join_dir}/${join_id}_Missing.tsv")

    def mem_arg = (params.mem_mode) ? "--mem_mode" : ""

    def delete_cmd = (params.overwrite) ? "rm -f $site_parquet $code_parquet $missing_tsv" : ":"

    """
    $delete_cmd &&
    python $score_site_script --out_dir $join_dir --join_id $join_id --bases $base_file --scaffold $scaffold_file $mem_arg &
    echo -n "${join_id},${join_dir}"
    """
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

    def join_directory = file(join_dir)

    if(!join_directory.isDirectory()){
        error "${join_directory} does not exist..."
    }

    """
    cd $join_dir

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
            echo "ERROR: Expected exactly one file ending with \$suf in $join_dir" >&2
            exit 1
        fi
        fname="\${files[0]}"
        prefix="\${fname%\$suf}"
        prefixes+=( "\$prefix" )
    done

    first="\${prefixes[0]}"
    for p in "\${prefixes[@]}"; do
        if [ "\$p" != "\$first" ]; then
            echo "ERROR: Prefix mismatch in $join_dir" >&2
            exit 1
        fi
    done

    echo -n "\$first,$join_dir"    
    """
}
