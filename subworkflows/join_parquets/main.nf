#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer

workflow joinCalledBases{
    
    take:
    called_bases_data
    pangenome_info
    output_dir
    join_id

    emit:
    joined_data

    main:
    
    head_output_dir = CHECK_OUTPUT_DIR(output_dir,join_id) | splitCsv | collect | flatten | collate(2) | first
    prep_join = pangenome_info.combine(head_output_dir)
    combined = called_bases_data.combine(prep_join)
    
    called_base_file = SAVE_CALLED_BASE_FILE(combined) | first


    joined_data = JOIN_CALLED_BASES(called_base_file,prep_join) | splitCsv
}

process CHECK_OUTPUT_DIR{
    
    executor = "local"
    cpus 1

    input:
    val(output_dir)
    val(join_id)

    output:
    stdout

    script:
    def join_dir = file("${output_dir}/${join_id}")
    def delete_cmd = (params.overwrite) ? "rm -rf $join_dir" : ":"

    """
    $delete_cmd &&
    echo -n $output_dir,$join_id
    """
}

process SAVE_CALLED_BASE_FILE{
    
    executor = "local"
    cpus 1

    input:
    tuple val(sample_id),val(called_base_path),val(pg_name),val(pg_fasta),val(head_output_dir),val(join_id)

    output:
    stdout

    script:
    def join_dir = file("${head_output_dir}/${join_id}")
    def called_base_file = file("${join_dir}/${join_id}_Called_Bases.txt")

    if(params.validate && !params.overwrite && file(join_dir).isDirectory()){
        error "Running in validation mode without --overwrite set, but ${join_dir} exists..."
    }

    """
    mkdir -p $join_dir &&
    echo "$called_base_path" >> $called_base_file
    echo -n $called_base_file 
    """
}

process JOIN_CALLED_BASES{
    cpus cpu

    input:
    val(called_base_file)
    tuple val(pg_name),val(pg_fasta),val(head_output_dir),val(join_id)

    output:
    stdout

    script:

    def join_script = file("${projectDir}/bin/join_parquets.py")
    def base_file = file(called_base_file)
    def join_dir = file("${head_output_dir}/${join_id}")

    """
    python $join_script -b $base_file -n $join_id -o $join_dir &&
    echo -n "$join_id,$join_dir"
    """
}

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
      "_Bases.parquet"
      "_Called_Bases.txt"
      "_Codes.parquet"
      "_Missing.tsv"
      "_Scaffold.parquet"
      "_Sites.parquet"
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
