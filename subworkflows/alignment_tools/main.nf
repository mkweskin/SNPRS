#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer

workflow getAlignment{
    
    take:
    input_data

    emit:
    alignment_file

    main:

    alignment_file = GET_ALIGNMENT(input_data)
}

process GET_ALIGNMENT{
    executor = "local"
    cpus 1

    input:
    tuple val(input_id),val(input_dir)

    output:
    stdout

    script:

    def base2align_script = file("${projectDir}/bin/helper_scripts/bases2align.py")
    def align_file = file("${input_dir}/${input_id}_aln.fasta")

    def align_exists = align_file.isFile()

    def base_file = file("${input_dir}/${input_id}_Bases.parquet")

    if(!base_file.exists()){
        error "${base_file} does not exist..."
    }

    delete_cmd = (params.overwrite) ? "rm -f $align_file" : ":"

    def generate_cmd
    if(params.overwrite || !align_exists){
        generate_cmd = "python $base2align_script $base_file $align_file" 
    } else{
        generate_cmd = ":"
    }

    """
    $delete_cmd &&
    $generate_cmd &&
    echo -n "${align_file}"
    """
}