#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer

workflow getAlignment{
    
    take:
    filtered_data

    emit:
    alignment_file

    main:

    pre_alignment_file = GET_ALIGNMENT(filtered_data) | collect | collate(1)
    alignment_file = pre_alignment_file | checkStop
}

process GET_ALIGNMENT{
    
    executor = "local"
    cpus 1

    input:
    tuple val(filter_id),val(filter_dir)

    output:
    stdout

    script:

    def base2align_script = file("${projectDir}/bin/helper_scripts/bases2align.py")

    def base_file = file("${filter_dir}/${filter_id}_Bases.parquet")

    if(!base_file.exists()){
        error "${base_file} does not exist..."
    }

    def align_file = file("${filter_dir}/${filter_id}_aln.fasta")

    def delete_cmd = (params.overwrite) ? "rm -f $align_file" 
    : """
if [ -e "$align_file" ] ; then
    echo "❌ Error: ${align_file} already exists! Use --overwrite to replace." >&2
    exit 1
fi"""    

    """
    $delete_cmd &&
    python $base2align_script $base_file $align_file &&
    echo -n "${align_file}"
    """
}

workflow checkStop{
    take:
    pre_alignment_file

    emit:
    alignment_file

    main:

    alignment_file = (params.alignment) ? Channel.empty() : pre_alignment_file

}
