#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = (params.validate || params.fast) ? 1 : params.cpus as Integer

workflow generateTree{
    
    take:
    input_data

    emit:
    tree_file

    main:

    tree_file = GENERATE_TREE(input_data)
}

process GENERATE_TREE{

    cpus cpu

    input:
    tuple val(input_id),val(input_dir),val(align_file)

    output:
    stdout

    script:

    def alignment_file = file(align_file)
    def tree_type = ((params.fast || params.validate)) ? "fasttree" : "iqtree"

    def tree_dir = ("${tree_type}" == "fasttree") ? file("${input_dir}/FastTree") : file("${input_dir}/iqTree")
    def tree_file = ("${tree_type}" == "fasttree") ? file("${tree_dir}/${input_id}_FastTree.nwk") : file("${tree_dir}/${input_id}.treefile")
    def log_file = ("${tree_type}" == "fasttree") ? file("${tree_dir}/out_${input_id}_FastTree") : file("${tree_dir}/out_${input_id}_iqTree")

    def bb = params.bb as Integer
    def iqtree_model = (params.gtr) ? "-m GTR+G" : "-m MFP+MERGE"

    def tree_cmd = ("${tree_type}" == "fasttree") ? "fasttree -nt -gtr -out $tree_file $alignment_file &> $log_file" : "iqtree -nt AUTO $iqtree_model -bb $bb -s $alignment_file --prefix $input_id --keep-ident &> $log_file"

    def delete_cmd = (params.overwrite) ? "rm -rf $tree_dir" : ":"

    if(!params.overwrite && tree_dir.isDirectory()){
        error "${tree_dir} exists but --overwrite is not set..."
    }

    """
    $delete_cmd &&
    mkdir $tree_dir &&
    cd $tree_dir &&
    $tree_cmd &&
    echo -n "${tree_file}"
    """
}