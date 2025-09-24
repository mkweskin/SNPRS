#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = (params.validate || params.fast) ? 1 : params.cpus as Integer

workflow generateTree{
    
    take:
    filtered_data

    emit:
    tree_data

    main:

    tree_type = (params.validate || params.fast) ? "fasttree" : "iqtree"

    tree_data = filtered_data
    .map{it-> tuple(it[0],it[1],it[3],"${tree_type}")} 
    | GENERATE_TREE
    | splitCsv
}

process GENERATE_TREE{

    cpus cpu

    input:
    tuple val(filter_id),val(filter_directory),val(alignment),val(tree_type)

    output:
    stdout

    script:

    def tree_dir = ("${tree_type}" == "fasttree") ? file("${filter_directory}/FastTree") : file("${filter_directory}/iqTree")
    
    def tree_file = ("${tree_type}" == "fasttree") ? file("${tree_dir}/${filter_id}_FastTree.nwk") : file("${tree_dir}/${filter_id}.treefile")
    def log_file = ("${tree_type}" == "fasttree") ? file("${tree_dir}/out_${filter_id}_FastTree") : file("${tree_dir}/out_${filter_id}_iqTree")

    def bb = params.bb as Integer
    def iqtree_model = (params.gtr) ? "-m GTR+G" : "-m MFP+MERGE"

    def tree_cmd = ("${tree_type}" == "fasttree") ? "fasttree -nt -gtr -out $tree_file $alignment &> $log_file" : "iqtree -nt AUTO $iqtree_model -bb $bb -s $alignment --prefix $filter_id --keep-ident &> $log_file"

    def delete_cmd = (params.overwrite) ? "rm -rf $tree_dir" : ":"

    if(!params.overwrite && tree_dir.isDirectory()){
        error "${tree_dir} exists but --overwrite is not set..."
    }

    """
    $delete_cmd &&
    mkdir $tree_dir &&
    cd $tree_dir &&
    $tree_cmd &&
    echo -n "${filter_id},${tree_file}"
    """
}