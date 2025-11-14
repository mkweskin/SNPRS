#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = (!params.iqtree) ? 1 : params.cpus as Integer

workflow generateTree{
    
    take:
    input_data

    emit:
    tree_file

    main:

    tree_file = input_data | GENERATE_TREE | splitCsv | collect | collate(1) | checkStop
}

workflow checkStop{
    take:
    pre_tree

    emit:
    tree_file

    main:

    tree_file = (params.tree) ? Channel.empty() : pre_tree

}


process GENERATE_TREE{

    cpus cpu

    input:
    tuple val(input_id),val(input_dir),val(align_file)

    output:
    stdout

    script:

    def alignment_file = file(align_file)
    def tree_type = (params.iqtree) ? "iqtree" : "fasttree"

    def tree_dir = ("${tree_type}" == "fasttree") ? file("${input_dir}/FastTree") : file("${input_dir}/iqTree")
    def tree_file = ("${tree_type}" == "fasttree") ? file("${tree_dir}/${input_id}_FastTree.nwk") : file("${tree_dir}/${input_id}.treefile")
    def log_file = ("${tree_type}" == "fasttree") ? file("${tree_dir}/out_${input_id}_FastTree") : file("${tree_dir}/out_${input_id}_iqTree")

    def bb = params.bb as Integer
    def iqtree_model = (params.gtr) ? "-m GTR+G" : "-m MFP+MERGE"

    def tree_cmd = ("${tree_type}" == "fasttree") ? "VeryFastTree -threads $cpu -nt -gtr gamma -out $tree_file $alignment_file &> $log_file" : "iqtree -nt AUTO $iqtree_model -bb $bb -s $alignment_file --prefix $input_id --keep-ident &> $log_file"

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

workflow fetchTree{

    take:
    tree_path

    emit:
    tree_file

    main:
    tree_file = FETCH_TREE(tree_path)
}

process FETCH_TREE{

    cpus 1
    executor = "local"

    input:
    val(tree_path)

    output:
    stdout

    script:

    def fetch_tree_script = file("${projectDir}/bin/helper_scripts/fetchTree.py")

    """
    python $fetch_tree_script $tree_path &&
    echo -n $tree_path
    """
}

workflow makeSplitTable{

    take:
    tree_file

    emit:
    tree_split_file

    main:
    tree_split_file = tree_file | map{it->it[0]} | MAKE_SPLIT_TABLE
}

process MAKE_SPLIT_TABLE{

    input:
    val(tree_file)
    
    output:
    stdout

    script:

    def get_split_script = file("${projectDir}/bin/tree2splits.py")

    def tree_path = file(tree_file)
    def tree_dir = tree_path.getParent()
    def tree_name = tree_path.getBaseName()
    def output_file = file("${tree_dir}/${tree_name}_TreeSplits.csv")

    def delete_cmd = (params.overwrite) 
    ? "rm -f $output_file"
    : """
if [ -e "$output_file" ]; then
    echo "❌ Error: ${output_file} already exists! Use --overwrite to replace." >&2
    exit 1
fi"""

    """
    $delete_cmd &&
    python $get_split_script $tree_path $output_file &&
    echo -n $output_file
    """
}

workflow makeSNPGroups{
    
    take:
    tree_split
    
    emit:
    group_file

    main:

    group_file = MAKE_SNP_GROUPS(tree_split)
}

process MAKE_SNP_GROUPS{

    input:
    tuple val(tree_file),val(split_file)

    output:
    stdout
    
    script:

    def get_group_script = file("${projectDir}/bin/splits2groups.py")

    def tree_path = file(tree_file)
    def split_path = file(split_file)

    def split_dir = split_path.getParent()
    def split_name = split_path.getBaseName().replaceFirst(/_TreeSplits$/, '')
  
    def output_file = file("${split_dir}/${split_name}_Monophyletic_Groups.csv")

    def mono_arg = (params.mono) ? "--mono" : ""

    def delete_cmd = (params.overwrite) 
    ? "rm -f $output_file"
    : """
if [ -e "$output_file" ]; then
    echo "❌ Error: ${output_file} already exists! Use --overwrite to replace." >&2
    exit 1
fi"""

    """
    $delete_cmd &&
    python $get_group_script --tree $tree_path --splits $split_path --out $output_file $mono_arg &&
    echo -n $output_file
    """
}
