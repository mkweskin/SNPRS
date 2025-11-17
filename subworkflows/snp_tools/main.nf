#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer

workflow generateSNPs{
    take:
    snp_group_info

    emit:
    snp_data

    main:

    snp_data = GENERATE_SNPS(snp_group_info)

}

process GENERATE_SNPS{

    input:
    tuple val(data_id),val(data_dir),val(snp_id),val(snp_dir),val(tree_file),val(group_file)

    output:
    stdout

    script:

    def generate_snp_script = file("${projectDir}/bin/generateSNPs.py")
    def tree_script = file("${projectDir}/bin/helper_scripts/groups2tree.py")

    if(params.snp_dir){
        error "Can't run GENERATE_SNPS if providing a --snp_dir"
    }

    def scaffold_parquet = file("${data_dir}/${data_id}_Scaffold.parquet")
    def base_parquet = file("${data_dir}/${data_id}_Bases.parquet")
    def snp_directory = file("${snp_dir}/${snp_id}")
    def group_tree = file("${snp_directory}/${snp_id}_SNP_Groups.nwk")

    def delete_cmd = (params.overwrite)
    ? "rm -rf $snp_directory" 
    : """
if [ -d "$snp_directory" ] ; then
    echo "❌ Error: ${snp_directory} already exists! Use --overwrite to replace." >&2
    exit 1
fi""" 

    """
    $delete_cmd &&
    mkdir -p $snp_directory &&
    python $tree_script $group_file $group_tree &&
    python $generate_snp_script --out $snp_directory --snp_id $snp_id --bases $base_parquet --scaffold $scaffold_parquet --groups $group_file --tree $tree_file &&
    echo -n "${snp_id},${snp_dir}"
    """
}

/*
workflow getFixedSites{

    take:
    called_base_file
    fixed_id

    emit:
    fixed_data

    main:

    fixed_data = GET_FIXED_SITES(called_base_file,fixed_id) | splitCsv
}

process GET_FIXED_SITES{

    input:
    val(called_base_file)
    val(fixed_id)

    output:
    stdout

    script:

    def get_fixed_script = file("${projectDir}/bin/getFixedSites.py")
    
    def called_base_path = file(called_base_file)
    
    def output_directory = file("${params.final_fixed_directory)}/${fixed_id}")
        
    def missing_arg = (params.missing) ? "--missing ${params.missing}" : ""
    def gap_arg = (params.no_gaps) ? "--no_gaps" : ""

    def delete_cmd = (params.overwrite) ? "rm -rf $output_directory" : 
    """
if [[ -d "${output_directory}" ]]; then
    echo "Error: Directory ${output_directory} already exists. Use --overwrite to replace it." >&2
    exit 1
fi"""  

    """
    $delete_cmd &&
    mkdir -p ${params.final_fixed_directory} &&
    mkdir $output_directory &&
    python $get_fixed_script --called_bases $called_base_path --fixed_id $fixed_id --out $output_directory $missing_arg $gap_arg &&
    echo -n "${fixed_id},${output_directory}"
    """
}
*/