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

    if(params.snp_dir){
        error "Can't run GENERATE_SNPS if providing a --snp_dir"
    }

    def base_parquet = file("${data_dir}/${data_id}_Bases.parquet")
    def snp_directory = file("${snp_dir}/${snp_id}")

    def output_json = file("${snp_directory}/${snp_id}.json")
    def output_comparisons = file("${snp_directory}/${snp_id}_Comparisons.csv")
    def row_numbers = file("${snp_directory}/${snp_id}_Row_Numbers.txt")
    def snp_parquet = file("${snp_directory}/${snp_id}_SNPs.parquet")
    
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
    python $generate_snp_script --out $snp_directory --snp_id $snp_id --bases $base_parquet --groups $group_file --tree $tree_file &&
    echo -n "${snp_id},${snp_dir}"
    """
}
