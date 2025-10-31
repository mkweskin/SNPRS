#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer
sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu

classified_directory = file(params.final_classified_directory)

workflow classifySample{

    take:
    classifier_data

    emit:
    classified_data

    main:

    classified_data = CLASSIFY_SAMPLE(classifier_data) | splitCsv
}

process CLASSIFY_SAMPLE{

    cpus sample_cpu

    input:
    tuple val(sample_id),val(called_bases_file),val(snp_directory),val(snp_id)

    output:
    stdout

    script:

    def snp_path = file("${snp_directory}/${snp_id}")
    def classified_dir = file("${classified_directory}/${snp_id}")
    def output_parquet = file("${classified_dir}/${sample_id}_Classified.parquet")
    def output_csv = file("${classified_dir}/${sample_id}_Classified.csv")

    def delete_cmd = (params.overwrite) ? "rm -rf $output_parquet $output_csv"
    : """
if [ -e "$output_parquet" ] || [ -e "$output_csv" ] ; then
    echo "❌ Error: Classification files for ${sample_id} already exist! Use --overwrite to replace." >&2
    exit 1
fi"""

    """
    mkdir -p $classified_dir &&
    $delete_cmd &&
    # classify &&
    echo -n "${sample_id},${output_parquet},${output_csv}"
    """
    

}