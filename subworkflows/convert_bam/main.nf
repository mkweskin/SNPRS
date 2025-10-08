#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer
sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu

///// Convert BAM files to Parquet /////
workflow bamToParquet{

    take:
    bam_data
    pangenome_info
    mapping_directory
    
    emit:
    raw_parquet_data

    main:
    raw_parquet_data = BAM_TO_PARQUET(bam_data,pangenome_info,mapping_directory) | splitCsv()
}

process BAM_TO_PARQUET{
    
    cpus sample_cpu

    tag "BAM2PQ_${sample_id}"

    input:
    tuple val(sample_id),val(sample_bam)
    tuple val(pg_name),val(fasta_path)
    val(mapping_directory)

    output:
    stdout

    script:

    def bam_convert_script = file("${projectDir}/bin/bam2parquet.py")

    def output_directory = file("${mapping_directory}/Raw_Parquet")
    def pileup_directory = file("${mapping_directory}/Pileups")
    def output_file = file("${output_directory}/${sample_id}_Raw.parquet")
    
    def pileup_file = file("${pileup_directory}/${sample_id}.pileup")
    def filter_bam = file("${output_directory}/${sample_id}_filter.bam")

    def mapq = params.mapq as Integer
    def baseq = params.baseq as Integer

    def delete_cmd = (params.overwrite)
    ? "rm -f $output_file $pileup_file"
    : """
if [ -e "$output_file" ] || [ -e "$pileup_file" ]; then
    echo "❌ Error: Output or pileup file already exists — use --overwrite to replace." >&2
    exit 1
fi"""    

    def pileup_cmd ="""samtools view -@ ${sample_cpu} -q ${mapq} -h -F 3844 -o ${filter_bam} ${sample_bam}; samtools index -@ ${sample_cpu} ${filter_bam}; samtools mpileup -Q ${baseq} --no-output-ends --no-output-del -f ${fasta_path} ${filter_bam} | awk 'NF >= 4 && (\$4 > 0 || \$4 != "")' > ${pileup_file};rm -f ${filter_bam}"""

    """
    mkdir -p ${output_directory} &&
    mkdir -p ${pileup_directory} &&
    $delete_cmd &&
    $pileup_cmd &&
    python ${bam_convert_script} --bam ${sample_bam} --pileup ${pileup_file} --fasta ${fasta_path} --parquet ${output_file} --mapq ${mapq} --baseq ${baseq} &&
    echo -n "${sample_id},${output_file}"
    """
}

///// Fetch Raw Parquets /////
workflow fetchRawParquet{
    take:
    input_raw_parquet

    emit:
    raw_parquet_data

    main:
    raw_parquet_data = FETCH_RAW_PARQUET(input_raw_parquet) | splitCsv()
}

process FETCH_RAW_PARQUET{
    executor = "local"
    cpus = 1

    input:
    val(input_raw_parquet)

    output:
    stdout

    script:

    def fetch_raw_parquet_script = file("${projectDir}/bin/fetchRawParquet.py")
    def full_parquet = file("${input_raw_parquet}")
    """
    python ${fetch_raw_parquet_script} -p ${full_parquet}
    """
}
