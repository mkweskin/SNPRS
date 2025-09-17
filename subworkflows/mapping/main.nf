#! /usr/bin/env nextflow
nextflow.enable.dsl=2

def cpu = params.cpus as Integer

///// Map reads and generate BAMS /////
workflow mapReads{

    take:
    read_data
    pangenome_info
    mapping_directory

    emit:
    bam_data

    main:
    
    mapping_info = pangenome_info
    .map{it -> tuple(it[0],it[1],"${read_data}","${mapping_directory}")}


    mapping_reads = FETCH_MAP_READS(mapping_info) 
    | splitCsv()

    bam_data = MAP_READS(mapping_reads) | splitCsv()
}

process FETCH_MAP_READS{

    executor = 'local'
    cpus = 1

    input:
    tuple val(pg_name),val(pg_fasta),val(read_data),val(mapping_directory)

    output:
    stdout

    script:

    def fetchMapScript = file("${projectDir}/bin/fetchMappingReads.py")

    def full_read = file("${read_data}")
    def fasta_file = file("${pg_fasta}")
    def fasta_dir = fasta_file.getParent()
    def bbmap_ref = file("${fasta_dir}/ref")

    def map_dir = file("${mapping_directory}/${pg_name}")
    def existing_ref = file("${map_dir}/ref")
    def bam_dir = file ("${map_dir}/BAM")
    def parquet_dir = file ("${map_dir}/Raw_Parquet")
    def base_call_dir = file ("${map_dir}/Base_Calls")

    """
    mkdir -p $map_dir &&
    mkdir -p $bam_dir &&
    mkdir -p $parquet_dir &&
    mkdir -p $base_call_dir &&
    cd $map_dir

    if [[ ! -d "${existing_ref}" ]]; then
        cp -as ${bbmap_ref} .
    fi

    python ${fetchMapScript} -d ${full_read} -e $params.map_ext -f $params.map_forward -r $params.map_reverse -m $map_dir
    """
}

process MAP_READS{

    cpus cpu

    input:
    tuple val(sample_id),val(forward),val(reverse),val(mapping_directory)

    output:
    stdout

    script:

    def map_dir = file("${mapping_directory}")
    def bam_dir = file ("${map_dir}/BAM")
    def bam_file = file("${bam_dir}/${sample_id}.bam")

    def mapping_cmd = reverse
    ? "bbmap.sh in=${forward} in2=${reverse} ambiguous=toss mappedonly=t out=stdout.sam | samtools view -b - | samtools sort -@ ${params.cpus} -o ${bam_file} -"
    : "bbmap.sh in=${forward} ambiguous=toss mappedonly=t out=stdout.sam | samtools view -b - | samtools sort -@ ${params.cpus} -o ${bam_file} -"

    """
    if [[ -f "${bam_file}" ]]; then
        echo "Error: BAM file exists at ${bam_file}" >&2
        exit 1
    fi

    cd $map_dir &&
    $mapping_cmd &&
    echo -n "${sample_id},${bam_file}"
    """
}

///// Fetch existing BAM files /////

workflow fetchBAM{

    take:
    input_bam_files

    emit:
    bam_files
    
    main:
    bam_files = FETCH_BAM(input_bam_files) | splitCsv()
}

process FETCH_BAM{
    executor = "local"
    cpus 1

    input:
    val(input_bam_files)

    output:
    stdout

    script:

    def fetchBAMScript = file("${projectDir}/bin/fetchBAMs.py")
    full_bam = file("${input_bam_files}")
    """
    python ${fetchBAMScript} -b ${full_bam}
    """
}

