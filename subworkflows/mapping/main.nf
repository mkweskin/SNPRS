#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer
sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu

///// Map reads and generate BAMS /////
workflow mapReads{

    take:
    read_data
    pangenome_info
    mapping_directory

    emit:
    bam_data

    main:

    mapping_reads = FETCH_MAP_READS("${read_data}")
    | splitCsv

    if(file("${mapping_directory}/ref").isDirectory()){
        bbmap_ref = "${mapping_directory}/ref"
    } else{
        bbmap_ref = pangenome_info.map{it->tuple(it[2],"${mapping_directory}")}
        | BBMAP_INDEX
        | collect
        | map{it->it[0]}
    }

    bam_data = MAP_READS(mapping_reads,bbmap_ref,mapping_directory) 
    | splitCsv
}

process FETCH_MAP_READS{

    executor = 'local'
    cpus = 1

    input:
    val(read_data)

    output:
    stdout

    script:

    def fetchMapScript = file("${projectDir}/bin/fetchMappingReads.py")
 
    def ext
    def forward
    def reverse

    if(params.validate){
        ext = "${params.pg_ext}"
        forward = "${params.pg_forward}"
        reverse = "${params.pg_reverse}"
    } else{
        ext = "${params.map_ext}"
        forward = "${params.map_forward}"
        reverse = "${params.map_reverse}"
    }

    """
    python ${fetchMapScript} -d ${read_data} -e $ext -f $forward -r $reverse
    """
}

process BBMAP_INDEX{
    
    cpus cpu

    input:
    tuple val(genome_file),val(mapping_directory)

    output:
    stdout

    script:

    def ref_directory = file("${mapping_directory}/ref")

    def bam_dir = file ("${mapping_directory}/BAMs")
    def parquet_dir = file ("${mapping_directory}/Raw_Parquet")
    def base_call_dir = file ("${mapping_directory}/Base_Calls")

    def fasta_file = file("${genome_file}")
    
    """
    TOTAL_MEM_MB=\$(free -m | awk '/^Mem:/{print \$2}')
    XMX_MB=\$((TOTAL_MEM_MB * 70 / 100))
    XMX_ARG="-Xmx\${XMX_MB}m"

    cd ${mapping_directory} &&
    mkdir -p $bam_dir &&
    mkdir -p $parquet_dir &&
    mkdir -p $base_call_dir &&
    bbmap.sh threads=${cpu} ref=${fasta_file} \$XMX_ARG &&
    echo -n $ref_directory
    """
}

process MAP_READS{

    tag "Map_${sample_id}"

    cpus sample_cpu

    input:
    tuple val(sample_id),val(forward),val(reverse)
    val(bbmap_ref)
    val(mapping_directory)

    output:
    stdout

    script:

    def bam_dir = file ("${mapping_directory}/BAMs")
    def bam_file = file("${bam_dir}/${sample_id}.bam")

    def raw_bam_file = file("${bam_dir}/${sample_id}_raw.bam")
    def sort_file = file("${bam_dir}/${sample_id}_sort.bam")
    def mate_file = file("${bam_dir}/${sample_id}_mate.bam")
    def dup_file = file("${bam_dir}/${sample_id}_dup.bam")

    def delete_cmd = (params.overwrite)
    ? "rm -f $bam_file $raw_bam_file $sort_file $mate_file $dup_file" 
    : """
if [ -e "$bam_file" ] || [ -e "$raw_bam_file" ] || [ -e "$sort_file" ] || [ -e "$mate_file" ] || [ -e "$dup_file" ]; then
    echo "❌ Error: BAM files or intermediates already exist! Use --overwrite to replace." >&2
    exit 1
fi"""    

    def mapping_cmd

    if(params.mem_mode){
        mapping_cmd = reverse
        ? """
bbmap.sh threads=${sample_cpu} in=${forward} in2=${reverse} ambiguous=toss mappedonly=t out=${raw_bam_file} && 
samtools sort -n -@ ${sample_cpu} -o ${sort_file} ${raw_bam_file} && rm -f ${raw_bam_file} && 
samtools fixmate -@ ${sample_cpu} -m ${sort_file} ${mate_file} && rm -f ${sort_file} &&
samtools sort -@ ${sample_cpu} -o ${sort_file} ${mate_file} && rm -f ${mate_file} &&
samtools markdup -@ ${sample_cpu} ${sort_file} ${dup_file} && rm -f ${sort_file} &&
samtools sort -@ ${sample_cpu} -o ${bam_file} ${dup_file} && rm -f ${dup_file} &&
samtools index -@ ${sample_cpu} ${bam_file}"""
        : """
bbmap.sh threads=${sample_cpu} in=${forward} ambiguous=toss mappedonly=t out=${raw_bam_file} && 
samtools sort -@ ${sample_cpu} -o ${bam_file} ${raw_bam_file} && rm -f ${raw_bam_file} && samtools index -@ ${sample_cpu} ${bam_file}"""
    } else{
        mapping_cmd = reverse ?
    """
bbmap.sh threads=${sample_cpu} in=${forward} in2=${reverse} ambiguous=toss mappedonly=t out=stdout.bam | \
samtools sort -n -@ ${sample_cpu} -T ${sample_id}_tmp - | \
samtools fixmate -@ ${sample_cpu} -m - - | \
samtools sort -@ ${sample_cpu} -T ${sample_id}_tmp - | \
samtools markdup -@ ${sample_cpu} - - | \
samtools sort -@ ${sample_cpu} -o ${bam_file} - && samtools index -@ ${sample_cpu} ${bam_file}""" :
    """
bbmap.sh threads=${sample_cpu} in=${forward} ambiguous=toss mappedonly=t out=stdout.bam| \
samtools sort -@ ${sample_cpu} -o ${bam_file} - && samtools index -@ ${sample_cpu} ${bam_file}"""
    }

    """
    cd $mapping_directory &&
    $delete_cmd &&
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

