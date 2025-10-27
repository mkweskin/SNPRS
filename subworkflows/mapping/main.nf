#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer
sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu

///// Map reads and generate BAMS /////
workflow mapReads{

    take:
    read_data
    genome_info
    mapping_directory

    emit:
    bam_data

    main:

    mapping_reads = genome_info
    .map{it -> tuple(it[0],read_data,mapping_directory)}
    | FETCH_MAP_READS
    | splitCsv

    if(file("${mapping_directory}/ref").isDirectory()){
        bbmap_ref = "${mapping_directory}/ref"
    } else{
        bbmap_ref = genome_info.map{it->tuple(it[2],"${mapping_directory}")}
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
    tuple val(genome_name),val(read_data),val(mapping_directory)

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

workflow mapSRA{
    
    take:
    sra_file
    sra_directory
    genome_info
    mapping_directory

    emit:
    sra_bam_data

    main:

    sra_ids = Channel.fromPath(sra_file).splitText().map{it.trim()}.filter{it}

    sra_reads = sra_ids.combine(Channel.of([sra_directory, params.map_forward, params.map_reverse, params.map_ext])) | STREAM_SRA | splitCsv

    if(file("${mapping_directory}/ref").isDirectory()){
        bbmap_ref = "${mapping_directory}/ref"
    } else{
        bbmap_ref = genome_info.map{it->tuple(it[2],"${mapping_directory}")}
        | BBMAP_INDEX
        | collect
        | map{it->it[0]}
    }

    sra_bam_data = MAP_READS(sra_reads,bbmap_ref,mapping_directory) | splitCsv
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
TOTAL_MEM_MB=\$(free -m | awk '/^Mem:/{print \$2}')
XMX_MB=\$((TOTAL_MEM_MB * 70 / 100))
XMX_ARG="-Xmx\${XMX_MB}m"

bbmap.sh threads=${sample_cpu} in=${forward} in2=${reverse} ambiguous=toss mappedonly=t out=${raw_bam_file} \$XMX_ARG && 
samtools sort -n -@ ${sample_cpu} -o ${sort_file} ${raw_bam_file} && rm -f ${raw_bam_file} && 
samtools fixmate -@ ${sample_cpu} -m ${sort_file} ${mate_file} && rm -f ${sort_file} &&
samtools sort -@ ${sample_cpu} -o ${sort_file} ${mate_file} && rm -f ${mate_file} &&
samtools markdup -@ ${sample_cpu} ${sort_file} ${dup_file} && rm -f ${sort_file} &&
samtools sort -@ ${sample_cpu} -o ${bam_file} ${dup_file} && rm -f ${dup_file} &&
samtools index -@ ${sample_cpu} ${bam_file}"""
        : """
TOTAL_MEM_MB=\$(free -m | awk '/^Mem:/{print \$2}')
XMX_MB=\$((TOTAL_MEM_MB * 70 / 100))
XMX_ARG="-Xmx\${XMX_MB}m"
bbmap.sh threads=${sample_cpu} in=${forward} ambiguous=toss mappedonly=t out=${raw_bam_file} \$XMX_ARG && 
samtools sort -@ ${sample_cpu} -o ${bam_file} ${raw_bam_file} && rm -f ${raw_bam_file} && samtools index -@ ${sample_cpu} ${bam_file}"""
    } else{
        mapping_cmd = reverse ?
    """
TOTAL_MEM_MB=\$(free -m | awk '/^Mem:/{print \$2}')
XMX_MB=\$((TOTAL_MEM_MB * 70 / 100))
XMX_ARG="-Xmx\${XMX_MB}m"
bbmap.sh threads=${sample_cpu} in=${forward} in2=${reverse} ambiguous=toss mappedonly=t out=stdout.bam \$XMX_ARG | \
samtools sort -n -@ ${sample_cpu} -T ${sample_id}_tmp - | \
samtools fixmate -@ ${sample_cpu} -m - - | \
samtools sort -@ ${sample_cpu} -T ${sample_id}_tmp - | \
samtools markdup -@ ${sample_cpu} - - | \
samtools sort -@ ${sample_cpu} -o ${bam_file} - && samtools index -@ ${sample_cpu} ${bam_file}""" :
    """
TOTAL_MEM_MB=\$(free -m | awk '/^Mem:/{print \$2}')
XMX_MB=\$((TOTAL_MEM_MB * 70 / 100))
XMX_ARG="-Xmx\${XMX_MB}m"
bbmap.sh threads=${sample_cpu} in=${forward} ambiguous=toss mappedonly=t out=stdout.bam \$XMX_ARG | \
samtools sort -@ ${sample_cpu} -o ${bam_file} - && samtools index -@ ${sample_cpu} ${bam_file}"""
    }

    """
    cd $mapping_directory &&
    mkdir -p $bam_dir &&
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

process STREAM_SRA {

    tag "Fetch_${srr_id}"

    cpus sample_cpu

    input:
    tuple val(srr_id),val(out_dir),val(forward),val(reverse),val(ext)

    output:
    stdout

    script:

    def out_directory = file("${out_dir}")
    def sra_log_directory = file("${out_directory}/logs")

    def forward_out = file("${out_directory}/${srr_id}${forward}")
    def reverse_out = file("${out_directory}/${srr_id}${reverse}")
    def se_out = file("${out_directory}/${srr_id}${ext}")

    def log_file = file("${sra_log_directory}/out_${srr_id}_Trim")

    def ow_arg = (params.overwrite) ? "overwrite=t": ""

    """
    mkdir -p $out_directory &&
    mkdir -p $sra_log_directory &&

    layout=\$(vdb-dump -R1 -C READ_LEN -f tab $srr_id | awk '{if(NF>1) print "PE"; else print "SE"}') &&

    if [[ "\$layout" == "PE" ]]; then
        fasterq-dump --split-spot --stdout --threads ${sample_cpu} $srr_id | bbduk.sh int=f in=stdin.fq out=${forward_out} out2=${reverse_out} ref=adapters ktrim=r k=23 mink=11 hdist=1 tbo threads=${sample_cpu} $ow_arg &> $log_file &&
        echo -n $srr_id,$forward_out,$reverse_out

    else
        fasterq-dump --stdout --threads ${sample_cpu} $srr_id | bbduk.sh in=stdin.fq out=${se_out} ref=adapters ktrim=r k=23 mink=11 hdist=1 threads=${sample_cpu} $ow_arg &> $log_file &&
        echo -n $srr_id,$se_out,
    fi
    """
}