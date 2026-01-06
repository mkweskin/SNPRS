#! /usr/bin/env nextflow
nextflow.enable.dsl=2

def cpu = params.cpus as Integer
def sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu

def read_ext = params.read_ext
def forward = params.forward
def reverse = params.reverse

def mapping_directory = file("${params.out}/Mapping")
def bbmap_dir = file("${mapping_directory}/ref")
def bam_dir = file("${mapping_directory}/BAMs")

///// Fetch Mapping Reads /////

workflow fetchMapReads{
    
    take:
    map_reads

    emit:
    read_data

    main:

    read_data = FETCH_MAP_READS(map_reads) | splitCsv
}

process FETCH_MAP_READS{

    executor = 'local'
    cpus = 1

    input:
    val(map_reads)

    output:
    stdout

    script:

    def fetchMapScript = file("${projectDir}/bin/fetchMappingReads.py")
    """
    python ${fetchMapScript} -d ${map_reads} -e $read_ext -f $forward -r $reverse
    """
}

///// Map reads and generate BAMS /////

workflow mapReads{

    take:
    mapping_data

    emit:
    bam_data

    main:
  
    reference_fasta = mapping_data
    .first()
    .map{it->it[3].toString()}

    bbmap_ref = (bbmap_dir.isDirectory()) 
    ? Channel.from(bbmap_dir)
    : BBMAP_INDEX(reference_fasta) | collect | map{it->it[0].toString()}
       
    bam_data = mapping_data
    .map{it->tuple(it[0],it[1],it[2])}
    .combine(bbmap_ref)
    | MAP_READS 
    | splitCsv
}


process BBMAP_INDEX{
    
    cpus cpu

    input:
    val(fasta_file)

    output:
    stdout

    script:
  
    """
    TOTAL_MEM_MB=\$(free -m | awk '/^Mem:/{print \$2}')
    XMX_MB=\$((TOTAL_MEM_MB * 70 / 100))
    XMX_ARG="-Xmx\${XMX_MB}m"

    mkdir -p $mapping_directory &&
    mkdir -p $bam_dir &&
    cd $mapping_directory &&
    bbmap.sh threads=${cpu} ref=${fasta_file} \$XMX_ARG &&
    echo -n $bbmap_dir
    """
}

process MAP_READS{

    tag "Map_${sample_id}"

    cpus sample_cpu

    input:
    tuple val(sample_id),val(forward_read),val(reverse_read),val(bbmap_ref)

    output:
    stdout

    script:

    bam_file = file("${bam_dir}/${sample_id}.bam")
    raw_sam_file = file("${bam_dir}/${sample_id}_raw.sam")

    def slow_arg
    if(params.vslow){
        slow_arg = "vslow=t"
    } else if(params.slow){
        slow_arg = "slow=t"
    } else{
        slow_arg = ""
    }

    delete_cmd = (params.overwrite)
    ? "rm -f $bam_file $raw_sam_file" 
    : """
if [ -e "$bam_file" ] || [ -e "$raw_sam_file" ] ; then
    echo "❌ Error: BAM files or intermediates already exist! Use --overwrite to replace." >&2
    exit 1
fi"""    


    mapping_cmd = (reverse_read) ?
    """
TOTAL_MEM_MB=\$(free -m | awk '/^Mem:/{print \$2}')
XMX_MB=\$((TOTAL_MEM_MB * 70 / 100))
XMX_ARG="-Xmx\${XMX_MB}m"
bbwrap.sh $slow_arg threads=${sample_cpu} in=${forward_read},${reverse_read} ambiguous=toss mappedonly=t maxindel=99 strictmaxindel=t append=t out=${raw_sam_file} \$XMX_ARG &&
samtools view -Su -@ ${sample_cpu} -F 4 ${raw_sam_file} | \
samtools sort -@ ${sample_cpu} - -o ${bam_file} && samtools index -@ ${sample_cpu} ${bam_file} && rm -f ${raw_sam_file}""" 
:
    """
TOTAL_MEM_MB=\$(free -m | awk '/^Mem:/{print \$2}')
XMX_MB=\$((TOTAL_MEM_MB * 70 / 100))
XMX_ARG="-Xmx\${XMX_MB}m"
bbmap.sh $slow_arg threads=${sample_cpu} in=${forward_read} ambiguous=toss mappedonly=t maxindel=99 strictmaxindel=t out=stdout.bam \$XMX_ARG | \
samtools sort -@ ${sample_cpu} -o ${bam_file} - && samtools index -@ ${sample_cpu} ${bam_file}"""

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
    """
    python ${fetchBAMScript} -b ${input_bam_files}
    """
}

