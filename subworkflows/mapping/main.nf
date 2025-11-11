#! /usr/bin/env nextflow
nextflow.enable.dsl=2

cpu = params.cpus as Integer
sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu

genome_directory = file(params.final_genome_directory)
mapping_directory = file(params.final_mapping_directory)
sra_directory = file(params.final_sra_directory)
bbmap_ref = file("${mapping_directory}/ref")

read_ext = params.read_ext
forward = params.forward
reverse = params.reverse

///// Fetch Mapping Reads /////

workflow fetchMapReads{
    
    take:
    read_dir

    emit:
    read_data

    main:

    read_data = FETCH_MAP_READS(read_dir) | splitCsv
}

process FETCH_MAP_READS{

    executor = 'local'
    cpus = 1

    input:
    val(read_dir)

    output:
    stdout

    script:

    def fetchMapScript = file("${projectDir}/bin/fetchMappingReads.py")
    full_read = file("${read_dir}")
    """
    python ${fetchMapScript} -d ${full_read} -e $read_ext -f $forward -r $reverse
    """
}

///// Download SRA Data /////
workflow fetchSRAReads{
    
    take:
    sra_file

    emit:
    read_data

    main:

    sra_ids = Channel.fromPath(sra_file).splitText().map{it.trim()}.filter{it}
    read_data = sra_ids | STREAM_SRA | splitCsv
}

process STREAM_SRA{

    tag "Fetch_${srr_id}"

    cpus sample_cpu

    input:
    val(srr_id)

    output:
    stdout

    script:

    def sra_log_directory = file("${sra_directory}/logs")

    def safe_ext = read_ext.startsWith('.') ? read_ext : ".${read_ext}"
    def forward_out = file("${sra_directory}/${srr_id}${forward}")
    def reverse_out = file("${sra_directory}/${srr_id}${reverse}")
    def se_out = file("${sra_directory}/${srr_id}${safe_ext}")

    def log_file = file("${sra_log_directory}/out_${srr_id}_Trim")

    def ow_arg = (params.overwrite) ? "overwrite=t" : ""
    def delete_cmd = (params.overwrite) ? "rm -f $forward_out $reverse_out $se_out" 
    : """
if [ -e "$forward_out" ] || [ -e "$reverse_out" ] || [ -e "$se_out" ] ; then
    echo "❌ Error: SRA read files already exist! Use --overwrite to replace." >&2
    exit 1
fi"""    

    """
    mkdir -p $sra_directory &&
    mkdir -p $sra_log_directory &&
    rm -rf $log_file &&
    $delete_cmd &&

    layout=\$(vdb-dump -R1 -C READ_LEN -f tab $srr_id | awk '{if(NF>1) print "PE"; else print "SE"}') &&
    
    {
        if [[ "\$layout" == "PE" ]]; then
            fasterq-dump --split-spot --stdout --threads ${sample_cpu} $srr_id | bbduk.sh int=f in=stdin.fq out=${forward_out} out2=${reverse_out} ref=adapters ktrim=r k=23 mink=11 hdist=1 tbo threads=${sample_cpu} $ow_arg &> $log_file &&
            echo -n $srr_id,$forward_out,$reverse_out

        else
            fasterq-dump --stdout --threads ${sample_cpu} $srr_id | bbduk.sh in=stdin.fq out=${se_out} ref=adapters ktrim=r k=23 mink=11 hdist=1 threads=${sample_cpu} $ow_arg &> $log_file &&
            echo -n $srr_id,$se_out,
        fi
    } || {
        echo -n ""
    }

    """
}

///// Map reads and generate BAMS /////
workflow mapReads{

    take:
    mapping_data

    emit:
    bam_data

    main:

    reference_fasta = mapping_data.first().map{it->it[3]}
    bbmap_ref = (bbmap_ref.isDirectory()) ? bbmap_ref : BBMAP_INDEX(reference_fasta) | collect | map{it->it[0]}
    
    read_data = mapping_data.map{it->tuple(it[0],it[1],it[2])}
    
    bam_data = MAP_READS(read_data,bbmap_ref) | splitCsv
}


process BBMAP_INDEX{
    
    cpus cpu

    input:
    val(genome_file)

    output:
    stdout

    script:

    def ref_directory = file("${mapping_directory}/ref")
    def fasta_file = file("${genome_file}")
    
    """
    TOTAL_MEM_MB=\$(free -m | awk '/^Mem:/{print \$2}')
    XMX_MB=\$((TOTAL_MEM_MB * 70 / 100))
    XMX_ARG="-Xmx\${XMX_MB}m"

    mkdir -p $mapping_directory &&
    cd $mapping_directory &&
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

    output:
    stdout

    script:

    def bam_dir = file ("${mapping_directory}/BAMs")
    def bam_file = file("${bam_dir}/${sample_id}.bam")

    def raw_sam_file = file("${bam_dir}/${sample_id}_raw.sam")

    def slow_arg
    if(params.vslow){
        slow_arg = "vslow=t"
    } else if(params.slow){
        slow_arg = "slow=t"
    } else{
        slow_arg = ""
    }

    def delete_cmd = (params.overwrite)
    ? "rm -f $bam_file $raw_bam_file $sort_file $mate_file $dup_file" 
    : """
if [ -e "$bam_file" ] || [ -e "$raw_sam_file" ] ; then
    echo "❌ Error: BAM files or intermediates already exist! Use --overwrite to replace." >&2
    exit 1
fi"""    


    def mapping_cmd = (reverse) ?
    """
TOTAL_MEM_MB=\$(free -m | awk '/^Mem:/{print \$2}')
XMX_MB=\$((TOTAL_MEM_MB * 70 / 100))
XMX_ARG="-Xmx\${XMX_MB}m"
bbwrap.sh $slow_arg threads=${sample_cpu} in=${forward},${reverse} ambiguous=toss mappedonly=t maxindel=99 strictmaxindel=t append=t out=${raw_sam_file} \$XMX_ARG &&
samtools view -Su -@ ${sample_cpu} -F 4 ${raw_sam_file} | \
samtools sort -@ ${sample_cpu} - -o ${bam_file} && samtools index -@ ${sample_cpu} ${bam_file} && rm -f ${raw_sam_file}""" 
:
    """
TOTAL_MEM_MB=\$(free -m | awk '/^Mem:/{print \$2}')
XMX_MB=\$((TOTAL_MEM_MB * 70 / 100))
XMX_ARG="-Xmx\${XMX_MB}m"
bbmap.sh $slow_arg threads=${sample_cpu} in=${forward} ambiguous=toss mappedonly=t maxindel=99 strictmaxindel=t out=stdout.bam \$XMX_ARG | \
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
    full_bam = file("${input_bam_files}")
    """
    python ${fetchBAMScript} -b ${full_bam}
    """
}

