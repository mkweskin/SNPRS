#! /usr/bin/env nextflow
nextflow.enable.dsl=2

params.sample_cpus = false
sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : Runtime.runtime.availableProcessors() as Integer

params.overwrite = false

// Get file with SRA accessions
params.sra_file = false

// Defaults
params.f = "_1.fastq.gz"
params.r = "_2.fastq.gz"
params.e = "fastq.gz"


workflow{

    if(!params.sra_file){
        error "Must provide SRR IDs via --sra_file"
    }

    sra_file = file(params.sra_file)
    if(!sra_file.exists()){
        error "${sra_file} does not exist...."
    }

    output_directory = (params.sra_out) ? file(params.sra_out) : file("SRA_Downloads")

    sra_ids = Channel.fromPath(params.sra_file).splitText().map{it.trim()}.filter{it}

    sra_ids.combine(Channel.of([output_directory, params.f, params.r, params.e])) | STREAM_SRA
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
        fasterq-dump --split-spot --stdout --threads ${sample_cpus} $srr_id | bbduk.sh int=f in=stdin.fq out=${forward_out} out2=${reverse_out} ref=adapters ktrim=r k=23 mink=11 hdist=1 tbo threads=${sample_cpus} $ow_arg &> $log_file &&
        echo -n $srr_id,$forward_out,$reverse_out

    else
        fasterq-dump --stdout --threads ${sample_cpus} $srr_id | bbduk.sh in=stdin.fq out=${se_out} ref=adapters ktrim=r k=23 mink=11 hdist=1 threads=${sample_cpus} $ow_arg &> $log_file &&
        echo -n $srr_id,$se_out,
    fi
    """
}
