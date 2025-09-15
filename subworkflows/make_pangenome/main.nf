#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// Only check args if reads are provided

if("${params.pg_reads}" != ""){
    
    read_directory = file("${params.pg_reads}")
    
    if(!read_directory.isDirectory()){
        error "Pangenome read directory (--pg_reads) ${params.pg_reads} does not exist..."
    }

    if(params.pg_name == ""){
        pg_name = "SNPRS_${params.timestamp}"
        params.pg_name = "${pg_name}"
    } else{
        pg_name = "${params.pg_name}"
    }

    if(params.pg_out == ""){
        output_directory = file("${pg_name}")
        params.pg_out = "${output_directory}"
    }
    else{
        output_directory = file("${params.pg_out}")
    }

    if(!output_directory.getParent().isDirectory()){
        error "Parent directory for output is not a valid directory [${output_directory.getParent()}]..."
    }

    if("${params.size}" == ""){
        error "Must provide approximate genome size in bp via --size"
    } else{
        size = params.size as Integer
    }

    coverage = params.coverage as Integer
    out_prop = params.out_prop as Float
    kmer = params.kmer as Integer
    cpu = params.cpus as Integer
    node = params.nodes as Integer

    ray_cores = cpu * node
}

workflow makePangenome{

    take:
    output_directory
    pg_name
    read_directory

    emit:
    pangenome_info

    main:

    input_pangenome_reads = FETCH_PG_READS(output_directory,pg_name,read_directory) 
    | splitCsv()


    base_count_file = input_pangenome_reads
    .map{it->tuple(it[0],it[3],it[4],"${output_directory}","${pg_name}")}
    | COUNT_BASES
    | collect
    | flatten
    | first

    subset_guide = CALCULATE_SUBSETS(base_count_file,output_directory,pg_name) |
    splitCsv() |
    branch{it ->
        link: it[5].toString() == "Link"
            return(tuple(it[0],it[1],it[2],it[3]))
        sample: true
            return(tuple(it[0],it[1],it[2],it[3],it[4]))
    }

    subset_folder_1 = SUBSET_READS(subset_guide.sample,output_directory,pg_name)
    subset_folder_2 = LINK_READS(subset_guide.link,output_directory,pg_name)

    subset_folder = subset_folder_1
    .concat(subset_folder_2)
    | collect
    | flatten
    | first
    
    ray_assembly = ASSEMBLE_PANGENOME(subset_folder,output_directory,pg_name)

    pangenome_info = PROCESS_RAY(ray_assembly,output_directory,pg_name) 
    | splitCsv
    | collect
    | flatten
    | collate(2)
}

process FETCH_PG_READS{

    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(output_directory)
    val(pg_name)
    val(read_directory)

    output:
    stdout

    script:

    def fetchPGScript = file("${projectDir}/bin/fetchPangenomeReads.py")
    
    def pangenome_directory = file("${output_directory}/${pg_name}")
    def pangenome_prep_directory = file("${pangenome_directory}/Prep_${pg_name}")
    def pangenome_subset_directory = file("${pangenome_prep_directory}/Subset_Reads")

    def group_file = file("${pangenome_prep_directory}/Read_Groups.csv")

    def full_out = file("${output_directory}")
    def full_read = file("${read_directory}")
    """
    mkdir -p $full_out &&
    mkdir $pangenome_directory &&
    mkdir $pangenome_prep_directory &&
    mkdir $pangenome_subset_directory &&
    python ${fetchPGScript} -d ${full_read} -e $params.pg_ext -f $params.pg_forward -r $params.pg_reverse -o $group_file 
    """
}

process COUNT_BASES {

    label 'basicTools'
    tag "Count_${sample_id}"

    cpus = 1
    memory = '4 GB'

    input:
    tuple val(sample_id), val(forward), val(reverse), val(output_directory), val(pg_name)

    output:
    stdout

    script:
    
    def pangenome_directory = file("${output_directory}/${pg_name}")
    def pangenome_prep_directory = file("${pangenome_directory}/Prep_${pg_name}")
    def base_count_file = file("${pangenome_prep_directory}/Read_Counts.csv")

    def reformat_cmd = reverse 
        ? "reformat.sh in=${forward} in2=${reverse} 2>&1" 
        : "reformat.sh in=${forward} 2>&1"

    """
    output=\$(${reformat_cmd})
    read_count_line=\$(echo "\$output" | grep "Output:")
    read_count=\$(echo "\$read_count_line" | awk '{print \$2}')
    base_count=\$(echo "\$read_count_line" | awk '{print \$5}')
    echo -e "${sample_id},\$read_count,\$base_count,${forward},${reverse}" >> ${base_count_file}
    echo -n ${base_count_file}
    """
}

process CALCULATE_SUBSETS{

    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(base_count_file)
    val(output_directory)
    val(pg_name)

    output:
    stdout

    script:

    def calculate_sub_script = file("${projectDir}/bin/calculateSubset.py")
    def pangenome_directory = file("${output_directory}/${pg_name}")
    def pangenome_prep_directory = file("${pangenome_directory}/Prep_${pg_name}")
    def pangenome_subset_directory = file("${pangenome_prep_directory}/Subset_Reads")
    def group_file = file("${pangenome_prep_directory}/Read_Groups.csv")
    
    """
    python ${calculate_sub_script} -b ${base_count_file} -g ${group_file} -s ${size} -c ${coverage} -o ${pangenome_subset_directory} -p ${out_prop}
    """
}

process SUBSET_READS {

    label 'basicTools'
    tag "Subset_${sample_id}"

    cpus 1
    executor="slurm"
    memory '4 GB'

    input:
    tuple val(sample_id), val(subsample_id),val(forward),val(reverse),val(allocated)
    val(output_directory)
    val(pg_name)

    output:
    stdout

    script:

    def pangenome_directory = file("${output_directory}/${pg_name}")
    def pangenome_prep_directory = file("${pangenome_directory}/Prep_${pg_name}")
    def pangenome_subset_directory = file("${pangenome_prep_directory}/Subset_Reads")

    def out1 = "${pangenome_subset_directory}/${subsample_id}_GenomeReads_1.fq.gz"
    def out2 = "${pangenome_subset_directory}/${subsample_id}_GenomeReads_2.fq.gz"
    def outs = "${pangenome_subset_directory}/${subsample_id}_GenomeReads_SE.fq.gz"

    def log_file = "${pangenome_subset_directory}/out_Subsample_${sample_id}"

    def reformat_cmd = reverse
        ? "reformat.sh in=${forward} in2=${reverse} out=${out1} out2=${out2} outs=${outs} samplebasestarget=${allocated} &> ${log_file}"
        : "reformat.sh in=${forward} out=${outs} samplebasestarget=${allocated} &> ${log_file}"

    """
    ${reformat_cmd}

    has_content () {
        zcat "\$1" 2>/dev/null | head -c 1 | grep -q .
    }
    for f in ${out1} ${out2} ${outs}; do
        if [[ -f "\$f" ]] && ! has_content "\$f"; then
            rm -f "\$f"
        fi
    done

    echo -n "${pangenome_subset_directory}"
    """
}

process LINK_READS {
    cpus 1
    executor = "local"

    input:
    tuple val(sample_id), val(subsample_id),val(forward),val(reverse)
    val(output_directory)
    val(pg_name)

    output:
    stdout

    script:   

    def pangenome_directory = file("${output_directory}/${pg_name}")
    def pangenome_prep_directory = file("${pangenome_directory}/Prep_${pg_name}")
    def pangenome_subset_directory = file("${pangenome_prep_directory}/Subset_Reads")

    def out1 = forward.endsWith(".gz")
    ? "${pangenome_subset_directory}/${subsample_id}_GenomeReads_1.fq.gz"
    : "${pangenome_subset_directory}/${subsample_id}_GenomeReads_1.fastq"   
        
    def out2 = forward.endsWith(".gz")
    ? "${pangenome_subset_directory}/${subsample_id}_GenomeReads_2.fq.gz"
    : "${pangenome_subset_directory}/${subsample_id}_GenomeReads_2.fastq"

    def outs = forward.endsWith(".gz")
    ? "${pangenome_subset_directory}/${subsample_id}_GenomeReads_SE.fq.gz"
    : "${pangenome_subset_directory}/${subsample_id}_GenomeReads_SE.fastq"

    def link_cmd = reverse
        ? "ln -s ${forward} ${out1}; ln -s ${reverse} ${out2}"
        : "ln -s ${forward} ${outs}"

    """
    ${link_cmd}
    echo -n "${pangenome_subset_directory}"
    """
}

process ASSEMBLE_PANGENOME {

    label 'assemblePangenome'
    tag "Assemble_Pangenome"

    executor = 'slurm'
    clusterOptions = "--nodes=${node} --ntasks-per-node=${cpu}"

    input:
    val(subset_folder)
    val(output_directory)
    val(pg_name)

    output:
    stdout

    script:

    def pangenome_directory = file("${output_directory}/${pg_name}")
    def pangenome_prep_directory = file("${pangenome_directory}/Prep_${pg_name}")
    
    def assembly_directory = file("${pangenome_prep_directory}/Ray_${pg_name}")
    def ray_log = file("${pangenome_prep_directory}/out_Ray_${pg_name}")
    def load_ray_module = params.ray_module == "" ? "" : "module load -s ${params.ray_module}"

    """
    $load_ray_module
    mpirun -np ${ray_cores} Ray -k ${kmer} -detect-sequence-files ${subset_folder} -o ${assembly_directory}  &> ${ray_log} &&
    echo -n "${assembly_directory}/Contigs.fasta"
    """
}

process PROCESS_RAY{
    cpus 1
    executor = "local"

    input:
    val(ray_assembly)
    val(output_directory)
    val(pg_name)

    output:
    stdout

    script:
    def pangenome_directory = file("${output_directory}/${pg_name}")

    def pangenome_file = file("${pangenome_directory}/${pg_name}.fasta")
    def stats_file = file("${pangenome_directory}/${pg_name}_BBStats")
    """
    rename.sh in=${ray_assembly} out=${pangenome_file} prefix=SNPRS addprefix=t trd=t
    stats.sh ${pangenome_file} &> ${stats_file}
    echo -n "${pg_name},${pangenome_file}"
    """
}

workflow{
    if (params.pg_out && params.pg_name && params.pg_reads) {
        pangenome_data = makePangenome(params.pg_out, params.pg_name,params.pg_reads)
    }
}
