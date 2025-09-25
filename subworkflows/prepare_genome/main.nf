#! /usr/bin/env nextflow
nextflow.enable.dsl=2

coverage = params.coverage as Integer
out_prop = params.out_prop as Float
kmer = params.kmer as Integer

cpu = params.cpus as Integer
node = params.nodes as Integer
ray_cpu = (params.ray_cpus) ? params.ray_cpus as Integer : cpu
ray_cores = ray_cpu * node

///// Create a SNPRS pangenome from reads /////
workflow assembleGenome{

    take:
    pangnome_directory
    pg_name
    pg_read_data

    emit:
    return_pangenome

    main:

    if(!params.size){
        error "Cannot assemble pangenome without --size estimate (genome size in basepairs)"
    }

    input_pangenome_reads = FETCH_PG_READS(pangnome_directory,pg_name,pg_read_data) | splitCsv

    base_count_file = input_pangenome_reads
    .map{it->tuple(it[0],it[3],it[4],"${pangnome_directory}","${pg_name}")}
    | COUNT_BASES
    | collect
    | flatten
    | first

    subset_guide = CALCULATE_SUBSETS(base_count_file,pangnome_directory,pg_name) | splitCsv
    | branch{it ->
        link: it[5].toString() == "Link"
            return(tuple(it[0],it[1],it[2],it[3]))
        sample: true
            return(tuple(it[0],it[1],it[2],it[3],it[4]))
    }

    subset_folder_1 = SUBSET_READS(subset_guide.sample,pangnome_directory,pg_name)
    subset_folder_2 = LINK_READS(subset_guide.link,pangnome_directory,pg_name)

    subset_folder = subset_folder_1
    .concat(subset_folder_2)
    | collect
    | flatten
    | first
    
    ray_assembly = ASSEMBLE_PANGENOME(subset_folder,pangnome_directory,pg_name)

    return_pangenome = PROCESS_RAY(ray_assembly,pangnome_directory,pg_name) | splitCsv | collect | flatten | collate(2)
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
    
    if(!params.overwrite && pangenome_directory.isDirectory()){
        error "${pangenome_directory} exists and --overwrite is not set..."
    }

    def pangenome_prep_directory = file("${pangenome_directory}/Prep_${pg_name}")
    def pangenome_subset_directory = file("${pangenome_prep_directory}/Subset_Reads")
    def group_file = file("${pangenome_prep_directory}/Read_Groups.csv")
    
    def validation_directory = file("${pangenome_directory}/Validation")
    def validation_reads_directory = file("${validation_directory}/Reads")

    def delete_cmd = (params.overwrite) ? "rm -rf $pangenome_directory" : ":"

    """
    mkdir -p $output_directory &&
    $delete_cmd &&
    mkdir $pangenome_directory &&
    mkdir $pangenome_prep_directory &&
    mkdir $pangenome_subset_directory &&
    mkdir $validation_directory &&
    mkdir $validation_reads_directory &&
    python ${fetchPGScript} --read_dir $read_directory --val_dir $validation_reads_directory --ext $params.pg_ext --forward $params.pg_forward --reverse $params.pg_reverse --group $group_file 
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

    def size = params.size as Integer
    
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
    clusterOptions = "--nodes=${node} --ntasks-per-node=${ray_cpu} --exclusive"
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

    cpus cpu

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
    cd ${pangenome_directory}
    bbmap.sh ref=${pangenome_file}
    samtools faidx ${pangenome_file}
    echo -n "${pg_name},${pangenome_file}"
    """
}

///// Get genome based on FASTA, and index if necessary /////
workflow prepareGenome{

    take:
    fasta_file

    emit:
    return_pangenome
    
    
    main:
        
    if (!fasta_file.isFile()) {
        error "Assembly provided by --fasta (${fasta_file}) does not exist..."
    }

    // Check if ref and .fai exist in the directory
    index_check = CHECK_INDEX(fasta_file) | splitCsv()

    branched = index_check.branch { it ->
        needs_indexing: (it[0].toString() == "BBMap_Absent" || it[1].toString() == "SAM_Absent")
        already_indexed: true
    }

    // Create ref and fai
    reindexed = branched.needs_indexing
    .map { fasta_file }
    | INDEX_FASTA
    | splitCsv()

    // ref and fai already exist
    existing = branched.already_indexed
    .map { tuple(it[2], it[3]) }

    return_pangenome = reindexed.concat(existing).collect().flatten().collate(2)
}

process CHECK_INDEX{
    executor = "local"
    cpus 1

    input:
    val(fasta_path)

    output:
    stdout

    script:
    
    def fasta_file = file("${fasta_path}") 
    def fasta_dir = fasta_file.getParent()
    def fasta_name = fasta_file.getName()
    def fasta_basename = fasta_file.getBaseName()
    
    def bbmap_ref = file("${fasta_dir}/ref")
    def sam_idx = file("${fasta_dir}/${fasta_name}.fai")

    def bbmap_check = bbmap_ref.exists()
    ? "BBMap_Present"
    : "BBMap_Absent"

    def sam_check = sam_idx.exists()
    ? "SAM_Present"
    : "SAM_Absent"

    """
    echo -n "${bbmap_check},${sam_check},${fasta_basename},${fasta_file}"
    """
}

process INDEX_FASTA{

    cpus cpu
    
    input:
    val(fasta_path)

    output:
    stdout

    script:
    
    def fasta_file = file("${fasta_path}") 
    def fasta_dir = fasta_file.getParent()
    def fasta_name = fasta_file.getName()
    def fasta_basename = fasta_file.getBaseName()
    
    def bbmap_ref = file("${fasta_dir}/ref")
    def sam_idx = file("${fasta_dir}/${fasta_name}.fai")

    def bbmap_cmd = bbmap_ref.exists()
    ? ":"
    : "bbmap.sh ref=${fasta_file}"

    def sam_cmd = sam_idx.exists()
    ? ":"
    : "samtools faidx ${fasta_file}"

    """
    cd ${fasta_dir} &&
    ${bbmap_cmd} &&
    ${sam_cmd} &&
    echo -n "${fasta_basename},${fasta_file}"
    """
}

///// Get genome based on pg_name (Look in SNPRS_Pangenomes) /////
workflow checkSNPRSGenome{

    take:
    pangenome_directory
    pg_name

    emit:
    pangenome_info

    main:
    
    genome_dir = file("${pangenome_directory}/${pg_name}")

    if(genome_dir.isDirectory()){
        pangenome_info = CHECK_GENOME(genome_dir) | splitCsv()
    } else{
        error "Directory ${genome_dir} does not exist..."
    }
}

process CHECK_GENOME{

    executor = "local"
    cpus 1

    input:
    val(genome_dir)

    output:
    stdout

    script:
    
    def genome_dir = file("${genome_dir}") 
    def bbmap_ref = file("${genome_dir}/ref")

    """
    cd ${genome_dir}

    get_fasta_from_fai() {
        local folder="\$1"
        local fai_files=("\$folder"/*.fai)
        if [[ ! -e "\${fai_files[0]}" ]]; then
            echo "Error: No .fai file found in \$folder" >&2
            return 1
        fi

        if (( \${#fai_files[@]} != 1 )); then
            echo "Error: Expected exactly 1 .fai file in \$folder, found \${#fai_files[@]}" >&2
            return 1
        fi

        local fasta="\${fai_files[0]%.fai}"
        echo -n "\$fasta"
    }

    FASTA=\$(get_fasta_from_fai "${genome_dir}")
    FASTA_NAME=\$(basename "\${FASTA%.*}")

    if [[ ! -d "${bbmap_ref}" ]]; then
        echo "Error: Directory ${bbmap_ref} does not exist" >&2
        exit 1
    elif [[ ! -f "\$FASTA" ]]; then
        echo "Error: FAI base file \$FASTA does not exist" >&2
        exit 1
    fi

    echo -n "\$FASTA_NAME,\$FASTA"
    """
}