#! /usr/bin/env nextflow
nextflow.enable.dsl=2

kmer = params.kmer as Integer

cpu = params.cpus as Integer
sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu

node = params.nodes as Integer
ray_cores = cpu * node

///// Create a SNPRS pangenome from reads /////
workflow assembleGenome{

    take:
    genome_directory
    genome_name
    pg_read_data

    emit:
    return_pangenome

    main:

    if(!params.size){
        error "Cannot assemble pangenome without --size estimate (genome size in basepairs)"
    }

    input_pangenome_reads = FETCH_PG_READS(genome_directory,genome_name,pg_read_data) | splitCsv

    base_count_file = input_pangenome_reads
    .map{it->tuple(it[0],it[3],it[4],"${genome_directory}","${genome_name}")}
    | COUNT_BASES | collect | map { it[0] }

    subset_guide = CALCULATE_SUBSETS(base_count_file,genome_directory,genome_name) | splitCsv
    | branch{it ->
        link: it[5].toString() == "Link"
            return(tuple(it[0],it[1],it[2],it[3]))
        sample: true
            return(tuple(it[0],it[1],it[2],it[3],it[4]))
    }

    subset_folder_1 = SUBSET_READS(subset_guide.sample,genome_directory,genome_name)
    subset_folder_2 = LINK_READS(subset_guide.link,genome_directory,genome_name)

    subset_folder = subset_folder_1
    .concat(subset_folder_2)  
    | collect | map { it[0] }

    ray_assembly = ASSEMBLE_PANGENOME(subset_folder,genome_directory,genome_name)

    return_pangenome = PROCESS_RAY(ray_assembly,genome_directory,genome_name) | splitCsv | collect | flatten | collate(3)
}

process FETCH_PG_READS{

    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(genome_directory)
    val(genome_name)
    val(read_directory)

    output:
    stdout

    script:

    def fetchPGScript = file("${projectDir}/bin/fetchPangenomeReads.py")

    def genome_directory = file("${genome_directory}")
    def genome_prep_directory = file("${genome_directory}/Prep_${genome_name}")
    
    if(!params.overwrite && genome_prep_directory.isDirectory()){
        error "${genome_prep_directory} exists and --overwrite is not set..."
    }

    def subset_directory = file("${genome_prep_directory}/Subset_Reads")
    def group_file = file("${genome_prep_directory}/Read_Groups.csv")

    def delete_cmd = (params.overwrite) ? "rm -rf $genome_prep_directory" : ":"
    
    """
    mkdir -p $genome_directory &&
    $delete_cmd &&
    mkdir $genome_prep_directory &&
    mkdir $subset_directory &&
    python ${fetchPGScript} --read_dir $read_directory --ext $params.pg_ext --forward $params.pg_forward --reverse $params.pg_reverse --group $group_file 
    """
}

process COUNT_BASES {

    tag "Count_${sample_id}"

    cpus sample_cpu

    input:
    tuple val(sample_id), val(forward), val(reverse), val(genome_directory), val(genome_name)

    output:
    stdout

    script:
    
    def genome_directory = file("${genome_directory}")
    def genome_prep_directory = file("${genome_directory}/Prep_${genome_name}")

    def base_count_file = file("${genome_prep_directory}/Read_Counts.csv")
    
    def stats_cmd = reverse 
    ? "seqkit stats -j 1 -a -T ${forward} ${reverse}" 
    : "seqkit stats -j 1 -a -T ${forward}"

    """
    output=\$(${stats_cmd})
    read_count=\$(echo "\$output" | awk -F'\\t' 'NR>1 {sum+=\$4} END{print sum}')
    base_count=\$(echo "\$output" | awk -F'\\t' 'NR>1 {sum+=\$5} END{print sum}')
    echo -e "${sample_id},\$read_count,\$base_count,${forward},${reverse}" >> ${base_count_file}
    echo -n ${base_count_file}
    """
}

process CALCULATE_SUBSETS{

    executor = 'local'
    cpus = 1

    input:
    val(base_count_file)
    val(genome_directory)
    val(genome_name)

    output:
    stdout

    script:

    def calculate_sub_script = file("${projectDir}/bin/calculateSubset.py")
    
    def genome_directory = file("${genome_directory}")
    def genome_prep_directory = file("${genome_directory}/Prep_${genome_name}")
    def subset_directory = file("${genome_prep_directory}/Subset_Reads")

    def group_file = file("${genome_prep_directory}/Read_Groups.csv")

    def size = params.size as Integer
    def out_prop = params.out_prop as Float
    def coverage = params.coverage as Integer
    
    """
    python ${calculate_sub_script} -b ${base_count_file} -g ${group_file} -s ${size} -c ${coverage} -o ${subset_directory} -p ${out_prop}
    """
}

process SUBSET_READS {

    tag "Subset_${sample_id}"

    cpus sample_cpu

    input:
    tuple val(sample_id), val(subsample_id),val(forward),val(reverse),val(allocated)
    val(genome_directory)
    val(genome_name)

    output:
    stdout

    script:

    def forward_ext = "${params.pg_forward}"
    def reverse_ext = "${params.pg_reverse}"
    def file_ext = "${params.pg_ext}"

    def genome_directory = file("${genome_directory}")
    def genome_prep_directory = file("${genome_directory}/Prep_${genome_name}")
    def subset_directory = file("${genome_prep_directory}/Subset_Reads")

    def genome_read_link_directory = file("${genome_directory}/Pangenome_Read_Links")

    def out1 = "${subset_directory}/${subsample_id}_GenomeReads${forward_ext}"
    def out2 = "${subset_directory}/${subsample_id}_GenomeReads${reverse_ext}"
    def outs = "${subset_directory}/${subsample_id}_GenomeReads${file_ext}"

    def log_file = "${subset_directory}/out_Subsample_${sample_id}"

    def reformat_cmd = reverse
        ? "reformat.sh in=${forward} in2=${reverse} out=${out1} out2=${out2} outs=${outs} samplebasestarget=${allocated} &> ${log_file}"
        : "reformat.sh in=${forward} out=${outs} samplebasestarget=${allocated} &> ${log_file}"

    def link_cmd = reverse
        ? "ln -s $forward ${genome_read_link_directory}/${sample_id}${forward_ext}; ln -s $reverse ${genome_read_link_directory}/${sample_id}${reverse_ext}"
        : "ln -s $forward ${genome_read_link_directory}/${sample_id}${file_ext}"

    """
    mkdir -p $genome_read_link_directory
    $link_cmd
    $reformat_cmd

    has_content () {
        zcat "\$1" 2>/dev/null | head -c 1 | grep -q .
    }
    for f in ${out1} ${out2} ${outs}; do
        if [[ -f "\$f" ]] && ! has_content "\$f"; then
            rm -f "\$f"
        fi
    done

    echo -n "${subset_directory}"
    """
}

process LINK_READS {
    cpus 1
    executor = "local"

    input:
    tuple val(sample_id), val(subsample_id),val(forward),val(reverse)
    val(genome_directory)
    val(genome_name)

    output:
    stdout

    script:   

    def forward_ext = "${params.pg_forward}"
    def reverse_ext = "${params.pg_reverse}"
    def file_ext = "${params.pg_ext}"

    def genome_directory = file("${genome_directory}")
    def genome_prep_directory = file("${genome_directory}/Prep_${genome_name}")
    def subset_directory = file("${genome_prep_directory}/Subset_Reads")

    def genome_read_link_directory = file("${genome_directory}/Pangenome_Read_Links")

    def out1 = "${subset_directory}/${subsample_id}_GenomeReads${forward_ext}"
    def out2 = "${subset_directory}/${subsample_id}_GenomeReads${reverse_ext}"
    def outs = "${subset_directory}/${subsample_id}_GenomeReads_SE${file_ext}"

    def link_cmd_1 = reverse
        ? "ln -s ${forward} ${out1}; ln -s ${reverse} ${out2}"
        : "ln -s ${forward} ${outs}"

    def link_cmd_2 = reverse
        ? "ln -s $forward ${genome_read_link_directory}/${sample_id}${forward_ext}; ln -s $reverse ${genome_read_link_directory}/${sample_id}${reverse_ext}"
        : "ln -s $forward ${genome_read_link_directory}/${sample_id}${file_ext}"

    """
    $link_cmd_1 &&
    $link_cmd_2 &&
    echo -n "${subset_directory}"
    """
}

process ASSEMBLE_PANGENOME {

    tag "Assemble_Pangenome"

    clusterOptions = "--nodes=${node} --ntasks-per-node=${cpu} --exclusive"
    
    input:
    val(subset_folder)
    val(genome_directory)
    val(genome_name)

    output:
    stdout

    script:

    def genome_directory = file("${genome_directory}")
    def genome_prep_directory = file("${genome_directory}/Prep_${genome_name}")
    def subset_directory = file(subset_folder)
    
    def assembly_directory = file("${genome_prep_directory}/Ray_${genome_name}")
    def ray_log = file("${genome_prep_directory}/out_Ray_${genome_name}")
    def load_ray_module = (params.ray_module) ? "module load -s ${params.ray_module}" : ":"

    """
    $load_ray_module
    mpirun --use-hwthread-cpus -np ${ray_cores} Ray -k ${kmer} -detect-sequence-files ${subset_folder} -o ${assembly_directory}  &> ${ray_log} &&
    echo -n "${assembly_directory}/Contigs.fasta"
    """
}

process PROCESS_RAY{

    cpus 1

    input:
    val(ray_assembly)
    val(genome_directory)
    val(genome_name)

    output:
    stdout

    script:
    def index_script = file("${projectDir}/bin/contig_idx.py")

    def genome_directory = file("${genome_directory}")
    def genome_prep_directory = file("${genome_directory}/Prep_${genome_name}")

    def genome_file = file("${genome_directory}/${genome_name}.fasta")
    def stats_file = file("${genome_directory}/${genome_name}_BBStats")

    def min_contig = params.min_contig as Integer

    """
    rename.sh in=${ray_assembly} out=${genome_file} prefix=SNPRS addprefix=t trd=t minscaf=${min_contig}
    stats.sh ${genome_file} &> ${stats_file}
    cd ${genome_directory}
    samtools faidx ${genome_file}
    python $index_script --fasta $genome_file --make_parquet
    echo -n "${genome_name},${genome_directory},${genome_file}"
    """
}


///// Get genome based on FASTA, and index if necessary /////
workflow useFASTA{

    take:
    genome_directory
    genome_name
    fasta_file

    emit:
    return_pangenome
    
    
    main:
        
    return_pangenome = USE_FASTA(genome_directory,genome_name,fasta_file) | splitCsv | collect | flatten | collate(3)
}

process USE_FASTA{

    cpus 1
    
    input:
    val(genome_directory)
    val(genome_name)
    val(fasta_path)

    output:
    stdout

    script:
    
    def index_script = file("${projectDir}/bin/contig_idx.py")

    def genome_directory = file("${genome_directory}")

    def fasta_file = file("${fasta_path}")
    def genome_file = file("${genome_directory}/${genome_name}.fasta")

    def link_cmd = ("${fasta_file}" == "${genome_file}") ? ":" : "ln -s $fasta_file $genome_file"

    def sam_idx = file("${genome_file}.fai")
    def index_parquet = file("${genome_directory}/${genome_name}.parquet")

    def sam_cmd = sam_idx.exists()
    ? ":"
    : "samtools faidx ${genome_file}"

    def idx_cmd = index_parquet.exists()
    ? ":"
    : "python $index_script --fasta $genome_file --make_parquet"
    

    """
    $link_cmd &&
    cd ${genome_directory} &&
    ${sam_cmd} &&
    ${idx_cmd} &&
    echo -n "${genome_name},${genome_directory},${genome_file}"
    """
}

///// Get genome based on genome_name (Look in SNPRS_Pangenomes) /////
workflow checkGenomeDir{

    take:
    genome_directory

    emit:
    pangenome_info

    main:
    
    genome_dir = file(genome_directory)

    if(!genome_dir.isDirectory()){
        error "Directory ${genome_dir} does not exist..."
    }

    pangenome_info = CHECK_GENOME_DIR(genome_dir) | splitCsv | collect | flatten | collate(3)
}

process CHECK_GENOME_DIR{

    cpus 1

    input:
    val(genome_directory)

    output:
    stdout

    script:
    
    def genome_dir = file(genome_directory)

    """
    cd ${genome_dir}

    get_fasta_from_fai() {
        local folder="\$1"
        local fai_files=( "\$folder"/*.fai )

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
    PARQUET_FILE="\${FASTA_NAME}.parquet"

    if [[ ! -f "\$PARQUET_FILE" ]]; then
        echo "Error: Parquet index file \$PARQUET_FILE does not exist in \$(pwd)" >&2
        exit 1
    elif [[ ! -f "\$FASTA" ]]; then
        echo "Error: FASTA file \$FASTA does not exist" >&2
        exit 1
    fi

    echo -n "\$FASTA_NAME,${genome_dir},\$FASTA"
    """
}