#! /usr/bin/env nextflow
nextflow.enable.dsl=2

def kmer = params.kmer as Integer
def cpu = params.cpus as Integer
def sample_cpu = (params.sample_cpus) ? params.sample_cpus as Integer : cpu
def count_cpu = (cpu >= 2) ? 2 : cpu

def node = params.nodes as Integer
def ray_cores = cpu * node

def read_ext = params.read_ext
def forward = params.forward
def reverse = params.reverse

def size = (params.size) ? params.size as Integer : 0
def out_prop = params.out_prop as Float
def coverage = params.coverage as Integer
def min_contig = params.min_contig as Integer

///// Create a SNPRS pangenome from reads /////

workflow assembleGenome{

    take:

    pg_reads
    genome_directory
    genome_name

    emit:
    return_pangenome

    main:
 
    // Sample_ID, Group_0, Group_1, Forward, Reverse
    input_pangenome_reads = FETCH_PG_READS(pg_reads,genome_directory,genome_name) | splitCsv

    // Get base counts
    base_count_file = (params.manual_counts) ? Channel.fromPath(params.manual_counts) | collect | map { it[0] } : input_pangenome_reads.map{it-> tuple(it[0],it[3],it[4],genome_directory,genome_name)} | COUNT_BASES  | collect | map { it[0] }

    // Subset or link reads for pangenome assembly
    subset_guide = CALCULATE_SUBSETS(base_count_file,genome_directory,genome_name) | splitCsv
    | branch{it ->
        link: it[5].toString() == "Link"
            return(tuple(it[0],it[1],it[2],it[3],it[6]))
        sample: true
            return(tuple(it[0],it[1],it[2],it[3],it[4],it[6]))
    }

    subset_folder_1 = subset_guide.sample | SUBSET_READS
    subset_folder_2 = subset_guide.link | LINK_READS

    subset_folder = subset_folder_1.concat(subset_folder_2) | collect | map { it[0] }

    // Assemble pangenome
    ray_fasta = ASSEMBLE_PANGENOME(subset_folder,genome_directory,genome_name)

    // Index pangenome
    return_pangenome = PROCESS_RAY(ray_fasta,genome_directory,genome_name) | splitCsv | collect | flatten | collate(3) 
}

process FETCH_PG_READS{

    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(pg_read_data)
    val(genome_dir)
    val(genome_name)

    output:
    stdout

    script:

    fetchPGScript = file("${projectDir}/bin/fetchPangenomeReads.py")

    pangenome_read_link_directory = file("${genome_dir}/Pangenome_Read_Links")
    genome_prep_directory = file("${genome_dir}/Prep_${genome_name}")

    subset_directory = file("${genome_prep_directory}/Subset_Reads")
    group_file = file("${genome_prep_directory}/Read_Groups.csv")

    delete_cmd = (params.overwrite) ? "rm -rf $genome_dir"
    : """
if [ -d "$genome_dir" ] ; then
    echo "❌ Error: $genome_dir already exists! Use --overwrite to replace." >&2
    exit 1
fi"""    

    
    """
    $delete_cmd &&
    mkdir $genome_dir &&
    mkdir $genome_prep_directory &&
    mkdir $subset_directory &&
    mkdir $pangenome_read_link_directory &&
    python ${fetchPGScript} --read_dir $pg_read_data --ext $read_ext --forward $forward --reverse $reverse --group $group_file --link_dir $pangenome_read_link_directory
    """
}

process COUNT_BASES {

    tag "Count_${sample_id}"

    cpus count_cpu

    input:
    tuple val(sample_id), val(forward_read), val(reverse_read), val(genome_dir), val(genome_name)

    output:
    stdout

    script:

    genome_prep_directory = file("${genome_dir}/Prep_${genome_name}")
    base_count_file = file("${genome_prep_directory}/Read_Counts.csv")
    
    stats_cmd = reverse_read 
    ? "seqkit stats -j $count_cpu -a -T ${forward_read} ${reverse_read}" 
    : "seqkit stats -j $count_cpu -a -T ${forward_read}"

    """
    output=\$(${stats_cmd})
    read_count=\$(echo "\$output" | awk -F'\\t' 'NR>1 {sum+=\$4} END{print sum}')
    base_count=\$(echo "\$output" | awk -F'\\t' 'NR>1 {sum+=\$5} END{print sum}')
    echo -e "${sample_id},\$read_count,\$base_count,${forward_read},${reverse_read}" >> ${base_count_file}
    echo -n ${base_count_file}
    """
}

process CALCULATE_SUBSETS{

    executor = 'local'
    cpus = 1

    input:
    val(base_count_file)
    val(genome_dir)
    val(genome_name)

    output:
    stdout

    script:

    calculate_sub_script = file("${projectDir}/bin/calculateSubset.py")
    
    genome_prep_directory = file("${genome_dir}/Prep_${genome_name}")
    subset_directory = file("${genome_prep_directory}/Subset_Reads")
    group_file = file("${genome_prep_directory}/Read_Groups.csv")

    data_args = (params.manual_counts) ? "-m ${file(params.manual_counts)}" : "-b ${base_count_file} -g ${group_file}"
    
    """
    python ${calculate_sub_script} $data_args -s ${size} -c ${coverage} -o ${subset_directory} -p ${out_prop}
    """
}

process SUBSET_READS {

    tag "Subset_${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), val(subsample_id),val(forward_read),val(reverse_read),val(allocated),val(subset_dir)

    output:
    stdout

    script:

    subset_directory = file(subset_dir)

    safe_ext = read_ext.startsWith('.') ? read_ext : ".${read_ext}"

    out1 = "${subset_directory}/${subsample_id}_GenomeReads${forward}"
    out2 = "${subset_directory}/${subsample_id}_GenomeReads${reverse}"
    outs = "${subset_directory}/${subsample_id}_GenomeReads${safe_ext}"

    log_file = "${subset_directory}/out_Subsample_${sample_id}"

    reformat_cmd = reverse_read
        ? "reformat.sh in=${forward_read} in2=${reverse_read} out=${out1} out2=${out2} outs=${outs} samplebasestarget=${allocated} &> ${log_file}"
        : "reformat.sh in=${forward_read} out=${outs} samplebasestarget=${allocated} &> ${log_file}"

    """
    $reformat_cmd

    has_content() {
        local file="\$1"

        if [[ ! -f "\$file" ]]; then
            echo "Error: File '\$file' not found" >&2
            return 1
        fi

        if [[ "\$file" == *.gz ]]; then
            zcat -- "\$file" 2>/dev/null | head -c 1 | grep -q .
        else
            cat -- "\$file" 2>/dev/null | head -c 1 | grep -q .
        fi
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
    tuple val(sample_id), val(subsample_id),val(forward_read),val(reverse_read),val(subset_dir)

    output:
    stdout

    script:   

    subset_directory = file(subset_dir)

    safe_ext = read_ext.startsWith('.') ? read_ext : ".${read_ext}"

    out1 = "${subset_directory}/${subsample_id}_GenomeReads${forward}"
    out2 = "${subset_directory}/${subsample_id}_GenomeReads${reverse}"
    outs = "${subset_directory}/${subsample_id}_GenomeReads${safe_ext}"

    link_cmd = reverse_read
        ? "ln -s ${forward_read} ${out1}; ln -s ${reverse_read} ${out2}"
        : "ln -s ${forward_read} ${outs}"

    """
    $link_cmd &&
    echo -n "${subset_directory}"
    """
}

process ASSEMBLE_PANGENOME {

    tag "Assemble_Pangenome"

    clusterOptions = "--nodes=${node} --ntasks-per-node=${cpu} --exclusive"
    
    input:
    val(subset_folder)
    val(genome_dir)
    val(genome_name)

    output:
    stdout

    script:
    
    genome_prep_directory = file("${genome_dir}/Prep_${genome_name}")
    assembly_directory = file("${genome_prep_directory}/Ray_${genome_name}")
    ray_log = file("${genome_prep_directory}/out_Ray_${genome_name}")
    load_ray_module = (params.ray_module) ? "module load -s ${params.ray_module}" : ":"

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
    val(genome_dir)
    val(genome_name)

    output:
    stdout

    script:
    index_script = file("${projectDir}/bin/contig_idx.py")

    genome_file = file("${genome_dir}/${genome_name}.fasta")
    stats_file = file("${genome_dir}/${genome_name}_BBStats")

    """
    rename.sh in=${ray_assembly} out=${genome_file} prefix=SNPRS addprefix=t trd=t minscaf=${min_contig}
    stats.sh ${genome_file} &> ${stats_file}
    cd ${genome_dir}
    samtools faidx ${genome_file}
    python $index_script --fasta $genome_file --make_parquet
    echo -n "${genome_name},${genome_dir},${genome_file}"
    """
}








///// Get genome based on FASTA, and index if necessary /////
workflow useFASTA{

    take:
    fasta_file
    genome_directory
    genome_name

    emit:
    processed_fasta
    
    main:
        
    processed_fasta = USE_FASTA(fasta_file,genome_directory,genome_name) | splitCsv | collect | flatten | collate(3)
}

process USE_FASTA{

    cpus 1
    
    input:
    val(fasta_file)
    val(genome_dir)
    val(genome_name)

    output:
    stdout

    script:
    
    index_script = file("${projectDir}/bin/contig_idx.py")

    fasta_parent = file(fasta_file).getParent()
    
    if(fasta_parent == genome_dir){
        error "Cannot use --fasta if already in the genome directory (move it out and SNPRS will link it)"
    }

    genome_file = file("${genome_dir}/${genome_name}.fasta")
    index_parquet = file("${genome_dir}/${genome_name}.parquet")
    sam_idx = file("${genome_file}.fai")
    stats_file = file("${genome_dir}/${genome_name}_BBStats")

    delete_cmd = (params.overwrite) ? "rm -rf $genome_dir"
    : """
if [ -d "$genome_dir" ] ; then
    echo "❌ Error: $genome_dir already exists! Use --overwrite to replace." >&2
    exit 1
fi"""    

    """
    $delete_cmd &&
    mkdir $genome_dir &&
    cd ${genome_dir} &&
    reformat.sh in=${fasta_file} out=${genome_file} minlength=${min_contig} &&
    stats.sh ${genome_file} &> ${stats_file} &&
    samtools faidx $genome_file &&
    python $index_script --fasta $genome_file --make_parquet &&
    echo -n "${genome_name},${genome_dir},${genome_file}"
    """
}

///// Check folder for all SNPRS genome components /////
workflow checkGenomeDir{

    take:
    genome_dir

    emit:
    pangenome_info

    main:
    
    check_dir = file(genome_dir)
    
    pangenome_info = (check_dir.isDirectory()) ? CHECK_GENOME_DIR(check_dir) | splitCsv | collect | flatten | collate(3) : Channel.empty()
}

process CHECK_GENOME_DIR{

    cpus 1

    input:
    val(genome_dir)

    output:
    stdout

    script:
    
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