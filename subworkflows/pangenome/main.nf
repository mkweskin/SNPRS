// Subworkflow for generating or fetching SNPRS Pangenome
pg_name = "${params.pg_name}" != "" ? "${params.pg_name}" : "${new java.util.Date().getTime()}"

base_pangenome_directory = file("${params.snprs_directory}/SNPRS_Pangenomes")
processed_pangenome_directory = file("${base_pangenome_directory}/${pg_name}")

pangenome_prep_directory = file("${base_pangenome_directory}/Prep_${pg_name}") 
pangenome_subset_directory = file("${pangenome_prep_directory}/Subset_Reads")
pangenome_conf_file = file("${pangenome_prep_directory}/Ray.conf")
ray_directory = file("${pangenome_prep_directory}/Ray_${pg_name}")
ref_path = "${params.ref_path}" == "" ? file("${ray_directory}/Contigs.fasta") : file("${params.ref_path}") // Fix for Ray

cpu = params.cores as Integer
node = params.nodes as Integer
ray_kmer = params.ray_kmer as Integer
ray_cores = cpu * node

workflow makePangenome{
        
        emit:
        pagenome_directory

        main:

        // Fetch reads and generate pangenome directory
        input_pangenome_reads = fetchPGReads(file("${params.pg_reads}")) 
        | splitCsv

        // Get read counts and base counts
        pangenome_read_info = input_pangenome_reads
        | map{it -> tuple(it[0],it[2])}
        | getReadInfo
        | splitCsv

        // Generate subsetting values
        pagenome_directory = input_pangenome_reads.join(pangenome_read_info,by:0)
        | map{it -> it.join(",")}
        | collect
        | getSubsetValues
        | splitCsv
        | collect
        | flatten
        | collate(6)
        | map{it -> tuple(it[0],it[2],it[5])}
        | subsetReads
        | collect
        | assemblePangenome
        | processPangenome
}

process fetchPGReads{

    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val read_info

    output:
    stdout

    script:

    fetchPGScript = file("${projectDir}/bin/fetchPangenomeReads.py")

    """
    # Error if new_pangenome_directory exists
    if [ -d ${pangenome_prep_directory} ]; then
        echo "Pangenome directory already exists...exiting..."
        exit 1
    fi
    mkdir -p $base_pangenome_directory &&
    mkdir $pangenome_prep_directory &&
    mkdir $pangenome_subset_directory &&
    python ${fetchPGScript} --read_dir ${read_info} --read_filetype $params.readext --forward_suffix $params.forward --reverse_suffix $params.reverse
    """
}

process getReadInfo {

    cpus = 1
    memory '4 GB'

    input:
    tuple val(sample_name), val(read_location)

    output:
    stdout

    script:
    def forward_reverse = read_location.contains(";") ? read_location.split(";") : [read_location]
    def reformat_cmd = forward_reverse.size() == 2 ? 
        "reformat.sh in=${forward_reverse[0]} in2=${forward_reverse[1]} 2>&1" :
        "reformat.sh in=${forward_reverse[0]} 2>&1"
        
    """
    output=\$(${reformat_cmd})
    read_count_line=\$(echo "\$output" | grep "Output:")
    read_count=\$(echo "\$read_count_line" | awk '{print \$2}')
    base_count=\$(echo "\$read_count_line" | awk '{print \$5}')
    echo "${sample_name},\$read_count,\$base_count"
    """
}

process getSubsetValues{

    cpus = 1
    executor = 'local'

    input:
    val(pangenome_read_info)

    output:
    stdout

    script:
    pg_info_file = file("${pangenome_prep_directory}/Pangenome_Read_Info.csv")
    read_subset_script = file("${projectDir}/bin/readSubsetter.py")

    genome_size = params.pg_size.toInteger()
    genome_coverage = params.pg_coverage.toInteger()

    pg_header = "Sample_ID,Pangenome_Group,Read_Location,Read_Count,Base_Count"
    pg_info = pangenome_read_info.join("\n")
    
    """
    echo $pg_header > $pg_info_file &&
    echo -e "${pg_info}" >> $pg_info_file &&
    python $read_subset_script -g $genome_size -c $genome_coverage -r $pg_info_file -o $pangenome_subset_directory
    """
}

process subsetReads{
    
        cpus = 1
        memory '4 GB'
    
        input:
        tuple val(sample_name), val(read_location), val(subset_count)
    
        output:
        stdout
    
        script:
        def forward_reverse = read_location.contains(";") ? read_location.split(";") : [read_location]
        def reformat_cmd = forward_reverse.size() == 2 ? 
            "reformat.sh in=${forward_reverse[0]} in2=${forward_reverse[1]} out=${pangenome_subset_directory}/${sample_name}_GenomeReads_1.fq.gz out2=${pangenome_subset_directory}/${sample_name}_GenomeReads_2.fq.gz samplebasestarget=${subset_count} &> ${pangenome_subset_directory}/out_Subsample_${sample_name}" :
            "reformat.sh in=${forward_reverse[0]} out=${pangenome_subset_directory}/${sample_name}_GenomeReads.fq.gz samplebasestarget=${subset_count} &> ${pangenome_subset_directory}/out_Subsample_${sample_name}"
        def subset_count = subset_count.toInteger()
        def echo_command = read_location.contains(";") ? "-p ${pangenome_subset_directory}/${sample_name}_GenomeReads_1.fq.gz ${pangenome_subset_directory}/${sample_name}_GenomeReads_2.fq.gz" : "-s ${pangenome_subset_directory}/${sample_name}_GenomeReads.fq.gz"
            
        """
        $reformat_cmd
        read_count_line=\$(cat ${pangenome_subset_directory}/out_Subsample_${sample_name} | grep "Output:")
        read_count=\$(echo "\$read_count_line" | awk '{print \$2}')
        base_count=\$(echo "\$read_count_line" | awk '{print \$5}')
        echo -n $echo_command
        """
}

process assemblePangenome{

    //clusterOptions = "--nodes=${node} --ntasks-per-node=${cpu}"
    executor = 'local'

    input:
    val(subset_reads)

    output:
    stdout

    script:
    subset_read_locs = subset_reads.join("\n")

    """
    echo "-k ${ray_kmer}" > ${pangenome_conf_file} &&
    echo "-o ${ray_directory}" >> ${pangenome_conf_file} &&
    echo -e "${subset_read_locs}" >> ${pangenome_conf_file}
    # mpirun -np ${ray_cores} Ray ${pangenome_conf_file} &> ${pangenome_prep_directory}/out_Ray_${params.pg_name}
    echo -n $ref_path
    """
}

process processPangenome{
    
    input:
    val(contig_file)

    output:
    stdout

    script:

    genome_script=file("${projectDir}/bin/Genome_SiteLengths.py")
    new_fasta = file("${processed_pangenome_directory}/contigs.fa")
    log_dir = file("${processed_pangenome_directory}/logs")
    
    """
    # Error if new_pangenome_directory exists
    if [ -d ${processed_pangenome_directory} ]; then
        echo "Pangenome directory already exists...exiting..."
        exit 1
    fi
    mkdir -p $base_pangenome_directory &&
    mkdir $processed_pangenome_directory &&
    mkdir $log_dir &&
    rename.sh in=${contig_file} out=${new_fasta} prefix=SNPRS addprefix=t trd=t &> ${log_dir}/out_Rename &&
    bowtie2-build ${new_fasta} ${processed_pangenome_directory}/contigs -p $cpu &> ${log_dir}/out_Bowtie2 &&
    bbmap.sh ref=${new_fasta} path=${processed_pangenome_directory} &> ${log_dir}/out_BBMap &&
    samtools faidx ${new_fasta} &> ${log_dir}/out_Samtools &&
    python $genome_script $processed_pangenome_directory &> ${log_dir}/out_Genome_SiteLengths &&
    stats.sh in=${new_fasta} &> ${processed_pangenome_directory}/BBStats &&
    echo -n $processed_pangenome_directory
    """
}

workflow fetchPangenome{

    emit:
    pangenome_directory

    main:

    if("${params.pg_path}" != ""){
        pangenome_directory = validatePangenome(file("${params.pg_path}"))
    }
    else if("${params.ref_path}" != ""){
        pangenome_directory = processPangenome(file("${params.ref_path}"))
    } 
    else if("${params.pg_name}" != ""){
        pangenome_directory = validatePangenome(processed_pangenome_directory)
    }
    else{
        error "No pangenome information provided...exiting..."
    }
}

process validatePangenome {

    executor = 'local'
    cpu = 1
    maxForks = 1

    input:
    val pangenome_directory

    output:
    stdout

    script:
    """
    files=(
        'BBStats' 'contigs.1.bt2' 'contigs.2.bt2' 'contigs.3.bt2' 'contigs.4.bt2'
        'contigs.fa' 'contigs.fa.fai' 'contigs_LocList' 'contigs.rev.1.bt2'
        'contigs.rev.2.bt2' 'contigs_SeqLength.tsv' 'logs' 'ref'
    )

    for file in "\${files[@]}"; do
        if [ ! -e "${pangenome_directory}/\${file}" ]; then
            echo "Error: \${file} is missing in ${pangenome_directory}" >&2
            exit 1
        fi
    done

    echo -n ${pangenome_directory}
    """
}