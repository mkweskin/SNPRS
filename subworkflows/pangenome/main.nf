// Subworkflow for generating SNPRS Pangenome
pangenome_directory = file(params.pangenome_directory)
new_pangenome_directory = file("${pangenome_directory}/${params.pg_name}")
pangenome_subset_directory = file("${new_pangenome_directory}/Subset_Reads")

cpu = params.cores as Integer
node = params.nodes as Integer
ray_cores = cpu * node

workflow makePangenome{
        
        take:
        pangenome_read_info

        emit:
        subset_reads

        main:

        // Fetch reads and generate pangenome directory
        input_pangenome_reads = fetchPGReads(pangenome_read_info) 
        | splitCsv

        // Get read counts and base counts
        pangenome_read_info = input_pangenome_reads
        | map{it -> tuple(it[0],it[2])}
        | getReadInfo
        | splitCsv

        // Generate subsetting values
        subset_reads = input_pangenome_reads.join(pangenome_read_info,by:0)
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
}

process fetchPGReads{

    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val dir // Directory containing read files

    output:
    stdout

    script:

    fetchPGScript = file("${projectDir}/bin/fetchPangenomeReads.py")

    """
    # Error if new_pangenome_directory exists
    if [ -d ${new_pangenome_directory} ]; then
        echo "Pangenome directory already exists...exiting..."
        exit 1
    fi
    python ${fetchPGScript} --read_dir ${dir} --read_filetype $params.readext --forward_suffix $params.forward --reverse_suffix $params.reverse && 
    mkdir -p $pangenome_directory &&
    mkdir $new_pangenome_directory &&
    mkdir $pangenome_subset_directory
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
    pg_info_file = file("${new_pangenome_directory}/Pangenome_Read_Info.csv")
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
            
        """
        $reformat_cmd
        read_count_line=\$(cat ${pangenome_subset_directory}/out_Subsample_${sample_name} | grep "Output:")
        read_count=\$(echo "\$read_count_line" | awk '{print \$2}')
        base_count=\$(echo "\$read_count_line" | awk '{print \$5}')
        echo "${sample_name},\$read_count,\$base_count"
        """
}

process assemblePangenome{

    clusterOptions = "--nodes=${node} --ntasks-per-node=${cpu}"

    input:
    val(subset_reads)

    output:
    stdout

    script:

    """
    echo "Contents of \$(pwd):"
    ls -lh
    echo "Contents of ${pangenome_subset_directory}:"
    ls -lh ${pangenome_subset_directory}
    mpirun -v -np ${ray_cores} Ray -k 31 -detect-sequence-files ${pangenome_subset_directory} -o ${new_pangenome_directory}/Ray_${params.pg_name}
    """
}























workflow fetchPangenome{

    take:
    pangenome_directory

    emit:
    pangenome_info

    main:
    
    print("Fetching pangenome...")
    pangenome_info = "Temp"
}
