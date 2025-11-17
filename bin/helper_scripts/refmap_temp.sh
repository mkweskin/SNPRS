# Map contigs
bowtie2 -p 40 -x Yeast -f -U SNPRS_1762186342821.fasta -S Yeast.sam
samtools view -Su -@ 40 -F 4 Yeast.sam | samtools sort -@ 40 - -o Yeast_Temp.bam
samtools view -@ 40 -H Yeast_Temp.bam > Yeast_Header.sam
samtools view -@ 40 Yeast_Temp.bam | grep -v "XS:" | cat Yeast_Header.sam - | samtools view -@ 40 -b - > Yeast.bam
samtools view Yeast.bam | awk 'BEGIN {OFS = "\t"} { print $1, $3, $4, $2, $6}' > Yeast_MapData.tsv
cut -f1 Yeast_MapData.tsv | sort > Uniquely_Mapping_Contigs
python genome_mapper.py Yeast_MapData.tsv