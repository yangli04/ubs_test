
SAMPLE=["XY-1_S10","XY-2_S11","XY-3_S12","XY-4_S13","XY-5_S14","XY-6_S15"]
ADDSAMPLE=['XY-1_S45','XY-2_S46','XY-3_S47','XY-4_S48','XY-5_S49','XY-6_S50']
#SAMPLE_ALL=SAMPLE+ADDSAMPLE
rule all:
	input:
#		sorted_sam=expand("map_sncRNA/sorted_sam/{sample}.sorted.sam", sample=SAMPLE_ALL),
#		rib_table=expand("map_sncRNA/table/{sample}_hisat3n_table.tsv", sample=SAMPLE),
#		m1="map_sncRNA/calculated_rate/methylation_level_filtered.tsv",	
#		m2="map_sncRNA/calculated_rate/methylation_level_notfiltered.tsv"
		sncrna_aligned=expand("map_sncRNA/{sample}.ribosomal.sam",sample=ADDSAMPLE),
                fq1=expand("sncRNA_depleted/{sample}_R1.fastq",sample=ADDSAMPLE),
                summary=expand("map_sncRNA/{sample}.sncRNA.summary",sample=ADDSAMPLE)

# # PE
# rule trim:
# 	input:
# 		input1="raw_data/{sample}_R1_001.fastq",
# 		input2="raw_data/{sample}_R2_001.fastq"
# 	output:
# 		output1="trimmed/trimmed_{sample}_R1.fastq",
# 		output2="trimmed/trimmed_{sample}_R2.fastq",
# 		tooshort1="trimmed/tooshort/{sample}.tooshort.1.fq",
# 		tooshort2="trimmed/tooshort/{sample}.tooshort.2.fq",
# 		log_file="trimmed/{sample}.cutadapt.log"
# 	shell:
# 		"""cutadapt -j 10 -m 20 -q 15 -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA;anywhere;" -A "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;anywhere;" -u 3 -U 3 -u -3 -U -3 --nextseq-trim=15 -o {output.output1} -p {output.output2} --too-short-output {output.tooshort1} --too-short-paired-output {output.tooshort2} {input.input1} {input.input2} > {output.log_file}"""
# 
# rule map_sncRNA:
# 	input:
# 		input1="trimmed/trimmed_{sample}_R1.fastq",
# 		input2="trimmed/trimmed_{sample}_R2.fastq"
# 	output:
# 		sncrna_aligned="map_sncRNA/{sample}.ribosomal.sam",
# 		fq1="sncRNA_depleted/{sample}_R1.fastq",
# 		fq2="sncRNA_depleted/{sample}_R2.fastq",
# 		summary="map_sncRNA/{sample}.sncRNA.summary"
# 	shell:
# 		"hisat-3n -q -x ~/data/reference/Homo_sapiens/hisat_18S_28SrRNA_CT/GRCh38_rRNA --summary-file {output.summary} --new-summary -1 {input.input1} -2 {input.input2} -S {output.sncrna_aligned} --base-change C,T -p 2 --un-conc sncRNA_depleted/{wildcards.sample}_R%.fastq"
# # PE end 
# 
# 
# SE
rule trim_SE:
	input:
		input1="raw_data/{sample}_R1_001.fastq",
	output:
		output1="trimmed/trimmed_{sample}_R1.fastq",
		tooshort1="trimmed/tooshort/{sample}.tooshort.1.fq",
		log_file="trimmed/{sample}.cutadapt.log"
	shell:
		"""cutadapt -j 10 -m 20 -q 15 -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA;anywhere;" -u 3 -u -3 --nextseq-trim=15 -o {output.output1} --too-short-output {output.tooshort1} {input.input1} > {output.log_file}"""

rule map_sncRNA_SE:
	input:
		input1="trimmed/trimmed_{sample}_R1.fastq",
	output:
		sncrna_aligned="map_sncRNA/{sample}.ribosomal.sam",
		fq1="sncRNA_depleted/{sample}_R1.fastq",
		summary="map_sncRNA/{sample}.sncRNA.summary"
	shell:
            "hisat-3n -q -x ~/data/reference/Homo_sapiens/hisat_18S_28SrRNA_CT/GRCh38_rRNA --summary-file {output.summary} --new-summary -U {input.input1} -S {output.sncrna_aligned} --base-change C,T -p 2 --un sncRNA_depleted/{wildcards.sample}.fastq"

# SE end

# 
# #merge SE and PE
# rule merge:
#		 input:
#				 "map_sncRNA/{sample}.ribosomal.sam"
# 
# 
# 
# rule flag_sort_depth:
# 	input:
# 		"map_sncRNA/{sample}.ribosomal.sam"
# 	output:
# 		sorted_sam="map_sncRNA/sorted_sam/{sample}.sorted.sam",
# 		flagstat="map_sncRNA/flagstat/{sample}.rRNA.flagstat",
# 		depth="map_sncRNA/depth/{sample}.rRNA.depth"
# 	shell:
# 		"""
# 		samtools flagstat {input} > {output.flagstat}
# 		samtools sort -O sam -o {output.sorted_sam} -@ 2 {input}
# 		samtools depth -a -o {output.depth} -@ 2 {output.sorted_sam}
# 		""" 
#				 
# rule hisat_table_sncRNA:
#		 input:
#				 "map_sncRNA/sorted_sam/{sample}.sorted.sam"
#		 output:
#				 "map_sncRNA/table/{sample}_hisat3n_table.tsv"
#		 shell:
#				 "hisat-3n-table -p 16 --alignments {input} --ref /home/yliuchicago/data/reference/Homo_sapiens/hisat_18S_28SrRNA_CT/Homo_sapiens.GRCh38.28S.18S.fa --output-name {output} --base-change C,T"
# 
# rule calculate_methylation_rate:
#	 input:
#		 expand(
#			 "map_sncRNA/table/{sample}_hisat3n_table.tsv",
#			 sample=SAMPLE
#		 ),
#	 output:
#		 "map_sncRNA/calculated_rate/methylation_level_notfiltered.tsv",
#	 threads: 4
#	 resources:
#		 mem_mb=8000,
#	 shell:
#		 """
#		 (
#		 echo -e "Sample\\tChrom\\tPos\\tStrand\\tUnconverted\\tDepth\\tRatio"
#		 for file in {input}; do
#			 sample=`basename $file | cut -d. -f1`
#			 cat $file | awk -v samplename="$sample" 'BEGIN{{FS="\\t"; OFS="\\t"}} NR > 1 {{ print samplename,$1,$2,$3,$7,$7+$5,$7/($7+$5) }}'
#		 done
#		 ) > {output}
#		 """
# 
# 
# 
# 
# rule hisat2_3n_filtering:
#		 input:
#				 "map_sncRNA/sorted_sam/{sample}.sorted.sam"
#		 output:
#				 "map_sncRNA/sorted_sam/{sample}.sorted.filtered.sam"
#		 shell:
#				 """samtools view -@ 6 -e "[Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" {input} -O SAM -o {output}"""
# 
#				 
# rule hisat_table_filtered_sncRNA:
#		 input:
#				 "map_sncRNA/sorted_sam/{sample}.sorted.filtered.sam"
#		 output:
#				 "map_sncRNA/table/{sample}_filtered_hisat3n_table.tsv"
#		 shell:
#				 "hisat-3n-table -p 16 --alignments {input} --ref /home/yliuchicago/data/reference/Homo_sapiens/hisat_18S_28SrRNA_CT/Homo_sapiens.GRCh38.28S.18S.fa --output-name {output} --base-change C,T"
#				 
# 
# 
# rule calculate_methylation_rate_filtered:
#	 input:
#		 expand(
#			 "map_sncRNA/table/{sample}_filtered_hisat3n_table.tsv", sample=SAMPLE
#		 )
#	 output:
#		 "map_sncRNA/calculated_rate/methylation_level_filtered.tsv",
#	 threads: 4
#	 resources:
#		 mem_mb=8000,
#	 shell:
#		 """
#		 (
#		 echo -e "Sample\\tChrom\\tPos\\tStrand\\tUnconverted\\tDepth\\tRatio"
#		 for file in {input}; do
#			 sample=`basename $file | cut -d. -f1`
#			 cat $file | awk -v samplename="$sample" 'BEGIN{{FS="\\t"; OFS="\\t"}} NR > 1 {{ print samplename,$1,$2,$3,$7,$7+$5,$7/($7+$5) }}'
#		 done
#		 ) > {output}
#		 """
# 



















#rule map_genome:
#	input:
#		input1="sncRNA_depleted/{sample}_R1.fastq",
#		input2="sncRNA_depleted/{sample}_R2.fastq"
#	output:
#		genome_aligned="genome/{sample}.genome.sam",
#		summary="genome/{sample}.genome.summary",
#		fq1="sncRNA_and_Genome_depleted/{sample}_contamination_R1.fastq",
#		fq2="sncRNA_and_Genome_depleted/{sample}_contamination_R2.fastq"
#	shell:
#		"hisat-3n -q -x /home/yliuchicago/data/reference/Homo_sapiens/hisat_3n_CT/GRCh38_CT -1 {input.input1} -2 {input.input2}  --summary-file {output.summary} --new-summary --base-change C,T -p 2 -S {output.genome_aligned} --un-conc sncRNA_and_Genome_depleted/{wildcards.sample}_contamination_R%.fastq"
#
#
#rule flag_sort_depth_genome:
#	input:
#		"genome/{sample}.genome.sam"
#	output:
#		sorted_sam="genome/sorted_sam/{sample}.genome.sorted.sam",
#		flagstat="genome/flagstat/{sample}.genome.flagstat",
#		depth="genome/depth/{sample}.genome.depth"
#	shell:
#		"""
#		samtools flagstat {input} > {output.flagstat}
#		samtools sort -O sam -o {output.sorted_sam} -@ 2 {input}
#		samtools depth -a -o {output.depth} -@ 2 {output.sorted_sam}
#		""" 
#
#rule hisat_table:
#	input: 
#		"genome/sorted_sam/{sample}.genome.sorted.sam"
#	output:
#		"genome_table/{sample}_hisat3n_table.tsv"
#	shell:
#		"hisat-3n-table -p 16 --alignments {input} --ref /home/yliuchicago/data/reference/Homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa --output-name {output} --base-change C,T"
