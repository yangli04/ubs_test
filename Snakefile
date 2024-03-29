SAMPLE=["HL-A1", "HL-A2", "HL-A3", "HL-A4","HL-A5","HL-A6","HL-A7"]
# start from S1 as paired end sequencing. S2 means single-end sequencing. 
BATCH=["S1"]

wildcard_constraints:
    sample="HL-A[0-9]"

rule all:
    input:
        sorted_sam=expand("map_sncRNA/sorted_sam/{sample}.sorted.sam",sample=SAMPLE),
        rib_table=expand("map_sncRNA/table/{sample}_hisat3n_table.tsv",sample=SAMPLE),
        m1="map_sncRNA/calculated_rate/methylation_level_filtered.tsv",
        m2="map_sncRNA/calculated_rate/methylation_level_notfiltered.tsv",
        sncrna_aligned=expand("map_sncRNA/{sample}_{batch}.ribosomal.sam",sample=SAMPLE,batch=BATCH),
        summary=expand("map_sncRNA/{sample}_{batch}.sncRNA.summary",sample=SAMPLE,batch=BATCH),
        sorted_sam_genome=expand("genome_mapping/{sample}.genome.sorted.sam",sample=SAMPLE),
        flagstat=expand("genome_mapping/flagstat/{sample}.genome.flagstat",sample=SAMPLE),
        depth=expand("genome_mapping/depth/{sample}.genome.depth",sample=SAMPLE),
        hisat_table=expand("genome_mapping/genome_table/{sample}_hisat3n_table.tsv",sample=SAMPLE),
        methylation_table_nonfiltered="genome_mapping/calculated_rate/methylation_level_notfiltered_genome.tsv",
        methylation_table_filtered="genome_mapping/calculated_rate/methylation_level_filtered_genome.tsv",
# PE
rule trim_PE:
    input:
        input1="raw_data/{sample}_S1_R1_001.fastq.gz",
        input2="raw_data/{sample}_S1_R2_001.fastq.gz",
    output:
        output1="trimmed_PE/trimmed_{sample}_R1.fastq.gz",
        output2="trimmed_PE/trimmed_{sample}_R2.fastq.gz",
        tooshort1="trimmed_PE/tooshort/{sample}.tooshort.1.fq.gz",
        tooshort2="trimmed_PE/tooshort/{sample}.tooshort.2.fq.gz",
        log_file="trimmed_PE/{sample}.cutadapt.log",
    shell:
        """cutadapt -j 10 -m 20 -q 15 -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA;anywhere;" -A "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;anywhere;" -u 3 -U 3 -u -3 -U -3 --nextseq-trim=15 -o {output.output1} -p {output.output2} --too-short-output {output.tooshort1} --too-short-paired-output {output.tooshort2} {input.input1} {input.input2} > {output.log_file}"""


rule map_sncRNA_PE:
    input:
        input1="trimmed_PE/trimmed_{sample}_R1.fastq.gz",
        input2="trimmed_PE/trimmed_{sample}_R2.fastq.gz",
    output:
        sncrna_aligned=temp("map_sncRNA/{sample}_S1.ribosomal.sam"),
        fq1="sncRNA_depleted/{sample}_S1_R1.fastq.gz",
        fq2="sncRNA_depleted/{sample}_S1_R2.fastq.gz",
        summary="map_sncRNA/{sample}_S1.sncRNA.summary",
    params:
        un="sncRNA_depleted/{sample}_S1_R%.fastq.gz",
    shell:
        "hisat-3n -q -x ~/data/reference/Homo_sapiens/hisat_18S_28SrRNA_CT/GRCh38_rRNA --summary-file {output.summary} --new-summary -1 {input.input1} -2 {input.input2} -S {output.sncrna_aligned} --base-change C,T -p 2 --un-conc-gz {params.un}"
# PE end 


# SE
rule trim_SE:
    input:
        input1="raw_data/{sample}_S2_R1_001.fastq.gz",
    output:
        output1="trimmed_SE/trimmed_{sample}_R1.fastq.gz",
        tooshort1="trimmed_SE/tooshort/{sample}.tooshort.1.fq.gz",
        log_file="trimmed_SE/{sample}.cutadapt.log",
    shell:
        """cutadapt -j 10 -m 20 -q 15 -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA;anywhere;" -u 3 -u -3 --nextseq-trim=15 -o {output.output1} --too-short-output {output.tooshort1} {input.input1} > {output.log_file}"""

rule map_sncRNA_SE:
    input:
        input1="trimmed_SE/trimmed_{sample}_R1.fastq.gz",
    output:
        sncrna_aligned=temp("map_sncRNA/{sample}_S2.ribosomal.sam"),
        fq1="sncRNA_depleted/{sample}_S2_R1.fastq.gz",
        summary="map_sncRNA/{sample}_S2.sncRNA.summary",
    shell:
        "hisat-3n -q -x ~/data/reference/Homo_sapiens/hisat_18S_28SrRNA_CT/GRCh38_rRNA --summary-file {output.summary} --new-summary -U {input.input1} -S {output.sncrna_aligned} --base-change C,T -p 2 --un-gz {output.fq1}"

# SE end

#merge SE and PE

# use S2 for SE and S1 for PE, since when S1 and S2 exists at the same time, usually we think SE is the added sequence
rule sort_sam_SE:
    input:
        s2="map_sncRNA/{sample}_S2.ribosomal.sam",
    output:
        S2=temp("map_sncRNA/{sample}_S2_sorted.sam"),
    shell:
        """
        samtools sort -O sam -o {output.S2} -@ 2 {input.s2}
        """
rule sort_sam_PE:
    input:
        s1="map_sncRNA/{sample}_S1.ribosomal.sam",
    output:
        S1=temp("map_sncRNA/{sample}_S1_sorted.sam"),
    shell:
        """
        samtools sort -O sam -o {output.S1} -@ 2 {input.s1}
        """

# Here merge the PE reads to SE reads if they both exists, or just rename the sam file to {sample}.ribosomal.sam 
rule merge:
    input:
        "map_sncRNA/{sample}_S1_sorted.sam" if "S1" in BATCH else [],
        "map_sncRNA/{sample}_S2_sorted.sam" if "S2" in BATCH else [],
    output:
                "map_sncRNA/{sample}.ribosomal.sam",
    run:
        if len(BATCH) == 2:
            shell(
                """
            samtools merge -@ 4 -o {output} {input}
                """
            )
        else:
            shell(
                """
            mv {input} {output}
                """
            )
rule flag_sort_depth:
    input:
        "map_sncRNA/{sample}.ribosomal.sam",
    output:
        sorted_sam="map_sncRNA/sorted_sam/{sample}.sorted.sam",
        flagstat="map_sncRNA/flagstat/{sample}.rRNA.flagstat",
        depth="map_sncRNA/depth/{sample}.rRNA.depth",
    shell:
        """
        samtools flagstat {input} > {output.flagstat}
        samtools sort -O sam -o {output.sorted_sam} -@ 2 {input}
        samtools depth -a -o {output.depth} -@ 2 {output.sorted_sam}
        """
# (for i in $(ls ./); do cat $i | awk -v samplename="${i%.rRNA.depth}" 'BEGIN{{OFS="\t"}} {{ print samplename,$0 }}'; done ) > output.tsv


rule hisat_table_sncRNA:
    input:
        "map_sncRNA/sorted_sam/{sample}.sorted.sam",
    output:
        "map_sncRNA/table/{sample}_hisat3n_table.tsv",
    shell:
        """
        hisat-3n-table -p 16 --alignments {input} --ref /home/yliuchicago/data/reference/Homo_sapiens/hisat_18S_28SrRNA_CT/Homo_sapiens.GRCh38.28S.18S.fa --output-name {output} --base-change C,T
        """

rule calculate_methylation_rate:
    input:
        expand(
        "map_sncRNA/table/{sample}_hisat3n_table.tsv", sample=SAMPLE
         ),
    output:
        "map_sncRNA/calculated_rate/methylation_level_notfiltered.tsv",
    threads: 4
    resources:
        mem_mb=8000,
    shell:
        """
         (
         echo -e "Sample\\tChrom\\tPos\\tStrand\\tUnconverted\\tDepth\\tRatio"
         for file in {input}; do
             sample=`basename $file | cut -d. -f1`
             cat $file | awk -v samplename="$sample" 'BEGIN{{FS="\\t"; OFS="\\t"}} NR > 1 {{ print samplename,$1,$2,$3,$7,$7+$5,$7/($7+$5) }}'
         done
         ) > {output}
         """

rule hisat2_3n_filtering:
    input:
        "map_sncRNA/sorted_sam/{sample}.sorted.sam",
    output:
        "map_sncRNA/sorted_sam/{sample}.sorted.filtered.sam",
    shell:
        """samtools view -@ 6 -e "[Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" {input} -O SAM -o {output}"""


rule hisat_table_filtered_sncRNA:
    input:
        "map_sncRNA/sorted_sam/{sample}.sorted.filtered.sam",
    output:
        "map_sncRNA/table/{sample}_filtered_hisat3n_table.tsv",
    shell:
        "hisat-3n-table -p 16 --alignments {input} --ref /home/yliuchicago/data/reference/Homo_sapiens/hisat_18S_28SrRNA_CT/Homo_sapiens.GRCh38.28S.18S.fa --output-name {output} --base-change C,T"


rule calculate_methylation_rate_filtered:
    input:
         expand(
             "map_sncRNA/table/{sample}_filtered_hisat3n_table.tsv", sample=SAMPLE
         )
    output:
        "map_sncRNA/calculated_rate/methylation_level_filtered.tsv",
    threads: 4
    resources:
        mem_mb=8000,
    shell:
        """
         (
         echo -e "Sample\\tChrom\\tPos\\tStrand\\tUnconverted\\tDepth\\tRatio"
         for file in {input}; do
             sample=`basename $file | cut -d. -f1`
             cat $file | awk -v samplename="$sample" 'BEGIN{{FS="\\t"; OFS="\\t"}} NR > 1 {{ print samplename,$1,$2,$3,$7,$7+$5,$7/($7+$5) }}'
         done
         ) > {output}
         """

rule map_genome_PE:
    input:
        input1="sncRNA_depleted/{sample}_S1_R1.fastq.gz",
        input2="sncRNA_depleted/{sample}_S1_R2.fastq.gz",
    output:
        genome_aligned=temp("mapping_PE/{sample}_S1.genome.sam"),
        summary=temp("mapping_SE/{sample}_S1.genome.summary"),
        fq2="sncRNA_and_Genome_depleted/{sample}_S1_contamination_R2.fastq.gz",
        fq1="sncRNA_and_Genome_depleted/{sample}_S1_contamination_R1.fastq.gz",
    params:
        un="sncRNA_and_Genome_depleted/{sample}_S1_contamination_R%.fastq.gz",
    shell:
        "hisat-3n -q -x /home/yliuchicago/data/reference/Homo_sapiens/hisat_3n_CT/GRCh38_CT -1 {input.input1} -2 {input.input2}  --summary-file {output.summary} --new-summary --base-change C,T -p 2 -S {output.genome_aligned} --un-conc-gz {params.un}"


rule map_genome_SE:
    input:
        fq1="sncRNA_depleted/{sample}_S2_R1.fastq.gz",
    output:
        genome_aligned=temp("mapping_SE/{sample}_S2.genome.sam"),
        summary=temp("mapping_SE/{sample}_S2.genome.summary"),
        fq1="sncRNA_and_Genome_depleted/{sample}_S2_contamination_R1.fastq.gz",
    shell:
        "hisat-3n -q -x /home/yliuchicago/data/reference/Homo_sapiens/hisat_3n_CT/GRCh38_CT --summary-file {output.summary} --new-summary -U {input.fq1} -S {output.genome_aligned} --base-change C,T -p 2 --un-gz {output.fq1}"

rule sort_sam_genome_SE:
    input:
        s2="mapping_SE/{sample}_S2.genome.sam",
    output:
        S2=temp("mapping_SE/{sample}_S2.genome.sorted.sam"),
    shell:
        """
        samtools sort -O sam -o {output.S2} -@ 2 {input.s2}
        """
rule sort_sam_genome_PE:
    input:
        s1="mapping_PE/{sample}_S1.genome.sam",
    output:
        S1=temp("mapping_PE/{sample}_S1.genome.sorted.sam"),
    shell:
        """
        samtools sort -O sam -o {output.S1} -@ 2 {input.s1}
        """

rule merge_genome:
    input:
        "mapping_SE/{sample}_S2.genome.sorted.sam" if "S2" in BATCH else [],
        "mapping_PE/{sample}_S1.genome.sorted.sam" if "S1" in BATCH else [],
    output:
        temp("genome_mapping/{sample}.genome.sam"),
    run:
        if len(BATCH) == 2:
            shell(
                """
            samtools merge -@ 4 -o {output} {input}
                """
            )
        else:
            shell(
                """
            mv {input} {output}
                """
            )

rule flag_sort_depth_genome:
    input:
        "genome_mapping/{sample}.genome.sam",
    output:
        sorted_sam="genome_mapping/{sample}.genome.sorted.sam",
        flagstat="genome_mapping/flagstat/{sample}.genome.flagstat",
        depth="genome_mapping/depth/{sample}.genome.depth",
    shell:
        """
        samtools flagstat {input} > {output.flagstat}
        samtools sort -O sam -o {output.sorted_sam} -@ 2 {input}
        samtools depth -a -o {output.depth} -@ 2 {output.sorted_sam}
        """

rule hisat_table:
    input:
        "genome_mapping/{sample}.genome.sorted.sam",
    output:
        "genome_mapping/genome_table/{sample}_hisat3n_table.tsv",
    shell:
        "hisat-3n-table -p 16 --alignments {input} --ref /home/yliuchicago/data/reference/Homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa --output-name {output} --base-change C,T"


rule calculate_methylation_rate_genome:
    input:
        expand(
        "genome_mapping/genome_table/{sample}_hisat3n_table.tsv", sample=SAMPLE
         ),
    output:
        "genome_mapping/calculated_rate/methylation_level_notfiltered_genome.tsv",
    threads: 4
    resources:
        mem_mb=8000,
    shell:
        """
         (
         echo -e "Sample\\tChrom\\tPos\\tStrand\\tUnconverted\\tDepth\\tRatio"
         for file in {input}; do
             sample=`basename $file | cut -d. -f1`
             cat $file | awk -v samplename="$sample" 'BEGIN{{FS="\\t"; OFS="\\t"}} NR > 1 {{ print samplename,$1,$2,$3,$7,$7+$5,$7/($7+$5) }}'
         done
         ) > {output}
         """

rule hisat2_3n_filtering_genome:
    input:
        "genome_mapping/{sample}.genome.sorted.sam",
    output:
        "genome_mapping/{sample}.genome.sorted.filtered.sam",
    shell:
        """samtools view -@ 6 -e "[Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" {input} -O SAM -o {output}"""


rule hisat_table_filtered_genome:
    input:
        "genome_mapping/{sample}.genome.sorted.filtered.sam",
    output:
        "genome_mapping/genome_table/{sample}_filtered_hisat3n_table.tsv",
    shell:
        "hisat-3n-table -p 16 --alignments {input} --ref /home/yliuchicago/data/reference/Homo_sapiens/hisat_18S_28SrRNA_CT/Homo_sapiens.GRCh38.28S.18S.fa --output-name {output} --base-change C,T"


rule calculate_methylation_rate_filtered_genome:
    input:
         expand(
        "genome_mapping/genome_table/{sample}_filtered_hisat3n_table.tsv", sample=SAMPLE
         ),
    output:
        "genome_mapping/calculated_rate/methylation_level_filtered_genome.tsv",
    threads: 4
    resources:
        mem_mb=8000,
    shell:
        """
         (
         echo -e "Sample\\tChrom\\tPos\\tStrand\\tUnconverted\\tDepth\\tRatio"
         for file in {input}; do
             sample=`basename $file | cut -d. -f1`
             cat $file | awk -v samplename="$sample" 'BEGIN{{FS="\\t"; OFS="\\t"}} NR > 1 {{ print samplename,$1,$2,$3,$7,$7+$5,$7/($7+$5) }}'
         done
         ) > {output}
         """

