# rule merge_and_mark_duplicates:
#     """
#     Mark all duplicate reads with samblaster.
#     """
#     input:
#         aligned_reads = rules.bwameth_mapping.output.aligned_reads
#     output:
#         merged_and_marked_aligned_reads = os.path.join("results", "aligned_reads", "{sample}", "aligned_reads_marked_duplicates.bam")
#     params:

#     threads: 8
#     log:
#         "logs/{sample}/{sample}_merge_and_mark_duplicates.log"
#     benchmark:
#         "benchmarks/{sample}/{sample}_merge_and_mark_duplicates.txt"
#     resources:
#         mem_gb = 8
#     conda:
#         "../envs/env.yaml"
#     shell:
#         """
#         # tmpdir=$(mktemp -d)

#         # samtools cat -b {input.aligned_reads} \
#         samtools view -h {input.aligned_reads} \
#         | samblaster 2> {log} \
#         | sambamba view \
#             -t {threads} \
#             -l 0 \
#             -S \
#             -f bam \
#             /dev/stdin \
#         | sambamba sort \
#             # --tmpdir=${{tmpdir}} \
#             -t {threads} \
#             -m {resources.mem_gb}GB \
#             -o {output.merged_and_marked_aligned_reads} \
#             /dev/stdin
#         """

# rule sum_nonconverted_reads:
#     """
#     """
#     input:
#         nonconverted_reads = rules.bwameth_mapping.log.mark_nonconverted_reads_log
#     output:
#         nonconverted_counts = os.path.join("results", "nonconverted_counts", "{sample}", "nonconverted_counts.tsv")
#     shell:
#         """
#         files=(*.tsv)	
#         paste *.tsv | awk -v numFiles=${#files[@]} -v OFS='\t' '	
#         {	
#         row = sep = ""	
#         for(i=1; i < NF/numFiles; ++i) { row = row sep $i; sep = OFS }	
#         sum = $(NF/numFiles) # last header col. / (1st) data col. to sum	
#         for(i=2; i<=numFiles; ++i) sum += $(NF/numFiles * i) # add other cols.	
#         printf "%s%s%s\\n", row, OFS, sum	
#         }' > tmp-counts.tsv	
#         awk '{print "!{library}\t" $0}' tmp-counts.tsv > !{library}-nonconverted-counts.tsv	
#         """

# rule combine_nonconversion:
#     """
#     """
#     shell:
#         """
#         cat *.tsv > combined-nonconverted.tsv
#         """

# rule select_human_reads:
#     """
#     Filter reads only from standard chromosomes.
#     """
#     input:
#         md_file = rules.merge_and_mark_duplicates.output.merged_and_marked_aligned_reads
#     output:
#         human_reads = os.path.join("results", "human_reads", "{sample}", "human_reads.bam")
#     threads: 8
#     log:
#         "logs/{sample}/{sample}_select_human_reads.log"
#     benchmark:
#         "benchmarks/{sample}/{sample}_select_human_reads.txt"
#     conda:
#         "../envs/env.yaml"
#     shell:
#         """
#         sambamba view \
#             -t {threads} \
#             -l 0 \
#             -f bam \
#             {input.md_file} \
#             chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
#         > {output.human_reads}
#         """

# rule human_insert_size:
#     """
#     Compute the insert size on the human-filtered reads.
#     """
#     input:
#         human_reads = rules.select_human_reads.output.human_reads
#     output:
#         human_insert_sizes = os.path.join("results", "human_insert_size", "{sample}", "{sample}.insertsize_metrics"),
#         insert_size_histogram = os.path.join("results", "human_insert_size", "{sample}", "{sample}.insertsize_metrics.pdf")
#     threads: 4
#     log:
#         "logs/{sample}/{sample}_human_insertsize.log"
#     benchmark:
#         "benchmarks/{sample}/{sample}_human_insertsize.txt"
#     conda:
#         "../envs/env.yaml"
#     shell:
#         """
#         picard -Xmx16g \
#             CollectInsertSizeMetrics \
#             VALIDATION_STRINGENCY=LENIENT \
#             I={input.human_reads} \
#             O={output.human_insert_sizes} \
#             MINIMUM_PCT=0.0001 \
#             HISTOGRAM_FILE={output.insert_size_histogram}
#         """

# rule human_gc_bias:
#     """
#     Compute GC bias on human-filtered reads.
#     """
#     input:
#         human_reads = rules.select_human_reads.output.human_reads,
#         genome_fasta = rules.get_genome.output.genome_fasta
#     output:
#         human_gc_bias = os.path.join("results", "human_gc_bias", "{sample}", "{sample}.gc_metrics"),
#         summary_metrics = os.path.join("results", "human_gc_bias", "{sample}", "{sample}.gc_summary_metrics"),
#         plot = os.path.join("results", "human_gc_bias", "{sample}", "{sample}.gc.pdf")
#     threads: 1
#     log:
#         "logs/{sample}/{sample}_human_gc_bias.log"
#     benchmark:
#         "benchmarks/{sample}/{sample}_human_gc_bias.txt"
#     conda:
#         "../envs/env.yaml"
#     shell:
#         """
#         picard -Xmx4g \
#             CollectGcBiasMetrics \
#             IS_BISULFITE_SEQUENCED=true \
#             VALIDATION_STRINGENCY=LENIENT \
#             I={input.human_reads} \
#             O={output.human_gc_bias} \
#             S={output.summary_metrics} \
#             CHART={output.plot} \
#             R={input.genome_fasta}
#         """

# rule methylDackel_mbias:
#     """
#     Compute methylation bias.
#     """
#     input:
#         md_file = rules.merge_and_mark_duplicates.output.merged_and_marked_aligned_reads,
#         genome_fasta = rules.get_genome.output.genome_fasta
#     output:
#         combined_mbias = os.path.join("results", "methylDackel", "mbias", "{sample}_combined_mbias.tsv"),
#         chn_inclusion_options = os.path.join("results", "methylDackel", "mbias", "{sample}_chn_inclusion_options.tsv"),
#         cpg_inclusion_options = os.path.join("results", "methylDackel", "mbias", "{sample}_cpg_inclusion_options.tsv")
#     params:
#         chn_file_prefix = os.path.join("results", "methylDackel", "mbias", "{sample}_chn"),
#         cpg_file_prefix = os.path.join("results", "methylDackel", "mbias", "{sample}_cpg")
#     threads: 8
#     log:
#         "logs/{sample}/{sample}_methylDackel_mbias.log"
#     benchmark:
#         "benchmarks/{sample}/{sample}_methylDackel_mbias.txt"
#     conda:
#         "../envs/env.yaml"
#     shell:
#         """
#         echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)" > {output.combined_mbias}
#         chrs=(`samtools view -H {input.md_file} | grep @SQ | cut -f 2 | sed 's/SN://'| grep -v _random | grep -v chrUn | sed 's/|/\\|/'`)

#         for chr in ${{chrs[*]}}; do
#             for context in CHH CHG CpG; do
#                 arg=''
#                 if [ $context = 'CHH' ]; then
#                     arg='--CHH --noCpG'
#                 elif [ $context = 'CHG' ]; then
#                     arg='--CHG --noCpG'
#                 fi
#                 # need two calls to add columns containing the counts without filtering duplicate reads (for rrEM-seq where start/end is constrained)
#                 join -t $'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
#                 <( \
#                     MethylDackel mbias --noSVG $arg -@ {threads} -r $chr {inupt.genome_fasta} {input.md_file} | \
#                     tail -n +2 | awk '{{print $1"-"$2"-"$3"\t"$0}}' | sort -k 1b,1 \
#                 ) \
#                 <( \
#                     MethylDackel mbias --noSVG --keepDupes -F 2816 $arg -@ {threads} -r $chr {input.genome_fasta} {input.md_file} | \
#                     tail -n +2 | awk '{{print $1"-"$2"-"$3"\t"$0}}' | sort -k 1b,1 \
#                 ) \
#                 | sed "s/^/${{chr}}\t${{context}}\t/" \
#                 >> {output.combined_mbias}
#             done
#         done
#         # makes the svg files for trimming checks
#         MethylDackel mbias \
#             -@ {threads} \
#             --noCpG \
#             --CHH \
#             --CHG \
#             -r ${{chrs[0]}} \
#             {input.genome_fasta} \
#             {input.md_file} \
#             {params.chn_file_prefix} \
#         2> {output.chn_inclusion_options}
        
#         for f in *chn*.svg; do 
#             sed -i "s/Strand<\\/text>/Strand $f ${{chrs[0]}} CHN <\\/text>/" $f
#         done;

#         MethylDackel mbias \
#             -@ {threads} \
#             -r ${{chrs[0]}} \
#             {input.genome_fasta} \
#             {input.md_file} \
#             {params.cpg_file_prefix} \
#         2> {output.cpg_file_prefix}

#         for f in *cpg*.svg; do 
#             sed -i "s/Strand<\\/text>/Strand $f ${{chrs[0]}} CpG<\\/text>/" $f
#         done;
#         """

# rule methylDackel_extract:
#     """
#     Compute per-base CpG, CHH, and CHG metrics.
#     """
#     input:
#         md_file = rules.merge_and_mark_duplicates.output.merged_and_marked_aligned_reads,
#         genome_fasta = rules.get_genome.output.genome_fasta,
#         mbias_chn_inclusion_options = rules.methylDackel_mbias.output.chn_inclusion_options,
#         mbias_cpg_inclusion_options = rules.methylDackel_mbias.output.cpg_inclusion_options
#     output:
#         chg_file = os.path.join("results", "methylDackel", "extract", "{sample}_CHG.methylKit.gz"),
#         chh_file = os.path.join("results", "methylDackel", "extract", "{sample}_CHH.methylKit.gz"),
#         cpg_file = os.path.join("results", "methylDackel", "extract", "{sample}_CpG.methylKit.gz")
#     params:
#         methyldackel_extract_prefix = os.path.join("results", "methylDackel", "{sample}")
#     threads: 8
#     log:
#         "logs/{sample}/{sample}_methylDackel_extract.log"
#     benchmark:
#         "benchmarks/{sample}/{sample}_methylDackel_extract.txt"
#     conda:
#         "../envs/env.yaml"
#     shell:
#         """
#         chn_inclusion_options=$(cat {input.mbias_chn_inclusion_options} | sed 's/Suggested inclusion options: //')
#         cpg_inclusion_options-$(cat {input.mbias_cpg_inclusion_options} | sed 's/Suggested inclusion options: //')

#         MethylDackel extract \
#             --methylKit \
#             $chn_inclusion_options \
#             -@ {threads} \
#             --CHH --CHG --noCpG \
#             -o {params.methyldackel_extract_prefix} \
#             {input.genome_fasta} \
#             {input.md_file};

#         MethylDackel extract \
#             --methylKit \
#             $cpg_inclusion_options \
#             -@ {threads} \
#             -o {params.methyldackel_extract_prefix} \
#             {input.genome_fasta} \
#             {input.md_file};
        
#         pigz -p {threads} *.methylKit
#         """

rule picard_gc_bias:
    """
    Compute GC bias on all aligned reads.
    """
    input:
        md_file = rules.merge_and_mark_duplicates.output.merged_and_marked_aligned_reads,
        genome_fasta = rules.get_genome.output.genome_fasta
    output:
        picard_gc_bias = os.path.join("results", "picard_gc_bias", "{sample}", "{sample}.gc_metrics"),
        summary_metrics = os.path.join("results", "picard_gc_bias", "{sample}", "{sample}.gc_summary_metrics"),
        plot = os.path.join("results", "picard_gc_bias", "{sample}", "{sample}.gc.pdf")
    threads: 1
    log:
        "logs/{sample}/{sample}_picard_gc_bias.log"
    benchmark:
        "benchmarks/{sample}/{sample}_picard_gc_bias.txt"
    conda:
        "../envs/env.yaml"
    shell:
        """
        picard \
            -Xmx4g \
            CollectGcBiasMetrics \
            IS_BISULFITE_SEQUENCED=true \
            VALIDATION_STRINGENCY=LENIENT \
            I={input.md_file} \
            O={output.picard_gc_bias} \
            S={output.summary_metrics} \
            CHART={output.plot} \
            R={input.genome_fasta}
        """

rule picard_stats:
    """
    Compute insert size metrics on all aligned reads.
    """
    input:
        md_file = rules.merge_and_mark_duplicates.output.merged_and_marked_aligned_reads
    output:
        insertsize_metrics = os.path.join("results", "picard_stats", "{sample}", "{sample}.insertsize_metrics"),
        insert_size_histogram = os.path.join("results", "picard_stats", "{sample}", "{sample}.insertsize_metrics.pdf")
    threads: 4
    log:
        "logs/{sample}/{sample}_picard_stats.log"
    benchmark:
        "benchmarks/{sample}/{sample}_picard_stats.txt"
    conda:
        "../envs/env.yaml"
    shell:
        """
        picard \
            -Xmx16g \
            CollectInsertSizeMetrics \
            VALIDATION_STRINGENCY=LENIENT  \
            I={input.md_file} \
            O={output.insertsize_metrics} \
            MINIMUM_PCT=0 \
            HISTOGRAM_FILE={output.insert_size_histogram}
        """

rule samtools_stats:
    """
    Collect statistics from BAM files and output to text format.
    """
    input:
        md_file = rules.merge_and_mark_duplicates.output.merged_and_marked_aligned_reads
    output:
        samstat_file = os.path.join("results", "samstats", "{sample}", "{sample}.samstat")
    threads: 2
    log:
    benchmark:
    conda:
        "../envs/env.yaml"
    shell:
        """
        samtools stats -@{threads} {input.md_file} > {output.samstat_file}
        """

rule samtools_flagstats:
    """
    Compute flag and index stats.
    """
    input:
        md_file = rules.merge_and_mark_duplicates.output.merged_and_marked_aligned_reads
    output:
        flagstat_file = os.path.join("results", "samtools_flagstats", "{sample}", "{sample}.flagstat"),
        idxstat_file = os.path.join("results", "samtools_flagstats", "{sample}", "{sample}.idxstat")
    threads: 2
    shell:
        """
        samtools flagstat -@{threads} {input.md_file} > {output.flagstat_file}
        samtools idxstats {input.md_file} > {output.idxstat_file}
        """
    

rule goleft_indexcov:
    """
    Compute genome coverage.
    """
    input:
        md_file = rules.merge_and_mark_duplicates.output.merged_and_marked_aligned_reads
    output:
        goleft_ped = os.path.join("results", "goleft_indexcov", "{sample}", "indexcov.ped"),
        goleft_roc = os.path.join("results", "goleft_indexcov", "{sample}", "indexcov.roc")
    params:
        directory = os.path.join("results", "goleft_indexcov", "{sample}")
    threads: 1
    log:
        "logs/goleft_indexcov/goleft_indexcov.log"
    benchmark:
        "benchmarks/goleft_indexcov/goleft_indexcov.txt"
    conda:
        "../envs/env.yaml"
    shell:
        """
        goleft indexcov --directory {params.directory} {input.md_file}
        """