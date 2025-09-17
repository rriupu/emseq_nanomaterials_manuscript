import pandas as pd

rule methylation_extraction:
    """
    Use bismark to extract methylation calls.
    """
    input:
        genome_fasta = rules.get_genome.output.genome_fasta,
        genome_indexing_flag_file = rules.genome_indexing.output,
        aligned_deduplicated_reads = rules.bismark_mapping.output.aligned_deduplicated_reads
    output:
        bedgraph = os.path.join("results", "methylation_extraction", "{sample}", "{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"),
        cytosine_report = os.path.join("results", "methylation_extraction", "{sample}", "{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.gz"),
        CHH_OT_methylation = os.path.join("results", "methylation_extraction", "{sample}", "CHH_OT_{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.txt.gz"),
        CHH_OB_methylation = os.path.join("results", "methylation_extraction", "{sample}", "CHH_OB_{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.txt.gz"),
        CHG_OT_methylation = os.path.join("results", "methylation_extraction", "{sample}", "CHG_OT_{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.txt.gz"),
        CHG_OB_methylation = os.path.join("results", "methylation_extraction", "{sample}", "CHG_OB_{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.txt.gz"),
        CpG_OT_methylation = os.path.join("results", "methylation_extraction", "{sample}", "CpG_OT_{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.txt.gz"),
        CpG_OB_methylation = os.path.join("results", "methylation_extraction", "{sample}", "CpG_OB_{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.txt.gz"),
        coverage_report = os.path.join("results", "methylation_extraction", "{sample}", "{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")
    log:
        os.path.join("logs", "{sample}", "{sample}_methylation_extraction.log")
    benchmark:
        os.path.join("benchmarks", "{sample}", "{sample}_methylation_extraction.txt")
    threads: 10
    resources:
        mem_gb = 20
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        bismark_methylation_extractor \
            -p \
            --gzip \
            --bedGraph \
            --buffer_size {resources.mem_gb}G \
            --cytosine_report \
            --parallel {threads} \
            --genome_folder $(realpath $(dirname {input.genome_fasta})) \
            --output_dir $(dirname {output.bedgraph}) \
            {input.aligned_deduplicated_reads} \
        2> {log}
        """

rule bismark_nucleotide_coverage_report:
    """
    
    """
    input:
        aligned_reads = rules.mapping.output.aligned_deduplicated_reads,
        genome_fasta = rules.get_genome.output.genome_fasta
    output:
        coverage_report = os.path.join(config["output_dir"], "aligned_reads", "{sample}", "{sample}.nucleotide_stats.txt")
    log:
        os.path.join("logs", "{sample}", "{sample}_bismark_nucleotide_coverage_report.log")
    benchmark:
        os.path.join("benchmarks", "{sample}", "{sample}_bismark_nucleotide_coverage_report.txt")
    threads: 1
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        bam2nuc \
            --dir $(dirname {input.aligned_reads}) \
            --genome_folder $(realpath $(dirname {input.genome_fasta})) \
            {input.aligned_reads} \
        2> {log}
        """

# rule bismark_processing_report:
#     """
#     """
#     input:
#         nucleotide_report = rules.bismark_nucleotide_coverage_report.output.coverage_report
#     output:
#         processing_report = os.path.join(config["output_dir"], "bismark_reports", "{sample}", "{sample}_processing_report.")
#     params:
#         alignment_report = os.path.join(config["output_dir"], "aligned_reads", "{sample}", "{sample}_R1_001_val_1_bismark_bt2_PE_report.txt"),
#         dedup_report = os.path.join(config["output_dir"], "aligned_reads", "{sample}", "{sample}_R1_001_val_1_bismark_bt2_pe.deduplication_report.txt"),
#         splitting_report = os.path.join(config["output_dir"], "methylation_extraction", "{sample}", "{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated_splitting_report.txt"),
#         mbias_report = os.path.join(config["output_dir"], "methylation_extraction", "{sample}", "{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.M-bias.txt")
#     log:
#         os.path.join(config["logs_dir"], "{sample}", "{sample}_bismark_processing_report.log")
#     benchmark:
#         os.path.join(config["benchmarks_dir"], "{sample}", "{sample}_bismark_processing_report.txt")
#     threads: 1 
#     conda:
#         "../envs/mapping.yaml"
#     shell:
#         """
#         bismark2report \
#             --output $(basename {output.processing_report}) \
#             --dir $(dirname {input.coverage_report}) \
#             --alignment_report {params.alignment_report} \
#             --dedup_report {params.dedup_report} \
#             --splitting_report {params.splitting_report} \
#             --mbias_report {params.mbias_report} \
#             --nucleotide_report {input.nucleotide_report} \
#         2> {log}
#         """

# rule bismark_summary_report:
#     """
#     """
#     input:
#         expand(rules.bismark_processing_report.output.processing_report, sample = SAMPLES)
#     output:
#         summary_report = os.path.join(config["output_dir"], "bismark_reports", "bismark_summary_report.html")
#     log:
#         os.path.join(config["logs_dir"], "bismark_summary_report", "bismark_summary_report.log")
#     benchmark:
#         os.path.join(config["benchmarks_dir"], "bismark_summary_report", "bismark_summary_report.txt")
#     threads: 1
#     conda:
#         "../envs/mapping.yaml"
#     shell:
#         """
#         bismark2summary \
#             -o {output.summary_report}
#         """