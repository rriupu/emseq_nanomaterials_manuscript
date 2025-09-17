import os

configfile: "config.yaml"

SAMPLES = glob_wildcards(os.path.join(config["raw_reads_dir"], "{sample}_R1_001.fastq.gz")).sample

include: "rules/qc.smk"
include: "rules/trimming.smk"
include: "rules/qc_after_trimming.smk"
include: "rules/mapping.smk"
include: "rules/methylation_extraction.smk"
# include: "rules/downstream_analyses.smk"

rule all:
    input:
        rules.multiqc.output,
        rules.multiqc_after_trim.output,
        expand(rules.bismark_nucleotide_coverage_report.output.coverage_report, sample = SAMPLES)
        