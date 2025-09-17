import os

rule get_genome:
    """
    Download the genome fasta sequence.
    """
    output:
        genome_fasta = os.path.join(
            config["data_dir"],
            "genome",
            "genome.fa.gz"
        )
    params:
        genome_url = config["genome_url"]
    threads: 1
    log:
        os.path.join(config["logs_dir"], "genome_download", "download.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "genome_download", "download.txt")
    shell:
        """
        wget -O {output.genome_fasta} {params.genome_url} 2> {log}
        """

rule genome_indexing:
    """
    Index the genome with bismark.
    """
    input:
        rules.get_genome.output.genome_fasta
    output:
        touch(os.path.join(config["data_dir"], "genome", "indexing_done"))
    threads: 10
    log:
        os.path.join(config["logs_dir"], "genome_indexing", "indexing.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "genome_indexing", "indexing.log")
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        bismark_genome_preparation \
            --verbose \
            --parallel {threads} \
            $(dirname {input}) \
        2> {log}
        """

rule bismark_mapping:
    """
    Map trimmed reads using bismark followed by deduplication.
    """
    input:
        genome_indexing_flag_file = rules.genome_indexing.output,
        genome_fasta = rules.get_genome.output.genome_fasta,
        trimmed_forward_read = rules.trimming.output.trimmed_forward_read,
        trimmed_reverse_read = rules.trimming.output.trimmed_reverse_read,
    output:
        aligned_deduplicated_reads = os.path.join(config["output_dir"], "aligned_reads", "{sample}", "{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.bam")
    params:
        unfiltered_bam = os.path.join(config["output_dir"], "aligned_reads", "{sample}", "{sample}_R1_001_val_1_bismark_bt2_pe.bam")
    threads: 10
    log:
        bismark_log = os.path.join(config["logs_dir"], "{sample}", "{sample}_mapping.log"),
        deduplication_log = os.path.join(config["logs_dir"], "{sample}", "{sample}_mark_nonconverted_reads.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "{sample}", "{sample}_mapping.log")
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        bismark \
            --genome $(dirname {input.genome_fasta}) \
            --parallel {threads} \
            --output_dir $(dirname {output.aligned_deduplicated_reads}) \
            -1 {input.trimmed_forward_read} \
            -2 {input.trimmed_reverse_read} \
        2> {log.bismark_log}

        deduplicate_bismark \
            -p \
            --output_dir $(dirname {output.aligned_deduplicated_reads}) \
            --bam \
            {params.unfiltered_bam} \
        2> {log.deduplication_log}

        rm {params.unfiltered_bam}
        """

