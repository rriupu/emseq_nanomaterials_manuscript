rule trimming:
    """
    Adapter trimming and quality filtering. Outputs compressed interleaved and 
    pre-processed reads.
    """
    input:
        forward_read = os.path.join(config["raw_reads_dir"], "{sample}_R1_001.fastq.gz"),
        reverse_read = os.path.join(config["raw_reads_dir"], "{sample}_R2_001.fastq.gz")
    output:
        trimmed_forward_read = os.path.join(config["output_dir"], "trimmed_reads", "{sample}", "{sample}_R1_001_val_1.fq.gz"),
        trimmed_reverse_read = os.path.join(config["output_dir"], "trimmed_reads", "{sample}", "{sample}_R2_001_val_2.fq.gz")
    threads: 8
    log:
        os.path.join(config["logs_dir"], "{sample}", "{sample}_trimming.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "{sample}", "{sample}_trimming.txt")
    conda:
        "../envs/trimming.yaml"
    shell:
        """
        trim_galore \
            --nextseq 20 \
            --clip_R1 10 \
            --clip_R2 10 \
            --three_prime_clip_R1 10 \
            --three_prime_clip_R2 10 \
            --cores {threads} \
            --paired \
            --output_dir $(dirname {output.trimmed_forward_read}) \
            {input.forward_read} \
            {input.reverse_read} \
        2> {log}
        """