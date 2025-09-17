rule fastqc_after_trim:
    """
    Run FastQC on all samples.
    """
    input:
        trimmed_reads = rules.trimming.output.trimmed_forward_read
    output:
        html = os.path.join(config["output_dir"], "qc_reports", "after_trimming", "{sample}", "{sample}_R1_001_val_1_fastqc.html"),
        zipfile = os.path.join(config["output_dir"], "qc_reports", "after_trimming", "{sample}", "{sample}_R1_001_val_1_fastqc.zip")
    params:
        wd = os.path.join(config["output_dir"], "qc_reports", "after_trimming", "{sample}/")
    threads: 4
    log:
        os.path.join(config["logs_dir"], "{sample}", "{sample}_fastqc_after_trimming.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "{sample}", "{sample}_fastqc_after_trimming.txt")
    container:
        "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
    shell:
        """
        mkdir -p {params.wd} 2>> {log}
        
        fastqc \
            --threads {threads} \
            --outdir {params.wd} \
            --dir {params.wd} \
            {input.trimmed_reads} \
        &>> {log}
        """

rule multiqc_after_trim:
    """
    Run mulqiQC on fastqc reports.
    """
    input:
        expand(rules.fastqc_after_trim.output.html, sample = SAMPLES)
    output:
        report = os.path.join(config["output_dir"], "qc_reports", "after_trimming", "multiqc_report.html")
    params:
        wd = os.path.join(config["output_dir"], "qc_reports", "after_trimming/")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "multiqc", "after_trimming.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "multiqc", "after_trimming.txt")
    container:
        "docker://multiqc/multiqc:v1.29"
    shell:
        """
        cd {params.wd}

        multiqc .
        """