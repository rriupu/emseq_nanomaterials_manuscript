rule methylKit:
    """
    """
    input:
        expand(rules.methylation_extraction.output.cytosine_report, sample = SAMPLES)
    output:
        sample_correlation_plot = os.path.join(config["output_dir"], "methylKit", "data_exploration", "sample_correlation.pdf")
    params:
        rscript = os.path.join("scripts", "methylkit.R"),
        input_dir = os.path.join(config["output_dir"], "methylation_extraction"),
        output_dir = os.path.join(config["output_dir"], "methylKit", "data_exploration")
    log:
        os.path.join(config["logs_dir"], "{sample}", "{sample}_methylKit_data_exploration.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "{sample}", "{sample}_methylKit_data_exploration.txt")
    threads: 10
    container:
        "docker://rriu/methylkit:1.33.3"
    shell:
        """
        Rscript {params.rscript} \
            --input_dir {params.input_dir} \
            --output_dir {params.output_dir} \
            --threads {threads}
        """