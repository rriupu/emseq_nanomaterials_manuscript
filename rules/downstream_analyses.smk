rule methylKit_data_exploration:
    """
    """
    input:
        expand(rules.methylation_extraction.output.cytosine_report, sample = SAMPLES)
    output:
        sample_correlation_plot = os.path.join(config["output_dir"], "methylKit", "data_exploration", "sample_correlation.pdf"),
        pca_plot = os.path.join(config["output_dir"], "methylKit", "data_exploration", "PCA.pdf"),
        clustering_plot = os.path.join(config["output_dir"], "methylKit", "data_exploration", "clustering.pdf")
    params:
        rscript = os.path.join(config["scripts_dir"], "methylkit_data_exploration.R"),
        input_dir = os.path.join(config["output_dir"], "methylation_extraction"),
        output_dir = os.path.join(config["output_dir"], "methylKit", "data_exploration")
    log:
        os.path.join(config["logs_dir"], "{sample}", "{sample}_methylKit_data_exploration.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "{sample}", "{sample}_methylKit_data_exploration.txt")
    threads: 10
    conda:
        "../envs/methylkit.yaml"
    shell:
        """
        Rscript {params.rscript} \
            --input_dir {params.input_dir} \
            --output_dir {params.output_dir} \
            --threads {threads}
        """
    
rule methylKit_differential_analysis:
    """
    """
    input:
        expand(rules.methylation_extraction.output.cytosine_report, sample = SAMPLES)
    output:
        os.path.join()
    params:
        rscript = os.path.join(config["scripts_dir"], "methylkit_differential_analysis.R"),
        input_dir = os.path.join(config["output_dir"], "methylation_extraction"),
        output_dir = os.path.join(config["output_dir"], "methylKit", "data_exploration")
    log:
        os.path.join(config["logs_dir"], "{sample}", "{sample}_methylKit_data_exploration.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "{sample}", "{sample}_methylKit_data_exploration.txt")
    threads: 1
    conda:
        "../envs/methylkit.yaml"
    shell:
        """
        Rscript {params.rscript} \
            --input_dir {params.input_dir} \
            --output_dir {params.output_dir}
        """