rule methylKit_data_exploration:
    """
    """
    input:
        expand(rules.methylation_extraction.output.cytosine_report, sample = SAMPLES)
    output:
        os.path.join()
    params:
        rscript = os.path.join(config["scripts_dir"], "methylkit_data_exploration.R")
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

        """
    