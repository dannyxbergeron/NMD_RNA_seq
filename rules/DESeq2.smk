rule DESeq2_genes:
    """ Differential expression for the different conditions """
    input:
        counts = 'results/kallisto/{counts}.tsv',
        samples = "data/design.tsv"
    output:
        results = directory("results/DESeq2/{counts}"),
    params:
        count_type = #TODO
    log:
        "logs/DESeq2/{counts}.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/DESeq2.R"

# rule DESeq2_transcripts:
#     """ Differential expression for the different conditions """
#     input:
#         counts = 'results/kallisto/{counts}.tsv',
#         samples = "data/design.tsv"
#     output:
#         results = directory("results/DESeq2/{counts}"),
#     log:
#         "logs/DESeq2/{counts}.log"
#     conda:
#         "../envs/R.yaml"
#     script:
#         "../scripts/DESeq2.R"
