############ I don't nees that for now #########################
# rule create_coco_annotation:
#     input:
#         gtf = config["path"]['annotation']
#     output:
#         "coco_annotation.tok"
#     conda:
#         "../envs/coco.yaml"
#     shell:
#         "coco correct_annotation {input.gtf} && touch {output}"

rule coco_cb:
    input:
        bam = rules.star_alignReads.output.bam,
        chrLength = "data/references/star_index/chrNameLength.txt"
    output:
        bedgraph = "results/CoCo/{id}.bedgraph"
    conda:
        "../envs/coco.yaml"
    threads:
        32
    shell:
        "coco cb -u -t {threads} -c 2500000 {input.bam} "
        "{output.bedgraph} {input.chrLength}"


rule keep_primary_chr:
    input:
        bedgraphs = expand("results/CoCo/{id}.bedgraph", id= simple_id),
        chrLength = "data/references/star_index/chrNameLength.txt",
    output:
        clean_bg = expand("results/CoCo/clean_{id}.bedgraph", id= simple_id),
        new_chrLength = "data/references/star_index/chrNameLength_modif.txt"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/clean_bg.py"

rule sort_bg:
    input:
        bg = "results/CoCo/clean_{id}.bedgraph"
    output:
        bg = "results/CoCo/sorted_clean_{id}.bedgraph"
    shell:
        "sort -k1,1 -k2,2n {input.bg} | sed 's/chrM/chrMT/g' > {output.bg}"

rule convert_bw:
    input:
        clean_bg = "results/CoCo/sorted_clean_{id}.bedgraph",
        new_chrLength = "data/references/star_index/chrNameLength_modif.txt"
    output:
        bw = "results/CoCo/bigwig/{id}.bw"
    conda:
        "../envs/bedgraphtobigwig.yaml"
    shell:
        "bedGraphToBigWig {input.clean_bg} {input.new_chrLength} {output.bw}"
