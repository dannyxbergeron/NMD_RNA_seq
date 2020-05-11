rule genomCov:
    input:
        bam = rules.star_alignReads.output.bam
    output:
        bedgraph = "results/genomCov/{id}.bedgraph"
    conda:
        "../envs/coco.yaml"
    shell:
         'TmpScale=$(bc <<< "scale=6;1000000/$(samtools view -f 0 -c {input.bam})") && '
         'bedtools genomecov -bg -ibam {input.bam} -scale $TmpScale '
         '| sort -k1,1 -k2,2n > {output.bedgraph}'


rule keep_primary_chr:
    input:
        bedgraphs = expand("results/genomCov/{id}.bedgraph", id= simple_id),
        chrLength = "data/references/star_index/chrNameLength.txt",
    output:
        clean_bg = expand("results/genomCov/clean_{id}.bedgraph", id= simple_id),
        new_chrLength = "data/references/star_index/chrNameLength_modif.txt"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/clean_bg.py"

rule convert_bw:
    input:
        clean_bg = "results/genomCov/clean_{id}.bedgraph",
        new_chrLength = "data/references/star_index/chrNameLength_modif.txt"
    output:
        bw = "results/genomCov/bigwig/{id}.bw"
    conda:
        "../envs/bedgraphtobigwig.yaml"
    shell:
        "bedGraphToBigWig {input.clean_bg} {input.new_chrLength} {output.bw} > convert_bw.log"
