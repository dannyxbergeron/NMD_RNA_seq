configfile: "config.json"

simple_id = list(config['datasets_test'].keys()) #CHANGED
counts = ['est_counts', 'transcript_est_counts']

rule all:
    input:
        tpm = "results/kallisto/tpm.tsv",
        est_counts = "results/kallisto/est_counts.tsv",
        start_out = expand("logs/STAR/{id}.log", id=simple_id),
        transcript_tpm = "results/kallisto/transcript_tpm.tsv",
        transcript_est_counts = "results/kallisto/transcript_est_counts.tsv",
        bedgraph = expand("results/CoCo/{id}.bedgraph", id=simple_id),
        clean_bg = expand("results/CoCo/bigwig/{id}.bw", id=simple_id),
        DESeq2_genes = "logs/DESeq2/genes.log",
        DESeq2_transcripts = "logs/DESeq2/transcripts.log",
        rename = "logs/DESeq2/rename.tok"



# rule download_genome:
#     """ Downloads the genome from Ensembl FTP servers """
#     output:
#         genome = config['path_test']['genome']
#     params:
#         link = config['download']['genome']
#     shell:
#         "wget --quiet -O {output.genome}.gz {params.link} && "
#         "gzip -d {output.genome}.gz "


rule create_transcriptome:
    """ Uses gffread to generate a transcriptome """
    input:
        genome = config['path_test']['genome'],
        gtf = config['path_test']['annotation']
    output:
        seqs = config['path_test']['transcriptome']
    conda:
        "envs/gffread.yaml"
    shell:
        "gffread {input.gtf} -g {input.genome} -w {output.seqs}"


rule generate_transcriptID_geneName:
    """
    Generating a two-column text file containing the gene -> transcript
    relationship
    """
    input:
        gtf = config['path_test']['annotation']
    output:
        map = config['path_test']['gene_name']
    conda:
        "envs/python.yaml"
    script:
        "scripts/generate_transcriptID_geneName.py"


rule trimming:
    """ Trims the FASTQ files using Trimmomatic """
    input:
        fq1 = "data/test_reads/{id}_1.fastq",
        fq2 = "data/test_reads/{id}_2.fastq"
    output:
        fq1 = "data/trimmed/{id}_1.fastq.gz",
        fq2 = "data/trimmed/{id}_2.fastq.gz",
        unpaired_fq1 = "data/trimmed/{id}_1.unpaired.fastq.gz",
        unpaired_fq2 = "data/trimmed/{id}_2.unpaired.fastq.gz"
    params:
        options = [
            "ILLUMINACLIP:data/adapters.fa:2:30:10", "LEADING:5",
            "TRAILING:5", "MINLEN:45"
        ]
    log:
        "logs/trimmomatic/{id}.log"
    threads:
        1 #32 CHANGED
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "-phred33 "
        "{input.fq1} {input.fq2} "
        "{output.fq1} {output.unpaired_fq1}  "
        "{output.fq2} {output.unpaired_fq2} "
        "{params.options} "
        "&> {log}"


rule qc:
    """ Assess the FASTQ quality using FastQC """
    input:
        fq1 = rules.trimming.output.fq1,
        fq2 = rules.trimming.output.fq2,
        unpaired_fq1 = rules.trimming.output.unpaired_fq1,
        unpaired_fq2 = rules.trimming.output.unpaired_fq2,
    output:
        fq1_out = "data/qc/{id}_1_fastqc.html"
    params:
        out_dir = "data/qc"
    log:
        "logs/fastqc/{id}.log"
    threads:
         1 #32 CHANGED
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc "
        "--outdir {params.out_dir} "
        "--format fastq "
        "--threads {threads} "
        "{input.fq1} {input.fq2} "
        "{input.unpaired_fq1} {input.unpaired_fq2} "
        "&> {log}"


rule kallisto_index:
    """ Generates the transcriptome index for Kallisto """
    input:
        qc = expand("data/qc/{id}_1_fastqc.html",
                    id=simple_id),
        transcriptome = rules.create_transcriptome.output.seqs
    output:
        idx = config['path_test']['kallisto_index']
    params:
        kmer = "31"
    log:
        "logs/kallisto/index.log"
    conda:
        "envs/kallisto.yaml"
    shell:
        "kallisto index "
        "--index={output.idx} "
        "--kmer-size={params.kmer} "
        "{input.transcriptome} "
        "&> {log}"


rule kallisto_quant:
    """ Generates counts using Kallisto pseudo-alignment """
    input:
        idx = rules.kallisto_index.output.idx,
        fq1 = rules.trimming.output.fq1,
        fq2 = rules.trimming.output.fq2
    output:
        quant = "results/kallisto/{id}/abundance.tsv",
        h5 = "results/kallisto/{id}/abundance.h5",
    params:
        bootstrap = "50",
        outdir = "results/kallisto/{id}"
    log:
        "logs/kallisto/{id}.log"
    threads:
        1 #32 CHANGED
    conda:
        "envs/kallisto.yaml"
    shell:
        "kallisto quant "
        "--bias "
        "--index={input.idx} "
        "--output-dir={params.outdir} "
        "--bootstrap-samples={params.bootstrap} "
        "--threads={threads} "
        "{input.fq1} {input.fq2} "
        "&> {log}"


rule combine_gene_quantification:
    """
    Custom Python script to collect and format Kallisto results for further
    processing.
    """
    input:
        datasets = expand(
            "results/kallisto/{id}/abundance.tsv",
            id=config['datasets_test'].keys()
        ),
        map = rules.generate_transcriptID_geneName.output.map
    output:
        tpm = "results/kallisto/tpm.tsv",
        est_counts = "results/kallisto/est_counts.tsv",
        transcript_tpm = "results/kallisto/transcript_tpm.tsv",
        transcript_est_counts = "results/kallisto/transcript_est_counts.tsv"
    conda:
        "envs/python.yaml"
    script:
        "scripts/combine_gene_quantification.py"


rule star_index:
    """ Generates the genome index for STAR """
    input:
        fasta = config["path_test"]["genome"],
        gtf = config["path_test"]['annotation']
    output:
        chrNameLength = "data/test_reference/star_index/chrNameLength.txt"
    params:
        dir = config['path_test']['star_index']
    log:
        "logs/STAR/index.log"
    conda:
        "envs/star.yaml"
    threads:
        1 #8 CHANGED
    shell:
        "mkdir -p {params.dir} && "
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 99"
        "--genomeSAindexNbases 12 " #CHANGED to remove !!
        "&> {log}"


rule star_alignReads:
    """ Generates a bam file using STAR """
    input:
        idx = rules.star_index.output,
        fq1 = rules.trimming.output.fq1,
        fq2 = rules.trimming.output.fq2
    output:
        bam = "results/STAR/{id}/Aligned.sortedByCoord.out.bam"
    params:
        index = config['path_test']['star_index'],
        output_dir = "results/STAR/{id}/"
    log:
        "logs/STAR/{id}.log"
    threads:
        1 #32 CHANGED
    conda:
        "envs/star.yaml"
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {params.index} "
        "--readFilesIn {input.fq1} {input.fq2}  "
        "--runThreadN {threads} "
        "--readFilesCommand zcat "
        "--outReadsUnmapped Fastx "
        "--outFilterType BySJout "
        "--outStd Log "
        "--outSAMunmapped None "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.output_dir} "
        "--outFilterScoreMinOverLread 0.3 "
        "--outFilterMatchNminOverLread 0.3 "
        "--outFilterMultimapNmax 100 "
        "--winAnchorMultimapNmax 100 "
        "--alignEndsProtrude 5 ConcordantPair "
        "&> {log}"

# include coco
include: "rules/coco.smk"

# include DESeq
include: "rules/DESeq2.smk"
