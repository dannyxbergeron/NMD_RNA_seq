configfile: "config.json"

original_name = list(config['datasets_colombo'].values())
simple_id = list(config['datasets_colombo'].keys())
counts = ['est_counts', 'transcript_est_counts']

rule all:
    input:
        tpm = "results/kallisto/tpm.tsv",
        est_counts = "results/kallisto/est_counts.tsv",
        start_out = expand("logs/STAR/{id}.log", id=simple_id),
        transcript_tpm = "results/kallisto/transcript_tpm.tsv",
        transcript_est_counts = "results/kallisto/transcript_est_counts.tsv",
        # DESeq2_genes = "logs/DESeq2/genes.log",
        # DESeq2_transcripts = "logs/DESeq2/transcripts.log",
        rename = "logs/DESeq2/rename.tok",
        bw = expand("results/genomCov/bigwig/{id}.bw", id=simple_id)


rule rename_files:
    input:
        fastq = expand("data/reads/{original}_{pair}.fastq",
                    original=original_name, pair=[1, 2])
    output:
        new_name = expand("data/reads/{id}_{pair}.fastq",
                    id=simple_id, pair=[1, 2])
    run:
        for id, original in config['datasets_colombo'].items():
            for num in [1, 2]:
                old = "data/reads/{}_{}.fastq".format(original, num)
                new_ = "data/reads/{}_{}.fastq".format(id, num)

                print(old, new_)
                # os.rename(old, new_)


rule download_genome:
    """ Downloads the genome from Ensembl FTP servers """
    output:
        genome = config['path']['genome']
    params:
        link = config['download']['genome']
    shell:
        "wget --quiet -O {output.genome}.gz {params.link} && "
        "gzip -d {output.genome}.gz "


rule create_transcriptome:
    """ Uses gffread to generate a transcriptome """
    input:
        genome = config['path']['genome'],
        gtf = config['path']['annotation']
    output:
        seqs = config['path']['transcriptome']
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
        gtf = config['path']['annotation']
    output:
        map = config['path']['gene_name']
    conda:
        "envs/python.yaml"
    script:
        "scripts/generate_transcriptID_geneName.py"


rule trimming:
    """ Trims the FASTQ files using Trimmomatic """
    input:
        fq = "data/reads/{id}.fastq",
    output:
        fq = "data/trimmed/{id}.fastq.gz",
    params:
        options = [
            "ILLUMINACLIP:adapters.fa:2:30:10", "LEADING:5",
            "TRAILING:5", "MINLEN:45"
        ]
    log:
        "logs/trimmomatic/{id}.log"
    threads:
        32
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "-phred33 "
        "{input.fq} "
        "{output.fq} "
        "{params.options} "
        "&> {log}"


rule qc:
    """ Assess the FASTQ quality using FastQC """
    input:
        fq = rules.trimming.output.fq
    output:
        fq_out = "data/qc/{id}_fastqc.html"
    params:
        out_dir = "data/qc"
    log:
        "logs/fastqc/{id}.log"
    threads:
         32
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc "
        "--outdir {params.out_dir} "
        "--format fastq "
        "--threads {threads} "
        "{input.fq} "
        "&> {log}"


rule kallisto_index:
    """ Generates the transcriptome index for Kallisto """
    input:
        qc = expand("data/qc/{id}_fastqc.html",
                    id=simple_id),
        transcriptome = rules.create_transcriptome.output.seqs
    output:
        idx = config['path']['kallisto_index']
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
        fq = rules.trimming.output.fq,
    output:
        quant = "results/kallisto/{id}/abundance.tsv",
        h5 = "results/kallisto/{id}/abundance.h5",
    params:
        bootstrap = "50",
        outdir = "results/kallisto/{id}"
    log:
        "logs/kallisto/{id}.log"
    threads:
        32
    conda:
        "envs/kallisto.yaml"
    shell:
        "kallisto quant "
        "--bias "
        "--index={input.idx} "
        "--output-dir={params.outdir} "
        "--bootstrap-samples={params.bootstrap} "
        "--threads={threads} "
        "--single -l 200 -s 20 "
        "{input.fq} "
        "&> {log}"


rule combine_gene_quantification:
    """
    Custom Python script to collect and format Kallisto results for further
    processing.
    """
    input:
        datasets_colombo = expand(
            "results/kallisto/{id}/abundance.tsv",
            id=config['datasets_colombo'].keys()
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
        fasta = config["path"]["genome"],
        gtf = config["path"]['annotation']
    output:
        chrNameLength = "data/references/star_index/chrNameLength.txt"
    params:
        dir = config['path']['star_index']
    log:
        "logs/STAR/index.log"
    conda:
        "envs/star.yaml"
    threads:
        32
    shell:
        "mkdir -p {params.dir} && "
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 99"
        "&> {log}"


rule star_alignReads:
    """ Generates a bam file using STAR """
    input:
        idx = rules.star_index.output,
        fq = rules.trimming.output.fq,
    output:
        bam = "results/STAR/{id}/Aligned.sortedByCoord.out.bam"
    params:
        index = config['path']['star_index'],
        output_dir = "results/STAR/{id}/"
    log:
        "logs/STAR/{id}.log"
    threads:
        32
    conda:
        "envs/star.yaml"
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {params.index} "
        "--readFilesIn {input.fq} "
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

# include DESeq
# include: "rules/DESeq2.smk"

# include genomeCov -bg
include: "rules/genomCov.smk"
