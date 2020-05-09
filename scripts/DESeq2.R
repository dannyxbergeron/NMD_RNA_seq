log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

#Loading count matrix from file, converting to integer
counts <- read.table(
	snakemake@input[["counts"]], header=TRUE,
	row.names=snakemake@params[["count_type"]], check.names=FALSE
)

counts <- as.matrix(counts)
mode(counts) <- "integer"


# Loading samples information from file.
all_samples <- read.table(
    snakemake@input[["samples"]], header=TRUE,
    row.names="sample", check.names=FALSE
)


dir.create(snakemake@output[["results"]], showWarnings=FALSE)
# Looping through the samples
for (cnd1 in c("A")) {
	for (cnd2 in c("B")) {

		exp <- sprintf("%s-%s", cnd1, cnd2)
		base <- sprintf("%s", cnd1)
		other <- sprintf("%s", cnd2)

		# Slicing
		samples <- subset(all_samples, condition==base | condition==other)
		count <- counts[,c(row.names(samples))]


		# Calculating DESeq2
		dds <- DESeqDataSetFromMatrix(
	    		countData=count,
	    		colData=samples,
	    		design= ~condition
		)


		dds$condition <- relevel(dds$condition, ref=base)
		dds <- DESeq(dds)
		results = results(dds, contrast=c("condition", cnd1, cnd2))

	    	# Writing results to file
		fname <- paste(
			snakemake@output[["results"]],
			paste(exp, "csv", sep='.'), sep='/'
		)

		write.csv(
		        as.data.frame(results),
		        file=fname,
		        quote=FALSE
		)
	}
}
