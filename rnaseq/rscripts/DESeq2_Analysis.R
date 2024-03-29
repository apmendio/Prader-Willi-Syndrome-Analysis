Author = "Aron Judd P. Mendiola"

# Raw text files can be found in the cluster directory path below
# "/share/lasallelab/Aron/RNA_seq_analysis/Analysis_March252021/GeneCounts/ReverseCounts"

# setwd = set working directory; equivalent to the Linux "cd".
# the R equivalent to the Linux pwd is getwd() = get working directory.
setwd("/Users/aron/Desktop/LaSalle Lab/Analysis/Deseq2/ReverseCounts")

# separated male and female samples into different folders for individual analysis of each timepoint
setwd("/Users/aron/Desktop/LaSalle Lab/Analysis/Deseq2/ReverseCountsMF/Male")

setwd("/Users/aron/Desktop/LaSalle Lab/Analysis/Deseq2/ReverseCountsMF/Female")
# load package DESeq2 (all functions)
library(DESeq2)

# read in the sample sheet
# header = TRUE: the first row is the "header", i.e. it contains the column names.
# sep = "\t": the columns/fields are separated with tabs.
sampletable <- read.table("sample_sheet.txt", header=T, sep="\t")

# sample analysis of male and female samples
malesampletable <- read.table("sample_sheet_male.txt", header=T, sep="\t")
  
femalesampletable <- read.table("sample_sheet_female.txt", header=T, sep="\t")

# add the sample names as row names (it is needed for some of the DESeq functions)
rownames(sampletable) <- sampletable$SampleName

# sample analysis of male and female samples
rownames(malesampletable) <- malesampletable$SampleName

rownames(femalesampletable) <- femalesampletable$SampleName
# display the first 6 rows
head(sampletable)

head(malesampletable)

head(femalesampletable)
# check the number of rows and the number of columns
nrow(sampletable) # if this is not 6, please raise your hand !
ncol(sampletable) # if this is not 4, also raise your hand !

nrow(malesampletable) 
ncol(malesampletable)

nrow(femalesampletable)
ncol(femalesampletable)

# only return file names with a given pattern
dir(pattern="counts.txt")

# save the results to a variable
files <- dir(pattern="counts.txt")

counts <- c()
for( i in seq_along(files) ){
  x <- read.table(file=files[i], sep="\t", header=F, as.is=T)
  counts <- cbind(counts, x[,2])
}

# sample analysis of male and female samples

mfiles <- dir(pattern="counts.txt")

mcounts <- c()
for( i in seq_along(mfiles) ){
  x <- read.table(file=files[i], sep="\t", header=F, as.is=T)
  mcounts <- cbind(mcounts, x[,2])
}

ffiles <- dir(pattern="counts.txt")

fcounts <- c()
for( i in seq_along(ffiles) ){
  x <- read.table(file=ffiles[i], sep="\t", header=F, as.is=T)
  fcounts <- cbind(fcounts, x[,2])
}

# set the row names
rownames(counts) <- x[,1]
# set the column names based on input file names, with pattern removed (if generated skip to line 54)
colnames(counts) <- sub("counts.txt","",files)

# sample analysis for male and female samples

# set the row names males
rownames(mcounts) <- x[,1]
# set the column names based on input file names, with pattern removed (if generated skip to line 54)
colnames(mcounts) <- sub("mcounts.txt","",files)

# set the row names females
rownames(fcounts) <- x[,1]
# set the column names based on input file names, with pattern removed (if generated skip to line 54)
colnames(fcounts) <- sub("fcounts.txt","",files)

# Option 1 that reads in a matrix (we will not do it here):

# first read in the matrix (only required if matrix file was generated in terminal)
# count_matrix <- read.delim("counts", header=T, sep="\t", row.names=1)

# then create the DESeq object
# countData is the matrix containing the counts
# sampletable is the sample sheet / metadata we created
# design is how we wish to model the data: what we want to measure here is the difference between the treatment times
se_star_matrix <- DESeqDataSetFromMatrix(countData = counts,
                                         colData = sampletable,
                                         design = ~ Timepoint)

pe_star_matrix <- DESeqDataSetFromMatrix(countData = counts,
                                         colData = sampletable,
                                         design = ~ Timepoint + Sex)

se_star_matrix_update <- DESeqDataSetFromMatrix(countData = counts,
                                         colData = sampletable,
                                         design = ~ Timepoint + Sex)

# sample analysis for male and female samples

male_star_matrix <- DESeqDataSetFromMatrix(countData = mcounts,
                                           colData = malesampletable,
                                           design = ~ Timepoint)

female_star_matrix <- DESeqDataSetFromMatrix(countData = fcounts,
                                             colData = femalesampletable,
                                             design = ~ Timepoint)
# Number of genes before filtering:
nrow(se_star_matrix)

nrow(pe_star_matrix)

nrow(se_star_matrix_update)

# sample analysis for male and female samples

nrow(male_star_matrix)

nrow(female_star_matrix)

# Filter
se_star_matrix <- se_star_matrix[rowSums(counts(se_star_matrix)) > 10, ]

pe_star_matrix <- pe_star_matrix[rowSums(counts(pe_star_matrix)) > 10, ]

se_star_matrix_update <- se_star_matrix_update[rowSums(counts(se_star_matrix_update)) > 10, ]

# Filter male and female samples

male_star_matrix_filtered <- male_star_matrix[rowSums(counts(male_star_matrix)) > 10, ]

female_star_matrix_filtered <- female_star_matrix[rowSums(counts(female_star_matrix)) > 10, ]

# Number of genes left after low-count filtering:
nrow(se_star_matrix)

nrow(pe_star_matrix)

nrow(se_star_matrix_update)

# Number of genes left after low-count filtering male and female samples:
nrow(male_star_matrix_filtered)

nrow(female_star_matrix_filtered)

# Run Deseq2
pe_star_matrix2 <- DESeq(pe_star_matrix)

se_star_matrix_update <- DESeq(se_star_matrix_update)

# Run Deseq2 of male and female samples
male_star_matrix2 <- DESeq(male_star_matrix_filtered)

female_star_matrix2 <- DESeq(female_star_matrix_filtered)

resultsNames(male_star_matrix2)

resultsNames(female_star_matrix2)
# compute normalized counts (log2 transformed); + 1 is a count added to avoid errors during the log2 transformation: log2(0) gives an infinite number, but log2(1) is 0.
# normalized = TRUE: divide the counts by the size factors calculated by the DESeq function
norm_counts <- log2(counts(se_star_matrix2, normalized = TRUE)+1)

# normalize DESeq counts of male and female samples
norm_mcounts <- log2(counts(male_star_matrix2, normalized = TRUE)+1)

norm_fcounts <- log2(counts(female_star_matrix2, normalized = TRUE)+1)

resultsNames(norm_mcounts)

resultsNames(norm_fcounts)

# Load the tximport package that we use to import Salmon counts
library(tximport)

# Read in the two-column data.frame linking transcript id (column 1) to gene id (column 2)
ensid <- read.table("ensembl_mm_100.tsv", 
                      sep="\t",
                      header=F)

# add the gene symbols
norm_counts_symbols <- merge(unique(tx2gene[,1:3]), data.frame(ID=rownames(norm_counts), norm_counts), by=1, all=F, verbose=T)

# add the gene symbols to male and female samples
norm_mcounts_symbols <- merge(unique(ensid[,1:3]), data.frame(ID=rownames(norm_mcounts), norm_mcounts), by=1, all=F, verbose=T)

norm_fcounts_symbols <- merge(unique(ensid[,1:3]), data.frame(ID=rownames(norm_fcounts), norm_fcounts), by=1, all=F, verbose=T)

# write normalized counts to text file
write.table(norm_mcounts_symbols, "normalized_mcounts.txt", quote=F, col.names=T, row.names=F, sep="\t")

write.table(norm_fcounts_symbols, "normalized_fcounts.txt", quote=F, col.names=T, row.names=F, sep="\t")

# Try with the vst transformation
vsd <- vst(se_star_matrix2)

vsd <- vst(pe_star_matrix2)

vsd2 <- vst(se_star_matrix_update)

# Try with the vst transformation male and female samples

vsdm <- vst(male_star_matrix2)
mmat <- assay(vsdm)
mmat <- limma::removeBatchEffect(mmat, vsdm$batch)
assay(vsdm) <- mmat
png("PCA_norm_male.png")
plotPCA(object = vsdm,
        intgroup = "Timepoint")
dev.off()

# also possible to perform custom transformation:
dds <- estimateSizeFactors(male_star_matrix2)
# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method.
png("PCA_norm_male2.png")
plotPCA( DESeqTransform( se ), intgroup = "Timepoint" )
dev.off()

vsdf <-vst(female_star_matrix2)


# load libraries pheatmap to create the heatmap plot
library(pheatmap)

# calculate between-sample distance matrix
sampleDistMatrix <- as.matrix(dist(t(assay(vsd))))

sampleDistMatrix2 <- as.matrix(dist(t(assay(vsd2))))

# calculate sample distance matrix between samples
msampleDisMatrix <- as.matrix(dist(t(assay(vsdm))))

# create figure in PNG format
png("sample_distance_heatmap_star3.png")
pheatmap(sampleDistMatrix2)
dev.off() 
# create figure in PNG format for male and female samples
png("male_counts_heatmap.png")
pheatmap(msampleDisMatrix)
dev.off() 

png("female_counts_heatmap.png")
pheatmap(msampleDisMatrix)
dev.off() 
# close PNG file after writing figure in it
dev.off() 

# create PCA plot

png("PCA_star7.png")
plotPCA(object = vsd,
        intgroup = "Sex", "Timepoint")
dev.off()

png("PCA_star6.png")
plotPCA(object = vsd,
        intgroup = "Timepoint")
dev.off()

png("PCA_star8.png")
plotPCA(object = vsd2,
        intgroup = "Timepoint")
dev.off()

# create PCA plot for male and female samples

png("PCA_male_Timepointdist.png")
plotPCA(object = vsdm,
        intgroup = "Timepoint")
dev.off()

png("PCA_female_Timepointdist.png")
plotPCA(object = vsdf,
        intgroup = "Timepoint")
dev.off()

# check results names: depends on what was modeled. Here it was the "Timepoint"
resultsNames(se_star_matrix2)

# extract results for t25 vs t0
# contrast: the column from the metadata that is used for the grouping of the samples (Time), then the baseline (t0) and the group compared to the baseline (t25) -> results will be as "t25 vs t0"
de <- results(object = se_star_matrix2, 
              name="Timepoint_ZT6_vs_ZT0")

ZT0_vs_ZT6de <- results(object = male_star_matrix2, 
              name="Timepoint_ZT6_vs_ZT0")

ZT0_vs_ZT3de <- results(object = male_star_matrix2, 
              name="Timepoint_ZT3_vs_ZT0")

ZT0_vs_ZT12de <- results(object = male_star_matrix2, 
              name="Timepoint_ZT12_vs_ZT0")


# processing the same results as above but including the log2FoldChange shrinkage
# useful for visualization and gene ranking
de_shrink <- lfcShrink(dds = se_star_matrix2,
                       coef="Timepoint_ZT6_vs_ZT0",
                       type="apeglm")

ZT0_vs_ZT6de_shrink <- lfcShrink(dds = se_star_matrix2,
                                 coef="Timepoint_ZT6_vs_ZT0",
                                 type="apeglm")

ZT0_vs_ZT3de_shrink <- lfcShrink(dds = se_star_matrix2,
                                 coef="Timepoint_ZT3_vs_ZT0",
                                 type="apeglm")

ZT0_vs_ZT12de_shrink <- lfcShrink(dds = se_star_matrix2,
                                  coef="Timepoint_ZT12_vs_ZT0",
                                  type="apeglm")

# check first rows of both results
head(de)
head(de_shrink)

head(ZT0_vs_ZT6de)
head(ZT0_vs_ZT6de_shrink)
head(ZT0_vs_ZT3de)
head(ZT0_vs_ZT3de_shrink)
head(ZT0_vs_ZT12de)
head(ZT0_vs_ZT12de_shrink)

# add the gene symbols to male and female samples
norm_trial_ZT0vZT6 <- merge(unique(ensid[,1:3]), data.frame(ID=rownames(ZT0_vs_ZT6de_shrink), ZT0_vs_ZT6de_shrink), by=1, all=F, verbose=T)

norm_fcounts_symbols <- merge(unique(ensid[,1:3]), data.frame(ID=rownames(norm_fcounts), norm_fcounts), by=1, all=F, verbose=T)

# write normalized counts to text file
write.table(norm_trial_ZT0vZT6, "Zt0vsZt6", quote=T, col.names=T, row.names=T, sep="\t")

