
library("edgeR")
library("NMF")
library("ggrepel")
library("ggplot2")
library("ggfortify")

########################## load the counts and sample annotation tables  ##################################

# read the counts and sample annotation tables
counts <- read.table("counts.txt",header=TRUE, row.names=1)
sampleAnnots <- read.table("sample-annotation.txt",row.names=1,header=TRUE)	
geneAnnots <- read.delim("gene-annotation.txt", header=TRUE)

################### check consistency of the input files and report problems if exist #####################

# get the sample names from the counts table and the sample annottaion table
samplesFromCounts <- colnames(counts)
samplesFromAnnots <- rownames(sampleAnnots)
genesFromCounts <- rownames(counts)
inconsistency = FALSE

# check if any name appears in only one of the files and not the other
differenceAc <- setdiff(samplesFromAnnots, samplesFromCounts)
differenceCa <- setdiff(samplesFromCounts, samplesFromAnnots)

# no differences, do nothing
if (length(differenceAc) + length(differenceCa) == 0)
	print ("Sample names in input files match.")
if (length(differenceAc) > 0) { #there is a difference, print informative message
	print("The following sample(s) are missing from the counts file, fix it and run again:")
	print (differenceAc)
	inconsistency = TRUE
}
if (length(differenceCa) > 0) { #there is a difference, print informative message
	print("The following sample(s) are missing from the sample annotation file, fix it and run again:")
	print (differenceCa)
	inconsistency = TRUE
}
if (inconsistency == TRUE) # abort
	quit(status=1)

################################# generate Count Per Million matrix ######################################

# generate an edgeR DGEList object that has the sample counts and anotations in it
group <- factor(c(sampleAnnots[,1]))
dgeList <- DGEList(counts=counts,group=group)

# cpm function from edgeR, prior.count is a starting value used to offset to prevent zero counts
cpmMatrix <- cpm(dgeList, normalized.lib.sizes=FALSE, log=FALSE, prior.count=0.25)
print("Generated CPM matrix.")

###################################### filter the cpm matrix ###########################################

numberOfSamples = ncol(counts)
numberOfGenes <- nrow(counts)

# for plotting of cpm per gene later
dgeListCpm <- DGEList(counts=cpmMatrix,group=group)
  
# run through the genes (rows), keep lines where at least 0.75 of the samples had at least 1 count.
threshold <- 0.75 * numberOfSamples
j <- 0
keep <- 0
newMatrix<-0
for (i in 1:numberOfGenes) {
	if (sum(counts[i,] >= 1) >= threshold) {
		j=j+1
		keep <- c(keep,i)
	}
}
newMatrix<-counts[keep,]
print(paste("Got the low expression genes filtered out, number of genes kept:", length(keep)))

########################## make a logCPM filtered matrix and save as a binary file #####################

filteredCounts <- newMatrix

## generate a dgeList object from the filtered matrix, use logCPM this time and 
filteredDgeList <- DGEList(counts=filteredCounts,group=group)

# cpm function from edgeR, prior.count is a starting value used to offset to prevent zero counts
filteredLogCpmMatrix <- cpm(filteredDgeList, normalized.lib.sizes=FALSE, log=TRUE, prior.count=0.25)
print("Generated filtered log CPM matrix.")

#save this into an RDS binary format:
saveRDS(filteredLogCpmMatrix, "filteredLogCpmMatrix.rds")
print("Saved the filtered LogCPM matrix as an rds binary file.")

################################### PCA and other plots, exclusion of outliers #############################

colors <- c("red", "blue") # red/"1"s are Lesion samples, blue/"2"s are Normal

# plot the raw counts per sample, color represents group (normal/lesion) to check for possible bias in coverage
densities=rep("-1", length(samplesFromCounts))
pdf("rawCountsPerSample.pdf")
barplot(colSums(dgeList$counts), las=2, main="Counts per sample", cex.axis=0.8, cex.names=0.5, col=colors[group], axisnames=FALSE, space=0.6, density=densities, border = NA, xlab="Samples", ylab="Counts")
dev.off()

# plot the gene expression levels per gene based on the cpm values
pdf("cpmPerGene.pdf")
barplot(rowSums(dgeListCpm$counts), las=2, main="CPM per gene", axisnames=FALSE, cex.axis=0.8, ylab="CPM", xlab="Genes")
dev.off()

# dendrogram that shows hierarchical clustering analysis, reflects the distance between samples based on logCPM data.
# transpose the filtered counnts matrix so each line describes an element on the plot and the columns are the variables
transposedFLCPM <- t(filteredLogCpmMatrix)
distance <- dist(transposedFLCPM)
clusters <- hclust(distance)
pdf("dendrogram.pdf")
plot(clusters, labels=group)
dev.off()

## PCA and detection of outliers:
# add the column with group information
transposedFLCPMwGroup <- cbind(transposedFLCPM, group)
pdf("PCA_FLCPM.pdf")
autoplot(prcomp(transposedFLCPM), data=transposedFLCPMwGroup, colour='group')
dev.off()

# add sample labels to identify the outlier lesion sample on the bottom right and one that seems mislabels as lesion
pdf("PCA_FLCPM_wLabels.pdf")
autoplot(prcomp(transposedFLCPM), data=transposedFLCPMwGroup, colour='group', label=TRUE, label.size=3)
dev.off()

# remove the outlier and the suspected mislabeled sample:
outlier <- "SRR1146078"
mislabeled <- "SRR1146216"
samplesToRemove <- c(outlier, mislabeled)
# vector of the samples we want to keep
samplesToKeep <- setdiff(samplesFromAnnots, samplesToRemove)

# subset of the samples and their annotations that we want to keep
sampleAnnotsToKeep <- as.matrix(sampleAnnots[samplesToKeep,]) 

# define a new 'group' factor that does not include the 2 samples that are removed
newGroup <- factor(c(sampleAnnotsToKeep[,1]))

#generate new objects for the downstream analysis
filteredCountsNoOutliers <- filteredCounts[,samplesToKeep]

## generate a dgeList object from the new matrix, use logCPM this time and 
filteredDgeListNoOutliers <- DGEList(counts=filteredCountsNoOutliers,group=newGroup)

# get logCPM
filteredLogCpmNoOuts <- cpm(filteredDgeListNoOutliers, normalized.lib.sizes=FALSE, log=TRUE, prior.count=0.25)
print("Removed outlier/mislabeles genes.")

################################# Differential expression analysis by edgeR #####################################

# Normalization is not needed for single factor analysis.
# The dispersion is estimated before testing for differential expression, using the filtered counts. 
# edgeR uses quantile-adjusted conditional maximum likelihood (qCML) method for experiments with a single factor.
# The design matrix defines the variables and levels considered in the model, in this case one factor:
design <- model.matrix(~newGroup)
filteredDgeListNoOutliers <- estimateDisp(filteredDgeListNoOutliers, design)

# exact test for the negative binomial distribution to compute exact p-values to assess differential expression
# between the normal and lesion samples.
# The function exactTest computes exact p-values by summing over all sums of counts that have a probability less
# than the probability under the null hypothesis of the observed sum of counts.
et <- exactTest(filteredDgeListNoOutliers)
topTags(et)

# modifying the row names of the differentially expressed genes table to be able to merge with the gene annotations
tempDE <- et$table
ENSEMBL <- rownames(tempDE)
rownames(tempDE) <- NULL
tempDEWithNames <- cbind(ENSEMBL, tempDE)
annotatedDE <- merge(tempDEWithNames, geneAnnots, by="ENSEMBL")

# sort the list of DEGs by ascending p-values and write to a tab-delimited text file
sorted <- annotatedDE[order(annotatedDE$PValue),]
write.table(sorted, "differntialExpression.txt", sep="\t")
print("Performed differential expression analysis, see 'differentialExpression.txt'")
