
library("edgeR")
library("NMF")
library("ggrepel")
library("ggplot2")
library("ggfortify")

########################## load the counts and sample annotation tables  ##################################

# read the counts and sample annotation tables
counts <- read.table("counts.txt",header=TRUE, row.names=1)
sampleAnnots <- read.table("sample-annotation.txt",row.names=1,header=TRUE)	

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

# this table will go through the filtration 
lowCountsNormal <- 0
lowCountsLesion <- 0
numberOfSamples = ncol(counts)
numberOfGenes <- nrow(counts)

# the factor assignment is determined by alphabetical order so lesion->"1", normal->"2"
lesionSamples <- with(dgeList, table(samples$group))[1]
normalSamples <- with(dgeList, table(samples$group))[2]

# a vector that keeps the gene IDs that were filtered out for later use
genesRemoved <- vector()
  
# run through the genes (rows, i) and the samples per gene (columns, j)
# whenever finding a gene woth low counts, keep the gene id for later use
# criteria for filtering out: cpm < 1 for 0.75 of the samples for each condition (normal, lesion)
i <- 1
while (i <= numberOfGenes) {
	for (j in 1:numberOfSamples) {
		temp <- (colnames(filteredCPM))[j]
		if ((sampleAnnots[c(temp),] == "normal") & (filteredCPM[i,j] < 1))
			lowCountsNormal = lowCountsNormal+1
		if ((sampleAnnots[c(temp),] == "lesional") & (filteredCPM[i,j] < 1))
			lowCountsLesion = lowCountsLesion+1
	}
	if ((lowCountsNormal > 0.75 * normalSamples) & (lowCountsLesion > 0.75 * lesionSamples)) {
		# this gene is lowly expressed and should be filtered out
		genesRemoved <- c(genesRemoved, genesFromCounts[i]) 
	}
	# move on to the next one
	lowCountsNormal <- 0
	lowCountsLesion <- 0
	i = i + 1	
}

print(paste("Got the list of low expression gene to filter out, number of genes to be excluded:", length(genesRemoved)))

########################## make a logCPM filtered matrix and save as a binary file #####################

# subtrsct the lines corresponding to genes that have low counts to generate a filtered logCPM matrix for the next steps
genesToKeep <- setdiff(genesFromCounts, genesRemoved)
filteredCounts <- counts[genesToKeep,]

## generate a dgeList object from the filtered matrix, use logCPM this time and 
filteredDgeList <- DGEList(counts=filteredCounts,group=group)

# cpm function from edgeR, prior.count is a starting value used to offset to prevent zero counts
filteredLogCpmMatrix <- cpm(filteredDgeList, normalized.lib.sizes=FALSE, log=TRUE, prior.count=0.25)
print("Generated filtered log CPM matrix.")

#save this into an RDS binary format:
saveRDS(filteredLogCpmMatrix, "filteredLogCpmMatrix.rds")
print("Saved the filtered LogCPM matrix as an rds binary file.")

