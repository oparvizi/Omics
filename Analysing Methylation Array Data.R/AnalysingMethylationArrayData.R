# A cross-package Bioconductor workflow for analysing methylation array data
# Source: https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("methylationArrayAnalysis")
#library(methylationArrayAnalysis)


##---------------------------------------------Differential methylation analysis
## Obtaining the data-----------------------------------------------------------
# set up a path to the data directory
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
# list the files
list.files(dataDirectory, recursive = TRUE)

## Loading the data-------------------------------------------------------------
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)

# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

# read in the sample sheet for the experiment
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet.csv")
targets

# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets)
rgSet

# give the samples descriptive names
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
rgSet

## Quality control--------------------------------------------------------------
# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)

# examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(2,1))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
       bg="white")

qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="qcReport.pdf")

# remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet

# remove poor quality samples from targets data
targets <- targets[keep,]
targets[,1:5]

# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)

## Normalisation----------------------------------------------------------------
# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet) 

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)

# visualise what the data looks like before and after normalisation
par(mfrow=c(2,1))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))

## Data exploration-------------------------------------------------------------
# MDS plots to look at largest sources of variation
par(mfrow=c(2,1))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
        col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.7)


# Examine higher dimensions to look at other sources of variation
par(mfrow=c(3,1))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal, 
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

# Filtering---------------------------------------------------------------------
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% 
                                                      c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

# exclude cross reactive probes 
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)

mSetSqFlt <- mSetSqFlt[keep,] 
mSetSqFlt

par(mfrow=c(2,1))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], cex=0.8)
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.65, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Source)])
legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

par(mfrow=c(3,1))
# Examine higher dimensions to look at other sources of variation
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Source)], dim=c(1,3))
legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Source)], dim=c(2,3))
legend("topright", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Source)], dim=c(3,4))
legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])

bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])

par(mfrow=c(2,1))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))

## Probe-wise differential methylation analysis---------------------------------
# this is the factor of interest
cellType <- factor(targets$Sample_Group)
# this is the individual effect that we need to account for
individual <- factor(targets$Sample_Source) 

# use the above to create a design matrix
design <- model.matrix(~0+cellType+individual, data=targets)
colnames(design) <- c(levels(cellType),levels(individual)[-1])

# fit the linear model 
fit <- lmFit(mVals, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(naive-rTreg,
                            naive-act_naive,
                            rTreg-act_rTreg,
                            act_naive-act_rTreg,
                            levels=design)
contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

# get the table of results for the first contrast (naive - rTreg)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs)

write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)

# plot the top 4 most significantly differentially methylated CpGs 
par(mfrow=c(4,1))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")
})

## Differential methylation analysis of regions---------------------------------
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "naive - rTreg", arraytype = "450K")


str(myAnnotation)


#endif /* NEWSTUFF */
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges

# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]

# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 2, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")

## Customising visualisations of methylation data-------------------------------
# indicate which genome is being used
gen <- "hg19"
# the index of the DMR that we will plot 
dmrIndex <- 1
# extract chromosome number and location from DMR results 
chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))
# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))

# CpG islands
islandHMM <- read.csv(paste0(dataDirectory,
                             "/model-based-cpg-islands-hg19-chr17.txt"),
                      sep="\t", stringsAsFactors=FALSE, header=FALSE)
head(islandHMM)

islandData <- GRanges(seqnames=Rle(islandHMM[,1]), 
                      ranges=IRanges(start=islandHMM[,2], end=islandHMM[,3]),
                      strand=Rle(strand(rep("*",nrow(islandHMM)))))
islandData

# DNAseI hypersensitive sites
dnase <- read.csv(paste0(dataDirectory,"/wgEncodeRegDnaseClusteredV3chr17.bed"),
                  sep="\t",stringsAsFactors=FALSE,header=FALSE)
head(dnase)

dnaseData <- GRanges(seqnames=dnase[,1],
                     ranges=IRanges(start=dnase[,2], end=dnase[,3]),
                     strand=Rle(rep("*",nrow(dnase))),
                     data=dnase[,5])
dnaseData

iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                    from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                    rstarts="exonStarts", rends="exonEnds", gene="name", 
                    symbol="name2", transcript="name", strand="strand", 
                    fill="darkblue",stacking="squish", name="RefSeq", 
                    showId=TRUE, geneSymbol=TRUE)

ann450kOrd <- ann450kSub[order(ann450kSub$chr,ann450kSub$pos),]
head(ann450kOrd)

bValsOrd <- bVals[match(ann450kOrd$Name,rownames(bVals)),]
head(bValsOrd)

# create genomic ranges object from methylation data
cpgData <- GRanges(seqnames=Rle(ann450kOrd$chr),
                   ranges=IRanges(start=ann450kOrd$pos, end=ann450kOrd$pos),
                   strand=Rle(rep("*",nrow(ann450kOrd))),
                   betas=bValsOrd)
# extract data on CpGs in DMR
cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])

# methylation data track
methTrack <- DataTrack(range=cpgData, groups=targets$Sample_Group,genome = gen,
                       chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                       type=c("a","p"), name="DNA Meth.\n(beta value)",
                       background.panel="white", legend=TRUE, cex.title=0.8,
                       cex.axis=0.8, cex.legend=0.8)
# CpG island track
islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.", 
                               chromosome=chrom,fill="darkgreen")

# DNaseI hypersensitive site data track
dnaseTrack <- DataTrack(range=dnaseData, genome=gen, name="DNAseI", 
                        type="gradient", chromosome=chrom)

# DMR position data track
dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                            chromosome=chrom,fill="darkred")

tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack, dnaseTrack,
               rTrack)
sizes <- c(2,2,5,2,2,2,3) # set up the relative sizes of the tracks
plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))


##-----------------------------------------------------------Additional analyses
## Gene ontology testing--------------------------------------------------------
# Get the significant CpG sites at less than 5% FDR
sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
# First 10 significant CpGs
sigCpGs[1:10]

# Total number of significant CpGs at 5% FDR
length(sigCpGs)

# Get all the CpG sites used in the analysis to form the background
all <- DMPs$Name
# Total number of CpG sites tested
length(all)

par(mfrow=c(1,1))
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE)


# Top 10 GO categories
topGSA(gst, number=10)

# load Broad human curated (C2) gene sets
load(paste(dataDirectory,"human_c2_v5.rdata",sep="/"))
# perform the gene set test(s)
gsa <- gsameth(sig.cpg=sigCpGs, all.cpg=all, collection=Hs.c2)

# top 10 gene sets
topGSA(gsa, number=10)

## Differential variability-----------------------------------------------------
load(file.path(dataDirectory,"ageData.RData"))

# calculate detection p-values
age.detP <- detectionP(age.rgSet)

# pre-process the data after excluding poor quality samples
age.mSetSq <- preprocessQuantile(age.rgSet)

# add sex information to targets information
age.targets$Sex <- getSex(age.mSetSq)$predictedSex

# ensure probes are in the same order in the mSetSq and detP objects
age.detP <- age.detP[match(featureNames(age.mSetSq),rownames(age.detP)),]
# remove poor quality probes
keep <- rowSums(age.detP < 0.01) == ncol(age.detP) 
age.mSetSqFlt <- age.mSetSq[keep,]

# remove probes with SNPs at CpG or single base extension (SBE) site
age.mSetSqFlt <- dropLociWithSnps(age.mSetSqFlt, snps = c("CpG", "SBE"))

# remove cross-reactive probes
keep <- !(featureNames(age.mSetSqFlt) %in% xReactiveProbes$TargetID)
age.mSetSqFlt <- age.mSetSqFlt[keep,] 


# tag sex chromosome probes for removal
keep <- !(featureNames(age.mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% 
                                                          c("chrX","chrY")])

age.pal <- brewer.pal(8,"Set1")
par(mfrow=c(2,1))
plotMDS(getM(age.mSetSqFlt), top=1000, gene.selection="common", 
        col=age.pal[factor(age.targets$Sample_Group)], labels=age.targets$Sex, 
        main="With Sex CHR Probes")
legend("topleft", legend=levels(factor(age.targets$Sample_Group)), 
       text.col=age.pal)

plotMDS(getM(age.mSetSqFlt[keep,]), top=1000, gene.selection="common", 
        col=age.pal[factor(age.targets$Sample_Group)], labels=age.targets$Sex, 
        main="Without Sex CHR Probes")
legend("top", legend=levels(factor(age.targets$Sample_Group)),
       text.col=age.pal)


# remove sex chromosome probes from data
age.mSetSqFlt <- age.mSetSqFlt[keep,]

# get M-values for analysis
age.mVals <- getM(age.mSetSqFlt)

design <- model.matrix(~factor(age.targets$Sample_Group)) 
# Fit the model for differential variability
# specifying the intercept and age as the grouping factor
fitvar <- varFit(age.mVals, design = design, coef = c(1,2))

# Summary of differential variability
summary(decideTests(fitvar))

topDV <- topVar(fitvar, coef=2)
# Top 10 differentially variable CpGs between old vs. newborns
topDV

# get beta values for ageing data
age.bVals <- getBeta(age.mSetSqFlt)
par(mfrow=c(4,1))
sapply(rownames(topDV)[1:4], function(cpg){
  plotCpg(age.bVals, cpg=cpg, pheno=age.targets$Sample_Group, 
          ylab = "Beta values")
})

## Cell type composition--------------------------------------------------------
# load sorted blood cell data package
library(FlowSorted.Blood.450k)
# ensure that the "Slide" column of the rgSet pheno data is numeric
# to avoid "estimateCellCounts" error
pData(age.rgSet)$Slide <- as.numeric(pData(age.rgSet)$Slide)
# estimate cell counts
cellCounts <- estimateCellCounts(age.rgSet)
# plot cell type proportions by age
par(mfrow=c(1,1))
a = cellCounts[age.targets$Sample_Group == "NewBorns",]
b = cellCounts[age.targets$Sample_Group == "OLD",]
boxplot(a, at=0:5*3 + 1, xlim=c(0, 18), ylim=range(a, b), xaxt="n", 
        col=age.pal[1], main="", ylab="Cell type proportion")
boxplot(b, at=0:5*3 + 2, xaxt="n", add=TRUE, col=age.pal[2])
axis(1, at=0:5*3 + 1.5, labels=colnames(a), tick=TRUE)
legend("topleft", legend=c("NewBorns","OLD"), fill=age.pal)


## Discussion-------------------------------------------------------------------
sessionInfo()













