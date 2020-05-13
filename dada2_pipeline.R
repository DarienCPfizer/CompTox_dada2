rm(list=ls())

setwd("C:\\Users\\CAPUND\\Desktop\\mb_5_8_20")

#source("https://bioconductor.org/biocLite.R")
#biocLite("dada2")

#biocLite("devtools")
library("devtools")
#devtools::install_github("benjjneb/dada2")
library(tibble)
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")


path <- "C:/Users/CAPUND/Desktop/mb_5_8_20"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2]) # bad sample
plotQualityProfile(fnRs[3:4]) # okay sample primer issue at tail

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(290,250),
                     maxN=0, 
                     maxEE=c(6,6), # previously turned off 
                     truncQ=2, 
                     rm.phix=TRUE,
                     compress=TRUE, 
                     multithread=FALSE) # On Windows set multithread=FALSE
head(out)


errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

dadaRs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[50]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa")
unname(head(taxa))

#### Start of plots in phyloseq ####
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(data.table)
library(reshape)

#import mapping file similar to mothur import
mapfile = "metadata.csv"
map <- read.csv(mapfile, row.names = c(1), sep=",", header=TRUE)
head(map)
map <- sample_data(map)

#importing changed taxa file
#taxfile = "B18-6D Pfizer study sample metadata.xlsx"
#taxa2 <- read.csv(taxfile, row.names = c(1), sep=",", header=TRUE)
#head(taxa2)
#taxa2 <- sample_data(taxa2)

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
ps
#merge in meta file
dada_merge <- merge_phyloseq(ps, map)
dada_merge

# http://joey711.github.io/phyloseq-demo/phyloseq-demo.html  
#alpha diversity
p =plot_richness(dada_merge, x="Sample_Project", measures=c("Shannon", "Simpson", "Observed"), color="Sample_Project") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_boxplot() + geom_jitter(position=position_dodge(0.8))

# reorder x-axis
newSTorder = c("BS", "Day0", "Day01", "Day03", 
               "Day07", "Day10", "Day14", "Day21",
               "14wks")
p$data$Sample_Project <- as.character(p$data$Sample_Project)
p$data$Sample_Project <- factor(p$data$Sample_Project, levels=newSTorder)
print(p)

#Beta diversity
#NMDS/Bray
ord.nmds.brayDada <- ordinate(dada_merge, method="NMDS", distance="bray")
plot_ordination(dada_merge, ord.nmds.brayDada, color="Sample_Project", 
                shape="Description", 
                title="Bray NMDS Site")+ geom_point(size = 2)+
                geom_polygon(aes(fill = Sample_Project, group = Description))

# test from elias input
#dada_merge$grp <- paste0(dada_merge$Sample_Project, dada_merge$Description, sep = "_")
plot_ordination(dada_merge, ord.nmds.brayDada, 
                title="Bray NMDS Site")+ geom_point(size = 2, aes(color = Sample_Project, fill = Description))+
  geom_polygon(aes(group = grp, fill = Description), alpha = 0.2)


#end of test
#unifrac or mds, unifrac typically needs tree
plot_ordination(dada_merge, ordinate(dada_merge, "MDS"), title="MDS Site", color = "Description") + geom_point(size = 2)

#### From Mothur scripting ####

### For rarefying data ###
sample_sums(dada_merge)
dada_merge <- rarefy_even_depth(dada_merge, sample.size = min(2000), # rarefy here (currently 2000)
                                rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
sample_sums(dada_merge)

ordASV.nmds.b <- ordinate(dada_merge, method="NMDS", distance="bray")

ordASV.nmds.b<-ordinate(dada_merge,"NMDS","unifrac")

ps.97.mm<-plot_ordination(dada_merge,ordASV.nmds.b,type="samples",color="Sample_Order",shape="Sample_Order")

ps.97.mm+scale_shape_manual(values=c(15,16,17,18,6,4,8))+geom_point(size=4)


#Export OTU table (all rows, only 3 columns)
glom <- tax_glom(seqtab.nochim)
otus <- otu_table(glom)
write.csv(Order_abundance, file='OrderSilva_abundance.csv')

#agglomerate taxa
glom <- tax_glom(dada_merge, taxrank = 'Order')
#create dataframe from phyloseq object
dat <- psmelt(glom)
#convert Phylum to a character vector from a factor
dat$Order <- as.character(dat$Order)
#aggregate
Order_abundance <- aggregate(Abundance~Sample+Order, dat, FUN=sum)
#reorganize the table so that each phylum is a column
library(reshape)
Order_abundance <- cast(Order_abundance, Sample ~ Order)
write.csv(Order_abundance, file='taxfilefix.csv')

### tax info at 100% bar chart###
physeq2 = filter_taxa(dada_merge, function(x) mean(x) > 0.1, TRUE)
physeq2
physeq3 = transform_sample_counts(physeq2, function(x) x / sum(x) )
physeq3

glom <- tax_glom(physeq3, taxrank = 'Class')
glom # should list # taxa as # phyla
data <- psmelt(glom) # create dataframe from phyloseq object
data$Class <- as.character(data$Class) #convert to character

#simple way to rename phyla with < 1% abundance
data$Class[data$Abundance < 0.01] <- "< 1% abund."

library(plyr)
medians <- ddply(data, ~Class, function(x) c(median=median(x$Abundance)))
remainder <- medians[medians$median <= 0.01,]$Class # removing phyla that is less than 1%

data[data$Class %in% remainder,]$Class <- "Phyla < 1% abund."
#rename phyla with < 1% relative abundance
data$Class[data$Abundance < 0.01] <- "Phyla < 1% abund."

#plot with condensed phyla into "< 1% abund" category
data$Sample_Project <- factor(data$Sample_Project,levels = c("BS", "Day0", "Day01", "Day03", 
                                         "Day07", "Day10", "Day14", "Day21",
                                         "14wks"))
data$Sample <- factor(data$Sample,levels = c("2018012858-1", # ordered time, crl, pfizer
"2018012858-2","2018012858-3","2018012858-4","2018012858-5","2018012858-6",
"2018012858-7","2018012858-8","2018016611-1","2018016611-2","2018016611-3",
"2018016611-4","2018016611-5","2018016611-6","2018016611-7","2018016611-8", 
"2018016640-1","2018016640-2","2018016640-3","2018016640-4","2018016640-5",
"2018016640-6","2018016640-7","2018016640-8","2018013076-1","2018013076-2",
"2018013076-3","2018013076-4","2018013076-5","2018013076-6","2018013076-7",
"2018013076-8", "2018018009-1","2018018009-2","2018018009-3","2018018009-4",
"2018018009-5","2018018009-6","2018018009-7","2018018009-8","2018016617-1",
"2018016617-2","2018016617-3","2018016617-4","2018016617-5","2018016617-6",
"2018016617-7","2018016617-8","2018018010-1","2018018010-2","2018018010-3",
"2018018010-4","2018018010-5","2018018010-6","2018018010-7","2018018010-8",
"2018016622-1","2018016622-2","2018016622-3","2018016622-4","2018016622-5",
"2018016622-6","2018016622-7","2018016622-8","2018018012-1","2018018012-2",
"2018018012-3","2018018012-4","2018018012-5","2018018012-6","2018018012-7",
"2018018012-8","2018016630-1","2018016630-2","2018016630-3","2018016630-4",
"2018016630-5","2018016630-6","2018016630-7","2018016630-8","2018018014-1",
"2018018014-2","2018018014-3","2018018014-4","2018018014-5","2018018014-6",
"2018018014-7","2018018014-8","2018016634-1","2018016634-2","2018016634-3",
"2018016634-4","2018016634-5","2018016634-6","2018016634-7","2018016634-8",
"2018018025-1","2018018025-2","2018018025-3","2018018025-4","2018018025-5",
"2018018025-6","2018018025-7","2018018025-8","2018016637-1","2018016637-2",
"2018016637-3","2018016637-4","2018016637-5","2018016637-6","2018016637-7",
"2018016637-8","2018018037-1","2018018037-2","2018018037-3","2018018037-4",
"2018018037-5","2018018037-6","2018018037-7","2018018037-8","2018037630-1",
"2018037630-2","2018037630-3","2018037630-4","2018037630-5","2018037630-6",
"2018037630-7","2018037630-8","2018037631-1","2018037631-2","2018037631-3",
"2018037631-4","2018037631-5","2018037631-6","2018037631-7","2018037631-8"
))

p <- ggplot(data=data, aes(x=Sample_Project, y=Abundance, fill=Class))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                               "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Description) # added site facet, ncol = 1 for vertical plot placement

#Tree for Unifrac
#source("https://bioconductor.org/biocLite.R")
#biocLite("DECIPHER")
library(DECIPHER)
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

#construct a neighbor-joining tree, and then fit a Generalized time-reversible with Gamma rate variation
#maximum likelihood tree using the neighbor-joining tree as a starting point.

library(phangorn)   #choose no when from compilation 
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

#Procrustes
pro <- procrustes(ord.nmds.bray97, ord.nmds.bray100)
summary(pro)
plot(pro)
plot(pro, kind=2)
residuals(pro)

protest(X = ord.nmds.bray100, Y = ord.nmds.brayDada, permutations = 999)
plot(protestTest)
