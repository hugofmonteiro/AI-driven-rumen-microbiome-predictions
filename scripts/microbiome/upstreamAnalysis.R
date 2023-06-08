#Loading the package
library(dada2); packageVersion("dada2")

#Loading fastq files
path <- "~/Fastq files/"
list.files(path)

# Fastq filenames have format: SAMPLENAME_R1_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles 
plotQualityProfile(fnFs[1:12])
plotQualityProfile(fnFs[60:72])
plotQualityProfile(fnFs[120:132])
plotQualityProfile(fnFs[180:192])
plotQualityProfile(fnFs[240:252])
plotQualityProfile(fnFs[300:312])
plotQualityProfile(fnFs[360:372])
plotQualityProfile(fnFs[420:432])

#Assigning file names to the filtered fastq files (Place filtered files in filtered/subdirectory)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

#Using standard filtering parameters
out <- filterAndTrim(fnFs, filtFs, truncLen=c(220), trimLeft=20,
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)

#Visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)

#Applying the core sample inference algorigth
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaFs[[1]] #Infererring true sequences

#Constructing a sequence table (ASV table)
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#frequency of chimeric sequences
sum(seqtab.nochim)/sum(seqtab)

#Checking the number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
sum(seqtab.nochim) 

#Taxonomy assignment
tax <- assignTaxonomy(seqtab.nochim, "~/Taxonomy Assignment/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
write.csv(tax, "Microbiome_tax_final_03Oct2022.csv")

#phyloseq
dim(seqtab.nochim) 
dim(tax)

rownames(seqtab.nochim)
samples.out <- rownames(seqtab.nochim)
samples.out # should be 454

# Load mapping file
map <- "ALL_complete_without_repeatedAnimals_wGroups.txt"

library(phyloseq); packageVersion("phyloseq")
library(microbiome)
library(devtools) # update all packages
library(ggplot2)
theme_set(theme_bw())

# Importing the mapping file as a sample metadata 
sample_metadata = import_qiime_sample_data (map) # here we import our mapping file
# look at this object in the "environment" -- are the variables mainly listed as "factors"?
colnames(sample_metadata)
head(sample_metadata)

# Make phyloseq object
# here we make our phyloseq ASV table, combining the reads and the metadata
mydata <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(sample_metadata), 
                   tax_table(tax)) # if you get an error, make sure that the sample names from seqtab.nochim (in environment) matches the first column of your .txt file, change .txt file accordingly
mydata_asv <- as.data.frame(otu_table(mydata)) # here we make ps a data.frame from our phyloseq object

# Item 2. Working with your phyloseq object "mydata"
summarize_phyloseq(mydata)
ps.r = transform_sample_counts(mydata, function(x) x/sum(x)) # this function transforms the data into a relative abundance (for each individual)
ps.r
ps.r_asv <- as.data.frame(otu_table(ps.r))


### TAXONOMY (Creating the relative abundances)
library(reshape2)
library(phyloseq)
ps_phylum <- tax_glom(mydata, "Phylum"); write.csv(otu_table(ps_phylum), "ps_phylum.csv")
ps_class <- tax_glom(mydata, "Class"); write.csv(otu_table(ps_class), "ps_class.csv")
ps_order <- tax_glom(mydata, "Order"); write.csv(otu_table(ps_order), "ps_order.csv")
ps_family <- tax_glom(mydata, "Family"); write.csv(otu_table(ps_family), "ps_family.csv")
ps_genus <- tax_glom(mydata, "Genus"); write.csv(otu_table(ps_genus), "ps_genus.csv")

ps_phylum.r <- transform_sample_counts(ps_phylum, function(x) x / sum(x)); write.csv(otu_table(ps_phylum.r), "ps_phylum.r.csv"); write.csv(tax_table(ps_phylum.r), "ps_phylum.r_tax.csv")
ps_class.r <- transform_sample_counts(ps_class, function(x) x / sum(x)); write.csv(otu_table(ps_class.r), "ps_class.r.csv"); write.csv(tax_table(ps_class.r), "ps_class.r_tax.csv")
ps_order.r <- transform_sample_counts(ps_order, function(x) x / sum(x)); write.csv(otu_table(ps_order.r), "ps_order.r.csv"); write.csv(tax_table(ps_order.r), "ps_order.r_tax.csv")
ps_family.r <- transform_sample_counts(ps_family, function(x) x / sum(x)); write.csv(otu_table(ps_family.r), "ps_family.r.csv"); write.csv(tax_table(ps_family.r), "ps_family.r_tax.csv")
ps_genus.r <- transform_sample_counts(ps_genus, function(x) x / sum(x)); write.csv(otu_table(ps_genus.r), "ps_genus.r.csv"); write.csv(tax_table(ps_genus.r), "ps_genus.r_tax.csv")

# CLR transformation on phyloseq object
phylum_clr<-microbiome::transform(ps_phylum, "clr")
class_clr<-microbiome::transform(ps_class, "clr")
order_clr<-microbiome::transform(ps_order, "clr")
family_clr<-microbiome::transform(ps_family, "clr")
genus_clr<-microbiome::transform(ps_genus, "clr")

write.csv(otu_table(phylum_clr), "phylum_clr.csv"); write.csv(tax_table(phylum_clr), "phylum_clr_tax.csv")
write.csv(otu_table(class_clr), "class_clr.csv"); write.csv(tax_table(class_clr), "class_clr_tax.csv")
write.csv(otu_table(order_clr), "order_clr.csv"); write.csv(tax_table(order_clr), "order_clr_tax.csv")
write.csv(otu_table(family_clr), "family_clr.csv"); write.csv(tax_table(family_clr), "family_clr_tax.csv")
write.csv(otu_table(genus_clr), "genus_clr.csv"); write.csv(tax_table(genus_clr), "genus_clr_tax.csv")

#### ALPHA-DIVERSITY Indexes
alpha <- microbiome::alpha
tab1 <- alpha(mydata, index = "all") # observed richness, diversity, evenness, dominance, rarity
head(tab1)
tab1$Sample_ID <- rownames(tab1) # this creates a column called Sample_ID using the rownames (which are the sample names)
write.csv(tab1, "alphadiversity_03Oct2022.csv")
