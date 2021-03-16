
########################################
### V4 processing for Abe Gonzalez
### Lou LaMartina, started March 15 2021
########################################


library(dada2)
library(decontam)
library(ggplot2)
library(phyloseq)


###################################
### prepare data for processing ###
###################################

# set working directory
setwd("~/Desktop/Lab/Projects/Gonzalez")


# set file paths
path <- "./cutadapt"
pathF <- "./cutadapt/fastqF"
pathR <- "./cutadapt/fastqR"


# set paths for filtered reads
filtered_pathF <- "./cutadapt/fastqF/Filtered"
filtered_pathR <- "./cutadapt/fastqR/Filtered"


# sort forward and reverse reads into file paths
fastqFs <- sort(list.files(pathF, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fastqRs <- sort(list.files(pathR, pattern = "_R2_001.fastq.gz", full.names = TRUE))


# extract file names
sample_names <- sapply(strsplit(basename(fastqFs), "_"), '[', 1)




################################
### Inspect sequence quality ###
################################

# visualize quality of reads
qualityF.plot <- plotQualityProfile(fastqFs[1:4]); qualityF.plot
qualityR.plot <- plotQualityProfile(fastqRs[1:4]); qualityR.plot


# check: if there are not the same number of F and R files, stop.
if(length(fastqFs) != length(fastqRs)) 
  stop("Forward and reverse files do not match.")


# save quality profiles
#ggsave("./Plots/qualityF.pdf", plot = qualityF.plot, device = "pdf", width = 12, height = 8, units = "in")
#ggsave("./Plots/qualityR.pdf", plot = qualityR.plot, device = "pdf", width = 12, height = 8, units = "in")




#######################
### Quality control ###
#######################

# give filtered files new names and paths
filteredFs <- file.path(filtered_pathF, paste0(sample_names, "_F_filt.fastq.gz"))
filteredRs <- file.path(filtered_pathR, paste0(sample_names, "_R_filt.fastq.gz"))


# filter based on quality and read length
filtered_out <- filterAndTrim(fastqFs, filteredFs, fastqRs, filteredRs, 
                              maxEE = 2, maxN = 0, truncQ = 10, #truncLen = 230,
                              rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE)


# inspect how many reads were filtered out of each sample
filtered_out
(1 - (filtered_out[,2] / filtered_out[,1])) * 100
cat(mean((1 - (filtered_out[,2] / filtered_out[,1])) * 100), "% removed\n")
# 230 truncate: 17.9638 % removed
# not truncate: 16.64507 % removed


# plot quality profiles of filtered reads
#filtF.plot <- plotQualityProfile(filteredFs[1:4]); filtF.plot
#filtR.plot <- plotQualityProfile(filteredRs[1:4]); filtR.plot


# set sample names to the ID only
names(filteredFs) <- sample_names
names(filteredRs) <- sample_names


# save quality profiles
#ggsave("./Plots/filt_qualityF.pdf", plot = filtF.plot, device = "pdf", width = 12, height = 8, units = "in")
#ggsave("./Plots/filt_qualityR.pdf", plot = filtR.plot, device = "pdf", width = 12, height = 8, units = "in")




############################
### Learning error rates ###
############################

# learn and visualize error rates of F reads
errorF <- learnErrors(filteredFs, multithread = TRUE)
errorF.plot <- plotErrors(errorF, nominalQ = TRUE, err_in = TRUE, err_out = TRUE); errorF.plot


# learn and visualize error rates of R reads
errorR <- learnErrors(filteredRs, multithread = TRUE)
errorR.plot <- plotErrors(errorR, nominalQ = TRUE, err_in = TRUE, err_out = TRUE); errorR.plot


# save error plots
#ggsave("./Plots/errorF.pdf", plot = errorF.plot, device = "pdf", width = 12, height = 8, units = "in")
#ggsave("./Plots/errorR.pdf", plot = errorR.plot, device = "pdf", width = 12, height = 8, units = "in")




################################
### Merging paired-end reads ###
################################

# create list of merged reads
mergers <- vector("list", length(sample_names))
names(mergers) <- sample_names


# sample inference and merging paired-end reads
for(i in sample_names) {
  cat("\nProcessing", i, "(", match(i, sample_names),"/", 
      length(sample_names), ") -", format(Sys.time(), "%I:%M %p"), "\n")
  derepF <- derepFastq(filteredFs[[i]])
  dadaF <- dada(derepF, err = errorF, multithread = TRUE)
  derepR <- derepFastq(filteredRs[[i]])
  dadaR <- dada(derepR, err = errorR, multithread = TRUE)
  merger <- mergePairs(dadaF, derepF, dadaR, derepR)
  mergers[[i]] <- merger
}


# remove dereps to save memory
rm(derepF, derepR)


# construct a sequence table: number of unique sequences (ASVs) in each sample
counts_table <- makeSequenceTable(mergers)
cat(ncol(counts_table), "ASVs in", nrow(counts_table), "samples\n")
# truncated: 7678 ASVs in 24 samples
# not truncated: 7193 ASVs in 24 samples




##################################
### Quality control: processed ###
##################################


########
### trim

# inspect sequence distribution
seq_distribution <- data.frame(table(nchar(getSequences(counts_table)))) 
colnames(seq_distribution) <- c("seqLength", "frequency")
seq_distribution$SeqLength <- as.numeric(as.character(seq_distribution$SeqLength))


# get stats
med_len <- median(nchar(getSequences(counts_table)))
min_len <- floor(median(nchar(getSequences(counts_table))) * 0.95)
max_len <- ceiling(median(nchar(getSequences(counts_table))) * 1.05)

lens <- data.frame(stat = c("median", "min", "max"),
                   x = c(med_len, min_len, max_len),
                   y = max(seq_distribution$Frequency) / 2)


# visualize and save
dist.plot <- 
  ggplot(seq_distribution, aes(x = seqLength, y = frequency)) +
  geom_bar(stat = "identity", width = 0.5, fill = "black") +
  geom_text(data = lens, aes(x = x, y = y, label = stat), vjust = 2, color = "red") +
  geom_vline(data = lens, aes(xintercept = x), linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = ceiling(quantile(seq_distribution$seqLength)),
                     labels = ceiling(quantile(seq_distribution$seqLength)))  +
  labs(x = "Sequence length (bp)", y = "Frequency")
dist.plot

#ggsave("./Plots/seq_distribution.pdf", plot = dist.plot, device = "pdf", width = 16, height = 4, units = "in")


# remove reads of non target length, 5% above and below the median 
counts_trimmed <- counts_table[ ,nchar(colnames(counts_table)) 
                                     %in% seq(min_len, max_len)]


cat(ncol(counts_trimmed), "ASVs in", nrow(counts_trimmed), "samples\n")
# with truncating:
# without adjusting max length: 6127 ASVs in 24 samples
# with adjusting max length: 6910 ASVs in 24 samples
# without truncating: 6098 ASVs in 24 samples



###################
### Remove chimeras

# removing chimeras with denovo screening
counts_noChim <- removeBimeraDenovo(counts_trimmed, method = "consensus",
                                            multithread = TRUE, verbose = TRUE)

# how many unique sequences were moved?
cat(ncol(counts_noChim), "out of", ncol(counts_trimmed), "ASVs passed (", 
    (1 - ncol(counts_noChim) / ncol(counts_trimmed)) * 100, "% removed )\n")
# 5690 out of 6098 ASVs passed ( 6.690718 % removed )


# what percentage of reads were identified as chimeras?
cat((1 - sum(counts_noChim) / sum(counts_table)) * 100, "% reads removed")
# 2.715508 % reads removed


# what is leftover?
cat(ncol(counts_noChim), "ASVs in", nrow(counts_noChim), "samples\n")
# 5690 ASVs in 24 samples




#######################
### Assign taxonomy ###
#######################

taxa_table <- assignTaxonomy(counts_noChim, 
                                      "~/Desktop/Lab/Projects/Misc/silva_nr_v132_train_set.fa.gz", 
                                      multithread = TRUE)


# good time to save
save.image("./RData/Gonzalez_dada2_env.RData")
 



##################################
### remove mock community ASVs ###
##################################

# create fasta
uniquesToFasta(counts_noChim, ids = paste0("nochim", sprintf("%06d", 1:ncol(counts_noChim))),
               fout = "./RData/Seqtabnochim.fasta")


# create taxonomy data frame
taxa_table.df <- data.frame(taxa_table)
taxa_table.df$sseqid <- paste0("nochim", sprintf("%06d", 1:ncol(counts_noChim)))
taxa_table.df$FASTA <- rownames(taxa_table.df)


# in terminal:
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# $ makeblastdb -dbtype nucl -in Seqtabnochim.fasta -input_type fasta
#
# $ blastn -db Seqtabnochim.fasta -query ~/Desktop/Lab/Projects/Misc/zymo_mock_all.fasta -task blastn -perc_identity 100 -outfmt 6 -out mock_align.txt
#
# $ echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n$(cat mock_align.txt)" > mock_align.txt
#
# $ sed 's/\t/,/g' mock_align.txt > mock_align.csv
#
# $ rm Seqtabnochim.fasta.n*
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# load
mock_align <- read.csv("RData/mock_align.csv")


# keep only best matches
mock_align <- subset(mock_align, bitscore >= 400)


# compare
mock_ASVs <- unique(merge(taxa_table.df, mock_align[1:2], by = "sseqid"))


# remove
counts_noMock <- counts_noChim[,! colnames(counts_noChim) %in% mock_ASVs$FASTA]
ncol(counts_noChim) - ncol(counts_noMock) == length(mock_ASVs$FASTA)


# save
mock_ASVs <- subset(taxa_table.df, FASTA %in% mock_ASVs$FASTA)
mock_ASVs$Source <- "Mock"




#######################################
### remove no template control ASVs ###
#######################################

sample_info <- data.frame(Sample_name = rownames(counts_noMock))
rownames(sample_info) <- sample_info$Sample_name


# create phyloseq object
counts_object <- phyloseq(otu_table(counts_noMock, taxa_are_rows = FALSE),
                           tax_table(taxa_table),
                           sample_data(sample_info))


# add logical variable to sample info stating if it's negative control
sample_data(counts_object)$NTC <- sample_data(counts_object)$Sample_name == "NTC"


# run check for contaminants using prevalence method and negative control
contam_prev <- isContaminant(counts_object, method = "prevalence", neg = "NTC")


# how many are contaminants?
table(contam_prev$contaminant)
# FALSE  TRUE 
# 5676     4 


# what taxa are they?
contam_prev$FASTA <- rownames(contam_prev)
NTC_ASVs <- subset(taxa_table.df, FASTA %in% subset(contam_prev, contaminant == TRUE)$FASTA)
NTC_ASVs$Source <- "NTC"



###############################
### remove nonspecific ASVs ###
###############################

# subset eukaryotes, mitochondria, chloroplasts
euk_ASVs <- subset(taxa_table.df, Kingdom == "Eukaryota")
chloro_ASVs <- subset(taxa_table.df, Order == "Chloroplast")
mito_ASVs <- subset(taxa_table.df, Family == "Mitochondria")


# combine them
euk_ASVs$Source <- "Eukaryota"
chloro_ASVs$Source <- "Chloroplasts"
mito_ASVs$Source <- "Mitochondria"
contaminants <- rbind(mock_ASVs, NTC_ASVs, euk_ASVs, mito_ASVs, chloro_ASVs)


# remove from data
counts_noContam <- counts_noMock[,! colnames(counts_noMock) %in% contaminants$FASTA]
ncol(counts_noChim) - ncol(counts_noContam) == length(unique(contaminants$FASTA))
counts_noContam <- counts_noContam[! rownames(counts_noContam) %in% c("NTC", "Mock"),]

taxa_noContam.df <- subset(taxa_table.df, ! FASTA %in% contaminants$FASTA)
taxa_noContam.df <- data.frame(FASTA = taxa_noContam.df$FASTA,
                               ASV = paste0("ASV", sprintf("%06d", 1:ncol(counts_noContam))), 
                               taxa_noContam.df[1:6])
identical(rownames(taxa_noContam.df), colnames(counts_noContam))




#########################
### Organize and save ### 
#########################

# save FASTA
uniquesToFasta(counts_noContam,
               ids = paste0(taxa_noContam.df$ASV, "__",
               taxa_noContam.df$Kingdom, "__",
               taxa_noContam.df$Phylum, "__",
               taxa_noContam.df$Class, "__",
               taxa_noContam.df$Order, "__",
               taxa_noContam.df$Family, "__",
               taxa_noContam.df$Genus),
               fout = "./RData/Gonzalez_v4_dada2.fasta")


# transpose
counts_noContam.df <- data.frame(t(counts_noContam))
counts_noContam.df <- data.frame(FASTA = rownames(counts_noContam.df), counts_noContam.df)
all <- merge(taxa_noContam.df, counts_noContam.df, by = "FASTA")
all <- all[order(all$ASV),]


# save
write.csv(all, "./RData/Gonzalez_v4_dada2.csv", row.names = FALSE, na = "")



# # # # # # # # # # #
# save R environment
save.image("./RData/Gonzalez_dada2_env.RData")



