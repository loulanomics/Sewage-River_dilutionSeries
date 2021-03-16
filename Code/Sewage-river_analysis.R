
########################################
### V4 analysis for Abe Gonzalez
### Lou LaMartina, started March 15 2021
########################################


library(ggplot2)
library(ape)
library(vegan)
library(reshape2)
library(phyloseq)
library(RColorBrewer)

setwd("~/Desktop/Lab/Projects/Gonzalez")


#################
### prep data ###
#################

# load counts and taxonomy
data <- read.csv("./RData/Gonzalez_v4_dada2.csv")
rownames(data) <- data$Sample_name


# load sample info
sample_info <- read.csv("./RData/Gonzalez_sample_info.csv")
rownames(sample_info) <- sample_info$Sample_name


# subset counts
counts <- data[c(2, 9:ncol(data))]
rownames(counts) <- counts$ASV
counts <- counts[-1]
counts <- counts[rowSums(counts) > 0, colSums(counts) > 0]
counts <- data.frame(t(counts))


# subset taxonomy
taxa <- data[1:8]
rownames(taxa) <- taxa$ASV
taxa <- subset(taxa, ASV %in% colnames(counts))


# convert to relative abundance
relabun <- counts / rowSums(counts)
relabun[1:5,1:5]




##################
### ordination ###
##################

# stat
PCoA.stat <- pcoa(vegdist(relabun, method = "bray"))


# extract values
PCoA.df <- data.frame(PCoA.stat$vectors[,1:2])
axis_scores <- PCoA.stat$values$Relative_eig[1:2] * 100


# add info
PCoA.df$Sample_name <- rownames(PCoA.df)
PCoA.df <- merge(PCoA.df, sample_info, by = "Sample_name")


# plot
pcoa.plot <- 
  ggplot(PCoA.df, aes(x = Axis.1, y = Axis.2, color = as.factor(Dilution))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", size = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80", size = 0.2) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        axis.title.x = element_text(size = 8, color = "black", face = "bold"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 8, color = "black", face = "bold"),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25)) +
  guides(color = guide_legend(keyheight = 0.4, keywidth = 0.4, units = "in", ncol = 1),
         shape = guide_legend(keyheight = 0.4, keywidth = 0.4, units = "in", ncol = 1)) +
  labs(x = "Axis 1\n(92.4%)", y = "Axis 2\n(2.12%)", color = "Dilution")

#ggsave("./Plots/PCoA.pdf", pcoa.plot, device = "pdf", width = 6, height = 3)




#############################################
### combine taxa to lowest classification ###
#############################################

classes <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")


# create list with ASVs that are unclassified at each level
unknownASVs.ls <- list()
for (i in 1:6){
  if (i == 1) {
    unknownASVs.ls[[i]] <- taxa[taxa[[2 + i]] == "", "ASV"]
  } else if (i > 1) {
    unknownASVs.ls[[i]] <- setdiff(taxa[taxa[[2 + i]] == "", "ASV"], taxa[taxa[[2 + i - 1]] == "", "ASV"])
  }
}
names(unknownASVs.ls) <- classes


# make their new names their lowest classification.
# for example, if an ASV is not classified down to genus, 
# its name will be "Family: " and whatever its family is
classes.ls <- list()
for (i in 1:6){
  if (i == 1) {
    classes.ls[[i]] <- data.frame(ASV = subset(taxa, ASV %in% unknownASVs.ls[[i]])$ASV,
                               Name = "Unclassified")    
  } else if (i > 1) {
    classes.ls[[i]] <- data.frame(ASV = subset(taxa, ASV %in% unknownASVs.ls[[i]])$ASV,
                       Name = paste0(classes[i - 1], ": ", subset(taxa, ASV %in% unknownASVs.ls[[i]])[[i + 1]]))
  }
}
names(classes.ls) <- classes


# combine into single data frame
classes.df <- dplyr::bind_rows(classes.ls)


# add this to taxa data
taxa <- merge(classes.df, taxa, by = "ASV", all.y = TRUE)
rownames(taxa) <- taxa$ASV


# for those classified down to genus, give them genus names
taxa$Name[is.na(taxa$Name)] <- paste0("Genus: ", taxa[is.na(taxa$Name), "Genus"])
length(unique(taxa$Name))


# aggregate by these names with phyloseq
object <- phyloseq(otu_table(as.matrix(relabun), taxa_are_rows = FALSE),
                   tax_table(as.matrix(data.frame(taxa[2]))))
glom <- speedyseq::tax_glom(object, "Name")


# extract data
relabun.glom <- data.frame(glom@otu_table@.Data)
tax.glom <- data.frame(glom@tax_table@.Data)
tax.glom$ASV <- rownames(tax.glom)


# get most abundant taxa
top11 <- names(sort(colSums(relabun.glom), decreasing = TRUE)[1:11])


# add to taxa df
tax.glom$Top <- tax.glom$Name
tax.glom$Top[! tax.glom$ASV %in% top11] <- "Other"


# add sample info
relabun.glom <- data.frame(Sample_name = rownames(relabun.glom), relabun.glom)
relabun.glom <- merge(sample_info, relabun.glom, by = "Sample_name")


# melt and combine
relabun.glom.m <- melt(relabun.glom, id.vars = colnames(sample_info), variable.name = "ASV", value.name = "Relabun")
relabun.glom.m <- relabun.glom.m[relabun.glom.m$Relabun > 0,]
relabun_tax.glom <- merge(relabun.glom.m, tax.glom, by = "ASV")


# plot
bars.plot <-
  ggplot(relabun_tax.glom, aes(x = as.factor(Dilution), y = Relabun, fill = Top, color = Top)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c(brewer.pal(11, "Paired"), "grey90")) +
  scale_color_manual(values = c(brewer.pal(11, "Paired"), "grey90"), guide = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        axis.title.x = element_text(size = 8, color = "black", face = "bold"),
        legend.title = element_text(size = 8, color = "black", face = "bold"),
        legend.text = element_text(size = 7, color = "black", face = "italic"),
        panel.border = element_blank(), 
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25)) +
  guides(fill = guide_legend(keyheight = 0.8, keywidth = 0.6, units = "in")) +
  labs(y = "Relative abundances of taxa", x = "Dilution", 
       fill = "Lowest\ntaxonomic\nclassification")
bars.plot

#ggsave("./Plots/bars.pdf", bars.plot, device = "pdf", width = 5, height = 4)

