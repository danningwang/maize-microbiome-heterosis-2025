library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(DESeq2)
library(vegan)
library(phyloseq)
library(dplyr)
library(reshape2)
library(biomformat)
library(qiime2R)
library(FSA)
library(lme4)
library(rcompanion)
library(viridis)

mytheme = theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size = 10),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 8)) 


setwd("~/heterosis")


#########################################################################
##                 import to phyloseq                                  ##
#########################################################################

# read asv biom table generated from QIIME pipeline
biom_table = read_biom("./data/ITS/feature-table.biom")
rm(biom_table)
feature.table = as.matrix(biom_data(biom_table))
dim(feature.table)                                             

# convert feature ID to ASV 
featureID2asv = cbind("FeatureID"=rownames(feature.table), "ASV"=paste("fASV", seq(1, nrow(feature.table)), sep = ""))
write.table(featureID2asv, "./intermediate_data/ITS/pre-analysis/fungi_feature_ID_to_ASV_file.txt", row.names = F, sep = "\t", col.names = T)

asv.table = feature.table
rownames(asv.table) = featureID2asv[,2]

# read taxonomy table generated from QIIME pipeline
taxonomy = read.table("./data/ITS/taxonomy.tsv", sep = "\t", header = T)
taxa.table = taxonomy %>%
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  mutate(Kingdom = gsub("k__", "", Kingdom),
         Phylum = gsub("p__", "", Phylum),
         Class = gsub("c__", "", Class),
         Order = gsub("o__", "", Order),
         Family = gsub("f__", "", Family),
         Genus = gsub("g__", "", Genus),
         Species = gsub("s__", "", Species))
table(rownames(feature.table) == taxa.table$Feature.ID)
rownames(taxa.table) = featureID2asv[, 2]
taxa.table = taxa.table[, -1]
taxa.table = as.matrix(taxa.table)

# read metadata file
md = read.csv("./data/16S/samples_metadata_new.csv", header = T)

# create phyloseq object
OTU.tab = otu_table(asv.table, taxa_are_rows = TRUE)
TAX.tab = tax_table(taxa.table)
md = md[match(colnames(OTU.tab), md$Sample.ID), ]
rownames(md) = md$Sample.ID
sample.tab = sample_data(md)
table(sample.tab$Sample.ID == colnames(OTU.tab))
ps = phyloseq(OTU.tab, TAX.tab, sample.tab)

# read phylogenetic tree file
rooted.tree = qza_to_phyloseq(tree = "./data/ITS/rooted-tree.qza")
table(rooted.tree$tip.label %in% rownames(feature.table))
# change tip label to ASV number
new_tiplabel = featureID2asv[match(rooted.tree$tip.label, featureID2asv[,1]), 1]
table(rooted.tree$tip.label == new_tiplabel)
new_tiplabel = featureID2asv[match(rooted.tree$tip.label, featureID2asv[,1]), 2]
rooted.tree$tip.label = new_tiplabel
ps.tree = merge_phyloseq(ps, phy_tree(rooted.tree))
ps.tree                                             
saveRDS(ps.tree, "./intermediate_data/ITS/pre-analysis/raw_fungi_asv_table_phyloseq.RDS")

###############################################################
##                  filter rare taxa                         ##
###############################################################

table(tax_table(ps.tree)[, "Phylum"])  # remove unidentified asv at phylum level
ps.tree.prefiltered = subset_taxa(ps.tree, Phylum!=' unidentified')

# function for filter taxa those express at least K reads in M*3 samples
expressKinM = function(K, M, otuTab){
  keep.taxa = c()
  for (i in 1:(ncol(otuTab)-1)) {
    print(i)
    if(sum(table(otuTab[otuTab[, i]>=K, ]$rep)>=3) >= M){
      keep.taxa = c(keep.taxa, colnames(otuTab)[i])
    }
  }
  return(keep.taxa)
}
# filter lowly-expressed taxa
tmp.otu = as.data.frame(cbind(t(otu_table(ps.tree.prefiltered)), "rep"=as.factor(rep(1:1056, each=3))))
keep.tax = expressKinM(10, 2, tmp.otu)  
ps.tree.filtered.abund = subset_taxa(ps.tree.prefiltered, taxa_names(ps.tree.prefiltered)%in%keep.tax) 


