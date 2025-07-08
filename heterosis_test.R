library(dplyr)
library(phyloseq)
library(FSA)
library(rcompanion)
library(vegan)
library(viridis)
library(ggplot2)
library(DESeq2)



bac.ps = readRDS("H:/Research Projects/heterosisProject/intermediate_data/16S/bac_asv_table_phyloseq_prefiltered.RDS")
# 31908 taxa and 3168 samples

# extract inbred line and hybrids samples
bac.ps.hyin = subset_samples(bac.ps, Germplasm %in% c('Hybrid', 'Inbred_line'))  #504 samples
bac.ps.hyin = filter_taxa(bac.ps.hyin, function(x) sum(x >= 10) >= 2, TRUE) 
table(tax_table(bac.ps.hyin)[, 'Phylum'])

# separate root samples and rhizosphere samples, do the following analysis separately
bac.ps.hyin.root = subset_samples(bac.ps.hyin, Compartment=='Root') 

# keep ASVs that have relative abundance > 0.05% at least in 20% samples
bac.ps.hyin.root.RA = transform_sample_counts(bac.ps.hyin.root, function(x) x/sum(x))
bac.ps.hyin.root.RA.abund = filter_taxa(bac.ps.hyin.root.RA, 
                                      function(x) sum(x > 0.0005) >= 0.2*nsamples(bac.ps.hyin.root.RA), TRUE)  
bac.ps.hyin.root.abund = subset_taxa(bac.ps.hyin.root, taxa_names(bac.ps.hyin.root)%in%taxa_names(bac.ps.hyin.root.RA.abund))

table(colSums(otu_table(bac.ps.hyin.root.abund)) >1000 )  # all true

##################################################
##               alpha diversity                ##
##################################################


# for root samples
pdf('rarecurve_root.pdf', width = 10, height = 10)
rarecurve(t(otu_table(bac.ps.hyin.root.abund)), step=100, label = F, col = "blue")
dev.off()

# rarefaction
bac.ps.hyin.root.abund.rare = rarefy_even_depth(bac.ps.hyin.root.abund, 1000, rngseed = 123555, replace = FALSE)

#Generate a data.frame with adiv measures
adiv.root = data.frame(
  "ObservedOTU" = phyloseq::estimate_richness(bac.ps.hyin.root.abund.rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(bac.ps.hyin.root.abund.rare, measures = "Shannon"),
  "ACE" = phyloseq::estimate_richness(bac.ps.hyin.root.abund.rare, measures = "ACE")
)
adiv.root = cbind(adiv.root, sample_data(bac.ps.hyin.root.abund.rare))


#################################################
##             Beta diversity                  ##
#################################################

# DESeq2 Normalization 

sample_data(bac.ps.hyin.root.abund)$batch = paste0(sample_data(bac.ps.hyin.root.abund)$Harvest.batch, 
                                                   sample_data(bac.ps.hyin.root.abund)$Treatment.batch)
bac.hyin.root.deseq = phyloseq_to_deseq2(bac.ps.hyin.root.abund, ~batch+Germplasm)
bac.hyin.root.deseq = estimateSizeFactors(bac.hyin.root.deseq, type="poscount")
bac.hyin.root.vst = varianceStabilizingTransformation(bac.hyin.root.deseq, blind=F)
bac.ps.hyin.root.abund.vst = bac.ps.hyin.root.abund
otu_table(bac.ps.hyin.root.abund.vst) = otu_table(assay(bac.hyin.root.vst), taxa_are_rows = T)

# Round negative values up to zeroes, to enable Bray-Curtis calculations
bac.ps.hyin.root.abund.vst = transformSampleCounts(bac.ps.hyin.root.abund.vst,function(x) ifelse(x<0,0,x)) 

# Constrained principal coordinate analysis
cap.bray.bac.root = ordinate(bac.ps.hyin.root.abund.vst, method = "CAP", distance = "bray", 
                           formula = ~Plant.Name+Condition(Harvest.batch:Treatment.batch+Block))
score.sites.root = scores(cap.bray.bac.root, display="sites",choice=c(1:5, 21:25))
beta.root = cbind.data.frame(score.sites.root, as.data.frame(sample_data(bac.ps.hyin.root.abund.vst)))

screeplot(cap.bray.bac.root, type="line")
eigenvals(cap.bray.bac.root) %>% summary() -> ev
percentage.bac = paste( c("CAP1", 'CAP2', 'CAP3', 'CAP4', "MDS1"), "(", paste( as.character(round(ev[2,c(1:5)]*100,2)), "%", ")", sep=""), sep = "" )

# Constrained principal coordinate analysis plot
temp = ggplot(data= beta.root, aes(x = CAP1, y = MDS1, color=Germplasm)) +
  geom_point(alpha=0.7, size=2) +
  #scale_color_viridis(discrete=TRUE, option = "plasma") +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values=c(3, 16, 17))+
  stat_ellipse(aes(color=Germplasm), linetype = 2) +
  mytheme +
  xlab(percentage.bac[1]) + ylab(percentage.bac[5]) +
  # scale_x_continuous(limits=c(-3,3), breaks = seq(-3, 3, 0.5)) +
   scale_y_continuous(limits=c(-1,1), breaks = seq(-1, 1, 0.5)) +
  coord_fixed() 


##  permANOVA test ##

adonis(t(as(otu_table(bac.ps.hyin.root.abund.vst),'matrix'))~Harvest.batch:Treatment.batch+Block+Treatment*Germplasm, 
        data = as(sample_data(bac.ps.hyin.root.abund.vst), "data.frame"), 
        strata = as(sample_data(bac.ps.hyin.root.abund.vst), "data.frame")$Treatment,
       parallel = 4, by="margin")

#############################################################
##                 heterosis test                          ##
#############################################################

## heterosis t-tests for each ASV in each individual crosses 

heterosis.test = function(hb.lev, vars.names, lmm.red ){
  
  # heterosis test for each hybrid-parents
  for (h in hb.lev){ 
    mat = strsplit(h,'x')[[1]][1]   # first parent is maternal line
    pat = strsplit(h,'x')[[1]][2]   # second parent is paternal line
    # do this for each ASV
    for (v in vars.names) { 
      
      # Subset data for the parents and their hybrids
      print(v)
      mat.data = lmm.red[lmm.red$Plant.Name==mat, v] 
      pat.data = lmm.red[lmm.red$Plant.Name==pat, v]
      hyb.data = lmm.red[lmm.red$Plant.Name==h, v]
      
      # Calculate mean and SE for hybrid, maternal, and paternal
      hyb.Mean = mean(hyb.data, na.rm=TRUE) 
      hyb.SEM = sd(hyb.data,na.rm=TRUE)/sqrt(length(hyb.data[!is.na(hyb.data)]))
      mat.Mean = mean(mat.data, na.rm=TRUE) 
      mat.SEM = sd(mat.data, na.rm=TRUE)/sqrt(length(mat.data[!is.na(mat.data)]))
      pat.Mean = mean(pat.data, na.rm=TRUE)
      pat.SEM = sd(pat.data, na.rm=TRUE)/sqrt(length(pat.data[!is.na(pat.data)]))
      
      # Calculate mid-parent value and identify extreme parent
      midparent = sum(mat.Mean, pat.Mean)/2   
      highparent = max(mat.Mean, pat.Mean) 
      lowparent = min(mat.Mean, pat.Mean) 
      
      # t-test for better-parent heterosis 
      if (hyb.Mean > highparent){
        BetterParent = highparent
        tt.BPH = t.test(hyb.data, mu=highparent, alternative='greater')
        # save results
        tt.results = data.frame('Maternal'=mat,'Paternal'=pat,'Hybrid'=h,
                                'MatMean'=mat.Mean,'PatMean'=pat.Mean,
                                'MatSE'=mat.SEM,'PatSE'=pat.SEM,'HybMean'=hyb.Mean,
                                'HybSE'=hyb.SEM,
                                'Midparent'=midparent,'BetterParent'=BetterParent, 'Trait'=v,
                                'MPH.t.stat'=NA,'MPH.t.df'=NA, 'MPH.p'=NA,
                                'BPH.t.stat'=tt.BPH$statistic,'BPH.t.df'=tt.BPH$parameter, 'BPH.p'=tt.BPH$p.value,
                                'MPH.pct'=NA, 'BPH.pct' = (hyb.Mean-BetterParent)/BetterParent) 
      }
      else if (hyb.Mean < lowparent){
        BetterParent = lowparent
        tt.BPH = t.test(hyb.data, mu=lowparent, alternative='less')
        # save results
        tt.results = data.frame('Maternal'=mat,'Paternal'=pat,'Hybrid'=h,
                                'MatMean'=mat.Mean,'PatMean'=pat.Mean,
                                'MatSE'=mat.SEM,'PatSE'=pat.SEM,'HybMean'=hyb.Mean,
                                'HybSE'=hyb.SEM,
                                'Midparent'=midparent,'BetterParent'=BetterParent, 'Trait'=v,
                                'MPH.t.stat'=NA,'MPH.t.df'=NA, 'MPH.p'=NA,
                                'BPH.t.stat'=tt.BPH$statistic,'BPH.t.df'=tt.BPH$parameter, 'BPH.p'=tt.BPH$p.value,
                                'MPH.pct'=NA, 'BPH.pct' = (hyb.Mean-BetterParent)/BetterParent) 
      }
      else if (hyb.Mean < highparent & hyb.Mean > lowparent) {
        # t-test for mid-parent heterosis
        tt.MPH = t.test(hyb.data, mu=midparent, alternative='two.sided') 
        # save results
        tt.results = data.frame('Maternal'=mat,'Paternal'=pat,'Hybrid'=h,
                                'MatMean'=mat.Mean,'PatMean'=pat.Mean,
                                'MatSE'=mat.SEM,'PatSE'=pat.SEM,'HybMean'=hyb.Mean,
                                'HybSE'=hyb.SEM,
                                'Midparent'=midparent,'BetterParent'=BetterParent, 'Trait'=v,
                                'MPH.t.stat'=tt.MPH$statistic,'MPH.t.df'=tt.MPH$parameter, 'MPH.p'=tt.MPH$p.value,
                                'BPH.t.stat'=NA,'BPH.t.df'=NA, 'BPH.p'=NA,
                                'MPH.pct'=(hyb.Mean-midparent)/midparent, 'BPH.pct' = NA) 
      }
      if (h==hb.lev[1] & v==vars.names[1]) { heterosis = tt.results} 
      else { heterosis = rbind(heterosis, tt.results) } 
    }
    
  }
  return(heterosis)
}



## Prepare data for heterosis test ##

bac.root.heterosis.df =  merge(beta.root[, 1:11], t(otu_table(bac.ps.hyin.root.abund.vst)), by.x='Sample.ID', by.y='row.names') %>%
  merge(., adiv.root[, c(1:3,5)], by='Sample.ID') %>%
  merge(., as(sample_data(bac.ps.hyin.root.abund.vst), 'data.frame'), by='Sample.ID')

vars.names.root = colnames(bac.root.heterosis.df)[c(2:179)]

# Remove treatment and batch effects using linear mixed model
for (var in vars.names.root){
  print(var)
  mod1 = lmer(as.formula(paste0(var,'~Treatment+(1|batch)+(1|Block)')),
               data=filter(bac.root.heterosis.df,!is.na(var)))
  mod1.res = resid(mod1)
  names(mod1.res) = bac.root.heterosis.df$Sample.ID
  residual.df = cbind.data.frame(mod1.res, names(mod1.res))
  colnames(residual.df) = c(var, 'Sample.ID')
  if (var == vars.names.root[1]){
    heterosis.res = residual.df
  }
  else{
    heterosis.res = merge(heterosis.res, residual.df, by='Sample.ID')
  }
}

heterosis.res = merge(heterosis.res, as(sample_data(bac.ps.hyin.root.abund.vst), 'data.frame')[, 1:4], by='Sample.ID')
  
# identify hybrids genotypes
hybrids = levels(factor(filter(bac.root.heterosis.df, Germplasm=='Hybrid')$Plant.Name))

het.res.root = heterosis.test(hybrids, vars.names.root, heterosis.res )  
  
# Adjust p-values to correct for multiple tests using BH method
heterosis.results.root = group_by(het.res.root, Hybrid) %>% mutate(BPH.padj = p.adjust(BPH.p, method='BH'),
                                                            MPH.padj = p.adjust(MPH.p, method='BH')) %>% 
  ungroup %>% as.data.frame %>%
  mutate(HetType=case_when(
    BPH.padj < 0.05 ~ 'BPH',
    MPH.padj < 0.05 ~ 'MPH',
    MPH.padj > 0.05 & BPH.padj > 0.05 ~ 'none')) %>% 
  mutate(HetType=factor(HetType,levels=c('none','MPH','BPH')))
# Keep only significant heterosis results
heterosis.sig.root = heterosis.results.root[heterosis.results.root$HetType != 'none', ]

# left join the two data frames
res.root = merge(heterosis.sig.root, as.data.frame(tax_table(bac.ps.hyin.root.abund)), all.x=T, by.x='Trait', by.y='row.names')
write.table(res.root, 'bac.root.heterosis.significant.res.txt', sep = '\t', quote = F, row.names = F)

