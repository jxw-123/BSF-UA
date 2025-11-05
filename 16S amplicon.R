# 加载必要的包
library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)




seqtab=ASV
sample.names=c("Cont1","Cont2","Cont3","UA1","UA2","UA3")

seqtab=seqtab[,-1]
rownames(seqtab)=sample.names
seqtab.nochim=as.matrix(seqtab)

samdf=data.frame(rownames(seqtab.nochim),c("Cont","Cont","Cont","UA","UA","UA"));
rownames(samdf)<- sample.names
colnames(samdf)<-c("id","group")
samdf;

taxa0=taxa
ta=taxa$...1

taxa0=as.matrix(taxa0)
rownames(taxa0)=ta


seqtab.nochim0<-seqtab.nochim
colnames(seqtab.nochim0)<-paste(rep("OTU",ncol(seqtab.nochim)),1:ncol(seqtab.nochim),sep="");
taxa0<-taxa0
rownames(taxa0)<-paste(rep("OTU",ncol(seqtab.nochim)),1:ncol(seqtab.nochim),sep="");

ps <- phyloseq(otu_table(seqtab.nochim0, taxa_are_rows=FALSE),sample_data(samdf),tax_table(taxa0));


bray_dist <- distance(ps, method = "bray")


pcoa_result <- ordinate(ps, method = "PCoA", distance = "bray")


plot_ordination(ps, pcoa_result, color = "group") +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +  # 添加95%置信椭圆
  scale_color_manual(values = c("Cont" = "blue", "UA" = "red")) +
  labs(title = "PCoA Plot (Bray-Curtis Distance)",
       color = "Group") +
  theme_bw()


estimate_richness(ps, measures =c("Chao1", "Shannon"));
plot_richness(ps);

plot_richness(ps,measures = "Chao1");
plot_richness(ps,measures = "Shannon");




ps_family <- tax_glom(ps, "Family", NArm = FALSE)


ps_family_rel <- transform_sample_counts(ps_family, function(x) x / sum(x))


top20_families <- names(sort(taxa_sums(ps_family_rel), decreasing = TRUE)[1:20])


ps_top20 <- prune_taxa(top20_families, ps_family_rel)


all_taxa <- taxa_names(ps_family_rel)
other_taxa <- setdiff(all_taxa, top20_families)


if(length(other_taxa) > 0) {

  other_otu <- as(otu_table(ps_family_rel)[, other_taxa, drop = FALSE], "matrix")
  

  other_sums <- rowSums(other_otu)
  

  other_tax <- matrix(rep("Other", 7), nrow = 1)
  colnames(other_tax) <- colnames(tax_table(ps_family_rel))
  rownames(other_tax) <- "Other"
  

  other_otu_table <- matrix(other_sums, ncol = 1)
  colnames(other_otu_table) <- "Other"
  rownames(other_otu_table) <- sample_names(ps_family_rel)
  

  merged_otu <- cbind(otu_table(ps_top20), other_otu_table)
  merged_tax <- rbind(tax_table(ps_top20), other_tax)
  

  ps_merged <- phyloseq(
    otu_table(merged_otu, taxa_are_rows = FALSE),
    sample_data(ps_family_rel),
    tax_table(merged_tax)
  )
} else {
  ps_merged <- ps_top20
}


plot_bar(ps_merged, x = "id", fill = "Family") + 
  facet_wrap(~group, scales = "free_x") +
  ylab("Relative Abundance (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))









