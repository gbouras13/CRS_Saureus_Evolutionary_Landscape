library(tidyverse)
library(gggenomes)
library(ggplot2)







#### 420 #####

path_dir <-  "../Structural_Locus_Deep_Dive"


# need to combine the fastas and index using samtools faids
# parse sequence length from fasta.fai file


sdrd_seqs <- read_fai(file.path(path_dir, "R_fastas/420_combo.fasta.fai")) 
                      
sdrd_seqs$seq_desc <- sdrd_seqs$seq_id


# read in teh combo gff
sdrd_genes <- read_gff(file.path(path_dir,"R_gffs/420_combo.gff"))              # per gene GC-content

# read in the linkds
# All-vs-all alignment | https://github.com/lh3/minimap2
# minimap2 -X -N 1 -p 0.1 -c 420_combo.fasta 420_combo.fasta > 420_combo.paf

sdrd_links <- read_paf(file.path(path_dir,"R_fastas/420_combo.paf"))
sdrd_links <- sdrd_links %>% 
  filter(start %in% c(1, 7149), end %in% c(15822, 4612))

# set the real endpoints of sdrC sdrD recombination based on igv
sdrd_links$end[2] <- 3207
sdrd_links$end2[2] <- 3207

sdrd_links$start[1] <- 7845
sdrd_links$start2[1] <- 3208



# get cluster of genes 
sdrd_cluster <- sdrd_genes %>%  dplyr::select(ID, gene)
sdrd_cluster$gene <- ifelse(is.na(sdrd_cluster$gene),  
                            sdrd_cluster$ID, 
                            sdrd_cluster$gene)

colnames(sdrd_cluster) <- c("feat_id", "cluster_id")


sdrd_cluster$cluster_id <- as.factor(sdrd_cluster$cluster_id)


# match with feat_id
sdrd_genes$feat_id <- sdrd_genes$ID


sdr_clust<- sdrd_genes %>%  filter(gene %in% c("sdrC", "sdrD", "sdrE") )

sdr_clust$end <- sdrd_genes%>%  filter(gene %in% c("sdrC", "sdrD", "sdrE") ) %>% pull(end)


sdrd_seqs$seq_id
sdrd_genes$seq_id
sdrd_links$seq_id


# finally draw plot
p <-  gggenomes(seqs=sdrd_seqs, genes=sdrd_genes, links=sdrd_links)

p+
  geom_seq()  +
  geom_gene() +
  geom_link()  

p2 <- p  %>%
  add_clusters(sdrd_cluster) %>% 
  add_feats(sdr_clust=sdr_clust) +
  geom_seq()  +
  geom_gene() +
  geom_link() + 
  geom_bin_label(size = 8) +
  geom_gene(aes(fill=cluster_id))+
  scale_fill_brewer("Gene", palette="Set3")+
  geom_feat(data=feats(sdr_clust), alpha=0.01, size=5, position="identity")


p2



png("Figures/420_sdrd_deletion_gggenomes.png",width=8.5, height=5, unit="in", res=600,type='cairo' )
p2+  theme_classic(base_size = 20)+
  labs(x = "Length (bp)", y= "")+
  theme(axis.line.y =element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y =  element_blank())

dev.off()









#### 4875 #####


# need to combine the fastas and index using samtools faids
# parse sequence length from fasta.fai file
blaZ_seqs <- read_fai(file.path(path_dir,"R_fastas/4875.combo.fasta.fai")) 

blaZ_seqs$seq_desc <- blaz_seqs$seq_id


# read in teh combo gff
blaZ_genes <- read_gff(file.path(path_dir,"R_gffs/4875_combo.gff"))              # per gene GC-content
blaZ_genes$gene <- ifelse(blaZ_genes$gene == "blaPC1", "blaZ", blaZ_genes$gene)

# read in the linkds
# All-vs-all alignment | https://github.com/lh3/minimap2
# minimap2 -X -N 1 -p 0.1 -c 4875.combo.fasta 4875.combo.fasta > 4875_combo.paf
blaZ_links <- read_paf(file.path(path_dir,"R_fastas/4875_combo.paf"))





# get cluster of genes 
blaZ_cluster <- blaZ_genes %>%  dplyr::select(ID, gene)
blaZ_cluster$gene <- ifelse(is.na(blaZ_cluster$gene),  
                            blaZ_cluster$ID, 
                            blaZ_cluster$gene)
# cluster blaZ, other mge 
'%ni%' <- Negate('%in%')
blaZ_cluster$gene <- ifelse(blaZ_cluster$gene  %ni% c("tra5", "resR","blaI", "blaR1", "blaZ"),
                            "Other Genes", blaZ_cluster$gene)





colnames(blaZ_cluster) <- c("feat_id", "cluster_id")


blaZ_cluster$cluster_id <- as.factor(blaZ_cluster$cluster_id)


# match with feat_id
blaZ_genes$feat_id <- blaZ_genes$ID


blaZ_clust<- blaZ_genes %>%  filter(gene %in% c("blaI", "blaR1", "blaZ") )

blaZ_clust$end <- blaZ_genes%>%  filter(gene %in% c("blaI", "blaR1", "blaZ") ) %>% pull(end)


blaZ_seqs$seq_id
blaZ_genes$seq_id
blaZ_links$seq_id


# finally draw plot
p <-  gggenomes(seqs=blaZ_seqs, genes=blaZ_genes, links=blaZ_links)

p+
  geom_seq()  +
  geom_gene() +
  geom_link()  

p2 <- p  %>%
  add_clusters(blaZ_cluster) %>% 
  add_feats(blaZ_clust=blaZ_clust) +
  geom_seq()  +
  geom_gene() +
  geom_link() + 
  geom_bin_label(size = 8) +
  geom_gene(aes(fill=cluster_id))+
  #scale_fill_brewer("Gene", palette="Set3")+
  geom_feat(data=feats(blaZ_clust), alpha=0.01, size=5, position="identity")


p2





### 10 highest and lowest

png("Figures/4875_blaZ_deletion_gggenomes.png",width=12, height=5, unit="in", res=600,type='cairo' )
p2+  theme_classic(base_size = 20)+
  labs(x = "Length (bp)", y= "")+
  theme(axis.line.y =element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y =  element_blank())

dev.off()

