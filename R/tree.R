library(ape)
library(phangorn)
library(ggtree)
library(tidyverse)

# load tree data and midpoint root 

tree <- read.tree("../poppunk/poppunk_viz/poppunk_viz_core_NJ.nwk")            # load in the tree
tree$edge.length <- pmax(tree$edge.length, 0.0)  # set any negative branch lengths to zero
tree <- midpoint(tree)                           # midpoint-root the tree


# load meta data and remove references 
meta <- read.csv(file = "../poppunk/poppunk_viz/poppunk_viz_microreact_clusters.csv" )
ref_ids <- meta %>%  filter(Status == "Reference") %>% pull(id)

tree <- drop.tip(tree, ref_ids)

meta <- meta %>%  filter(Status == "Query")


# load clinical metadata and merge 
meta1 <- read.csv(file = "../metadata/gess_time.csv")
meta2 <- read.csv(file = "../metadata/metadata_phylogentic_tree.csv")
merged_meta<- merge(x = meta1,
                    y = meta2,
                    by.x = "Sample", by.y = "Cnumber_ID",
                    all = TRUE)

meta<- merge(x = meta,
             y = merged_meta,
             by.x = "id", by.y = "Gess_plus_Time",
             all = TRUE)

rm(meta1, meta2, merged_meta)


# change tree tipname 

meta$tipname<- paste0( "Host:" ,meta$rid)

meta$tipname<- paste(meta$tipname, paste0("CI:", meta$Sample), sep = " | " )

oldtip<-data.frame(tree[["tip.label"]])


meta<-  meta %>% left_join(x = oldtip,
                           y = meta,
                           by = c("tree...tip.label..."= "id"))


tree[["tip.label"]]<- meta$tipname

meta<-meta %>%
  relocate(tipname, .before = 1)

rm(oldtip)


# make tree 

p <- ggtree(tree) %<+% meta + 
  geom_tippoint(aes(color=Time), size =2.5)+
  geom_tiplab(align=TRUE, linetype='dashed', linesize=.05, offset=.00001, colour= "black")+
  scale_color_manual(values = c("#ff005f","#520479"), name="Timepoint")+ 
  geom_treescale(x=0, y=65)+
  xlim(NA,0.0175)+
  ylim(NA,72)

p





########
##make colour pallete 
########


library(paletteer)

col1<- c(paletteer_d("colorBlindness::paletteMartin",15, type = "discrete"), paletteer_d("ggthemr::solarized",1))
col2<- paletteer_d("vapoRwave::hotlineBling",8, type = "discrete")
col2[8]<- "grey" 
col3<-paletteer_d("ggthemr::fresh", 2)
col4<-paletteer_d("ggthemes::hc_darkunica",2)


########################
### popunk & mlst  clusters
########################

########################
### POPUNK
########################

heatmapData=meta%>% 
  dplyr::select(Cluster_Cluster__autocolour)
rownames(heatmapData) <- meta$tipname

heatmapData$Cluster_Cluster__autocolour <- as.factor(heatmapData$Cluster_Cluster__autocolour)
colnames(heatmapData)[1] <- "PopPUNK" 


p2 <- gheatmap(p,
               heatmapData, 
               offset = 0.0015, 
               width = 0.05,
               colnames_position='top', 
               colnames_angle=45,
               colnames_offset_y = 0, 
               hjust=0, 
               font.size = 5,
               legend_title = "PopPUNK Cluster", )+
   scale_fill_manual(values=col1, name= "PopPUNK Cluster" ) 
  # theme_classic(base_size = 18)

p2






########################
# add CC
########################

CC=meta%>% 
  dplyr::select(CC)
rownames(CC) <- meta$tipname

colnames(CC)[1] <- "CC" 

CC$CC<- ifelse(CC$CC == "Not_Assigned","none assigned", CC$CC )



# need to do this for some reason
library(ggnewscale)


p3 <- p2 + new_scale_fill()

p3 <- gheatmap(p3,
               CC, 
               offset = 0.0021, 
               width = 0.05,
               colnames_position='top', 
               colnames_angle=45,
               colnames_offset_y = 0, 
               hjust=0, 
               font.size = 5,
               legend_title = "CC")+
  scale_fill_manual(values=col2, name="CC") 
  # theme_classic(base_size = 18)

p3

########################
# add phenotype
########################

crs_phen=meta%>% 
  dplyr::select(CRS_pheno) %>%  mutate_all(factor)

rownames(crs_phen) <- meta$tipname

colnames(crs_phen)[1] <- "phenotype" 



# need to do this for some reason
library(ggnewscale)


p4 <- p3 + new_scale_fill()

p4 <- gheatmap(p4,
               crs_phen, 
               offset = 0.0027, 
               width = 0.05,
               colnames_position='top', 
               colnames_angle=45,
               colnames_offset_y = 0, 
               hjust=0, 
               font.size = 5,
               legend_title = "CRS")+
  scale_fill_manual(values=col3, name="CRS") 
  # theme_classic(base_size = 18)

p4


########################
# add asthma
########################

asthma_stat=meta%>% 
  dplyr::select(asthma) %>%  mutate_all(factor)

rownames(asthma_stat) <- meta$tipname

colnames(asthma_stat)[1] <- "asthma" 

asthma_stat$asthma<- ifelse(asthma_stat$asthma == 1, "Asthmatic", "Non-asthmatic")



# need to do this for some reason
library(ggnewscale)


p5 <- p4 + new_scale_fill()


p5 <- gheatmap(p5,
               asthma_stat, 
               offset = 0.0033, 
               width = 0.05,
               colnames_position='top', 
               colnames_angle=45,
               colnames_offset_y = 0, 
               hjust=0, 
               font.size = 5,
               legend_title = "Asthma status")+
  scale_fill_manual(values=col4, name= "Asthma status")+
  theme_classic(base_size = 18) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())




p5



grDevices::dev.size("in")  

png(filename = paste0("Figures/","ggtree.png"), units = "in", width = 20, height = 12, res = 600,  type='cairo') 
p5
dev.off()


svglite:: svglite(filename =paste0("Figures/","ggtree.svg"), width = 20, height = 12)
p5
dev.off()


