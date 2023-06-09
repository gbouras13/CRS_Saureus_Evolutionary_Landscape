library(pheatmap)
library(viridis)
library(reshape2)
library(ggplot2)
library(hexbin)
library(mixOmics)


dir.create("Figures/")



# read in mash similarity matrix 
plas_mash_df <- read.csv(file = "../Plasmid_Snakemake_Out/RESULTS/mash_matrix.csv")

# covnert to matrix, assign rownames
plas_mash_mat <- as.matrix(plas_mash_df)

rownames(plas_mash_mat) <- colnames(plas_mash_df)

# mash similarity 
plas_mash_mat <- 1 - plas_mash_mat


# read in jaccard matrix 
jaccard_df <- read.csv(file = "../Plasmid_Snakemake_Out/RESULTS/jaccard_matrix.csv")

# covnert to matrix, assign rownames
jaccard_mat <- as.matrix(jaccard_df)
rownames(jaccard_mat) <- colnames(jaccard_mat)

# convert NAs to 0 (no/few genes)
jaccard_mat[is.na(jaccard_mat)] <- 0

# set diagonals to 1 
for (i in 1:length(jaccard_df$ARC_1_1)){
  
  jaccard_mat[i,i] <- 1
}


###### calculate mash vs jaccard 
# melt dfs
jaccard_long_df <- melt(jaccard_mat, varnames = c("P1", "P2"))
colnames(jaccard_long_df)[3] <- "Jaccard"
mash_long_df <- melt(plas_mash_mat, varnames = c("P1", "P2"))
colnames(mash_long_df)[3] <- "Mash"

# combine mash and jaccard
combo_long_df <- merge(jaccard_long_df, mash_long_df, by = c("P1", "P2"))

# remove diagonals 
combo_long_df <- combo_long_df %>%  filter(P1 != P2)

# get samples of plasmids for plot 
combo_long_df$P1 <- as.character(combo_long_df$P1)
combo_long_df$P2 <- as.character(combo_long_df$P2)
combo_long_df$P1_Sample <- substr(combo_long_df$P1,1,nchar(combo_long_df$P1)-2)
combo_long_df$P2_Sample <- substr(combo_long_df$P2,1,nchar(combo_long_df$P2)-2)



# Plot


jaccard_vs_mash_plot <- ggplot(data = combo_long_df, aes(x = Mash, y = Jaccard)) +
  stat_binhex() + 
  scale_fill_distiller(palette = "Spectral", limits= c(1,300), direction = -1) +
  xlim(c(0,1))+
  geom_hline(yintercept = 0.7, linetype='dashed', col = 'red', size=1) + 
  geom_vline(xintercept = 0.98, linetype='dashed', col = 'red', size=1) + 
  theme_bw(base_size = 16) + 
  labs( x = "Mash Similarity", y = "Jaccard Gene Similarity") 



png(filename =paste0("Figures/", "jaccard_vs_mash_plot.png"), units = "in", width = 5.77, height = 4, res = 1200,  type='cairo') 
jaccard_vs_mash_plot
dev.off()




#################
##################


combo_long_df %>%  
  filter(Jaccard > 0.7, Mash > 0.98)

summary(combo_long_df$Mash)


# read in ST  from the poppunk clusters csv
meta_orig_df <- read.csv(file = "../metadata/poppunk_mlst.csv" )
colnames(meta_orig_df)[1] <- "Sample"


meta_df <- meta_orig_df %>% dplyr::select(Sample, ST, CC)


# get P1 and P2 dataframes
meta_p1_df <- meta_df
meta_p2_df <- meta_df
meta_p1_df$P1_Sample <- meta_p1_df$Sample
meta_p1_df$Sample <- NULL
colnames(meta_p1_df)[1] <- "ST_P1"
meta_p2_df$P2_Sample <- meta_p2_df$Sample
meta_p2_df$Sample <- NULL
colnames(meta_p2_df)[1] <- "ST_P2"


# combine into final long dataframe
combo_long_df <- merge(combo_long_df, meta_p1_df, by = c("P1_Sample"))
combo_long_df <- merge(combo_long_df, meta_p2_df, by = c("P2_Sample"))


combo_long_df$Same_ST <- ifelse(combo_long_df$ST_P1 == combo_long_df$ST_P2,
                                "Yes", "No")



#### check increase of plasmids in same strain group 

test<-combo_long_df %>% filter( Same_ST=="Yes")

test<-data.frame(unique(test$P1))

test <- separate(test, "unique.test.P1.", into = c("new_column1", "new_column2"), sep = "_")

test <- test %>%
  group_by(new_column1) %>%
  slice_max(new_column2)

meta2 <- read.csv(file = "../metadata/metadata_phylogentic_tree.csv")
meta3 <- read.csv(file = "../metadata/gess_time.csv")


test<-left_join(x = test, y = meta2,
                           by = c("new_column1"= "Cnumber_ID" ))

test<-left_join(x = test, y = meta3,
                by = c("new_column1"= "Sample" ))

test<- test %>% filter(Same_Strain=="Yes")


####

plasmids <- data.frame(Plamsid = rownames(plas_mash_mat))

plasmids$Sample = substring(plasmids$Plamsid , 1, nchar(plasmids$Plamsid )-2)

plasmids <- plyr::join(plasmids, meta_df)


annotdf <- data.frame( CC = plasmids$CC) 


# annotdf$CC<- ifelse(annotdf$CC== "Not_Assigned","none assigned",annotdf$CC )

rownames(annotdf) <- plasmids$Plamsid



# mycolors <-  c("blue","seagreen3","darkgreen",
#                "brown","tan","red","orange")

# paletteer::paletteer_d("vapoRwave::hotlineBling",7, type = "discrete")


mycolors <-  c("#00B19DFF","grey" ,"#D98594FF",
               "#8B2E8BFF","#EB4762FF","#228BDCFF","#D96237FF")
#
# mycolors <- paletteer::paletteer_d("vapoRwave::hotlineBling",7, type = "discrete")
# 
# mycolors[7]<- "grey" 

names(mycolors) <- unique(plasmids$CC)

mycolors <- list(mycolors = mycolors)

names(mycolors) <- "CC"




# show heatmap 
h_all <- pheatmap(plas_mash_mat,
                  color = magma(100),
                  fontsize_col = 12,
                  fontsize_row = 12,
                  fontsize = 24,
                  annotation_col  = annotdf,
                  annotation_colors = mycolors
) 


h_all 





# jpeg(file = paste0("Figures/", "mash_heatmap_plasmids.jpeg"), 
#      units = "in", width = 5.77, height = 5.77,  bg = "transparent", res = 600)
# h_all
# dev.off()


svglite:: svglite(filename =paste0("Figures/","mash_heatmap_plasmids.svg"), width = 5.77*3, height = 4*3)
h_all
dev.off()

png(filename = paste0("Figures/","mash_heatmap_plasmids.png"), units = "in", width = 5.77*3, height = 4*3 , res = 600,  type='cairo') 
h_all
dev.off()






###########
## plasmid CARD 
##########





### extract matching plasmids for blaz with identical plasmids
same_plasmid_df <- combo_long_df %>%  filter(Mash >= 0.98, Jaccard >= 0.7 )
plasmid_table <- table(same_plasmid_df$P1)

plasmid_count_df <- as.data.frame(t(plasmid_table))
colnames(plasmid_count_df) <- c("A", "Plasmid", "Count")
plasmid_count_df$A   <- NULL
plasmid_count_df <- plasmid_count_df %>%  arrange(desc(Count))

# get major plasmid 
major_conserved_plasmid <- same_plasmid_df %>% 
  filter(P1 == "C21_1") %>% 
  pull(P2)

major_conserved_plasmid <- append(major_conserved_plasmid, "C21_1")


write.csv(major_conserved_plasmid, 
          file=paste0("../R_Output_CSVs/", "major_conserved_plasmid.csv" ),
          row.names = F)





