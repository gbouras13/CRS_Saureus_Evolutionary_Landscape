library(tidyverse)


dir.create("R_Output_CSVs/")


# Set the path to the main folder
main_folder <- "../Snakemake_Output/SNIPPY_PAIR"

# Get a list of subfolder paths in the main folder
subfolders <- list.files(main_folder, full.names = TRUE)

# # Loop through the subfolders and extract a file from each one
# for (subfolder in subfolders) {
#   file_path <- file.path(subfolder, "snps.gff")
#   file_contents <-rtracklayer::import(file_path)
#   file_contents <- as.data.frame(file_contents)
#   assign(paste(basename(subfolder),"snps",sep="_"), as.data.frame(file_contents)) 
#   rm(file_contents)
# }



datalist = vector("list")

# Loop through the subfolders and extract a file from each one
for (subfolder in subfolders) {
  file_path <- file.path(subfolder, "snps.gff")
  file_contents <-rtracklayer::import(file_path)
  file_contents <- as.data.frame(file_contents)
  file_contents$rid<- basename(subfolder)
  datalist[[subfolder]] <- file_contents
}


main_df= do.call(rbind, datalist)

snp_df<-main_df %>% group_by(rid) %>% tally()



# load meta data and remove references 
meta <- read.csv(file = "../poppunk/poppunk_viz/poppunk_viz_microreact_clusters.csv" )


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



cluster_SS1<- meta %>%  filter(timepoint== "T0" ) %>%  select(c(rid, CC)) %>% arrange(rid)

cluster_SS2<- meta %>%  filter(timepoint== "T1" ) %>%  select(c(rid, CC)) %>% arrange(rid)

cluster_SS3<-meta %>%  filter(timepoint== "T0" ) %>%  select(c(rid, Cluster_Cluster__autocolour)) %>% arrange(rid)

cluster_SS4<-meta %>%  filter(timepoint== "T1" ) %>%  select(c(rid, Cluster_Cluster__autocolour)) %>% arrange(rid)


cluster_SS<-cbind(cluster_SS1, cluster_SS2,cluster_SS3,cluster_SS4)

cluster_SS$cc_same<- ifelse(cluster_SS[,2]== cluster_SS[,4], 1,0)

cluster_SS$vlkc_same<- ifelse(cluster_SS[,6]== cluster_SS[,8], 1,0)

cluster_SS$same<- ifelse(cluster_SS[,9]==1 & cluster_SS[,10]==1, 1,0)

table(cluster_SS$same)

rm( cluster_SS1,cluster_SS2,cluster_SS3,cluster_SS4 )

cluster_SS<-cluster_SS[,c(1,11)]

cluster_SS$rid<- as.character( cluster_SS$rid)

cluster_SS<- left_join(cluster_SS, snp_df )

cluster_SS<- cluster_SS %>%  arrange(n)

cluster_SS$same<- as.factor(cluster_SS$same)

cluster_SS$same<- ifelse(cluster_SS$same == 1, "similair", "different")


library(ggforce)

snp_graph<-ggplot(cluster_SS) +
  geom_bar(stat = "identity", position = "dodge", alpha=0.2,
           aes(x = reorder(rid, n) , y = n, fill= same, colour=same)) +
  theme_classic(base_size = 16)+
  scale_fill_manual(values = c("#ff00f4", "#38FF13")) +
  scale_color_manual(values = c("#ff00f4", "#38FF13")) +
  geom_hline(yintercept=100,linetype="dashed", size=2)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Host ID")+
  ylab("SNVs divergence between pairs")+
  guides(fill=guide_legend(title="CC and VLKCs\nclustering"))+
  guides(colour=guide_legend(title="CC and VLKCs\nclustering"))+
  facet_zoom(ylim = c(0, 150), zoom.size = 1.1)



grDevices::dev.size("in")    


svglite:: svglite(filename =paste0("Figures/","snp.svg"), width = 20, height = 6)
snp_graph
dev.off()




##########
### snp genes 
##########

main_df<-main_df %>%
  separate(note, paste("",1:14), sep = " ", extra = "merge")

colnames(main_df) <- c("seqnames", "start", "end", "width", "strand", "source", "type", "score", "phase", "Type", "SNP","Position","NA", "Chromosome", "plus", 
                         "nucleotide_Acid_Position", "Amino_Acid_Position","Variant_Type", 
                         "Change_Coding_Sequence", "Change_Protein", 
                         "Genome_Tag",  "gene_Name", "Functional_Annotation","rid")


main_df$rid<- as.integer( main_df$rid)

meta_snps<- meta %>%  filter(Time=="T0") %>% select (rid, Same_Strain)

main_df<- left_join( main_df,meta_snps)

main_df<- main_df %>% filter( Same_Strain== "Yes") %>% filter( Chromosome=="CDS")

main_df_all<- main_df %>% filter( Same_Strain== "Yes") %>% 
  filter(Chromosome=="CDS") %>% 
  filter(Functional_Annotation!= "hypothetical protein") 


main_df_all_table <-main_df_all %>%  select(10, 14, 18, 20, 22, 23, 24 )

write.csv( main_df_all_table, "R_Output_CSVs/snp_table.csv")

main_df<- main_df %>% filter( Same_Strain== "Yes") %>% 
  filter(Chromosome=="CDS") %>% 
  filter(Functional_Annotation!= "hypothetical protein") %>% 
  filter(Variant_Type != "synonymous_variant")



main_df1<-main_df %>%  filter(gene_Name=="") 



# isolates<- c("c13.gff", "c52.gff", "c9.gff", "c22.gff", "c133.gff", "c241.gff")

dir_path<- "../CHROMOSOME_GFFS"

file_list <- list.files(dir_path, pattern = "gff$")


datalist = vector("list")

for (i in file_list) {
  file_paths <- file.path(dir_path, i) # Combine the directory path with the file name to create a file path
  df<- rtracklayer::import(file_paths)
  df<-as.data.frame(df)
  datalist[[i]] <- df
}


lookup_df= do.call(rbind, datalist)

lookup_df<-unnest(lookup_df,Dbxref)

lookup_df<- lookup_df%>% filter(product %in% main_df$Functional_Annotation)


lookup_df<- lookup_df %>%
  filter(grepl("^GO", Dbxref))


my_vector <- main_df$Functional_Annotation

values_not_in_df <- setdiff(unique(my_vector), unique(lookup_df$product))


GO_table<-  data.frame( table(lookup_df$Dbxref))



