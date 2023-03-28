library(tidyverse)

# Set the path to the folder that contains the TSV files
folder_path <-"../Snakemake_Output/ABRICATE"

# Set the pattern to match at the end of the TSV filenames
filename_pattern <- "_card.tsv"


# Get a list of TSV files that match the pattern
tsv_files <- list.files(folder_path, pattern = filename_pattern, full.names = TRUE)

tsv_files<- tsv_files[-69]


datalist = vector("list")

# Loop through the TSV files and do something with each one
for (tsv_file in tsv_files) {
  # Do something with the TSV file here
  file_path <- tsv_file
  card_abricate <- read_tsv(file_path)
  card_abricate  <- as.data.frame(card_abricate)
  card_abricate$CI<- basename(tsv_file)
  datalist[[tsv_file]] <- card_abricate 
}




CARD_df= do.call(rbind, datalist)

CARD_df$CI <- sub("_card.tsv*", "", CARD_df$CI)

### only keep coverage and identity above 90%

CARD_df<- CARD_df %>%   filter(`%COVERAGE` >90) %>% filter( `%IDENTITY`  > 90)

df<- CARD_df %>% select ( CI, GENE )

df$value<- 1


df <-df %>%  group_by( CI, GENE ) %>% summarise_all(sum)

df <- df %>%
  mutate(GENE = recode(GENE, 
                      "Staphylococcys_aureus_LmrS" = "lmrS",
                      "PC1_beta-lactamase_(blaZ)" = "PC1 beta-lactamase (blaZ)", 
                      "Staphylococcus_aureus_FosB" = "fosB",
                      "Staphylococcus_aureus_norA"= "norA",
                      "ErmA"= "ermA"))


df<- df%>% spread(GENE, value)

df[is.na(df)] <- 0




# load clinical metadata and merge 
meta1 <- read.csv(file = "../metadata/gess_time.csv")
meta2 <- read.csv(file = "../metadata/metadata_phylogentic_tree.csv")
meta<- merge(x = meta1,
                    y = meta2,
                    by.x = "Sample", by.y = "Cnumber_ID",
                    all = TRUE)


rm(meta1, meta2)


meta$rowname<- paste0( "Host:" ,meta$rid)

meta$rowname<- paste(meta$rowname, paste0("CI:", meta$Sample), sep = " | " )

meta<-meta %>% select("rowname" ,"Sample", "Gess_plus_Time", "GESS_ID", "Timepoint", "Same_Strain")

meta$Same_Strain<- ifelse(meta$Same_Strain== "Yes", 1, 0 )


df<- merge.data.frame(meta, df, by.x = "Gess_plus_Time"  , by.y ="CI" )


df<- df %>% unite("ref", Timepoint, Same_Strain , remove = FALSE)


df$ref<- factor (df$ref, levels = c( "T1_0","T0_0", "T1_1", "T0_1") )

df_calc_card<- df

temp<- df


mat <- data.matrix(df)

rownames(mat) = df$rowname


######heatmapt Card 

library(ComplexHeatmap)


colors = structure(c("#000000","#EDBB8AFF", "#EB4A40FF", "#6C2167FF"), names = c("0", "1", "2","3")) # black, red


ht1<- Heatmap(mat[,-c(1:7)], name = "gene",
              cluster_columns = TRUE,
              cluster_rows  = FALSE,
              col = colors,
              row_split = mat[, 5],
              row_gap = unit(c(1, 4, 1), "mm"),
              cluster_row_slices= FALSE ,
              rect_gp = gpar(col = "white", lwd = 0.75),
              column_names_rot = 45,
              column_names_side = "top",
              show_column_dend = FALSE, 
              show_row_dend =  FALSE,  
              column_names_gp = gpar(fontsize = 15),
              row_names_side = "right",
              row_names_gp = gpar(fontsize = 15),
              row_title = c("Different strain T0", "Different strain T1", "Same strain T0", "Same strain T1"),
              row_title_rot = 90, 
              row_title_gp = gpar(fontsize = 20), 
              show_parent_dend_line = FALSE, 
              heatmap_legend_param = list(title = "Gene copies",
                                          at = c( 0,1,2,3),
                                          labels = c("0", "1","2", "3")))


              
              
ht1<-draw(ht1, padding = unit(c(2,2, 2, 2), "mm"))


grDevices::dev.size("in")    


svglite:: svglite(filename =paste0("Figures/","CARD.svg"), width = 10, height = 22)
ht1
dev.off()




##############
###### vfdb 
#############


# Set the path to the folder that contains the TSV files
folder_path <-"../Snakemake_Output/ABRICATE"

# Set the pattern to match at the end of the TSV filenames
filename_pattern <- "_vfdb.tsv"


# Get a list of TSV files that match the pattern
tsv_files <- list.files(folder_path, pattern = filename_pattern, full.names = TRUE)

tsv_files<- tsv_files[-69]


datalist = vector("list")

# Loop through the TSV files and do something with each one
for (tsv_file in tsv_files) {
  # Do something with the TSV file here
  file_path <- tsv_file
  vfdb_abricate <- read_tsv(file_path)
  vfdb_abricate  <- as.data.frame(vfdb_abricate )
  vfdb_abricate $CI<- basename(tsv_file)
  datalist[[tsv_file]] <- vfdb_abricate 
}




vfdb_df= do.call(rbind, datalist)


vfdb_df$CI <- sub("_vfdb.tsv*", "", vfdb_df$CI)

### only keep coverage and identity above 90%

vfdb_df<- vfdb_df%>%   filter(`%COVERAGE` >90) %>% filter( `%IDENTITY`  > 90)

df<- vfdb_df %>% select ( CI, GENE )

df$value<- 1


df <-df %>%  group_by( CI, GENE ) %>% summarise_all(sum)



df<- df%>% spread(GENE, value)

df[is.na(df)] <- 0




# load clinical metadata and merge 
meta1 <- read.csv(file = "../metadata/gess_time.csv")
meta2 <- read.csv(file = "../metadata/metadata_phylogentic_tree.csv")
meta<- merge(x = meta1,
             y = meta2,
             by.x = "Sample", by.y = "Cnumber_ID",
             all = TRUE)


rm(meta1, meta2)


meta$rowname<- paste0( "Host:" ,meta$rid)

meta$rowname<- paste(meta$rowname, paste0("CI:", meta$Sample), sep = " | " )

meta<-meta %>% select("rowname" ,"Sample", "Gess_plus_Time", "GESS_ID", "Timepoint", "Same_Strain")

meta$Same_Strain<- ifelse(meta$Same_Strain== "Yes", 1, 0 )


df<- merge.data.frame(meta, df, by.x = "Gess_plus_Time"  , by.y ="CI" )


df<- df %>% unite("ref", Timepoint, Same_Strain , remove = FALSE)


df$ref<- factor (df$ref, levels = c( "T1_0","T0_0", "T1_1", "T0_1") )

df_calc_vfdb<- df

temp<- df


mat <- data.matrix(df)

rownames(mat) = df$rowname


######heatmapt Card 

library(ComplexHeatmap)


colors = structure(c("#000000","#EDBB8AFF", "#EB4A40FF", "#6C2167FF"), names = c("0", "1", "2","3")) # black, red


ht1<- Heatmap(mat[,-c(1:7)], name = "gene",
              cluster_columns = TRUE,
              cluster_rows  = FALSE,
              col = colors,
              row_split = mat[, 5],
              row_gap = unit(c(1, 4, 1), "mm"),
              cluster_row_slices= FALSE ,
              rect_gp = gpar(col = "white", lwd = 0.75),
              column_names_rot = 45,
              column_names_side = "top",
              show_column_dend = FALSE, 
              show_row_dend =  FALSE,  
              column_names_gp = gpar(fontsize = 15),
              row_names_side = "right",
              row_names_gp = gpar(fontsize = 15),
              row_title = c("Different strain T0", "Different strain T1", "Same strain T0", "Same strain T1"),
              row_title_rot = 90, 
              row_title_gp = gpar(fontsize = 20), 
              show_parent_dend_line = FALSE, 
              heatmap_legend_param = list(title = "Gene copies",
                                          at = c( 0,1,2,3),
                                          labels = c("0", "1","2", "3")))




ht1<-draw(ht1, padding = unit(c(2,2, 2, 2), "mm"))


grDevices::dev.size("in")    


svglite:: svglite(filename =paste0("Figures/","VFDB.svg"), width = 22, height = 22)
ht1
dev.off()


###### calculations of total genes present or absent




table_amr_isolates<- df_calc_card %>% select(8:27) %>% rowSums() %>%  data.frame()
table_amr_genes<- df_calc_card %>% select(8:27) %>% colSums() %>%  data.frame()

table_amr<- df_calc_card %>%
  group_by(Timepoint, Same_Strain) %>%
  select(-c(1:5)) %>%
  summarise_all(sum) %>% t() %>%  data.frame()

table_amr$num_df<- table_amr$X1==table_amr$X3

table_amr$num_df2<- table_amr$X2==table_amr$X4


table_vfdb_isolates<- df_calc_vfdb %>% select(8:84) %>% rowSums() %>%  data.frame()
table_vfdb_genes<- df_calc_vfdb %>% select(8:84) %>% colSums() %>%  data.frame()

table_vfdb<- df_calc_vfdb %>%
  group_by(Timepoint, Same_Strain) %>%
  select(-c(1:5)) %>%
  summarise_all(sum) %>% t() %>%  data.frame()


table_vfdb$num_df<- table_vfdb$X1==table_vfdb$X3

table_vfdb$num_df2<- table_vfdb$X2==table_vfdb$X4

