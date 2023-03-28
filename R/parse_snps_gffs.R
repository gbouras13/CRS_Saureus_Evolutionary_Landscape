library(tidyverse)

same_strain_gess <- c("420", "1415", "1676", "1992", "3344", 
                      "3997", "4681", "4875", "5047",
                      "5060", "5142", "5519", "5562", "5728")

# read in each gff 

myfiles <- c()


col_names <- c("contig", "program", "type",
               "start", "stop", "phase1", "phase2", "phase3",
               "description")

for (gess in same_strain_gess){
  
  tmp_df <- read.csv(file = paste0("../Snakemake_Output/SNIPPY_PAIR/", gess, "/snps.gff"), 
                     sep = "\t", 
                     skip = 0,
                     col.names = col_names)
  tmp_df$Gess <- gess
  tmp_df <- tmp_df %>%  relocate(Gess)
  myfiles[[gess]] <- tmp_df
  
}

# get total df
total_df <- bind_rows(myfiles)

# information about the variant
total_df$Coding_Region <- grepl("CDS", total_df$description)
# 136 variants in CDS in the 14 strains
table(total_df$Coding_Region)

total_df$synonymous_variant <- grepl("synonymous_variant", total_df$description)
total_df$stop_gained <- grepl("stop_gained", total_df$description)
total_df$missense_variant <- grepl("missense_variant", total_df$description)
total_df$frameshift_variant <- grepl("frameshift_variant", total_df$description)
total_df$inframe_variant <- grepl("inframe", total_df$description)

# get the coding df
coding_df <- total_df %>% 
  filter(Coding_Region)

# get the gene information
coding_df$gene <-  sapply(strsplit(as.character(coding_df$description), "_"), tail, 1)

# remove first 6 chars
coding_df$gene <- substr(coding_df$gene, 7, nchar(coding_df$gene))
# fix one with mistake due to underscore
coding_df$gene[25] <- "yoqJ UPF0398 protein SAOUHSC_01463"
# get genes that are mutated more than once 
table(coding_df$gene)[table(coding_df$gene) > 1]

dir.create("../R_Output_CSVs/")

coding_df$MSCRAMM <- ifelse(grepl("fibron", coding_df$description ),
                            TRUE, FALSE)

coding_df$MSCRAMM <- ifelse(grepl("sdrC", coding_df$description ),
                            TRUE, coding_df$MSCRAMM)

coding_df$MSCRAMM <- ifelse(grepl("isd", coding_df$description ),
                            TRUE, coding_df$MSCRAMM)

coding_df$MSCRAMM <- ifelse(grepl("Surface protein G", coding_df$description ),
                            TRUE, coding_df$MSCRAMM)





write.csv(coding_df, 
          file = paste0("../R_Output_CSVs/", "coding_snps.csv")
          )

write.csv(total_df, 
          file = paste0("../R_Output_CSVs/", "all_snps.csv")
)


# counts
table(total_df$Gess)
table(total_df$synonymous_variant) 
table(total_df$stop_gained )
table(total_df$missense_variant) 
table(total_df$frameshift_variant) 
table(total_df$inframe_variant )


# get genes that are mutated more than once 
table(coding_df$gene)[table(coding_df$gene) > 1]



