library(tidyverse)
# 4784 is a maybe
same_strain_gess <- c("420", "1415", "1676", "1992", "3344", 
                      "3997", "4681", "4875", "5047",
                      "5060", "5142", "5519", "5562", "5728")

#####
# query = T1 change 
####

# read in each gff 
myfiles <- c()


col_names <- c("contig", "program", "type",
               "start", "stop", "phase1", "phase2", "phase3",
               "description")

for (gess in same_strain_gess){
  

  
  tmp_df <- read.csv(file = paste0("../Snakemake_Output/NUCDIFF/", gess, "/results/", gess, "_query_struct.gff" ), 
                     sep = "\t", 
                     skip = 1,
                     col.names = col_names)
  # only if length is > 0
  if (length(tmp_df$contig) > 0){
  
  tmp_df$Gess <- gess
  tmp_df <- tmp_df %>%  relocate(Gess)
  
  myfiles[[gess]] <- tmp_df
  }
  
}


# get total df
T1_df <- bind_rows(myfiles)

# parse

ref_df$description

T1_df$type <-  sapply(strsplit(as.character(T1_df$description), ";"), "[", 2)
T1_df$type <-  sapply(strsplit(as.character(T1_df$type), "="), "[", 2)
T1_df$len <-  sapply(strsplit(as.character(T1_df$description), ";"), "[", 3)
T1_df$len <-  sapply(strsplit(as.character(T1_df$len), "="), "[", 2)
T1_df$len <- as.numeric(T1_df$len)

#####
# ref = T0 change 
####

# read in each gff 
myfiles2 <- c()


for (gess in same_strain_gess){
  
  
  
  tmp_df <- read.csv(file = paste0("../Snakemake_Output/NUCDIFF/", gess, "/results/", gess, "_ref_struct.gff" ), 
                     sep = "\t", 
                     skip = 1,
                     col.names = col_names)
  # only if length is > 0
  if (length(tmp_df$contig) > 0){
    
    tmp_df$Gess <- gess
    tmp_df <- tmp_df %>%  relocate(Gess)
    
    myfiles2[[gess]] <- tmp_df
  }
  
}


# get total df
T0_df <- bind_rows(myfiles2)
T0_df$type <-  sapply(strsplit(as.character(T0_df$description), ";"), "[", 2)
T0_df$type <-  sapply(strsplit(as.character(T0_df$type), "="), "[", 2)
T0_df$len <-  sapply(strsplit(as.character(T0_df$description), ";"), "[", 3)
T0_df$len <-  sapply(strsplit(as.character(T0_df$len), "="), "[", 2)
T0_df$len <- as.numeric(T0_df$len)


T0_df$description <- NULL
T0_df$phase1 <- NULL
T0_df$phase2 <- NULL
T0_df$phase3 <- NULL

T1_df$description <- NULL
T1_df$phase1 <- NULL
T1_df$phase2 <- NULL
T1_df$phase3 <- NULL

dir.create("../R_Output_CSVs/")



write.csv(T0_df, 
          file = paste0("../R_Output_CSVs/", "nucdiff_T0.csv"),
          row.names = FALSE
)

write.csv(T1_df, 
          file = paste0("../R_Output_CSVs/", "nucdiff_T1.csv"),
          row.names = FALSE
)