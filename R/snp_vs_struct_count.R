library(tidyverse)

nucdiff_T0 <- read.csv(file = paste0("../R_Output_CSVs/", "nucdiff_T0.csv"))
table(nucdiff_T0 %>%  filter(len > 100) %>% pull(Gess))


snps_df <- read.csv(file = "../R_Output_CSVs/all_snps.csv")

# number of SNPs vs structural variants
table(snps_df$Gess)
table(nucdiff_T0 %>%  filter(len > 100) %>% pull(Gess))
