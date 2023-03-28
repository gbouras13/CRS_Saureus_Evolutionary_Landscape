##### get list of plasmids with blaR1 ####

pres_abs_csv <- read.csv(file = "../Plasmid_Snakemake_Out/PANAROO/gene_presence_absence_roary.csv")

pres_abs_blaz <- pres_abs_csv %>% 
  filter(Gene %in% c("blaR1", "blaPC1~~~blaZ", "blaI")) %>% 
  dplyr::select(Gene, C100_1:C91_1)

# long DF
pres_abs_blaz_long <- pres_abs_blaz %>%  
  pivot_longer(cols = -Gene)
pres_abs_blaz_long$value <- ifelse(pres_abs_blaz_long$value != "", 
                                   1, 0)
# Wide DF
pres_abs_blaz_wide <- pres_abs_blaz_long %>%  
  pivot_wider(names_from = Gene,values_from = value)

pres_abs_blaz_wide$sum <- pres_abs_blaz_wide$`blaPC1~~~blaZ` + pres_abs_blaz_wide$blaR1 + pres_abs_blaz_wide$blaI

# all isolates with the full beta lactamase operon
pres_abs_blaz_final_list <- pres_abs_blaz_wide %>% 
  filter(sum == 3) %>% pull(name)

write.csv(pres_abs_blaz_final_list, file = "../R_Output_CSVs/all_blaz_plasmids.csv",
          row.names = F)




