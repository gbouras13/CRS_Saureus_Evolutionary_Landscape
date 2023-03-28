library(tidyverse)
'%ni' <- Negate('%in%')

# get C and Gess match 
ref_df <- read.csv("../metadata/gess_time.csv")

# merge in C numbers to biofilm sheet
biofilm_df <- read.csv("../metadata/biofilm_data.csv")
biofilm_df[1] <- NULL
biofilm_df$Gess_plus_Time <- paste0(biofilm_df$id, "_", biofilm_df$group)

# merge in the data 
biofilm_df <- merge(biofilm_df, ref_df)


# baseline and control
biofilm_baseline_control_df <- biofilm_df %>%  filter(Antibiotic == "augmentin",
                                     Concentration %in% c("baseline", "control") )


# read in blaZ plasmids 
pres_abs_blaz_final_list <- read.csv(file = "../R_Output_CSVs/all_blaz_plasmids.csv")

pres_abs_blaz_final_list$Sample <-  sapply(strsplit(as.character(pres_abs_blaz_final_list$x), "_"), "[", 1)



# copy number 
copy_df <- read.csv(file = "../metadata/plassembler_copy_number.csv")
copy_df <- copy_df %>% filter(Plasmid %in% pres_abs_blaz_final_list$x)
copy_df <- merge(copy_df, biofilm_baseline_control_df)

# blaz plasmids vs not 
biofilm_df$blaZ <- ifelse( biofilm_df$Sample %in% pres_abs_blaz_final_list$Sample,
                           TRUE, FALSE)

biofilm_baseline_control_df$blaZ <- ifelse( biofilm_baseline_control_df$Sample %in% pres_abs_blaz_final_list$Sample,
                           TRUE, FALSE)

# baseline
jpeg(file = paste0("../Figures/", "blaz_baseline_biofilm.jpeg"), 
     units = "in", width = 13, height = 10,  bg = "transparent", res = 600)
ggplot(biofilm_baseline_control_df %>% filter(Concentration == "baseline"), 
       aes(x = blaZ, y = fluorescence, 
           color = Timepoint)) +
  geom_boxplot() + geom_jitter()
dev.off()

jpeg(file = paste0("../Figures/", "blaz_control_biofilm.jpeg"), 
     units = "in", width = 13, height = 10,  bg = "transparent", res = 600)
ggplot(biofilm_baseline_control_df %>% filter(Concentration == "control"), 
       aes(x = blaZ, y = fluorescence, 
           color = Timepoint)) +
  geom_boxplot() + geom_jitter()
dev.off()

# big differences at baseline and control for this.
linear_model <- lm(fluorescence~blaZ*Timepoint+ Concentration, 
                   data = biofilm_baseline_control_df)
summary(linear_model)


# seems reasonable to conclude that there is a difference
biofilm_baseline_control_df %>% 
  group_by(Timepoint, blaZ,  Concentration) %>%  
  summarise(mean = mean(fluorescence),
            sd = sd(fluorescence),
            count = n())

##### Same vs Diff Strain #####

# baseline
jpeg(file = paste0("../Figures/", "same_strain_baseline_biofilm.jpeg"), 
     units = "in", width = 13, height = 10,  bg = "transparent", res = 600)
ggplot(biofilm_baseline_control_df %>% filter(Concentration == "baseline"), 
       aes(x = Same_Strain, y = fluorescence, 
           color = Timepoint)) +
  geom_boxplot() + geom_jitter()
dev.off()

jpeg(file = paste0("../Figures/", "same_strain_control_biofilm.jpeg"), 
     units = "in", width = 13, height = 10,  bg = "transparent", res = 600)
ggplot(biofilm_baseline_control_df %>% filter(Concentration == "control"), 
       aes(x = Same_Strain, y = fluorescence, 
           color = Timepoint)) +
  geom_boxplot() + geom_jitter()
dev.off()

# big differences at baseline and control for this.
linear_model <- lm(fluorescence~Same_Strain*Timepoint+ Concentration, 
                   data = biofilm_baseline_control_df)
summary(linear_model)


# seems reasonable to conclude that there is a difference
biofilm_baseline_control_df %>% 
  group_by(Timepoint, Same_Strain, Concentration) %>%  
  summarise(mean = mean(fluorescence),
            sd = sd(fluorescence),
            count = n())


# corr of 0.4

# correlation with copy number only for baseline not control
cor(copy_df %>% filter(Concentration == "baseline") %>%  
      pull(plasmid_copy_number_long), 
    copy_df %>% filter(Concentration == "baseline") %>%  
      pull(fluorescence), method = 'spearman')

cor(copy_df %>% filter(Concentration == "control") %>%  
      pull(plasmid_copy_number_long), 
    copy_df %>% filter(Concentration == "control") %>%  
      pull(fluorescence), method = 'spearman')

