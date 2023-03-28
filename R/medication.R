library(tidyverse)


med_df<-read.csv("med_df.csv" )
med_df$rid<- as.character(med_df$rid)




med_df$Medication <- factor(med_df$Medication, levels= c("Augmentin","Ciprofloxacin 500mg", "mupirocin 0.05% Nasal rinse",
                                                         "Tobramycin Nasal rinse", "Trimethoprim/Sulfamethoxazole", "Doxycycline 50mg","Clarithromycin 250mg",
                                                         "Clindamycin 150mg","Doxycycline 100mg","Moxifloxacin 400mg", "Follow-up"))



k<- med_df %>% filter (Same_Strain=="Yes") %>%  ggplot( aes(x = start, y=Medication, colour=Medication)) +
  geom_segment(aes(xend= end, yend= Medication), linewidth=10) +
  facet_grid(rid ~., drop = TRUE, scales = "free", space = "free") +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("#0bd3d3", "#fb733c", "#d0d0d0", "#f9d70b", "#f20505",
                                "#711c91", "#133e7c", "#ea00d9" , "#564d4d", "Black", 
                                "#9bff37"))

k<-k+ ggtitle("Same strain CIs") +
  xlab("Time (days)") + ylab("")+
  theme(legend.title=element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.y=element_blank())


k2<- med_df %>% filter(Same_Strain=="No") %>%  ggplot( aes(x = start, y=Medication, colour=Medication)) +
  geom_segment(aes(xend= end, yend= Medication), linewidth =10)+
  facet_grid(rid ~ ., drop = TRUE, scales = "free", space = "free")+
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("#0bd3d3", "#fb733c", "#d0d0d0", "#f9d70b", "#f20505",
                                "#711c91", "#133e7c", "#ea00d9" , 
                                "#9bff37"))

k2<-k2+ ggtitle("Different strain CIs") +
  xlab("Time (days)") + ylab("")+
  theme(legend.title=element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.y=element_blank())



png(filename =paste0("Figures/", "med.png"), units = "in", width = 5.77*3, height = 10.19*2, res = 600,  type='cairo') 
ggpubr::ggarrange(k,k2, nrow = 1, common.legend = TRUE)
dev.off()



med_df %>% filter(Medication != "Follow-up") %>% 
  group_by(Medication) %>% 
  summarise(count = n(),
            sum_value = sum(dur))

 
 med_df %>% filter(Medication != "Follow-up") %>% 
   group_by(Same_Strain) %>% 
   summarise(count = n(),
             mean_value = mean(dur), 
             sd= sd(dur),
             min_value = min(dur),
             max_value = max(dur),
             range_value = max_value - min_value)
 
 
 
 ######  antibiotics correlation with biomas and baseline viability 
 
 
 data_comb_redone <- read.csv("data_comb_redone.csv")
 
 
 df<- data_comb_redone %>% filter(antibiotic == "baseline")
 
 df_first<- df %>% filter(group=="T0") %>% select(c(1,2,3,5,6,7,8,9,10)) 
 
 df_second<- df %>% filter(group=="T1") %>% select(c(1,2,3,5,6,7,8,9,10))
 
 df_comb<- merge(df_first, df_second, by="id" )
 
 df_comb$dif_biomass<- df_comb$biomass.y- df_comb$biomass.x
 df_comb$dif_baseline<- df_comb$fluorescence.y- df_comb$fluorescence.x
 
 
 df_exposure<-med_df %>% filter(Medication != "Follow-up") %>% 
   group_by(rid) %>% 
   summarise(sum_exposure = sum(dur))
 
 
 df_comb<- merge( df_comb, df_exposure, by.x = "id", by.y = "rid")
 
 
 library(ggpubr)
 
 ggplot(df_comb, aes(x= sum_exposure, y= dif_biomass, color=same_strain.y ))+
   geom_point()+
   stat_cor(method = "spearman")
 
 ggplot(df_comb, aes(x= sum_exposure, y= dif_baseline, color=same_strain.y ))+
   geom_point()+
   stat_cor(method = "spearman")
 
 ggplot(df_comb, aes(x= sum_exposure, y= dif_biomass ))+
   geom_point()+
   stat_cor(method = "spearman")

df_med_table<-med_df %>% filter(Medication=="Follow-up")





####### tables 
meta1 <- read.csv(file = "../metadata/gess_time.csv")
meta2 <- read.csv(file = "../metadata/metadata_phylogentic_tree.csv")
merged_meta<- merge(x = meta1,
                    y = meta2,
                    by.x = "Sample", by.y = "Cnumber_ID",
                    all = TRUE)

merged_meta<- merged_meta %>% filter(timepoint=="T1")

df_med_table<- merge(df_med_table, merged_meta, by.x = "rid", by.y = "rid")

df_med_table %>% 
  group_by(CRS_pheno) %>% 
  summarise(count = n(),
            mean_value = mean(end_follow),
            mean_age = mean(age),
            sd= sd(end_follow),
            min_value = min(end_follow),
            max_value = max(end_follow),
            range_value = max_value - min_value)


df_med_table %>% 
  summarise(count = n(),
            mean_value = mean(end_follow),
            mean_age = mean(age),
            sd_age = sd(age),
            sd= sd(end_follow),
            min_value = min(end_follow),
            max_value = max(end_follow),
            range_value = max_value - min_value)

 
# med_df<-read.csv("med_df.csv" )
# med_df$rid<- as.character(med_df$rid)
# 
# 
# 
# test<- med_df %>% filter( rid== 4875 )
# 
# test2<-test %>% 
#   uncount(dur, .id = "copy") %>% 
#   mutate(row = start + copy - 1) %>% 
#   complete(rid, row = 1:max(end_follow)) %>% 
#   dplyr:: select(c(1,2,3,5))
#   
# 
# 
# rid<- unique(med_df$rid)
# 
# 
# datalist = vector("list")
# 
# # Loop through 
# 
# 
# for (i in rid) {
#   
#     df<- med_df %>% filter( rid == i )
#    
#     df2<-df %>% 
#     uncount(dur, .id = "copy") %>% 
#     mutate(row = start + copy - 1) %>% 
#     complete(rid, row = 1:max(end_follow)) %>% 
#     dplyr:: select(c(1,2,3,5))
#   
#   
#   datalist[[i]] <- df2
# }
# 
# main_df= do.call(rbind, datalist)
# 
# main_df<-main_df[, -3]
# 
# 
# dup_rows <- main_df %>% 
#   group_by(rid, row) %>% 
#   filter(n() > 1)
# 
# 
# 
# colnames(main_df) <- c("ID", "day", "medication")
# 
# 
# df<- main_df[c(1:23),]
# 
# df$ID<-c(1,1,1,1,1,1,1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4)
# 
# df$day<-c(1,2,3,4,5,5,6,1,2,3,1,2,2,3,4,5,1,1,2,2,3,3,4)
# 
# df$medication<-c(NA,"A",NA,"A","A","B","A",NA,"C",NA,"A","C","D","D",NA,"A","C","D","C","D","C","D",NA)
# 
# 
# write.csv(main_df, "df_SO")
# 
# 
# k<-df%>% ggplot() +
#   geom_line(aes(x=day, y= medication ,color= medication, group=medication),linewidth = 1,) +
#   facet_grid(ID~ ., drop = TRUE, scales = "free", space = "free")
# 
# 
# k<-main_df %>% ggplot() +
#   geom_line(aes(x=row, y= Medication ,color= Medication, group=Medication),linewidth = 1,) +
#   facet_grid(rid ~ ., drop = TRUE, scales = "free", space = "free")
# 
# 
# png(filename =paste0("Figures/", "med.png"), units = "in", width = 10, height = 10, res = 300,  type='cairo') 
# k
# dev.off()
# 
# 
# df <- data.frame( ID = 1, day_start = c(2, 4, 5, 6), day_end = c(2, 5, 5, 6), medication = c("A", "A", "B", "A") )
# 



 