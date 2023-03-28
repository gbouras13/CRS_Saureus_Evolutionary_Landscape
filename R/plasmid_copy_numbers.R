library(tidyverse)

# gets copy numbers
copy_df <- read.csv(file = "../metadata/plassembler_copy_number.csv")
meta_df <- read.csv(file = "../metadata/gess_time.csv")
# copy number DF
copy_df <- merge(copy_df, meta_df)


#blaz plasmids 
blaz_plasmids <- read.csv(file = "../R_Output_CSVs/all_blaz_plasmids.csv")


copy_df$'blaz plasmid'<- ifelse(copy_df$Plasmid %in% blaz_plasmids$x, "yes", "no")


#### same strain copy numbers 
same_strain_copy_df <- copy_df %>% filter(Same_Strain == "Yes")


# same_strain_copy_df <- copy_df %>% filter(Same_Strain == "No") 
# 
# same_strain_copy_df<- same_strain_copy_df [, c(1:3, 14:18) ]

# conserved plasmids - there are 8 conserved between T0 and T1 SS isolates
conserved_plasmids <- c("C80_1", "C208_1", "C67_1",
                        "C294_2", "C72_1", "C351_1",
                        "C133_1", "C179_1", "C45_1",
                        "C149_1", "C224_1", "C349_1",
                        "C222_1", "C333_1", "C241_1",
                        "C309_1")


conserved_df <- same_strain_copy_df %>% filter(Plasmid %in% conserved_plasmids)

conserved_df <- conserved_df %>%  arrange(Gess_plus_Time)


# long
T0 <- conserved_df %>%  filter(Timepoint == "T0") %>% 
  pull(plasmid_copy_number_long)

T1 <- conserved_df %>%  filter(Timepoint == "T1") %>% 
  pull(plasmid_copy_number_long)

plot_df <- data.frame(T0, T1)


# long
library(ggpubr)


plasmid_CN_long_reads<- ggpaired(plot_df,
         cond1 = "T0",
         cond2 = "T1",
         color = "condition",
         fill = "condition",
         color_palette= c("#ff9a00", "#c900ff"),
         palette =adjustcolor(c("#ff9a00", "#c900ff"), alpha.f = 0.3),
         xlab = "",
         ylab = "Plasmid Copy Number (long reads)",
         line.color = "gray",
         line.size = 1,
         point.size = 2,
         legend.title = "Timepoint") +
  stat_compare_means(paired = TRUE)+
  theme_classic(base_size = 16)
  


# wilcox test suggests trend to higher copy number

wilcox.test(plot_df$T0, plot_df$T1, paired = TRUE)


# short
T0 <- conserved_df %>%  filter(Timepoint == "T0") %>% 
  pull(plasmid_copy_number_short)
T1 <- conserved_df %>%  filter(Timepoint == "T1") %>% 
  pull(plasmid_copy_number_short)

plot_df <- data.frame(T0, T1)

plasmid_CN_short_reads<- ggpaired(plot_df,
                                 cond1 = "T0",
                                 cond2 = "T1",
                                 color = "condition",
                                 fill = "condition",
                                 color_palette= c("#ff9a00", "#c900ff"),
                                 palette = adjustcolor(c("#ff9a00", "#c900ff"), alpha.f = 0.3),
                                 xlab = "",
                                 ylab = "Plasmid Copy Number (short reads)",
                                 line.color = "gray",
                                 line.size = 1,
                                 point.size = 2,
                                 legend.title = "Timepoint") +
  stat_compare_means(paired = TRUE)+
  theme_classic(base_size = 16)



png(filename =paste0("Figures/", "same_strain_plasmid_copy_long.png"), units = "in", width = 5.77*1.5, height = 4*1.5, res = 600,  type='cairo') 
ggarrange(plasmid_CN_long_reads, plasmid_CN_short_reads, common.legend = TRUE )
dev.off()





############
# blaz 
##########

blaz_plasmids <- read.csv(file = "../R_Output_CSVs/all_blaz_plasmids.csv")

blaz_plasmid_df <- copy_df %>% filter(Plasmid %in% blaz_plasmids$x)


blaZ_plasmid_copy_short<- ggplot(blaz_plasmid_df,
       aes(x = Timepoint,
           y = plasmid_copy_number_short) ) +
  geom_boxplot()+ geom_jitter() + ylab("Copy Number Short")+
  stat_compare_means()


blaZ_plasmid_copy_long<-ggplot(blaz_plasmid_df,
       aes(x = Timepoint,
           y = plasmid_copy_number_long) ) +
  geom_boxplot() + geom_jitter() + ylab("Copy Number Long")+
  stat_compare_means()



ggarrange(blaZ_plasmid_copy_short, blaZ_plasmid_copy_long, common.legend = TRUE )


# some evidence higher
wilcox.test(plasmid_copy_number_long~Timepoint,
            data = blaz_plasmid_df)







### all plasmids short vs long 

# p<-ggplot(copy_df,aes(x = plasmid_copy_number_short,y = plasmid_copy_number_long, color=blaz))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)
#   
# 
# 
# p+stat_cor(method='spearman')


p <- ggscatter(copy_df, x = "plasmid_copy_number_short", 
                y = "plasmid_copy_number_long",
               xlab = "Plasmid Copy Number short reads",
               ylab = "Plasmid Copy Number long reads",
               size = 2)

p<-p + stat_cor(aes(label = ..r.label..),  label.x = 3, method = 'spearman', size = 10 )

p<-p+ theme_classic(base_size = 16)

p<-p+ geom_abline(linetype = "dashed", size= 1.5, color= "grey")






p2<- ggplot(copy_df, aes(x=length, y=plasmid_copy_number_long/plasmid_copy_number_short, colour=`blaz plasmid`))+
  geom_point(size=2)+
  theme_classic(base_size = 16)+
  xlab("Plasmid lenght (bp)")+
  ylab("Ratio plasmid copy number (long/short)")+
  scale_color_manual(values = c("#00ffb8", "#7c00ff"))+
  theme(legend.position = "top")
  

  
  
p2<-p2+ geom_hline(linetype = "dashed", size= 1.5, color= "grey", yintercept = 1)
  
       

png(filename =paste0("Figures/", "all_plasmids_short_vs_long_correlation.png"), units = "in", width = 5.77*1.5, height = 5.77*2, res = 600,  type='cairo') 

ggpubr::ggarrange(p, NULL, p2, nrow = 3, heights = c(1,0.1,1), labels = c("A","", "B"))

dev.off()


# cor of 0.64
cor(copy_df$plasmid_copy_number_long, copy_df$plasmid_copy_number_short,
    method = 'spearman')




#### copy numbers
long_div_short <- copy_df$plasmid_copy_number_long / copy_df$plasmid_copy_number_short

wilcox.test(long_div_short)

median(long_div_short)


long_div_short <- blaz_plasmid_df$plasmid_copy_number_long / blaz_plasmid_df$plasmid_copy_number_short

wilcox.test(long_div_short)

median(long_div_short)





######## diff strain 

#### same strain copy numbers 
diff_strain_copy_df <- copy_df %>% filter(Same_Strain == "No")
diff_strain_copy_df %>% 
  group_by(Timepoint) %>% 
  summarise(count = n())

### blaZ

blaz_plasmid_df %>% 
  group_by(Same_Strain, Timepoint) %>%  
  summarise(count = n())


### all 

copy_df %>% 
  group_by(Same_Strain, Timepoint) %>%  
  summarise(count = n())



