library(tidyverse)

'%nin%' <- Negate('%in%')

# Load data

data_comb_redone <- read.csv("data_comb_redone.csv")

data_comb_redone$concentration<- factor(data_comb_redone$concentration, levels = c("baseline","control","0.625", "1.25","2.5","5","10","20","40","80","160","320","640"))


biomass<-data_comb_redone %>% 
  filter(concentration %in% c("baseline"))

biomass<- biomass%>%  unite("ref2", timepoint,same_strain, sep = "_" ,na.rm = TRUE, remove = FALSE)


biomass2<- biomass%>% filter(id != "5562")

k<- ggpubr::ggpaired(biomass, x = "timepoint", 
                     y = "biomass", 
                     facet.by = "same_strain", 
                     color = "ref2",
                     fill = "ref2", 
                     color_palette=c("#00f9ff","#000000", "#030056","#0900ff"),
                     palette = adjustcolor(c("#00f9ff","#000000", "#030056","#0900ff"), alpha.f = 0.3),
                     line.color = "gray",
                     line.size = 1,
                     point.size = 3,
                     short.panel.labs = FALSE, 
                     id= "id", 
                     xlab = "Timepoint", 
                     ylab = "ABS (600nm)",
                     ggtheme =  theme_classic(base_size = 16), 
                     panel.labs = list( same_strain = c("Different strain pairs (n=20)", "Same strain pairs (n=14)")),  
                     panel.labs.background = list(fill = "grey", color = "black"))


k<- k + stat_compare_means(method = "wilcox.test",label.y = 3.5, label.x = 1.35, paired = TRUE, size =5)

k<- k + theme(legend.title = element_blank(), legend.position = "")

k


png(filename =paste0("Figures/", "biomass.png"), units = "in", width = 5.77, height = 5.77, res = 1200,  type='cairo') 
k
dev.off()




k2<- ggpubr::ggpaired(biomass, x = "timepoint", 
                     y = "fluorescence", 
                     facet.by = "same_strain", 
                     color = "ref2",
                     fill = "ref2", 
                     color_palette=c("#00f9ff","#000000", "#030056","#0900ff"),
                     palette = adjustcolor(c("#00f9ff","#000000", "#030056","#0900ff"), alpha.f = 0.3),
                     line.color = "gray",
                     line.size = 1,
                     point.size = 3,
                     short.panel.labs = FALSE, 
                     id= "id", 
                     xlab = "Timepoint", 
                     ylab = "fluorescence (rfu)",
                     ggtheme =  theme_classic(base_size = 16), 
                     panel.labs = list( same_strain = c("Different strain pairs (n=20)", "Same strain pairs (n=14)")),  
                     panel.labs.background = list(fill = "grey", color = "black"))



k2<- k2 + stat_compare_means(method = "wilcox.test",label.y = 3.5, label.x = 1.35, paired = TRUE, size =5)

k2<- k2 + theme(legend.title = element_blank(), legend.position = "")

k2




library(rstatix)

res.aov <-biomass %>% anova_test(fluorescence ~ ref2)

res.aov

stat.test<- t_test(biomass, fluorescence ~ ref2, p.adjust.method = "bonferroni")

stat.test

stat.test <-stat.test %>% add_xy_position(x ="ref2")

stat.test


k2<- ggpubr::ggboxplot(biomass, x = "ref2", 
                      y = "fluorescence",
                      color = "ref2",
                      fill = "ref2", 
                      color_palette=c("#00f9ff","#000000", "#030056","#0900ff"),
                      palette = adjustcolor(c("#00f9ff","#000000", "#030056","#0900ff"), alpha.f = 0.3),
                      line.color = "gray",
                      line.size = 1,
                      point.size = 3,
                      short.panel.labs = FALSE,
                      xlab = "Timepoint", 
                      ylab = "fluorescence (rfu)",
                      ggtheme =  theme_classic(base_size = 16), 
                      panel.labs = list( same_strain = c("Different strain pairs (n=20)", "Same strain pairs (n=14)")),  
                      panel.labs.background = list(fill = "grey", color = "black"))+
  scale_x_discrete(labels=c("T0 Different strain", "T0 Same strain","T1 Different strain", "T1 Same strain")) +
  stat_pvalue_manual(stat.test, hide.ns = FALSE) +
  labs(subtitle = get_test_label(res.aov, detailed = FALSE))+
  rremove("xlab")+ 
  rremove("legend")


k2<-k2+ rotate_x_text(45)


png(filename =paste0("Figures/", "baseline_viability.png"), units = "in", width = 5.77, height = 5.77, res = 1200,  type='cairo') 
k2
dev.off()
