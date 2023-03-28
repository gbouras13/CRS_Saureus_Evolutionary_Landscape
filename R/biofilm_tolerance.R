library(tidyverse)

'%nin%' <- Negate('%in%')

# Load data

data_comb_redone <- read.csv("data_comb_redone.csv")

data_comb_redone$concentration<- factor(data_comb_redone$concentration, levels = c("baseline","control","0.625", "1.25","2.5","5","10","20","40","80","160","320","640"))



library(plotrix)

#get mean line control# 

data_comb_redone %>% 
  filter(concentration %in% c("control"))  %>% summarise(
    mean = mean(fluorescence),
    sd = sd(fluorescence),
    se=std.error(fluorescence))


graph_all_antibiotics<- data_comb_redone %>% 
  filter(concentration %nin% c( "baseline", "control")) %>% 
    ggplot( aes(x=concentration,y=fluorescence))+
    stat_summary(geom = "bar",
                 fun = mean, color= "#5a09f9", fill= adjustcolor(c("#5a09f9"), alpha.f = 0.1)) +
    stat_summary(fun.data = mean_se,  
                 geom = "errorbar", color= "#5a09f9")+
    facet_wrap(.~ antibiotic,
               drop = FALSE,
               nrow = 1)+
    geom_hline(yintercept = 9564.583, colour= "red", size=2, linetype = "dashed" )+
    theme_classic(base_size = 16)+
    rotate_x_text (angle = 45)+
    labs(
      x = "Concentration (mg/L)",
      y = "Fluorescence (rfu)")
  
  

  
same_strain_tolerance<-    data_comb_redone %>%  filter(same_strain == "Yes") %>%
    filter(concentration %nin% c( "baseline", "control")) %>% 
    
    ggplot(aes(x=concentration, y=fluorescence, colour = timepoint, group = timepoint )) + 
    stat_summary(fun= mean, geom = "line", size=1.5) +
    stat_summary(fun.data = mean_se, geom = "errorbar", linewidth = 0.8)+
    labs(
      x = "Concentration (mg/L)",
      y = "Fluorescence (rfu)",
      tag = "")+
    theme_classic(base_size = 16) +
    scale_colour_manual(values = c("#7570B3","#E7298A"), labels = c("T0 same strain", "T1 same strain")) +
    rotate_x_text (angle = 45) +
    theme(legend.position="top", legend.title=element_blank())+
    facet_wrap(.~ antibiotic, drop = FALSE, nrow = 1)
  

svglite:: svglite(filename =paste0("Figures/","biofilm_tolerance.svg"), width = 5.77*3, height = 12)
ggpubr::ggarrange( graph_all_antibiotics, same_strain_tolerance, nrow = 2)
dev.off()
  
  
  library(lme4)
  library(lmerTest)
  library(predictmeans)
  library(multcomp)
  library(lsmeans)

  
  test<- data_comb_redone %>% filter( concentration %nin% c("baseline"))  #%>% filter(All.Courses <250 )

  model<- lmer(fluorescence~ group * same_strain  + concentration + antibiotic +(1|id), data= test)
  
  summary(model)
  
  drop1(model)
  
  res <- resid(model)
  
  residplot(model, level=1)
  
  
  
  ######get table of GLMM ####
  
  sjPlot::tab_model(
    model,
    show.r2 = TRUE,
    show.icc = FALSE,
    show.ci = FALSE,
    show.se = TRUE, 
    show.std = TRUE,
    show.re.var = TRUE,
    p.style = "scientific_stars",
    emph.p = TRUE,
    file = "model.doc")
  

  
  # # then take this html file and make .png file
  # webshot::webshot("plot.html", "plot.png")
  