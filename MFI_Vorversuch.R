##
# This R script was created for R version 3.4
#
##GLobal MFIs
#load libraries
library(ggplot2)
library(export)
library(ggpubr)
library(ggsci)

#set working directory
setwd(dir = "C:/Users/drmurke/Documents/Masterarbeit/Daten/MFI")

FLS <- read.csv("./FLS_MFI_Vorversuch.csv", header=T)

###################
#plots
global_MFI_CCL2.p<-ggplot(FLS, aes(x = Time, y = CCL2 , color = Dose , group = Dose)) + #label = ID für Labeling
  geom_jitter(size = 2.5, width=0.1) +
#  geom_line() +
 #facet_grid(~Time, scales = "free") +
  # facet_wrap(~Time) +
  labs(title = "FLS CCL2", y = "MFI", x = "Time") +
  theme_bw()+
  scale_fill_jama()+
  scale_color_jama()
global_MFI_CCL2.p

global_MFI_CCL3.p<-ggplot(FLS, aes(x = Time, y = CCL3 , color = Dose , group = Dose)) + #label = ID für Labeling
  geom_jitter(size = 2.5, width=0.1) +
  #  geom_line() +
  #facet_grid(~Time, scales = "free") +
  # facet_wrap(~Time) +
  labs(title = "FLS CCL3", y = "MFI", x = "Time") +
  theme_bw()+
  scale_fill_jama()+
  scale_color_jama()
global_MFI_CCL3.p


global_MFI_CD106.p<-ggplot(FLS, aes(x = Time, y = CD106 , color = Dose , group = Dose)) + #label = ID für Labeling
  geom_jitter(size = 2.5, width=0.1) +
  #  geom_line() +
  #facet_grid(~Time, scales = "free") +
  # facet_wrap(~Time) +
  labs(title = "FLS CD106", y = "MFI", x = "Time") +
  theme_bw()+
  scale_fill_jama()+
  scale_color_jama()
global_MFI_CD106.p

###Survival nCD4
nCD4 <- read.csv("./Survival_nCD4.csv", header=T)
surv_nCD4 <- ggplot(nCD4, aes(x = IL2, y= Survival, color = Dose)) +
  geom_point(size = 2.5) +
  labs(title = "Alive nCD4 cells ± IL2", y = "Alive nCD4 (%)", x = "") + 
  theme_bw()+
  scale_fill_jama()+
  scale_color_jama()



vorversuch_plot<-ggarrange(global_MFI_CD106.p, global_MFI_CCL3.p, surv_nCD4, common.legend = T, legend="bottom", nrow =1, ncol=3, labels = c("A","B","C"))
vorversuch_plot

export::graph2ppt(vorversuch_plot, file = "Vorversuch.pptx", width = 6.3, height=3, append =T)
