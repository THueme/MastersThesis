##
# This R script was created for R version 3.4
#
##GLobal MFIs
#load libraries
library(ggplot2)
library(ggbiplot)
library(export)
library(ggsci)
library(ggpubr)

#set working directory
setwd(dir = "C:/Users/drmurke/Documents/Masterarbeit/Daten/Survival")

FLS <- read.csv("./FLS_survival.csv", header=T)
FLS_nCD4 <- read.csv("./FLS+nCD4_survival.csv", header=T)
FLS_TH1 <- read.csv("./FLS+TH1_survival.csv", header=T)
FLS_TH17 <- read.csv("./FLS+TH17_survival.csv", header=T)
FLS_TH17TH1 <- read.csv("./FLS+TH17TH1_survival.csv", header=T)
FLS_all <- rbind (FLS, FLS_nCD4,FLS_TH1,FLS_TH17,FLS_TH17TH1)


##########################
##normalise data
#FLS
FLS_Ctrl <- subset(FLS, Dose == "D0")
FLS_D25 <- subset(FLS, Dose == "D25")
FLS_norm <-FLS_D25 #use D25 as template
# reduce matrixes to only MFI values
FLS_Ctrl <- FLS_Ctrl[6]
FLS_D25 <- FLS_D25[6]
#devide by Ctrl to normalise
FLS_norm[6]<- FLS_D25-FLS_Ctrl
#substract 1 for absolute change
#FLS_norm[8:18]<-FLS_norm[8:18] - 1
#replace any INF (when deviding by 0) with 0
FLS_norm[sapply(FLS_norm, is.infinite)] <- 0

#FLS_nCD4
#create subset matrixes
FLS_nCD4_Ctrl <- subset(FLS_nCD4, Dose == "D0")
FLS_nCD4_D25 <- subset(FLS_nCD4, Dose == "D25")
FLS_nCD4_norm <- FLS_nCD4_D25 #use D25 as template
# reduce matrixes to only MFI values
FLS_nCD4_Ctrl <- FLS_nCD4_Ctrl[6]
FLS_nCD4_D25 <- FLS_nCD4_D25[6]
#devide by Ctrl to normalise
FLS_nCD4_norm[6] <- FLS_nCD4_D25-FLS_nCD4_Ctrl
#substract 1 for absolute change
#FLS_nCD4_norm[6] <- FLS_nCD4_norm[6] - 1
#replace any INF (when deviding by 0) with 0
FLS_nCD4_norm[sapply(FLS_nCD4_norm, is.infinite)] <- 0

#FLS_TH1
FLS_TH1_Ctrl <- subset(FLS_TH1, Dose == "D0")
FLS_TH1_D25 <- subset(FLS_TH1, Dose == "D25")
FLS_TH1_norm <- FLS_TH1_D25
FLS_TH1_Ctrl <- FLS_TH1_Ctrl[6]
FLS_TH1_D25 <- FLS_TH1_D25[6]
FLS_TH1_norm[6] <- FLS_TH1_D25-FLS_TH1_Ctrl
#FLS_TH1_norm[6] <- FLS_TH1_norm[6] - 1
FLS_TH1_norm[sapply(FLS_TH1_norm, is.infinite)] <- 0

#FLS_TH17
FLS_TH17_Ctrl <- subset(FLS_TH17, Dose == "D0")
FLS_TH17_D25 <- subset(FLS_TH17, Dose == "D25")
FLS_TH17_norm <- FLS_TH17_D25
FLS_TH17_Ctrl <- FLS_TH17_Ctrl[6]
FLS_TH17_D25 <- FLS_TH17_D25[6]
FLS_TH17_norm[6] <- FLS_TH17_D25-FLS_TH17_Ctrl
#FLS_TH17_norm[6] <- FLS_TH17_norm[6] - 1
FLS_TH17_norm[sapply(FLS_TH17_norm, is.infinite)] <- 0

#FLS_TH17TH1
FLS_TH17TH1_Ctrl <- subset(FLS_TH17TH1, Dose == "D0")
FLS_TH17TH1_D25 <- subset(FLS_TH17TH1, Dose == "D25")
FLS_TH17TH1_norm <- FLS_TH17TH1_D25
FLS_TH17TH1_Ctrl <- FLS_TH17TH1_Ctrl[6]
FLS_TH17TH1_D25 <- FLS_TH17TH1_D25[6]
FLS_TH17TH1_norm[6] <- FLS_TH17TH1_D25-FLS_TH17TH1_Ctrl
#FLS_TH17TH1_norm[6] <- FLS_TH17TH1_norm[6] - 1
FLS_TH17TH1_norm[sapply(FLS_TH17TH1_norm, is.infinite)] <- 0

FLS_MFI_norm <- rbind(FLS_norm, FLS_nCD4_norm, FLS_TH1_norm, FLS_TH17_norm, FLS_TH17TH1_norm)

##########################
# plots
##
##plot D25 survival
surv_24h <- ggplot(subset(FLS_all, Dose == "D25" & Time == "24h"), aes(x = Group, y = Live , fill = Group)) + #label = ID für Labeling
#  ylim(30) +
  geom_point(aes(color=Group),size = 2.5) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, size = 0.1 ) +
  geom_pwc(ref.group = 1, label = "p.adj.signif", method = "t.test", label.size=3,
           hide.ns=T,vjust = 0.5, p.adjust.method = "fdr")  +
  labs(title = "24h" , x = "") +
  theme_bw() +
  scale_fill_jama()+
  scale_color_jama() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
surv_24h

surv_8h<-ggplot(subset(FLS_all, Dose == "D25" & Time == "08h"), aes(x = Group, y = Live , fill = Group)) + #label = ID für Labeling
  geom_point(aes(color=Group),size = 2.5) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, size = 0.1 ) +
  geom_pwc(ref.group = 1, label = "p.adj.signif", method = "t.test", label.size=3,
           hide.ns=T,vjust = 0.5, p.adjust.method = "fdr")  +
  labs(title = "8h", y = "Alive (%)" , x = "") +
  theme_bw() +
  scale_fill_jama()+
  scale_color_jama() +
  ylim = c(30, 115) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
surv_8h

surv<- ggarrange(surv_8h + rremove("ylab") + rremove("xlab"),
                 surv_24h + rremove("ylab") + rremove("xlab"),
                 ncol=2, nrow=1, common.legend = T, legend ="bottom",
                 labels = c("A", "B"))

surv<-annotate_figure(surv, left = textGrob("Live cells (%)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                      top = text_grob("Live FLS Cells", color="black",face="bold", size="16"))
surv
export::graph2ppt(surv, file = "survival.pptx", height = 2.25, width =4.2, append = T)

##plot normalised Data
global_MFI.p<-ggplot(FLS_MFI_norm, aes(x = Group, y = Live  , group = Group, fill = Group)) + #label = ID für Labeling
  geom_point(size = 2.5) +
  #  geom_line() +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, size = 0.1 ) +
  # facet_wrap(~Time) +
  labs(title = "Δ of live cells (%)", y = "ΔAlive (%)" , x = "") +
  theme_bw() +
  scale_fill_jama()+
  scale_color_jama()
global_MFI.p
export::graph2ppt(global_MFI.p, file = "survival.pptx", height = 2, width =2*2.5, append = T)
