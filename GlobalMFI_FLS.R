##
# This R script was created for R version 3.4
#
##GLobal MFIs
#load libraries
library(ggplot2)
library(ggbiplot)
library(export)
library(ggpubr)
library(ggsci)

#### functions
## MFI plot
mfi_plot<- function(time, param, data){
  ggplot(subset(FLS_MFI_norm, Time == time), aes(x = Group, y = eval(parse(text = paste(param))) , fill = Group)) +
    stat_summary(fun.y = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.5, size = 0.1 ) +
    geom_point(aes(color=Group), size = 2.5) +
    geom_pwc(ref.group = 1, label = "p.adj.signif", method = "t.test", label.size=3,
             hide.ns=T,vjust = 0.5, p.adjust.method = "fdr")  +
    geom_pwc(ref.group = 2, label = "p.adj.signif", method = "t.test", label.size=3,
             hide.ns=T,vjust = 0.5, bracket.nudge.y = 0.12, p.adjust.method = "fdr") +
    labs(title = param, y = "MFI change", x = "") +
    theme_bw()+
    scale_fill_jama()+
    scale_color_jama()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}



#set working directory
setwd(dir = "C:/Users/drmurke/Documents/Masterarbeit/Daten/MFI")

#import csv files
FLS <- read.csv("./FLS_FLS_MFI.csv", header=T)
FLS_nCD4 <- read.csv("./FLS+nCD4_FLS_MFI.csv", header=T)
FLS_TH1 <- read.csv("./FLS+TH1_FLS_MFI.csv", header=T)
FLS_TH17 <- read.csv("./FLS+TH17_FLS_MFI.csv", header=T)
FLS_TH17TH1 <- read.csv("./FLS+TH17TH1_FLS_MFI.csv", header=T)



##########################
##normalise data
#FLS
#create subset matrixes
FLS_Ctrl <- subset(FLS, Dose == "D0")
FLS_D25 <- subset(FLS, Dose == "D25")
FLS_norm <-FLS_D25 #use D25 as template
# reduce matrixes to only MFI values
FLS_Ctrl <- FLS_Ctrl[7:17]
FLS_D25 <- FLS_D25[7:17]
#devide by Ctrl to normalise
FLS_norm[7:17]<- FLS_D25/FLS_Ctrl
#replace any INF (when deviding by 0) with 0
FLS_norm[sapply(FLS_norm, is.infinite)] <- NaN

#FLS_nCD4
FLS_nCD4_Ctrl <- subset(FLS_nCD4, Dose == "D0")
FLS_nCD4_D25 <- subset(FLS_nCD4, Dose == "D25")
FLS_nCD4_norm <- FLS_nCD4_D25
FLS_nCD4_Ctrl <- FLS_nCD4_Ctrl[7:17]
FLS_nCD4_D25 <- FLS_nCD4_D25[7:17]
FLS_nCD4_norm[7:17] <- FLS_nCD4_D25/FLS_nCD4_Ctrl
FLS_nCD4_norm[sapply(FLS_nCD4_norm, is.infinite)] <- NaN

#FLS_TH1
FLS_TH1_Ctrl <- subset(FLS_TH1, Dose == "D0")
FLS_TH1_D25 <- subset(FLS_TH1, Dose == "D25")
FLS_TH1_norm <- FLS_TH1_D25
FLS_TH1_Ctrl <- FLS_TH1_Ctrl[7:17]
FLS_TH1_D25 <- FLS_TH1_D25[7:17]
FLS_TH1_norm[7:17] <- FLS_TH1_D25/FLS_TH1_Ctrl
FLS_TH1_norm[sapply(FLS_TH1_norm, is.infinite)] <- NaN

#FLS_TH17
FLS_TH17_Ctrl <- subset(FLS_TH17, Dose == "D0")
FLS_TH17_D25 <- subset(FLS_TH17, Dose == "D25")
FLS_TH17_norm <- FLS_TH17_D25
FLS_TH17_Ctrl <- FLS_TH17_Ctrl[7:17]
FLS_TH17_D25 <- FLS_TH17_D25[7:17]
FLS_TH17_norm[7:17] <- FLS_TH17_D25/FLS_TH17_Ctrl
FLS_TH17_norm[sapply(FLS_TH17_norm, is.infinite)] <- NaN

#FLS_TH17TH1
FLS_TH17TH1_Ctrl <- subset(FLS_TH17TH1, Dose == "D0")
FLS_TH17TH1_D25 <- subset(FLS_TH17TH1, Dose == "D25")
FLS_TH17TH1_norm <- FLS_TH17TH1_D25
FLS_TH17TH1_Ctrl <- FLS_TH17TH1_Ctrl[7:17]
FLS_TH17TH1_D25 <- FLS_TH17TH1_D25[7:17]
FLS_TH17TH1_norm[7:17] <- FLS_TH17TH1_D25/FLS_TH17TH1_Ctrl
FLS_TH17TH1_norm[sapply(FLS_TH17TH1_norm, is.infinite)] <- NaN

FLS_MFI_norm <- rbind(FLS_norm, FLS_nCD4_norm, FLS_TH1_norm, FLS_TH17_norm, FLS_TH17TH1_norm)



###########################
#                         #
#                         #
#        PCA Plots        #
#                         #
###########################
FLS_MFI_norm8h<-subset(FLS_MFI_norm, Time == "08h")
FLS_MFI_norm8h_extra<-cbind(FLS_MFI_norm8h[1:6],FLS_MFI_norm8h[8],FLS_MFI_norm8h[10:12],FLS_MFI_norm8h[16:17])
FLS_MFI_norm8h_intra<-cbind(FLS_MFI_norm8h[1:7],FLS_MFI_norm8h[9],FLS_MFI_norm8h[13:15])
FLS_MFI_norm8h<-subset(FLS_MFI_norm, Time == "08h")
FLS_MFI_norm24h<-subset(FLS_MFI_norm, Time == "24h")
FLS_MFI_norm24h_extra<-cbind(FLS_MFI_norm24h[1:6],FLS_MFI_norm24h[8],FLS_MFI_norm24h[10:12],FLS_MFI_norm24h[16:17])
FLS_MFI_norm24h_intra<-cbind(FLS_MFI_norm24h[1:7],FLS_MFI_norm24h[9],FLS_MFI_norm24h[13:15])
#PCA_intra
#8h
df2 <- FLS_MFI_norm8h_intra
df2.m <- df2[7:11] #matrix
df2.g <- df2[, 4] #group annotation

df2.pca <- prcomp(df2.m, center = TRUE, scale. = TRUE)

PCA_FLS_8h_intra <- ggbiplot(df2.pca, obs.scale = 1, var.scale = 2, groups = as.factor(df2.g),
                  var.axes = TRUE) +
  geom_point(aes(color = groups), size = 2) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+
  scale_fill_jama()+
  scale_color_jama()+
  labs(title = "PCA FLS cells 8h intracellular")

PCA_FLS_8h_intra


#24h
df2 <- FLS_MFI_norm24h_intra
df2.m <- df2[7:11] #matrix
df2.g <- df2[, 4] #group annotation

df2.pca <- prcomp(df2.m, center = TRUE, scale. = TRUE)

PCA_FLS_24h_intra <- ggbiplot(df2.pca, obs.scale = 1, var.scale = 2, groups = as.factor(df2.g),
                             var.axes = TRUE) +
  geom_point(aes(color = groups), size = 2) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+
  scale_fill_jama()+
  scale_color_jama()+
  labs(title = "PCA FLS cells 24h intracellular")

PCA_FLS_24h_intra

#### extrazellulär
#8h
df2 <- FLS_MFI_norm8h_extra
df2.m <- df2[7:12] #matrix
df2.g <- df2[, 4] #group annotation

df2.pca <- prcomp(df2.m, center = TRUE, scale. = TRUE)

PCA_FLS_8h_extra <- ggbiplot(df2.pca, obs.scale = 1, var.scale = 2, groups = as.factor(df2.g),
                             var.axes = TRUE) +
  geom_point(aes(color = groups), size = 2) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+
  scale_fill_jama()+
  scale_color_jama()+
  labs(title = "PCA FLS cells 8h extra cellular")

PCA_FLS_8h_extra


#24h
df2 <- FLS_MFI_norm24h_extra
df2.m <- df2[7:12] #matrix
df2.g <- df2[, 4] #group annotation

df2.pca <- prcomp(df2.m, center = TRUE, scale. = TRUE)

PCA_FLS_24h_extra <- ggbiplot(df2.pca, obs.scale = 1, var.scale = 2, groups = as.factor(df2.g),
                             var.axes = TRUE) +
  geom_point(aes(color = groups), size = 2) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+
  scale_fill_jama()+
  scale_color_jama()+
  labs(title = "PCA FLS cells 24h extra cellular")

PCA_FLS_24h_extra

##export graphs
export::graph2ppt(PCA_FLS_24h_extra, file = "PCA_FLS.pptx", append =T, width=6.3, height = 4.5)
export::graph2ppt(PCA_FLS_24h_intra, file = "PCA_FLS.pptx", append =T, width=6.3, height = 4.5)
export::graph2ppt(PCA_FLS_8h_extra, file = "PCA_FLS.pptx", append =T, width=6.3, height = 4.5)
export::graph2ppt(PCA_FLS_8h_intra, file = "PCA_FLS.pptx", append =T, width=6.3, height = 4.5)



###############################
#                             #
#                             #
#         Plot galore         #
#                             #
#                             #
###############################
#create plots 8h


onesamplesigns <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("####", "###", "##", "#", "ns"))



FLS.08h.CD54<-mfi_plot(time="08h",param="CD54")
FLS.08h.CD54

FLS.08h.CD90<-mfi_plot(time="08h",param="CD90")
FLS.08h.CD90

FLS.08h.CD106<-mfi_plot(time="08h",param="CD106")
FLS.08h.CD106

FLS.08h.CD80<-mfi_plot(time="08h",param="CD80")
FLS.08h.CD80

FLS.08h.CD86<-mfi_plot(time="08h",param="CD86")
FLS.08h.CD86

FLS.08h.MHCII<-mfi_plot(time="08h",param="MHCII")
FLS.08h.MHCII

#####################
#create plots 24h
test<-subset(FLS_MFI_norm, Time == "08h" & Group == "FLS")
compare_means(IL6 ~ 1, test, mu = 1, method = "t.test", symnum.args = onesamplesigns )

FLS.24h.CD54<-mfi_plot(time="24h",param="CD54")
FLS.24h.CD54

FLS.24h.CD90<-mfi_plot(time="24h",param="CD90")

FLS.24h.CD90

FLS.24h.CD106<-mfi_plot(time="24h",param="CD106")
FLS.24h.CD106

FLS.24h.CD80<-mfi_plot(time="24h",param="CD80")
FLS.24h.CD80

FLS.24h.CD86<-mfi_plot(time="24h",param="CD86")
FLS.24h.CD86

FLS.24h.MHCII<-mfi_plot(time="24h",param="MHCII")
FLS.24h.MHCII

########################
#Intrazellulär 8h


FLS.08h.CCL2<-mfi_plot(time="08h",param="CCL2")
FLS.08h.CCL2

FLS.08h.CCL3<-mfi_plot(time="08h",param="CCL3")
FLS.08h.CCL3

FLS.08h.IL10<-mfi_plot(time="08h",param="IL10")
FLS.08h.IL10

FLS.08h.IL6<-mfi_plot(time="08h",param="IL6")
FLS.08h.IL6


###########################
#intracell 24h
FLS.24h.CCL2<-mfi_plot(time="24h",param="CCL2")
FLS.24h.CCL2

FLS.24h.CCL3<-mfi_plot(time="24h",param="CCL3")
FLS.24h.CCL3

FLS.24h.IL10<-mfi_plot(time="24h",param="IL10")
FLS.24h.IL10

FLS.24h.IL6<-mfi_plot(time="24h",param="IL6")
FLS.24h.IL6

FLS.24h.TNFa<-mfi_plot(time="24h",param="TNFa")
FLS.24h.TNFa

###TNFa
##38 = ausreißer
FLS.08h.TNFa<-ggplot(subset(FLS_MFI_norm, Time == "08h" & ID != "38"), aes(x = Group, y = TNFa , fill = Group)) +
  geom_point(aes(color=Group), size = 2.5) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, size = 0.1 ) +
  geom_pwc(ref.group = 1, label = "p.adj.signif", method = "t.test", label.size=3,
           hide.ns=T,vjust = 0.5, p.adjust.method = "fdr")  +
  geom_pwc(ref.group = 2, label = "p.adj.signif", method = "t.test", label.size=3,
           hide.ns=T,vjust = 0.5, bracket.nudge.y = 0.12, p.adjust.method = "fdr") +
  labs(title = "TNFa", y = "MFI change", x = "") +
  theme_bw()+
  scale_fill_jama()+
  scale_color_jama()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
FLS.08h.TNFa






######################
#
#   Arrange and Export Graphs
#
######################
##FLS 24h extra cellular
FLS.24h.extracell<-ggarrange(FLS.24h.CD54 + rremove("ylab") + rremove("xlab"),
          FLS.24h.CD90 + rremove("ylab") + rremove("xlab"),
          FLS.24h.CD106 + rremove("ylab") + rremove("xlab"),
          FLS.24h.CD80 + rremove("ylab") + rremove("xlab"),
          FLS.24h.CD86 + rremove("ylab") + rremove("xlab"),
          FLS.24h.MHCII + rremove("ylab") + rremove("xlab"),
          ncol=3, nrow=2, common.legend =T, legend = "bottom",
          labels = c("A", "B", "C","D","E","F"))
FLS.24h.extracell<-annotate_figure(FLS.24h.extracell, left = textGrob("relative MFI change",rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                top = text_grob("FLS 24h infection", color="black",face="bold", size="16"))
FLS.24h.extracell
export::graph2ppt(FLS.24h.extracell, file="FLS.pptx" , width = 6.3, height=4.5, append =T)

##FLS 8h extra cellular
FLS.08h.extracell<-ggarrange(FLS.08h.CD54 + rremove("ylab") + rremove("xlab"),
                             FLS.08h.CD90 + rremove("ylab") + rremove("xlab"),
                             FLS.08h.CD106 + rremove("ylab") + rremove("xlab"),
                             FLS.08h.CD80 + rremove("ylab") + rremove("xlab"),
                             FLS.08h.CD86 + rremove("ylab") + rremove("xlab"),
                             FLS.08h.MHCII + rremove("ylab") + rremove("xlab"),
                             ncol=3, nrow=2, common.legend =T, legend = "bottom",
                             labels = c("A", "B", "C","D","E","F"))
FLS.08h.extracell<-annotate_figure(FLS.08h.extracell, #, left = textGrob("relative MFI change", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                top = text_grob("FLS 8h infection", color="black",face="bold", size="16"))


#FLS_extracell_all<-annotate_figure(FLS_extracell_all, left = textGrob("relative MFI change",
#                                                    rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
FLS.08h.extracell
export::graph2ppt(FLS.08h.extracell, file="FLS.pptx" , width = 6.3, height=4.5, append =T)

##FLS 8h intra cellular
FLS.08h.intracell<-ggarrange(FLS.08h.CCL2 + rremove("ylab") + rremove("xlab"),
                             FLS.08h.CCL3 + rremove("ylab") + rremove("xlab"),
                             FLS.08h.IL10 + rremove("ylab") + rremove("xlab"),
                             FLS.08h.IL6 + rremove("ylab") + rremove("xlab"),
                             FLS.08h.TNFa + rremove("ylab") + rremove("xlab"),
                             ncol=3, nrow=2, common.legend =T, legend = "bottom",
                             labels = c("A", "B", "C","D","E","F"))
FLS.08h.intracell<-annotate_figure(FLS.08h.intracell, left = textGrob("relative MFI change", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                   top = text_grob("FLS 8h infection", color="black",face="bold", size="16"))

FLS.08h.intracell

export::graph2ppt(FLS.08h.intracell, file="FLS.pptx" , width = 6.3, height=4.5, append =T)
##24h intra
##FLS 24h intra cellular
FLS.24h.intracell<-ggarrange(FLS.24h.CCL2 + rremove("ylab") + rremove("xlab"),
                             FLS.24h.CCL3 + rremove("ylab") + rremove("xlab"),
                             FLS.24h.IL10 + rremove("ylab") + rremove("xlab"),
                             FLS.24h.IL6 + rremove("ylab") + rremove("xlab"),
                             FLS.24h.TNFa + rremove("ylab") + rremove("xlab"),
                             ncol=3, nrow=2, common.legend =T, legend = "bottom",
                             labels = c("A", "B", "C","D","E","F"))
FLS.24h.intracell<-annotate_figure(FLS.24h.intracell, left = textGrob("relative MFI change", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                   top = text_grob("FLS 24h infection", color="black",face="bold", size="16"))

FLS.24h.intracell
export::graph2ppt(FLS.24h.intracell, file="FLS.pptx" , width = 6.3, height=4.5, append =T)



##PCA

FLS_PCA<-ggarrange(PCA_FLS_24h_extra,
                   PCA_FLS_8h_intra,
                   ncol=1, nrow=2, legend = "right",
                   common.legend =T, labels = c("A","B"))

export::graph2ppt(FLS_PCA,  file= "PCA_FLS", width = 6.3, height = 6, append =T)

