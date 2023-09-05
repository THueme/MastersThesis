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

#custom palette damit Gruppen gleich wie FLS
pal_fill<-pal_jama("default", alpha = 1)(5)
pal_fill
pal_fill<-pal_fill[2:5]

onesamplesigns <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("####", "###", "##", "#", "ns"))



#### functions
mfi_plot<- function(time, param){
  ggplot(subset(CD4_MFI_norm, Time == time), aes(x = Group, y = eval(parse(text = paste(param))) , fill = Group)) +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.5, size = 0.1 ) +  #add median
    geom_point(aes(color=Group), size = 2.5) +   # plot datapoints
    geom_pwc(ref.group = 1, label = "p.adj.signif", method = "t.test", label.size=3,
             hide.ns=T,vjust = 0.5, p.adjust.method = "fdr")  +
    geom_pwc(ref.group = 2, label = "p.adj.signif", method = "t.test", label.size=3,
             hide.ns=T,vjust = 0.5, bracket.nudge.y = 0.12, p.adjust.method = "fdr") + ##add group comparison
    labs(title = param, y = "MFI change", x = "") +
    theme_bw()+
    scale_color_manual(values = pal_fill)+
    scale_fill_manual(values = pal_fill) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
}



#set working directory
setwd(dir = "C:/Users/drmurke/Documents/Masterarbeit/Daten/MFI")

#import csv files
#variable names: Celltype_Group
CD4_nCD4 <- read.csv("./FLS+nCD4_CD4_MFI.csv", header=T)
CD4_TH1 <- read.csv("./FLS+TH1_CD4_MFI.csv", header=T)
CD4_TH17 <- read.csv("./FLS+TH17_CD4_MFI.csv", header=T)
CD4_TH17TH1 <- read.csv("./FLS+TH17TH1_CD4_MFI.csv", header=T)

CD4<- rbind(CD4_nCD4, CD4_TH1, CD4_TH17, CD4_TH17TH1)
#PCA



##########################
##normalise data
#Group FLS

#CD4_nCD4
#create subset matrixes
CD4_nCD4_Ctrl <- subset(CD4_nCD4, Dose == "D0")
CD4_nCD4_D25 <- subset(CD4_nCD4, Dose == "D25")
CD4_nCD4_norm <- CD4_nCD4_D25 #use D25 as template
# reduce matrixes to only MFI values
CD4_nCD4_Ctrl <- CD4_nCD4_Ctrl[7:17]
CD4_nCD4_D25 <- CD4_nCD4_D25[7:17]
#devide by Ctrl to normalise
CD4_nCD4_norm[7:17] <- CD4_nCD4_D25/CD4_nCD4_Ctrl 

#replace any INF (when deviding by 0) with 0
CD4_nCD4_norm[sapply(CD4_nCD4_norm, is.infinite)] <- NaN

#CD4_TH1
CD4_TH1_Ctrl <- subset(CD4_TH1, Dose == "D0")
CD4_TH1_D25 <- subset(CD4_TH1, Dose == "D25")
CD4_TH1_norm <- CD4_TH1_D25
CD4_TH1_Ctrl <- CD4_TH1_Ctrl[7:17]
CD4_TH1_D25 <- CD4_TH1_D25[7:17]
CD4_TH1_norm[7:17] <- CD4_TH1_D25/CD4_TH1_Ctrl
CD4_TH1_norm[sapply(CD4_TH1_norm, is.infinite)] <- NaN

#CD4_TH17
CD4_TH17_Ctrl <- subset(CD4_TH17, Dose == "D0")
CD4_TH17_D25 <- subset(CD4_TH17, Dose == "D25")
CD4_TH17_norm <- CD4_TH17_D25
CD4_TH17_Ctrl <- CD4_TH17_Ctrl[7:17]
CD4_TH17_D25 <- CD4_TH17_D25[7:17]
CD4_TH17_norm[7:17] <- CD4_TH17_D25/CD4_TH17_Ctrl
CD4_TH17_norm[sapply(CD4_TH17_norm, is.infinite)] <- NaN

#CD4_TH17TH1
CD4_TH17TH1_Ctrl <- subset(CD4_TH17TH1, Dose == "D0")
CD4_TH17TH1_D25 <- subset(CD4_TH17TH1, Dose == "D25")
CD4_TH17TH1_norm <- CD4_TH17TH1_D25
CD4_TH17TH1_Ctrl <- CD4_TH17TH1_Ctrl[7:17]
CD4_TH17TH1_D25 <- CD4_TH17TH1_D25[7:17]
CD4_TH17TH1_norm[7:17] <- CD4_TH17TH1_D25/CD4_TH17TH1_Ctrl
CD4_TH17TH1_norm[sapply(CD4_TH17TH1_norm, is.infinite)] <- NaN

CD4_MFI_norm <- rbind(CD4_nCD4_norm, CD4_TH1_norm, CD4_TH17_norm, CD4_TH17TH1_norm)







###############################
#                             #
#                             #
#             PCA             #
#                             #
#                             #
###############################


CD4_MFI_norm8h<-subset(CD4_MFI_norm, Time == "08h")
CD4_MFI_norm8h_extra<-cbind(CD4_MFI_norm8h[1:6],CD4_MFI_norm8h[9],CD4_MFI_norm8h[11],CD4_MFI_norm8h[16:17])
CD4_MFI_norm8h_intra<-cbind(CD4_MFI_norm8h[1:7],CD4_MFI_norm8h[10],CD4_MFI_norm8h[13], CD4_MFI_norm8h[15])

CD4_MFI_norm24h<-subset(CD4_MFI_norm, Time == "24h")
CD4_MFI_norm24h_extra<-cbind(CD4_MFI_norm24h[1:6],CD4_MFI_norm24h[9],CD4_MFI_norm24h[11],CD4_MFI_norm24h[16:17])
CD4_MFI_norm24h_intra<-cbind(CD4_MFI_norm24h[1:8],CD4_MFI_norm24h[10],CD4_MFI_norm24h[12], CD4_MFI_norm24h[15])
#PCA_intra
#8h
df2 <- CD4_MFI_norm8h_intra
df2.m <- df2[7:10] #matrix
df2.g <- df2[, 4] #group annotation

df2.pca <- prcomp(df2.m, center = TRUE, scale. = TRUE)

PCA_CD4_8h_intra <- ggbiplot(df2.pca, obs.scale = 1, var.scale = 2, groups = as.factor(df2.g),
                             var.axes = TRUE) +
  geom_point(aes(color = groups), size = 2) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+
  scale_fill_jama()+
  scale_color_jama()+
  labs(title = "PCA CD4 cells 8h intracellular")

PCA_CD4_8h_intra


#24h
df2 <- CD4_MFI_norm24h_intra
df2.m <- df2[7:10] #matrix
df2.g <- df2[, 4] #group annotation

df2.pca <- prcomp(df2.m, center = TRUE, scale. = TRUE)

PCA_CD4_24h_intra <- ggbiplot(df2.pca, obs.scale = 1, var.scale = 2, groups = as.factor(df2.g),
                              var.axes = TRUE) +
  geom_point(aes(color = groups), size = 2) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+
  scale_fill_jama()+
  scale_color_jama()+
  labs(title = "PCA CD4 cells 24h intracellular")

PCA_CD4_24h_intra

#### extrazellulär
#8h
CD4_MFI_norm8h_extra[1,10]<-1 #set inf to 1
df2 <- CD4_MFI_norm8h_extra
df2.m <- df2[7:10] #matrix
df2.g <- df2[, 4] #group annotation

df2.pca <- prcomp(df2.m, center = TRUE, scale. = TRUE)

PCA_CD4_8h_extra <- ggbiplot(df2.pca, obs.scale = 1, var.scale = 2, groups = as.factor(df2.g),
                             var.axes = TRUE) +
  geom_point(aes(color = groups), size = 2) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+
  scale_fill_jama()+
  scale_color_jama()+
  labs(title = "PCA CD4 cells 8h extra cellular")

PCA_CD4_8h_extra


#24h
df2 <- CD4_MFI_norm24h_extra
df2.m <- df2[7:10] #matrix
df2.g <- df2[, 4] #group annotation

df2.pca <- prcomp(df2.m, center = TRUE, scale. = TRUE)

PCA_CD4_24h_extra <- ggbiplot(df2.pca, obs.scale = 1, var.scale = 2, groups = as.factor(df2.g),
                              var.axes = TRUE) +
  geom_point(aes(color = groups), size = 2) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+
  scale_fill_jama()+
  scale_color_jama()+
  labs(title = "PCA CD4 cells 24h extra cellular")

PCA_CD4_24h_extra

##export graphs
export::graph2ppt(PCA_CD4_24h_extra, file = "PCA_CD4.pptx", append =T, width=6.3, height = 4.5)
export::graph2ppt(PCA_CD4_24h_intra, file = "PCA_CD4.pptx", append =T, width=6.3, height = 4.5)
export::graph2ppt(PCA_CD4_8h_extra, file = "PCA_CD4.pptx", append =T, width=6.3, height = 4.5)
export::graph2ppt(PCA_CD4_8h_intra, file = "PCA_CD4.pptx", append =T, width=6.3, height = 4.5)





###########
#PCA
df2 <- CD4_MFI_norm

df2.m <- df2[7:17] #matrix
df2.g <- df2[, 4] #group annotation

df2.pca <- prcomp(df2.m, center = TRUE, scale. = TRUE)

df2.p <- ggbiplot(df2.pca, obs.scale = 1, var.scale = 2, groups = as.factor(df2.g),
                  var.axes = TRUE) +
  geom_point(aes(color = groups), size = 2) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+
  labs(title = "PCA CD4 extracellular markers 24h")+
  scale_color_manual(values = pal_fill)+
  scale_fill_manual(values = pal_fill)

df2.p

export::graph2ppt(df2.p, file="PCA_CD4_norm.pptx", append =T, height = 4, width =6.3,)


###############################
#                             #
#                             #
#         Plot galore         #
#                             #
#                             #
###############################

###extrazellulär 8h
CD4.08h.CD25<-mfi_plot(time="08h", param = "CD25")
#CD4.08h.CD25

CD4.08h.CD80<-mfi_plot(time="08h", param = "CD80")
#CD4.08h.CD80

CD4.08h.CD86<-mfi_plot(time="08h", param = "CD86")
#CD4.08h.CD86

CD4.08h.MHCII<-mfi_plot(time="08h", param = "MHCII")
CD4.08h.MHCII


#extrazellulär 24h
CD4.24h.CD25<-mfi_plot(time="24h", param = "CD25")
#CD4.24h.CD25

CD4.24h.CD80<-mfi_plot(time="24h", param = "CD80")
#CD4.24h.CD80

CD4.24h.CD86<-mfi_plot(time="24h", param = "CD86")
#CD4.24h.CD86

CD4.24h.MHCII<-mfi_plot(time="24h", param = "MHCII")
#CD4.24h.MHCII


######
#Intrazellulär 8h
CD4.08h.CCL2<-mfi_plot(time="08h", param = "CCL2")
#CD4.08h.CCL2

CD4.08h.CCL3<-mfi_plot(time="08h", param = "CCL3")
#CD4.08h.CCL3

CD4.08h.IL10<-mfi_plot(time="08h", param = "IL10")
#CD4.08h.IL10

CD4.08h.IL6<-mfi_plot(time="08h", param = "IL6")
#CD4.08h.IL6

CD4.08h.TNFa<-mfi_plot(time="08h", param = "TNFa")
#CD4.08h.TNFa

CD4.08h.IL17<-mfi_plot(time="08h", param = "IL17")
#CD4.08h.IL17

CD4.08h.IFNg<-mfi_plot(time="08h", param = "IFNg")
#CD4.08h.IFNg



###intrazellulär 24h
CD4.24h.CCL2<-mfi_plot(time="24h", param = "CCL2")
#CD4.24h.CCL2

CD4.24h.CCL3<-mfi_plot(time="24h", param = "CCL3")
#CD4.24h.CCL3

CD4.24h.IL10<-mfi_plot(time="24h", param = "IL10")
#CD4.24h.IL10

CD4.24h.IL6<-mfi_plot(time="24h", param = "IL6")
#CD4.24h.IL6

CD4.24h.TNFa<-mfi_plot(time="24h", param = "TNFa")
#CD4.24h.TNFa

CD4.24h.IL17<-mfi_plot(time="24h", param = "IL17")
#CD4.24h.IL17

CD4.24h.IFNg<-mfi_plot(time="24h", param = "IFNg")
#CD4.24h.IFNg

test<-subset(CD4_MFI_norm, Time == "08h" & Group == "FLS_TH17TH1")
compare_means(IFNg ~ 1, test, mu = 1, method = "t.test", symnum.args = onesamplesigns )

######################
#
#   Arrange and Export Graphs
#
######################
##FLS 24h extra cellular
CD4.24h.extracell<-ggarrange(CD4.24h.CD25 + rremove("ylab") + rremove("xlab"),
                             CD4.24h.CD80 + rremove("ylab") + rremove("xlab"),
                             CD4.24h.CD86 + rremove("ylab") + rremove("xlab"),
                             CD4.24h.MHCII + rremove("ylab") + rremove("xlab"),
                             ncol=2, nrow=2, common.legend =T, legend = "bottom",
                             labels = c("A", "B", "C","D","E","F"))
CD4.24h.extracell<-annotate_figure(CD4.24h.extracell, #left = textGrob("relative MFI change",rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                   top = text_grob("CD4 24h extracellular markers", color="black",face="bold", size="16"),
                                   left = textGrob("relative MFI change", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
CD4.24h.extracell
export::graph2ppt(CD4.24h.extracell, file ="CD4", width = 4.2, height = 4.5, append =T)

##FLS 8h extra cellular
CD4.08h.extracell<-ggarrange(CD4.08h.CD25 + rremove("ylab") + rremove("xlab"),
                             CD4.08h.CD80 + rremove("ylab") + rremove("xlab"),
                             CD4.08h.CD86 + rremove("ylab") + rremove("xlab"),
                             CD4.08h.MHCII + rremove("ylab") + rremove("xlab"),
                             ncol=2, nrow=2, common.legend =T, legend = "bottom",
                             labels = c("A","B","C", "D", "E","F","G","L"))
CD4.08h.extracell<-annotate_figure(CD4.08h.extracell, #, left = textGrob("relative MFI change", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                   top = text_grob("CD4 8h extracellular markers", color="black",face="bold", size="16"),
                                   left = textGrob("relative MFI change", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
export::graph2ppt(CD4.08h.extracell, file ="CD4", width = 4.2, height = 4.5, append =T)


##FLS 8h intra cellular
CD4.08h.intracell<-ggarrange(CD4.08h.CCL2 + rremove("ylab") + rremove("xlab"),
                             CD4.08h.CCL3 + rremove("ylab") + rremove("xlab"),
                             CD4.08h.IL10 + rremove("ylab") + rremove("xlab"),
                             CD4.08h.IL6 + rremove("ylab") + rremove("xlab"),
                             CD4.08h.TNFa + rremove("ylab") + rremove("xlab"),
                             CD4.08h.IL17 + rremove("ylab") + rremove("xlab"),
                             CD4.08h.IFNg + rremove("ylab") + rremove("xlab"),
                             ncol=3, nrow=2, common.legend =T, legend = "bottom",
                             labels = c("A", "B", "C","D","E","F"))
CD4.08h.intracell<-annotate_figure(CD4.08h.intracell, left = textGrob("relative MFI change", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                   top = text_grob("CD4 8h intracellular markers", color="black",face="bold", size="16"))
export::graph2ppt(CD4.08h.intracell, file="CD4.pptx" , width = 6.3, height=4.5, append =T)

##8h intracellulär ohne CCL2, CCL3 und IL6
CD4.08h.intracell<-ggarrange(CD4.08h.IL10 + rremove("ylab") + rremove("xlab"),
                             CD4.08h.TNFa + rremove("ylab") + rremove("xlab"),
                             CD4.08h.IL17 + rremove("ylab") + rremove("xlab"),
                             CD4.08h.IFNg + rremove("ylab") + rremove("xlab"),
                             ncol=2, nrow=2, common.legend =T, legend = "bottom",
                             labels = c("A", "B", "C","D","E","F"))
CD4.08h.intracell<-annotate_figure(CD4.08h.intracell, left = textGrob("relative MFI change", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                   top = text_grob("CD4 8h intracellular markers", color="black",face="bold", size="16"))

export::graph2ppt(CD4.08h.intracell, file="CD4.pptx" , width = 4.2, height=4.5, append =T)

#CD4.08h.intracell
CD4.08h.intracell<-ggarrange(CD4.08h.CCL2 + rremove("ylab") + rremove("xlab"),
                             CD4.08h.CCL3 + rremove("ylab") + rremove("xlab"),
                             CD4.08h.IL6 + rremove("ylab") + rremove("xlab"),
                             ncol=2, nrow=2, common.legend =T, legend = "bottom",
                             labels = c("A", "B", "C","D","E","F"))
CD4.08h.intracell<-annotate_figure(CD4.08h.intracell, left = textGrob("relative MFI change", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                   top = text_grob("CD4 8h intracellular markers", color="black",face="bold", size="16"))
export::graph2ppt(CD4.08h.intracell, file="CD4.pptx" , width = 4.2, height=4.5, append =T)


#24h intrazellulär
CD4.24h.intracell<-ggarrange(CD4.24h.CCL2 + rremove("ylab") + rremove("xlab"),
                             CD4.24h.CCL3 + rremove("ylab") + rremove("xlab"),
                             CD4.24h.IL10 + rremove("ylab") + rremove("xlab"),
                             CD4.24h.IL6 + rremove("ylab") + rremove("xlab"),
                             CD4.24h.TNFa + rremove("ylab") + rremove("xlab"),
                             CD4.24h.IL17 + rremove("ylab") + rremove("xlab"),
                             ncol=3, nrow=2, common.legend =T, legend = "bottom",
                             labels = c("A", "B", "C","D","E","F"))
CD4.24h.intracell<-annotate_figure(CD4.24h.intracell, left = textGrob("relative MFI change", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                   top = text_grob("CD4 24h intracellular markers", color="black",face="bold", size="16"))

#CD4.24h.intracell
export::graph2ppt(CD4.24h.intracell, file="CD4.pptx" , width = 6.3, height=4.5, append =T)



