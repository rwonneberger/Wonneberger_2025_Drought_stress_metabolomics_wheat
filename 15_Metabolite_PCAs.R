
#setwd("/data/home/ronja.wonneberger/Projects/CresWheat/Results")
source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")



## Make PCAs

d1<-fread(paste0("Results/Normalized_T_MA_IS_Weight_norm_time_1_trt_all_tol_all.txt"), sep="\t")
d2<-fread(paste0("Results/Normalized_T_MA_IS_Weight_norm_time_2_trt_all_tol_all.txt"), sep="\t")
d3<-fread(paste0("Results/Normalized_T_MA_IS_Weight_norm_time_3_trt_all_tol_all.txt"), sep="\t")
d4<-fread(paste0("Results/Normalized_T_MA_IS_Weight_norm_time_4_trt_all_tol_all.txt"), sep="\t")

d1$trt<-ifelse(substr(d1$Sample , 10, 10)=="c", "control", "drought")
d2$trt<-ifelse(substr(d2$Sample , 10, 10)=="c", "control", "drought")
d3$trt<-ifelse(substr(d3$Sample , 10, 10)=="c", "control", "drought")
d4$trt<-ifelse(substr(d4$Sample , 10, 10)=="c", "control", "drought")


d1$genotype<-substr(d1$Sample , 1, 8)
d2$genotype<-substr(d2$Sample , 1, 8)
d3$genotype<-substr(d3$Sample , 1, 8)
d4$genotype<-substr(d4$Sample , 1, 8)

d1%<>%relocate(Sample, trt, genotype)
d2%<>%relocate(Sample, trt, genotype)
d3%<>%relocate(Sample, trt, genotype)
d4%<>%relocate(Sample, trt, genotype)

pca.d1<-prcomp(d1[,-c(1:3)])
pca.d2<-prcomp(d2[,-c(1:3)])
pca.d3<-prcomp(d3[,-c(1:3)])
pca.d4<-prcomp(d4[,-c(1:3)])

summ1 <- summary(pca.d1)
summ2 <- summary(pca.d2)
summ3 <- summary(pca.d3)
summ4 <- summary(pca.d4)

pca1_PC1<-round(summ1$importance[2,][1]*100, 2)
pca1_PC2<-round(summ1$importance[2,][2]*100, 2)

pca2_PC1<-round(summ2$importance[2,][1]*100, 2)
pca2_PC2<-round(summ2$importance[2,][2]*100, 2)

pca3_PC1<-round(summ3$importance[2,][1]*100, 2)
pca3_PC2<-round(summ3$importance[2,][2]*100, 2)

pca4_PC1<-round(summ4$importance[2,][1]*100, 2)
pca4_PC2<-round(summ4$importance[2,][2]*100, 2)

pca.d1<-as.data.frame(pca.d1$x)
pca.d2<-as.data.frame(pca.d2$x)
pca.d3<-as.data.frame(pca.d3$x)
pca.d4<-as.data.frame(pca.d4$x)

pca.d1$trt<-d1$trt
pca.d2$trt<-d2$trt
pca.d3$trt<-d3$trt
pca.d4$trt<-d4$trt

pca.d1$genotype<-d1$genotype
pca.d2$genotype<-d2$genotype
pca.d3$genotype<-d3$genotype
pca.d4$genotype<-d4$genotype

colors<-c('#f58231', '#ffe119', '#3cb44b', '#4363d8', '#dcbeff', '#800000', '#000075', '#e6194B', '#a9a9a9')


size=0.5
fontsize_text=7

p0<-ggplot(pca.d1, aes(x = PC1, y = PC2, color = trt)) +
  geom_point(size=size) +
  scale_color_manual(values = colors) +theme_bw()+
  guides(color = guide_legend(title = "Treatment"))+ xlab(paste0("PC1 (", pca1_PC1, "%)")) + ylab(paste0("PC2 (", pca1_PC2, "%)"))+ theme(legend.position="bottom")+  theme(legend.title = element_text(size=8),legend.text = element_text(size=7) ,text = element_text(size = fontsize_text),axis.text.x = element_text(size = fontsize_text), plot.title = element_text(size=fontsize_text)) 

p1<-ggplot(pca.d1, aes(x = PC1, y = PC2, color = trt)) +
  geom_point(size=size) +
  scale_color_manual(values = colors) +theme_bw()+
  xlab(paste0("PC1 (", pca1_PC1, "%)")) + ylab(paste0("PC2 (", pca1_PC2, "%)")) + theme(legend.position = "none")+ ggtitle("7 DAD") + theme(text = element_text(size = fontsize_text),axis.text.x = element_text(size = fontsize_text), plot.title = element_text(size=fontsize_text))

p2<-ggplot(pca.d2, aes(x = PC1, y = PC2, color = trt)) +
  geom_point(size=size) +
  scale_color_manual(values = colors) +theme_bw()+
  xlab(paste0("PC1 (", pca2_PC1, "%)")) + ylab(paste0("PC2 (", pca2_PC2, "%)")) + theme(legend.position = "none")+ ggtitle("13 DAD")+ theme(text = element_text(size = fontsize_text),axis.text.x = element_text(size = fontsize_text), plot.title = element_text(size=fontsize_text))

p3<-ggplot(pca.d3, aes(x = PC1, y = PC2, color = trt)) +
  geom_point(size=size) +
  scale_color_manual(values = colors) +theme_bw()+
  xlab(paste0("PC1 (", pca3_PC1, "%)")) + ylab(paste0("PC2 (", pca3_PC2, "%)")) + theme(legend.position = "none")+ ggtitle("22 DAD")+ theme(text = element_text(size = fontsize_text),axis.text.x = element_text(size = fontsize_text), plot.title = element_text(size=fontsize_text))

p4<-ggplot(pca.d4, aes(x = PC1, y = PC2, color = trt)) +
  geom_point(size=size) +
  scale_color_manual(values = colors) +theme_bw()+
  xlab(paste0("PC1 (", pca4_PC1, "%)")) + ylab(paste0("PC2 (", pca4_PC2, "%)")) + theme(legend.position = "none")+ ggtitle("28 DAD")+ theme(text = element_text(size = fontsize_text),axis.text.x = element_text(size = fontsize_text), plot.title = element_text(size=fontsize_text))


legend<-get_legend(p0)
p5<-as_ggplot(legend)



legend_y = 0.1
height = (1-legend_y)/2
pos1 = legend_y
pos2 = pos1 + height

tiff("Figs/Fig5.tiff", width=8, height=9.5, units="cm", compression = "lzw",res=600)
ggdraw()+
  draw_plot(p1, x = 0, y = pos2, width = 0.5, height = height)+
  draw_plot(p2, x = 0.5, y = pos2, width = 0.5, height = height)+
  draw_plot(p3, x = 0, y = pos1, width = 0.5, height = height) +
  draw_plot(p4, x = 0.5, y = pos1, width = 0.5, height = height)+
  
  draw_plot(legend, x = 0, y = 0, width = 1, height = legend_y)+
  draw_plot_label(label = c("a", "b", "c", "d"), size = 7,
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, pos2, pos2 ))
dev.off()







d1234<-fread(paste0("Results/Normalized_T_MA_IS_Weight_norm_time_1234_trt_all_tol_all.txt"), sep="\t")

d1234$trt<-ifelse(substr(d1234$Sample , 10, 10)=="c", "control", "drought")
d1234$time<-substr(d1234$Sample , 12, 12)
d1234$genotype<-substr(d1234$Sample , 1, 8)

d1234%<>%relocate(Sample, trt, time, genotype)

pca.d1234<-prcomp(d1234[,-c(1:4)])

summ1234 <- summary(pca.d1234)

pca1234_PC1<-round(summ1234$importance[2,][1]*100, 2)
pca1234_PC2<-round(summ1234$importance[2,][2]*100, 2)

pca.d1234<-as.data.frame(pca.d1234$x)

pca.d1234$trt<-d1234$trt
pca.d1234$time<-d1234$time
pca.d1234$genotype<-d1234$genotype



tiff("Figs/Fig4.tiff", width=8, height=10, units="cm", compression = "lzw",res=600)
ggplot(pca.d1234, aes(x = PC1, y = PC2, color = time, shape=trt)) +
  geom_point(size=0.5) +
  scale_color_manual(values = colors, labels=c('7 DAD', '13 DAD', '22 DAD', '28 DAD')) +theme_bw()+
  guides(shape = guide_legend(title = "Treatment", ncol=1,title.position="top"), color = guide_legend(title = "Timepoint", ncol=2, title.position="top")) + xlab(paste0("PC1 (", pca1234_PC1, "%)")) + ylab(paste0("PC2 (", pca1234_PC2, "%)"))+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7))+  
  theme(legend.title = element_text(size=7),legend.text = element_text(size=7))+ theme(legend.position="bottom")
dev.off()

