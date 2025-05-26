


source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

library(Rtsne)


colors<-c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')




size=0.5
fontsize_text=7

d2<-fread(paste0("Results/Normalized_T_MA_IS_Weight_norm_time_2_trt_all_tol_all.txt"), sep="\t")
d3<-fread(paste0("Results/Normalized_T_MA_IS_Weight_norm_time_3_trt_all_tol_all.txt"), sep="\t")
d4<-fread(paste0("Results/Normalized_T_MA_IS_Weight_norm_time_4_trt_all_tol_all.txt"), sep="\t")

# plot 1
Genotype <- substr(d1$Sample, start = 1, stop = 8)
Treatment<-substr(d1$Sample, start = 10, stop = 10)

Treatment<-gsub("c", "control", Treatment)
Treatment<-gsub("d", "drought", Treatment)

met_matrix_1 <- as.matrix(d1[,-1])

tsne_out_1 <- Rtsne(met_matrix_1)


tsne_plot_1 <- data.frame(x = tsne_out_1$Y[,1], 
                        y = tsne_out_1$Y[,2])
tsne_plot_1$Genotype<-Genotype
tsne_plot_1$Treatment<-Treatment

p1<-ggplot2::ggplot(tsne_plot_1) + geom_point(aes(x=x,y=y, color=Genotype, shape=Treatment))+ theme_bw()+
  scale_color_manual(values = colors) + xlab("Dim 1")+ ylab("Dim 2")+ theme(legend.position = "none")+ggtitle("7 DAD") + theme(legend.title = element_text(size=8),legend.text = element_text(size=7) ,text = element_text(size = fontsize_text),axis.text.x = element_text(size = fontsize_text), plot.title = element_text(size=fontsize_text)) 



# plot 2
Genotype <- substr(d2$Sample, start = 1, stop = 8)
Treatment<-substr(d2$Sample, start = 10, stop = 10)

Treatment<-gsub("c", "control", Treatment)
Treatment<-gsub("d", "drought", Treatment)

met_matrix_2 <- as.matrix(d2[,-1])

tsne_out_2 <- Rtsne(met_matrix_2)


tsne_plot_2 <- data.frame(x = tsne_out_2$Y[,1], 
                          y = tsne_out_2$Y[,2])
tsne_plot_2$Genotype<-Genotype
tsne_plot_2$Treatment<-Treatment

p2<-ggplot2::ggplot(tsne_plot_2) + geom_point(aes(x=x,y=y, color=Genotype, shape=Treatment))+theme_bw()+
  scale_color_manual(values = colors) + xlab("Dim 1")+ ylab("Dim 2")+ theme(legend.position = "none")+ggtitle("13 DAD")+ theme(legend.title = element_text(size=8),legend.text = element_text(size=7) ,text = element_text(size = fontsize_text),axis.text.x = element_text(size = fontsize_text), plot.title = element_text(size=fontsize_text)) 





# plot 3
Genotype <- substr(d3$Sample, start = 1, stop = 8)
Treatment<-substr(d3$Sample, start = 10, stop = 10)

Treatment<-gsub("c", "control", Treatment)
Treatment<-gsub("d", "drought", Treatment)

met_matrix_3 <- as.matrix(d3[,-1])

tsne_out_3 <- Rtsne(met_matrix_3)


tsne_plot_3 <- data.frame(x = tsne_out_3$Y[,1], 
                          y = tsne_out_3$Y[,2])
tsne_plot_3$Genotype<-Genotype
tsne_plot_3$Treatment<-Treatment

p3<-ggplot2::ggplot(tsne_plot_3) + geom_point(aes(x=x,y=y, color=Genotype, shape=Treatment))+theme_bw()+
  scale_color_manual(values = colors) + xlab("Dim 1")+ ylab("Dim 2")+ theme(legend.position = "none")+ggtitle("22 DAD")+ theme(legend.title = element_text(size=8),legend.text = element_text(size=7) ,text = element_text(size = fontsize_text),axis.text.x = element_text(size = fontsize_text), plot.title = element_text(size=fontsize_text)) 





# plot 4
Genotype <- substr(d4$Sample, start = 1, stop = 8)
Treatment<-substr(d4$Sample, start = 10, stop = 10)

Treatment<-gsub("c", "control", Treatment)
Treatment<-gsub("d", "drought", Treatment)

met_matrix_4 <- as.matrix(d4[,-1])

tsne_out_4 <- Rtsne(met_matrix_4)


tsne_plot_4 <- data.frame(x = tsne_out_4$Y[,1], 
                          y = tsne_out_4$Y[,2])
tsne_plot_4$Genotype<-Genotype
tsne_plot_4$Treatment<-Treatment

p4<-ggplot2::ggplot(tsne_plot_4) + geom_point(aes(x=x,y=y, color=Genotype, shape=Treatment))+theme_bw()+
  scale_color_manual(values = colors) + xlab("Dim 1")+ ylab("Dim 2")+ theme(legend.position = "none")+ggtitle("28 DAD")+ theme(legend.title = element_text(size=8),legend.text = element_text(size=7) ,text = element_text(size = fontsize_text),axis.text.x = element_text(size = fontsize_text), plot.title = element_text(size=fontsize_text)) 





# legend
p0<-ggplot2::ggplot(tsne_plot_4) + geom_point(aes(x=x,y=y, color=Genotype, shape=Treatment))+theme_bw()+
  scale_color_manual(values = colors) + xlab("Dim 1")+ xlab("Dim 2")+ theme(legend.position="bottom")


legend<-get_legend(p0)
p5<-as_ggplot(legend)



legend_y = 0.1
height = (1-legend_y)/2
pos1 = legend_y
pos2 = pos1 + height

tiff("Figs/AdditionalFileY.tiff", width=20, height=22, units="cm", compression = "lzw",res=600)
ggdraw()+
  draw_plot(p1, x = 0, y = pos2, width = 0.5, height = height)+
  draw_plot(p2, x = 0.5, y = pos2, width = 0.5, height = height)+
  draw_plot(p3, x = 0, y = pos1, width = 0.5, height = height) +
  draw_plot(p4, x = 0.5, y = pos1, width = 0.5, height = height)+
  
  draw_plot(legend, x = 0, y = 0, width = 1, height = legend_y)+
  draw_plot_label(label = c("a", "b", "c", "d"), size = 7,
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, pos2, pos2 ))
dev.off()




