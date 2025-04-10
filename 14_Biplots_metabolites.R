# Corr metabolites with each other

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")
library(corrplot)
library(factoextra)
library(multcompView)
library(psych)

########## Corr per genotype, either ctrl or drt ########## 

lookuptable<-data.frame(Day = c(1:4), DAD = c(7, 13, 22, 28))

plotlist<-list()
counter<-1

for (i in 1:4){
  for (j in c("control", "drought")){
    
      metab<-fread(paste0("Results/Normalized_T_MA_IS_Weight_norm_time_", i, "_trt_", j, "_tol_all.txt"))

    
    names(metab)[1]<-"Unique_name"
    
    metab$genotype<-substring(metab$Unique_name, 1, nchar(metab$Unique_name)-6)
    metab$day<-substring(metab$Unique_name, 12, nchar(metab$Unique_name)-2)
    metab$trt<-substring(metab$Unique_name, 10, nchar(metab$Unique_name)-4)
    
    metab%<>%filter(trt == substring(j, 1, 1), day == i)%>%relocate(genotype)
    
    metab$day<-metab$trt<-metab$Unique_name<-NULL
    
    metab_mean<-metab%>%
      group_by(genotype) %>%
      summarise_all("mean")%>%
      arrange(genotype)
    
    
    res.pca <- prcomp(metab_mean[,-1],  scale = T, center=T)
    rownames(res.pca$x)<-metab_mean$genotype
    
    plotlist[[counter]]<-fviz_pca_biplot(res.pca, label =c( "ind"), title=paste0(lookuptable[lookuptable$Day == i,]$DAD, " DAD ", j), repel = TRUE, alpha.var=0.3, labelsize = 3)+
      theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"))
    counter=counter+1
    
  }
}




tiff("Figs/AdditionalFile16.tiff", width = 17, height=20 ,units="cm", res=600, compression = "lzw")
ggarrange(plotlist[[1]], plotlist[[2]], plotlist[[3]], plotlist[[4]], plotlist[[5]], plotlist[[6]], plotlist[[7]], plotlist[[8]], ncol=2, nrow=4, labels = c("a", "b", "c", "d", "e", "f", "g", "h"))
dev.off()


