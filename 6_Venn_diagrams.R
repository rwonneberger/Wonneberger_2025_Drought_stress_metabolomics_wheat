# joint analysis of data across days etc, outputs from the MetaboAnalyst pipeline

setwd("/Results")
setwd("//storage-al.slu.se/home$/rawo0001/My Documents/Projects/CresWheat/Manuscripts/Metabolome_Drought_stress/Paper_analyses/")
library(data.table)
library("ggVennDiagram")

###################### Venn diagram of DAMs (ttest FDR 0.05) per day ###################### 

d1<-fread("/Results/T_MA_IS_Weight_norm_time_1_trt_all_tol_all_tt.txt")
d2<-fread("/Results/T_MA_IS_Weight_norm_time_2_trt_all_tol_all_tt.txt")
d3<-fread("/Results/T_MA_IS_Weight_norm_time_3_trt_all_tol_all_tt.txt")
d4<-fread("/Results/T_MA_IS_Weight_norm_time_4_trt_all_tol_all_tt.txt")

# Add a day
d1$Day<-"7 DAD"
d2$Day<-"13 DAD"
d3$Day<-"22 DAD"
d4$Day<-"28 DAD"


d1_up<-d1[d1$direction == "up", ]
d1_down<-d1[d1$direction == "down", ]

d2_up<-d2[d2$direction == "up", ]
d2_down<-d2[d2$direction == "down", ]

d3_up<-d3[d3$direction == "up", ]
d3_down<-d3[d3$direction == "down", ]

d4_up<-d4[d4$direction == "up", ]
d4_down<-d4[d4$direction == "down", ]



set_size=3
label_size=2.5

up_list<-list(d1_up$Metabolites, d2_up$Metabolites, d3_up$Metabolites, d4_up$Metabolites)
names(up_list)<-c("7 DAD", "13 DAD", "22 DAD", "28 DAD")

p1<-ggVennDiagram(up_list, label_alpha = 0, set_size = set_size, label_size = label_size)+
  ggplot2::scale_fill_distiller(palette = "Reds", direction=1) +  scale_color_manual(values = c("darkgrey", "darkgrey", "darkgrey", "darkgrey"))+ scale_x_continuous(expand = expansion(mult = .1))


down_list<-list(d1_down$Metabolites, d2_down$Metabolites, d3_down$Metabolites, d4_down$Metabolites)
names(down_list)<-c("7 DAD", "13 DAD", "22 DAD", "28 DAD")

p2<-ggVennDiagram(down_list, label_alpha = 0, set_size = set_size, label_size = label_size)+
  ggplot2::scale_fill_distiller(palette = "Blues", direction=1) +  scale_color_manual(values = c("darkgrey", "darkgrey", "darkgrey", "darkgrey"))+ scale_x_continuous(expand = expansion(mult = .1))


tiff(paste0("Figs/Fig6.tiff"), width=17, height=6, compression = "lzw", res=600, units="cm")
ggarrange(p1, p2, ncol=2, labels="auto")

dev.off()





# function from https://stackoverflow.com/questions/23559371/get-the-list-of-items-in-venn-diagram to get the intersects

getVennOverlap_up <- function(lsvenn = list(Day1 = d1_up$Metabolites, Day2=d2_up$Metabolites, Day3=d3_up$Metabolites, Day4=d4_up$Metabolites)
) {
  
  ItemsList <- gplots::venn(lsvenn, show.plot = FALSE)
  print(lengths(attributes(ItemsList)$intersections))
  return(attributes(ItemsList)$intersections)
}
getVennOverlap_up()


getVennOverlap_down <- function(lsvenn = list(Day1 = d1_down$Metabolites, Day2=d2_down$Metabolites, Day3=d3_down$Metabolites, Day4=d4_down$Metabolites)
) {
  
  ItemsList <- gplots::venn(lsvenn, show.plot = FALSE)
  print(lengths(attributes(ItemsList)$intersections))
  return(attributes(ItemsList)$intersections)
}
getVennOverlap_down()






## Rewatering

getVennOverlap <- function(lsvenn = list(Day3up = d3_up$Metabolites, Day4up=d4_up$Metabolites, Day3down=d3_down$Metabolites, Day4down=d4_down$Metabolites)
) {
  
  ItemsList <- gplots::venn(lsvenn, show.plot = FALSE)
  print(lengths(attributes(ItemsList)$intersections))
  return(attributes(ItemsList)$intersections)
}
getVennOverlap()




