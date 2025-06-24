# Run an ANOVA to identify traits most strongly affected by drought - on a pot basis so an interaction term cna be included


source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

library("factoextra")

# ANOVA of harvest and imaging traits
# Read in postharvest data
harvest<-fread("Data/Harvest_data_filtered2.txt")

# Read in imaging data from last day of drought
ctrl<-fread("Data/All_control_sd_NA_filt2.txt")%>%filter(timePoint == "2021-09-01")%>%dplyr::select(plotId, genotype, Treatment, biomass, height, compactness, width, hsv_h_mean, hsv_h_brown2green, hsv_h_yellow2green, hsv_h_red2green, lab_b_mean, lab_l_mean, lab_a_mean)
drought<-fread("Data/All_drought_sd_NA_filt2.txt")%>%filter(timePoint == "2021-09-01", plotId != "2132TB245")%>%dplyr::select(plotId, genotype, Treatment, biomass, height, compactness, width, hsv_h_mean, hsv_h_brown2green, hsv_h_yellow2green, hsv_h_red2green, lab_b_mean, lab_l_mean, lab_a_mean)
names(ctrl)[1:3]<-names(drought)[1:3]<-c("PlantID", "Genotype", "Treatment")

ctrl_drought<-rbind(ctrl, drought)

harvest<-left_join(harvest, ctrl_drought)

# Run ANOVA
output<-matrix(NA, nrow = length(names(harvest)[-c(1:3)]), ncol = 4)
j=1

for(i in names(harvest)[-c(1:3)]){
  two.way <- aov(as.formula(paste(i, "~ Genotype * Treatment")), data = harvest)
  
  output[j,1]<-i
  output[j,2]<-summary(two.way)[[1]][["Pr(>F)"]][1]
  output[j,3]<-summary(two.way)[[1]][["Pr(>F)"]][2]
  output[j,4]<-summary(two.way)[[1]][["Pr(>F)"]][3]
  j=j+1
}

output<-as.data.frame(output)
names(output)<-c("Trait", "p_Genotype", "p_Treatment", "p_Genotype_Treatment")
output$p_Genotype<-signif(as.numeric(output$p_Genotype), 3)
output$p_Treatment<-signif(as.numeric(output$p_Treatment), 3)
output$p_Genotype_Treatment<-signif(as.numeric(output$p_Genotype_Treatment), 3)
output%<>%arrange(p_Treatment,desc=T)

fwrite(output, "Results/Table1.txt", sep="\t")




# Loss of trait of imaging data per day
# Read in all imaging data (all days)
ctrl<-fread("Data/All_control_sd_NA_filt2.txt")%>%filter(Day != "day 1")%>%dplyr::select(plotId, genotype, Treatment, Day_Int_, Weight_A_g_ , biomass, height, compactness, width, hsv_h_mean, hsv_h_brown2green, hsv_h_yellow2green, hsv_h_red2green, lab_b_mean, lab_l_mean, lab_a_mean)
drought<-fread("Data/All_drought_sd_NA_filt2.txt")%>%filter(plotId != "2132TB245", Day != "day 1")%>%dplyr::select(plotId, genotype, Treatment, Day_Int_, Weight_A_g_ , biomass, height, compactness, width, hsv_h_mean, hsv_h_brown2green, hsv_h_yellow2green, hsv_h_red2green, lab_b_mean, lab_l_mean, lab_a_mean)
names(ctrl)[1:4]<-names(drought)[1:4]<-c("PlantID", "Genotype", "Treatment", "Day")


# Calculate PAW
all_weight<-2015.4
weight_PAW<-1151.17

ctrl$PAW<-(ctrl$Weight_A_g_-all_weight)/weight_PAW*100
drought$PAW<-(drought$Weight_A_g_-all_weight)/weight_PAW*100

#Remove lines with odd PAW values
ctrl%<>%filter(PAW > 0)
drought%<>%filter(PAW > 0)

df_cd<-rbind(ctrl, drought)
df_cd$Day<-df_cd$Day+ 20
df_cd$Weight_A_g_<-NULL

anova_df_cd<-df_cd


#ANOVA time course series

df_melt<-melt(df_cd, id.vars = c("PlantID", "Genotype", "Treatment", "Day"))
df_melt$value<-as.numeric(df_melt$value)

df_cd<-as.data.table(dcast(df_melt,  PlantID + Genotype + Treatment  ~ variable + Day, value.var = "value"))

df_cd<-df_cd[, -1]

## Which traits are most strongly affected by drought?
output<-matrix(NA, nrow = length(names(df_cd)[-c(1:3)]), ncol = 4)
j=1

for(i in names(df_cd)[-c(1:3)]){
  two.way <- aov(as.formula(paste(i, "~ Genotype * Treatment")), data = df_cd)
  
  output[j,1]<-i
  output[j,2]<-summary(two.way)[[1]][["Pr(>F)"]][1]
  output[j,3]<-summary(two.way)[[1]][["Pr(>F)"]][2]
  output[j,4]<-summary(two.way)[[1]][["Pr(>F)"]][3]
  j=j+1
}

output<-as.data.table(output)
names(output)<-c("Trait", "p_Genotype", "p_Treatment", "p_Genotype_Treatment")
output$p_Genotype<-signif(as.numeric(output$p_Genotype), 3)
output$p_Treatment<-signif(as.numeric(output$p_Treatment), 3)
output$p_Genotype_Treatment<-signif(as.numeric(output$p_Genotype_Treatment), 3)

names(output)[1]<-"Trait_Day"

fwrite(output, "Results/AdditionalFile11.txt", sep="\t")






# Loss of trait of imaging data per genotype and day

genotype_mean<-anova_df_cd[, .(biomass=mean(biomass, na.rm=T),
                           Height=mean(height, na.rm=T),
                           compactness=mean(compactness, na.rm=T),
                           width=mean(width, na.rm=T),
                           hsv_h_mean=mean(hsv_h_mean, na.rm=T),
                           hsv_h_brown2green=mean(hsv_h_brown2green, na.rm=T),
                           hsv_h_yellow2green=mean(hsv_h_yellow2green, na.rm=T),
                           hsv_h_red2green=mean(hsv_h_red2green, na.rm=T),
                           lab_b_mean=mean(lab_b_mean, na.rm=T),
                           lab_l_mean=mean(lab_l_mean, na.rm=T),
                           lab_a_mean=mean(lab_a_mean, na.rm=T),
                           PAW=mean(PAW, na.rm=T)
), by=list(Genotype, Treatment, Day)]

ct_mean<-genotype_mean%>%filter(Treatment == "control")
dt_mean<-genotype_mean%>%filter(Treatment == "drought")

table(dt_mean$Genotype == ct_mean$Genotype)


calc_loss<-round((1- (dt_mean[,-c(1:3, 15)]/ct_mean[,-c(1:3, 15)]))*100, 1)
genotype_mean<-cbind(ct_mean[,c(1,3)], calc_loss)
genotype_mean$DAS<-genotype_mean$Day-20
names(genotype_mean)[2]<-"DAD"
genotype_mean%<>%relocate(Genotype, DAS, DAD)%>%arrange(DAS, Genotype )

fwrite(genotype_mean, "Results/Loss_of_trait_per_day_geno.txt", sep="\t")


# Loss of trait of imaging data per day

genotype_mean_mean<-genotype_mean[,.(biomass=mean(biomass), 
                                     Height=mean(Height),
                                     compactness=mean(compactness), 
                                     width=mean(width), 
                                     hsv_h_mean=mean(hsv_h_mean), 
                                     hsv_h_brown2green=mean(hsv_h_brown2green), 
                                     hsv_h_yellow2green=mean(hsv_h_yellow2green), 
                                     hsv_h_red2green=mean(hsv_h_red2green), 
                                     lab_b_mean=mean(lab_b_mean), 
                                     lab_l_mean=mean(lab_l_mean), 
                                     lab_a_mean=mean(lab_a_mean) 
),
by=.(DAS, DAD)]

fwrite(genotype_mean_mean, "Results/AdditionalFile12.txt", sep="\t")




# Loss of trait of harvest data

df_in<-fread("Data/Harvest_data_filtered2.txt")

genotype_mean<-df_in[, .(PH_TopSpike=mean(PH_TopSpike, na.rm=T), 
                         PH_BaseSpike=mean(PH_BaseSpike, na.rm=T), 
                         PH_FlagLeaf=mean(PH_FlagLeaf, na.rm=T),
                         PH_TopNode=mean(PH_TopNode, na.rm=T),
                         #PeduncleLength=mean(PeduncleLength, na.rm=T),
                         #Length_SpikeBase_FlagLeaf=mean(Length_SpikeBase_FlagLeaf, na.rm=T),
                         SpikeCulmRatio=mean(SpikeCulmRatio, na.rm=T),
                         SpikeNumber=mean(SpikeNumber, na.rm=T),
                         #InfertileSpikes=mean(InfertileSpikes, na.rm=T),
                         PlantWeight=mean(PlantWeight, na.rm=T),
                         SpikeWeight=mean(SpikeWeight, na.rm=T),
                         #SpikePlantRatioWeight=mean(SpikePlantRatioWeight, na.rm=T),
                         #GrainPlantRatioWeight=mean(GrainPlantRatioWeight, na.rm=T),
                         SpikeLength=mean(SpikeLength, na.rm=T),
                         AwnLength=mean(AwnLength, na.rm=T),
                         #LengthLastInternode=mean(LengthLastInternode, na.rm=T),
                         NumberFertileSpikes=mean(NumberFertileSpikes, na.rm=T),
                         GrainNumber=mean(GrainNumber, na.rm=T),
                         GrainWeight=mean(GrainWeight, na.rm=T),
                         TGW=mean(TGW, na.rm=T),
                         GrainDensity=mean(GrainDensity, na.rm=T),
                         GrainArea=mean(GrainArea, na.rm=T),
                         GrainWidth=mean(GrainWidth, na.rm=T),
                         GrainLength=mean(GrainLength, na.rm=T),
                         GrainLengthWidth=mean(GrainLengthWidth, na.rm=T),
                         SpikeDensity=mean(SpikeDensity, na.rm=T),
                         RachisNodeNumber=mean(RachisNodeNumber, na.rm=T),
                         FertileGrains=mean(FertileGrains, na.rm=T),
                         InfertileGrains=mean(InfertileGrains, na.rm=T),
                         TillerNumber1=mean(TillerNumber1, na.rm=T),
                         TillerNumber2=mean(TillerNumber2, na.rm=T),
                         TillerNumberGain=mean(TillerNumberGain, na.rm=T)
), by=list(Genotype, Treatment)]


ct_mean<-genotype_mean%>%filter(Treatment == "control")
dt_mean<-genotype_mean%>%filter(Treatment == "drought")

table(dt_mean$Genotype == ct_mean$Genotype)


calc_loss<-(1- (dt_mean[,-c(1:2)]/ct_mean[,-c(1:2)]))*100
genotype_mean<-cbind(ct_mean[,1:2], calc_loss)

fwrite(genotype_mean, "Data/Harvest_data_filtered_loss_of_trait2.txt")





