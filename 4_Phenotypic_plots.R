# Plot phenotypic traits under drought and control and calculate loss of trait

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")
library(lme4)
library(multcompView)
library(psych)


df<-fread("Data/Harvest_data_filtered2.txt")


ctrl<-fread("Data/All_control_sd_NA_filt2.txt")%>%filter(timePoint == "2021-09-01")%>%dplyr::select(plotId, genotype, Treatment, biomass, height, compactness, width, hsv_h_mean, hsv_h_brown2green, hsv_h_yellow2green, hsv_h_red2green, lab_b_mean, lab_l_mean, lab_a_mean)
drought<-fread("Data/All_drought_sd_NA_filt2.txt")%>%filter(timePoint == "2021-09-01", plotId != "2132TB245")%>%dplyr::select(plotId, genotype, Treatment, biomass, height, compactness, width, hsv_h_mean, hsv_h_brown2green, hsv_h_yellow2green, hsv_h_red2green, lab_b_mean, lab_l_mean, lab_a_mean) # Remove 2132TB245 bc outlier
names(ctrl)[1:3]<-names(drought)[1:3]<-c("PlantID", "Genotype", "Treatment")

ctrl_drought<-rbind(ctrl, drought)

df<-left_join(df, ctrl_drought)


genotype_mean<-df[, .(PH_TopSpike=mean(PH_TopSpike, na.rm=T), 
                      PH_BaseSpike=mean(PH_BaseSpike, na.rm=T), 
                      PH_FlagLeaf=mean(PH_FlagLeaf, na.rm=T),
                      PH_TopNode=mean(PH_TopNode, na.rm=T),
                      #PeduncleLength=mean(PeduncleLength, na.rm=T),
                      #Length_SpikeBase_FlagLeaf=mean(Length_SpikeBase_FlagLeaf, na.rm=T),
                      SpikeCulmRatio=mean(SpikeCulmRatio, na.rm=T),
                      SpikeNumber=mean(SpikeNumber, na.rm=T),
                      InfertileSpikes=mean(InfertileSpikes, na.rm=T),
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
                      TillerNumberGain=mean(TillerNumberGain, na.rm=T),
                      biomass=mean(biomass, na.rm=T),
                      Height=mean(height, na.rm=T),
                      compactness=mean(compactness, na.rm=T),
                      width=mean(width, na.rm=T),
                      hsv_h_mean=mean(hsv_h_mean, na.rm=T),
                      hsv_h_brown2green=mean(hsv_h_brown2green, na.rm=T),
                      hsv_h_yellow2green=mean(hsv_h_yellow2green, na.rm=T),
                      hsv_h_red2green=mean(hsv_h_red2green, na.rm=T),
                      lab_b_mean=mean(lab_b_mean, na.rm=T),
                      lab_l_mean=mean(lab_l_mean, na.rm=T),
                      lab_a_mean=mean(lab_a_mean, na.rm=T)
), by=list(Genotype, Treatment)]

ct_mean<-genotype_mean%>%filter(Treatment == "control")
dt_mean<-genotype_mean%>%filter(Treatment == "drought")

table(dt_mean$Genotype == ct_mean$Genotype)


calc_loss<-round((1- (dt_mean[,-c(1:2)]/ct_mean[,-c(1:2)]))*100, 1)
genotype_mean<-cbind(ct_mean[,1], calc_loss)


loss_long<-melt(genotype_mean)
names(loss_long)[3]<-"Perc_loss"


df_melt<-melt(df, id.vars = c("PlantID", "Genotype", "Treatment"))
df_melt$value<-as.numeric(df_melt$value)

df_out<-data.frame()

for (i in unique(df_melt$variable)){
  for (j in unique(df_melt$Genotype)){
    
    sub_df<-df_melt%>%filter(variable == i, Genotype == j)
    
    sub_df$Treatment<-as.factor(sub_df$Treatment)
    
    
    #Runs Tukeys HSD test on all line/year combos to determine which are significantly different and assigns letters for comparison
    exp_tukey <- TukeyHSD(exp_aov <- aov(value ~ Treatment, data = sub_df))
    exp_tukey_result <- data.frame(exp_tukey$Treatment)
    
    
    #Create a matrix to hold p-values
    p_mat <- matrix(NA,length(levels(sub_df$Treatment)),length(levels(sub_df$Treatment)))
    
    
    #Extract p-values from Tukey HSD
    p_mat[lower.tri(p_mat)] <- rownames(exp_tukey_result)
    p_mat[lower.tri(p_mat)] <- round(exp_tukey_result$p.adj,3)
    p_mat[upper.tri(p_mat)] <- t(p_mat[lower.tri(p_mat)])
    rownames(p_mat) <- levels(sub_df$Treatment)
    colnames(p_mat) <- levels(sub_df$Treatment)
    
    #Assign letters to indicate if groups are significantly different (adjusted p-value <0.01)
    letters <- multcompLetters(p_mat,threshold=0.01)
    letters <- data.frame(letters$Letters)
    letters$Treatment <- row.names(letters)
    letters$Genotype<-j
    letters$variable<-i
    
    #Merge letters with trait data
    df_subset <- merge(df_melt,letters,by=c("Genotype", "Treatment", "variable"))
    df_out<-rbind(df_out, df_subset)
  }
}

df_out$value<-as.numeric(df_out$value)


# Select most interesting traits and plot them

df_filt_1<-df_out%>%filter(variable%in%c("TGW", "GrainWeight", "GrainArea", "GrainWidth", "GrainNumber", "NumberFertileSpikes", "TillerNumberGain", "biomass"))
df_filt_1[df_filt_1$variable == "biomass",]$value <-df_filt_1[df_filt_1$variable == "biomass",]$value/1000000

helper_df<-as.data.table(df_filt_1)[,.(maxval=max(value, na.rm=T)),by=.(variable)]

df_filt_1<-left_join(df_filt_1, helper_df)

df_filt_1<-merge(df_filt_1, loss_long, by=c("Genotype", "variable"))
df_filt_1$Perc_loss<-paste0(df_filt_1$Perc_loss, "%")
df_filt_1[df_filt_1$Treatment == "control",]$Perc_loss<-""

df_filt_1$variable<-as.character(df_filt_1$variable)
df_filt_1[df_filt_1$variable == "biomass",]$variable <- "Biomass at 22 DAD"


scaleFUN <- function(x) sprintf("%.3f", x)

tiff("Figs/Fig2.tiff", width = 17, height=20 ,compression = "lzw",units="cm", res=600)
ggplot(df_filt_1, aes(x=Treatment, y=value)) + 
  geom_violin()+ 
  geom_boxplot(width=0.3, fill="#D55E00") +
  geom_jitter(width=0.2, size=0.2)+ 
  facet_grid(variable~Genotype, scales="free_y") + 
  geom_text(aes(y=maxval*1.3,label=paste0(letters.Letters, "\n", Perc_loss)), size=2, vjust="inward") +   
  theme_classic()+ 
  theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA), strip.text = element_text(size = 5.3), text = element_text(size = 7),
        axis.text.x = element_text(size=7))+ 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
ylab("Phenotypic value")
dev.off()

