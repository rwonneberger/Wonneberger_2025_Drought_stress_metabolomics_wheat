# Ranking genotypes by performance under drought

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")
library(factoextra)


# Read in the datasets we need
df_in<-fread("Data/Harvest_data_filtered2.txt")

selected_day=3

if(selected_day == 1){exp_day <- "2021-08-17"}
if(selected_day == 2){exp_day <- "2021-08-23"}
if(selected_day == 3){exp_day <- "2021-09-01"}
if(selected_day == 4){exp_day <- "2021-09-07"}


drought<-fread("Data/All_drought_sd_NA_filt2.txt")%>%filter(timePoint == exp_day)%>%dplyr::select(plotId, genotype, Treatment, biomass, height, compactness, width, hsv_h_mean, hsv_h_brown2green, hsv_h_yellow2green, hsv_h_red2green, lab_b_mean, lab_l_mean, lab_a_mean)
control<-fread("Data//All_control_sd_NA_filt2.txt")%>%filter(timePoint == exp_day)%>%dplyr::select(plotId, genotype, Treatment, biomass, height, compactness, width, hsv_h_mean, hsv_h_brown2green, hsv_h_yellow2green, hsv_h_red2green, lab_b_mean, lab_l_mean, lab_a_mean)
names(drought)[1:3]<-names(control)[1:3]<-c("PlantID", "Genotype", "Treatment")

all_image<-rbind(drought, control)
df_in<-left_join(df_in, all_image)



genotype_mean<-df_in[, .(PH_TopSpike=mean(PH_TopSpike, na.rm=T), 
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



# Make a biplot with all  
res.pca <- prcomp(genotype_mean[, -c(1:2)],  scale = TRUE, center=T)
genotype_mean$Treatmentshort<-substr(genotype_mean$Treatment, 1, 1)
rownames(res.pca$x)<-paste0(genotype_mean$Genotype, "_", genotype_mean$Treatmentshort)

eig.val <- get_eigenvalue(res.pca)
eig.val[1,2]
eig.val[2,2]

tiff("Figs/Fig3.tiff", width=17, height=17,  compression = "lzw", res=600, units="cm")
fviz_pca_biplot(res.pca, alpha.var=0.3,invisible="quali",labelsize=3, label =c("var", "ind"), habillage=genotype_mean$Treatment, title="", repel=T)
dev.off()



# Rank product
dt_mean<-genotype_mean%>%filter(Treatment == "drought")


dt_mean%<>%arrange(desc(GrainNumber))
dt_mean$GrainNumber_Rank<-seq(1:12)
dt_mean%<>%arrange(desc(biomass))
dt_mean$Biomass_Rank<-seq(1:12)
dt_mean%<>%arrange(desc(GrainWeight))
dt_mean$GrainWeight_Rank<-seq(1:12)
dt_mean%<>%arrange(desc(TillerNumberGain))
dt_mean$TillerNumberGain_Rank<-seq(1:12)
dt_mean%<>%arrange(desc(SpikeNumber))
dt_mean$SpikeNumber_Rank<-seq(1:12)


dt_mean$Rank_product<-dt_mean$GrainNumber_Rank*dt_mean$Biomass_Rank*dt_mean$GrainWeight_Rank*dt_mean$TillerNumberGain_Rank*dt_mean$SpikeNumber_Rank
dt_mean%<>%arrange(Rank_product)%>%relocate(Genotype, GrainNumber_Rank, Biomass_Rank, GrainWeight_Rank, TillerNumberGain_Rank, SpikeNumber_Rank, Rank_product)

dt_mean_sub<-dt_mean%>%dplyr::select(Genotype, GrainNumber_Rank, Biomass_Rank, GrainWeight_Rank , TillerNumberGain_Rank, SpikeNumber_Rank, Rank_product )

# Make a biplot with the traits to see if we get the rank product results confirmed
df_biplot<-dt_mean%>%dplyr::select(Genotype, GrainNumber, biomass,GrainWeight ,TillerNumberGain, SpikeNumber )

res.pca <- prcomp(df_biplot[, -c(1)],  scale = TRUE, center=T)
rownames(res.pca$x)<-df_biplot$Genotype

p1<-fviz_pca_biplot(res.pca, label =c("var", "ind"), labelsize = 3, title="", repel=T)



# Make a biplot with loss of trait


genotype_mean$Treatmentshort<-NULL
dt_mean<-genotype_mean%>%filter(Treatment == "drought")
ct_mean<-genotype_mean%>%filter(Treatment == "control")

dt_mean$Genotype == ct_mean$Genotype

calc_loss<-(1- (dt_mean[,-c(1:2)]/ct_mean[,-c(1:2)]))*100
loss<-cbind(ct_mean[,1], calc_loss)

loss<-do.call(data.frame,lapply(loss, function(x) replace(x, is.infinite(x),0)))
names(loss)
loss1<-loss%>%dplyr::select(Genotype, GrainNumber, biomass, GrainWeight , TillerNumberGain, SpikeNumber )

res.pca <- prcomp(loss1[, -c(1)],  scale = T, center=T)
rownames(res.pca$x)<-loss$Genotype

p2<-fviz_pca_biplot(res.pca, label =c("var", "ind"),  labelsize = 3, title="", repel=T)

tiff("Figs/Fig10.tiff", width=17, height=8,  compression = "lzw", res=600, units="cm")
ggarrange(p1, p2, ncol=2, labels="auto")
dev.off()



