# PLSR exp2 - genotype basis - only the plots sampled for metabolites on a given day

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")
library(pls)
library(ggpubr)


selected_day<-3

# Selected phenotypic dataset: Choose between control, drought and loss_of_trait
dataset<-"control"

# selected metabolic dataset
dataset_metab<-"control"

# Read in the datasets we need
df_in<-fread("Harvest_data_filtered2.txt")

anno<-read.csv("Annotation.csv")[, c(2:5,9)]

metab<-fread(paste0("Results/Normalized_T_MA_IS_Weight_norm_time_", selected_day, "_trt_all_tol_all.txt"))


if (dataset != "loss_of_trait"){
  anno%<>%filter(day == selected_day, condition == dataset)
} else if (dataset == "loss_of_trait"){
  anno%<>%filter(day == selected_day)
}


if(selected_day == 1){exp_day <- "day 8"}
if(selected_day == 2){exp_day <- "day 14"}
if(selected_day == 3){exp_day <- "day 23"}
if(selected_day == 4){exp_day <- "day 29"}

ctrl<-fread("Data/All_control_sd_NA_filt2.txt")%>%filter(Day == exp_day)%>%dplyr::select(plotId, genotype, Treatment, biomass, height, compactness, width, hsv_h_mean, hsv_h_brown2green, hsv_h_yellow2green, hsv_h_red2green, lab_b_mean, lab_l_mean, lab_a_mean)
drought<-fread("Data/All_drought_sd_NA_filt2.txt")%>%filter(Day == exp_day, plotId != "2132TB245")%>%dplyr::select(plotId, genotype, Treatment, biomass, height, compactness, width, hsv_h_mean, hsv_h_brown2green, hsv_h_yellow2green, hsv_h_red2green, lab_b_mean, lab_l_mean, lab_a_mean) # Remove 2132TB245 because it's always been an outlier
names(ctrl)[1:3]<-names(drought)[1:3]<-c("PlantID", "Genotype", "Treatment")

ctrl_drought<-rbind(ctrl, drought)

df_in<-left_join(df_in, ctrl_drought)

samples_keep<-intersect(metab$Sample, anno$Unique_name)

anno%<>%filter(Unique_name%in%samples_keep)%>%arrange(Unique_name)
metab%<>%filter(Sample%in%samples_keep)%>%arrange(Sample)

table(anno$Unique_name == metab$Sample)

df<-merge(metab, anno, by.x="Sample", by.y= "Unique_name")
df$Sample<-df$name_geno<-df$condition<-df$day<-NULL
df%<>%arrange(name_pot)%>%relocate(name_pot)

pot_keep<-intersect(df$name_pot, df_in$PlantID)

df_in%<>%filter(PlantID%in%pot_keep)%>%arrange(PlantID)

table(df$name_pot == df_in$PlantID)


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

# Select trait of interest
trait<-names(genotype_mean)[-c(1:2)]

if (dataset== "control"){
  genotype_mean%<>%filter(Treatment == "control")
  
} else if (dataset== "drought"){
  genotype_mean%<>%filter(Treatment == "drought")
  
} else if (dataset == "loss_of_trait") {
  
  ct_mean<-genotype_mean%>%filter(Treatment == "control")
  dt_mean<-genotype_mean%>%filter(Treatment == "drought")
  
  table(dt_mean$Genotype == ct_mean$Genotype)
  
  
  calc_loss<-(1- (dt_mean[,-c(1:2)]/ct_mean[,-c(1:2)]))*100
  genotype_mean<-cbind(ct_mean[,1:2], calc_loss)
  
}



metab$Genotype<-substring(metab$Sample, 1, nchar(metab$Sample)-6)

metab$trt<-substring(metab$Sample, 10, nchar(metab$Sample)-4)

if (dataset_metab == "control"){
  metab%<>%filter(trt == "c")
} else if (dataset_metab == "drought"){
  metab%<>%filter(trt == "d")
}

metab$trt<-NULL

metab_mean<-metab[,-1]%>%
  group_by(Genotype) %>%
  summarise_all("mean")%>%
  arrange(Genotype)


genotype_mean%<>%arrange(Genotype)

table(metab_mean$Genotype == genotype_mean$Genotype)

## For correlation analysis
correlation_table<-merge(metab_mean, genotype_mean[,-2])

rownames(correlation_table)<-correlation_table$Genotype
correlation_table$Genotype<-NULL

last_metabolite<-(ncol(metab_mean)-1)

cor_mat<-cor(correlation_table[,1:last_metabolite],correlation_table[,(last_metabolite+1):ncol(correlation_table)] , use = "complete.obs", method="spearman")


# matrix of the p-value of the correlation
p.mat <- cor.mtest.selected(correlation_table, last_metabolite, "spearman")




sig_out<-data.frame()


traitlist<-names(genotype_mean)[-c(1,2)]

for(trait in traitlist){

sub_df<-genotype_mean%>%dplyr::select(all_of(trait))
meta_mat<-as.matrix(metab_mean[,-1])

pls_df<-data.frame(sub_df, I(meta_mat))
rownames(pls_df)<-genotype_mean$Genotype


######### Do the modeling

rank_product_df<-data.frame(Metabolites = colnames(meta_mat))
prediction_df<-data.frame()
q_vec<-c()
pearson_vec<-c()
j=1

## Randomly subdivide the dataset into 6 parts without replacement
randnum<-as.list(cvsegments(dim(pls_df)[1], 4))

for (i in randnum){
  randvec<-unlist(randnum)
  randvec<-randvec[!randvec%in%i]
  
  metabTrain<-pls_df[randvec,] # Make the training set
  metabTest<-pls_df[i,] # Make the test set
  
  ## Train the model on the training set and predict the test set
  out1 <- plsr(as.formula(paste(trait, "~ meta_mat")), center=T,scale=T,data = metabTrain, validation = "LOO")
  
  ## Select number of variables
  cv = RMSEP(out1)
  dims <- which.min(cv$val[estimate = "adjCV", , ]) - 1
  
  if (dims == 0){
    dims = 1
  }
  
  ncomp.onesigma <- selectNcomp(out1, method = "onesigma", plot = TRUE)
  
  out1 <- plsr(as.formula(paste(trait, "~ meta_mat")), ncomp = dims, center=T,scale=T,data = pls_df, validation = "none")
  
  ## Predict the test set
  prediction <- as.data.frame(predict(out1, ncomp = dims, newdata = metabTest))
  prediction$Plot<-rownames(prediction)
  names(prediction)[1]<-"predicted_value"
  prediction_df<-rbind(prediction_df, prediction )
  
  ## Calculate Q2
  pls_trait<-data.frame(pls_df[,1])
  pls_trait$Plot<-rownames(pls_df)
  
  pls_trait<-merge(pls_trait, prediction)
  
  ## Calc Q2 
  PRESS <- sum((pls_trait[,2]-pls_trait$predicted_value)^2, na.rm=T) # Predictive residual error sum of squares
  TSS <- sum((pls_trait[,2]-mean(pls_trait$predicted_value, na.rm=T))^2, na.rm=T) # Total sum of squares
  Q2<- 1 - PRESS/TSS
  Q2
  q_vec<-c(q_vec, Q2)
  
  ## Calc R2
  #correlation<-cor.test(pls_trait[,2], pls_trait$predicted_value, method="pearson")
  #correlation$estimate
  
  #pearson_vec<-c(pearson_vec, correlation$estimate)
  
  
  ## Rank the predictors for each model based on their absolute regression coefficient
  coef_df<-data.frame(coef(out1))
  names(coef_df)<-"coefficients"
  coef_df$Metabolites<-rownames(coef_df)
  coef_df%<>%arrange(desc(abs(coefficients)))
  coef_df$Rank<-seq(1:length(coef_df$Metabolites))/1000 #divide by 1000 to keep the numbers small, otherwise R will complain
  
  coef_df<-coef_df[, 2:3]
  names(coef_df)[2]<-paste0("Rank", j)
  
  ## Bind all together
  rank_product_df<-merge(rank_product_df, coef_df, by ="Metabolites")
  j=j+1
}

## Calc the rank-product

outrank<- rank_product_df[,-1] %>% 
  mutate(prod = reduce(., `*`))

outrank$Metabolites<-rank_product_df$Metabolites
outrank%<>%arrange(prod)

mean(q_vec)
mean(pearson_vec)

## plot original vs predicted values
pls_trait<-data.frame(pls_df[,1])
pls_trait$Plot<-rownames(pls_df)

pls_trait<-merge(pls_trait, prediction_df)

#plot(pls_trait[,2], pls_trait$predicted_value)



names(pls_trait)[2]<-"observed_value"


#Plot the top metabolite(s) against phenotype
#number_met<-20


output<-matrix(NA, nrow = dim(outrank)[1], ncol = 5)

for (i in 1:dim(outrank)[1]){
  
  
  met<-outrank[i, ncol(outrank)]
  
  
  table(metab_mean$Genotype == genotype_mean$Genotype)
  sub_plot<-data.frame(Metabolite = metab_mean%>%dplyr::select(all_of(met)), Trait = genotype_mean%>%dplyr::select(all_of(trait)))
  names(sub_plot)<-c("Metabolite", "Trait")
  
  output[i,1]<-met
  output[i,2]<-outrank[i, (ncol(outrank)-1)]
  output[i,3]<-round(cor(sub_plot$Metabolite, sub_plot$Trait, use = "complete.obs", method="spearman"), 10)
  output[i,4]<-round(cor.test(sub_plot$Metabolite, sub_plot$Trait, use = "complete.obs", method="spearman")$p.value, 10)
  
  
  
  
  Model <- aov(Trait ~Metabolite, data = sub_plot)
  Model$residuals
  Residuals <- Model$residuals
  SumResSquared <- sum(Residuals^2)
  TotalSumSquares <- sum((sub_plot$Trait - mean(sub_plot$Trait))^2)
  RSquared <- 1 - (SumResSquared/TotalSumSquares)
  RSquared
  
  output[i,5]<-round(RSquared, 2)
  
  

}
output<-as.data.frame(output)
names(output)<-c("Metabolite", "Rankprod", "Spearman corr", "Corr p value",  "R2")

output$FDR<-p.adjust(output$`Corr p value`, method="fdr")
output%<>%relocate(Metabolite, Rankprod, `Spearman corr`, `Corr p value`, FDR, R2)
output$Rank<-seq(1:length(output$Metabolite))


output$Rankprod<-as.numeric(output$Rankprod)
output$`Spearman corr`<-as.numeric(output$`Spearman corr`)
output$`Corr p value`<-as.numeric(output$`Corr p value`)
output$R2<-as.numeric(output$R2)

fwrite(output, paste0("PLS_Met_", trait, "_", dataset, "_", dataset_metab, "_Day", selected_day, ".txt"), sep="\t")


significant<-output%>%filter(`Corr p value` <= 0.05)


if (dim(significant)[1] != 0){
  significant$trait<-trait
  sig_out<-rbind(sig_out, significant)
}
}


sig_out$Sig<-"ns"
sig_out[sig_out$`Corr p value` <= 0.05,]$Sig<-"*"
sig_out[sig_out$`Corr p value` <= 0.01,]$Sig<-"**"
sig_out[sig_out$`Corr p value` <= 0.001,]$Sig<-"***"
sig_out[sig_out$`Corr p value` <= 0.0001,]$Sig<-"****"
sig_out[sig_out$`Corr p value` <= 0.00001,]$Sig<-"*****"


sig_out$Spearman<-paste0(round(as.numeric(sig_out$`Spearman corr`), 3), sig_out$Sig)
sig_out1<-sig_out

sig_out1$`FDR-adjusted P-value`<-"ns"
sig_out1[sig_out1$FDR <= 0.05,]$`FDR-adjusted P-value`<-"*"
sig_out1[sig_out1$FDR <= 0.01,]$`FDR-adjusted P-value`<-"**"
sig_out1[sig_out1$FDR <= 0.001,]$`FDR-adjusted P-value`<-"***"
sig_out1[sig_out1$FDR <= 0.0001,]$`FDR-adjusted P-value`<-"****"
sig_out1[sig_out1$FDR <= 0.00001,]$`FDR-adjusted P-value`<-"*****"
sig_out1<-sig_out1%>%arrange(`Corr p value`)%>%dplyr::select(Metabolite, trait, `Spearman corr`, `Corr p value`, Sig,FDR, `FDR-adjusted P-value`, R2, Rank)
names(sig_out1)<-c("Metabolite", "Trait", "Spearman rho","Corr p value", "P-value","FDR", "FDR-adjusted P-value", "R2", "Rank")

fwrite(sig_out1, paste0("PLSR_alltraits_mostsig_", dataset, "_", dataset_metab, "_metab_Day", selected_day, ".txt"), sep="\t")

sig_out_wide<-as.data.table(dcast(sig_out,  Metabolite ~ trait, value.var = "Spearman"))
sig_out_wide$count<-rowSums(!is.na(sig_out_wide[,-1]))

fwrite(sig_out_wide, paste0("PLSR_alltraits_", dataset, "_", dataset_metab, "_metab_Day", selected_day, ".txt"), sep="\t")
