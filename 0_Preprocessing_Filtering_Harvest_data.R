# Analysis of Marvin, tiller data in exp2

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")
library(lme4)
library(multcompView)
library(psych)

exp = "Exp2"

setwd("//storage-al.slu.se/home$/rawo0001/My Documents/Projects/CresWheat/Manuscripts/Metabolome_Drought_stress/Scripttest/")

df<-fread("Data/Harvest_data2.txt")


df_melt<-melt(df, id.vars = c("PlantID", "Genotype", "Treatment"))
df_melt$value<-as.numeric(df_melt$value)
df_filt<-df_melt%>%
  group_by( Genotype, Treatment, variable)%>%
  mutate(value = replace(value, abs(value - mean(value)) > 3*sd(value), NA))%>%
  mutate(num_na = sum(is.na(value))) %>% 
  filter(!any(num_na >=18)) %>%
  select(-num_na)%>%ungroup




df_in<-as.data.table(dcast(df_filt,  PlantID + Genotype + Treatment ~ variable, value.var = "value"))

df_in$PeduncleLength<-df_in$Length_SpikeBase_FlagLeaf<-df_in$SpikePlantRatioWeight<-df_in$GrainPlantRatioWeight<-df_in$LengthLastInternode<-NULL

fwrite(df_in, "Data/Harvest_data_filtered2.txt")
