# filtering of raw post-harvest data

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")
library(lme4)
library(multcompView)
library(psych)


# Experiment 1

df<-fread("Harvest_data1.txt") # Additional file 4 sheet 2

df%<>%dplyr::select(-PeduncleLength, -Length_SpikeBase_FlagLeaf, -SpikePlantRatioWeight)

#Select the genotypes of interest from the second exp
goi<-c("PPPW_003", "PPPW_004", "PPPW_007", "PPPW_011", "PPPW_012", "PPPW_016", "PPPW_017", "PPPW_023",
       "PPPW_025", "PPPW_033", "PPPW_034", "PPPW_046")

df_melt<-melt(df, id.vars = c("PlantID", "Genotype", "Treatment"))
df_melt$value<-as.numeric(df_melt$value)
df_filt<-df_melt%>%
  group_by("PlantID", "Genotype", "Treatment", "variable")%>%
  mutate(filt_val = replace(value, abs(value - mean(value)) > 2*sd(value), NA))


df_in<-as.data.table(dcast(df_filt,  PlantID + Genotype + Treatment ~ variable, value.var = "value"))

df_in$PeduncleLength<-df_in$Length_SpikeBase_FlagLeaf<-df_in$SpikePlantRatioWeight<-df_in$GrainPlantRatioWeight<-df_in$LengthLastInternode<-NULL

fwrite(df_in, "Harvest_data_filtered1.txt") 





# Experiment 2

df<-fread("Harvest_data2.txt") # Additional file 4 sheet 1

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

fwrite(df_in, "Harvest_data_filtered2.txt") 
