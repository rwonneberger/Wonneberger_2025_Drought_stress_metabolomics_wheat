source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")


#df <-  list.files(path = ".", pattern = "PLSR_alltraits_mostsig_") %>% map_df(~fread(.))

classes<-fread("Data/Metabolite_chemical_class.txt")[,c(1,3)]
names(classes)[1]<-"Metabolite"


df2<-fread("PLSR_alltraits_mostsig_control_control_metab_Day3.txt")

df5<-fread("PLSR_alltraits_mostsig_drought_drought_metab_Day3.txt")

df8<-fread("PLSR_alltraits_mostsig_loss_of_trait_control_metab_Day3.txt")

df11<-fread("PLSR_alltraits_mostsig_loss_of_trait_drought_metab_Day3.txt")



df2$Dataset<-"C_C_3"

df5$Dataset<-"D_D_3"

df8$Dataset<-"L_C_3"

df11$Dataset<-"L_D_3"



df1<-bind_rows(df2, df5, df8, df11) # Day 3


df2<-df1

df2$`Corr p value`<-as.character(df2$`Corr p value`)

df2[df2$`Corr p value` == "0",]$`Corr p value` <- "<0.00001"

df2$FDR<-as.character(df2$FDR)

df2[df2$FDR == "0",]$FDR <- "<0.00001"

df2%<>%filter(`P-value` != "*" , `P-value` != "**")

df2<-left_join(df2, classes)

df2$Direction<-ifelse(df2$`Spearman rho` < 0, "Neg", "Pos")

fwrite(df2, "Day3_corr_PLS_results.txt", sep="\t")



table(df2$Class, df2$Dataset, df2$Direction)
table(df2$Class, df2$Dataset)
table(df2$Trait, df2$Dataset)
table(df2$Class, df2$Dataset, df2$Trait)

table(df2$Metabolite, df2$Dataset)

table(df2$Trait)


unique(df2$Metabolite)
