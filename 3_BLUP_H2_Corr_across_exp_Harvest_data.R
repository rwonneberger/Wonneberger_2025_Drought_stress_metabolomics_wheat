# Calculate BLUEs, Spearman correlation and broad sense heritability (H2) for harvest traits

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")


exp1<-fread("Data/Harvest_data_filtered1.txt")
exp2<-fread("Data/Harvest_data_filtered2.txt")

traits_keep<-intersect(names(exp1), names(exp2))
exp1%<>%select(all_of(traits_keep))
exp2%<>%select(all_of(traits_keep))

exp1$exp<-"Exp1"
exp2$exp<-"Exp2"

exp1%<>%relocate(PlantID, Genotype, Treatment, exp)
exp2%<>%relocate(PlantID, Genotype, Treatment, exp)

exp1_c<-exp1%>%filter(Treatment == "control")
exp1_d<-exp1%>%filter(Treatment == "drought")

exp2_c<-exp2%>%filter(Treatment == "control")
exp2_d<-exp2%>%filter(Treatment == "drought")


df_d<-rbind(exp1_d,  exp2_d)
df_c<-rbind(exp1_c,  exp2_c)


goi<-c("PPPW_003", "PPPW_004", "PPPW_007", "PPPW_011", "PPPW_012", "PPPW_016", "PPPW_017", "PPPW_023",
       "PPPW_025", "PPPW_033", "PPPW_034", "PPPW_046")

df_sub_d<-df_d%>%filter(Genotype%in%goi)
df_sub_c<-df_c%>%filter(Genotype%in%goi)




# Blues/Blups
# Code modified from https://github.com/mighster/BLUPs_Heritability/blob/master/BLUP_Tutorial.R

trt ="drought"

if (trt == "control"){
  df<-as.data.frame(df_sub_c)
} else if (trt == "drought"){
  df<-as.data.frame(df_sub_d)
}


nloc = 2
nrep = 20

#create function to identify maximum value in a column
colMax <- function(data) sapply(data, max, na.rm = TRUE)
#create function to identify minimum value in a column
colMin <- function(data) sapply(data, min, na.rm = TRUE)
#create function to identify mean value in a column
colMean <- function(data) sapply(data, mean, na.rm = TRUE)
#create function to identify median value in a column
colMedian <- function(data) sapply(data, median, na.rm = TRUE)
#create function to identify standard dev value in a column
colStdev <- function(data) sapply(data, sd, na.rm = TRUE)



#Create empty data frame for BLUP output
DataOutput <- data.frame(matrix(vector(),12,1, dimnames=list(c(), c("Genotype"))))

#fill empty dataframe with 1-300 so that the cbind will work later on
DataOutput$Genotype <- unique(df[,2]) #fill in Entry numbers


#this empty dataframe is for variance components
DataVarComp <- data.frame()
DataVarCompOutput <- data.frame()
HeritabilityData <- data.frame()
DataOutputBLUEs<-data.frame()
DataOutputBLUPs<-data.frame()


#this empty dataframe is for dropped variance components
drops <- c("var1","var2","sdcor") 

colnames(df)
str(df[,5:ncol(df)])
#take only the columns with numerical data
colnum=c(5:(ncol(df)))


correlations<-data.frame(Trait = "",  Corr = "", Pval ="")


for (i in 1:18){  #this second loop runs through each TRAIT, one at a time
    x=colnum[i]  #set the current [i] column as x
    trait=colnames(df)[x] #sets the current column header as the trait name
    df1_sub <- df #make a copy that we can use for our analysis

    colnames(df1_sub)[x]="y"  #renames the trait variable as "y" for the model analysis below
    
    # Calculate correlations between experiments
    
    genotype_mean<-as.data.table(df1_sub)[, .(ymean=mean(y, na.rm=T)), by=list(Genotype, exp)]
    genotype_mean%<>%arrange(Genotype, exp)
    genotype_exp1<-genotype_mean%>%filter(exp == "Exp1")
    genotype_exp2<-genotype_mean%>%filter(exp == "Exp2")
    
    
    corr_df<-data.frame(Trait = trait , Corr = cor(genotype_exp1$ymean, genotype_exp2$ymean, method="spearman"), Pval =cor.test(genotype_exp1$ymean, genotype_exp2$ymean, method="spearman")$p.value)
    
    correlations<-rbind(correlations, corr_df)
    
    #cor(genotype_exp1$ymean, genotype_exp2$ymean, method="spearman")
    #cor.test(genotype_exp1$ymean, genotype_exp2$ymean, method="spearman")$p.value
    
    #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
    #Knowing variability components allows us to calculate Heritability.
    #Random-effects terms are distinguished by vertical bars or pipes (|) 
    model <- lmer(y ~ (1| Genotype) + (1|exp) + (1|Genotype:exp),df1_sub) #our random effects mixed model
    summary(model)
    
    model_blue <- lmer(y ~  Genotype + (1|exp) + (1|Genotype:exp),df1_sub)
    blues<-as.data.frame(fixef(model_blue))
    #blues$Genotype<-rownames(blues)
    #blues$Genotype[1]<-"PPPW_003"
    #blues$Genotype<-gsub("genotype", "", blues$Genotype)
    blues$BLUE<-blues$`fixef(model_blue)`[1]+blues$`fixef(model_blue)`
    blues$BLUE[1]<-blues$`fixef(model_blue)`[1]
    blues$`fixef(model_blue)`<-NULL
    colnames(blues) <- "BLUE" #rename the BLUP column by the trait in our loop
    #add columns to existing dataframe   
    blues$Trait<-trait
    DataOutputBLUEs <- rbind(DataOutputBLUEs,blues)
    
    blup = coef(model)$Genotype #coef extracts model coefficients from lmer objects returned by modeling functions
    hist(blup[,1]) #plot it out
    blup$Trait<-trait
    names(blup)[1]<-"BLUP"
    DataOutputBLUPs <- rbind(DataOutputBLUPs,blup) #ammends our dataframe with the new BLUP column for the new trait.
    
    
    varComp<-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Genotype
    #Modify variance component df by
    #deleting columns in variance component dataframe (we don't need it)
    varComp<-varComp[ , !(names(varComp) %in% drops)]
    #set the trait name from the loop
    varComp$Trait<-trait
    # Add the Day
    #add columns to existing dataframe
    DataVarComp <- rbind(DataVarComp,varComp) 
    
  }


#Reshape the Blues/Blups to something more readable
DataOutputBLUEs$Genotype<-rep(goi)
DataOutputBLUEOutput <- reshape(DataOutputBLUEs, idvar = c("Genotype"), timevar = "Trait", direction = "wide")

DataOutputBLUPs$Genotype<-rep(goi)
DataOutputBLUPOutput <- reshape(DataOutputBLUPs, idvar = c("Genotype"), timevar = "Trait", direction = "wide")

#reshape our variance components dataframe so that we can run the heritability script
DataVarCompOutput <- reshape(DataVarComp, idvar = c("Trait"), timevar = "grp", direction = "wide")

#the broad sense heritability script
#Taken from Gioia et al. 2017

HeritabilityData <- (DataVarCompOutput[,3]) / ((DataVarCompOutput[,3]) +  (DataVarCompOutput[,2]/nloc) + (((DataVarCompOutput[,5])) / (nrep*nloc)))

#summary statistics
DataColMax <- colMax(df[,colnum])
DataColMin <- colMin(df[,colnum])
DataColMean <- colMean(df[,colnum])
DataColMedian <- colMean(df[,colnum])
DataColStdev <- colStdev(df[,colnum])

#CVg <- ((DataVarCompOutput[,3])) / ((DataColMean)^2)
CVg <- (sqrt(DataVarCompOutput[,3])) / ((DataColMean))
#out <- LSD.test(model,"Entry", p.adj="bonferroni")

#bind the heritability to the variance component data frame

DataVarCompOutput <- cbind(DataVarCompOutput,HeritabilityData,CVg,DataColMin,DataColMax,DataColMean,DataColMedian,DataColStdev)


#output that beast
write.csv(correlations,paste0("Results/Corr_Exp1_Exp2_harvesttraits_", trt, ".csv"))

write.csv(DataVarCompOutput,paste0("Results/H2_Exp1_Exp2_harvesttraits_", trt, ".csv"))
fwrite(DataOutputBLUEOutput,paste0("Results/BLUEs_Exp1_Exp2_harvesttraits_", trt, ".txt"), sep="\t", row.names = F)
fwrite(DataOutputBLUPOutput,paste0("Results/BLUPs_Exp1_Exp2_harvesttraits_", trt, ".txt"), sep="\t", row.names = F)



#######################
# Plot corr between exp

corr_drought<-fread("Results/Corr_Exp1_Exp2_harvesttraits_drought.csv")[-1,]
corr_control<-fread("Results/Corr_Exp1_Exp2_harvesttraits_control.csv")[-1,]
corr_drought$Treatment<-"drought"
corr_control$Treatment<-"control"


corr<-rbind(corr_drought, corr_control)
#corr$Day<-as.numeric(sub("day ", "", corr$Day))
corr$V1<-NULL


###################################################################
# Plot h2

heri_drought<-fread("Results/H2_Exp1_Exp2_harvesttraits_drought.csv")
heri_control<-fread("Results/H2_Exp1_Exp2_harvesttraits_control.csv")
heri_drought$Treatment<-"drought"
heri_control$Treatment<-"control"

heri<-rbind(heri_drought, heri_control)
heri$V1<-NULL

all_df<-left_join(corr, heri)

all_df_order<-all_df%>%filter(Treatment == "control")
all_df_order%<>%
  arrange(desc(Corr)) %>%
  mutate(Trait = factor(Trait, levels = unique(Trait)))%>%select(Trait)

all_df$Trait<-factor(all_df$Trait, levels=all_df_order$Trait)

names(all_df)[c(2, 9)]<-c("Spearman correlation", "Broad sense heritability (H2)")

all_df%<>%select(Trait, Treatment, "Spearman correlation", "Broad sense heritability (H2)")%>%melt()

tiff(paste0("Figs/AdditionalFile9.tiff"), width=8, height=10, compression = "lzw", res=600, units="cm")
ggplot(all_df, aes(x=Trait, y=value, fill=Treatment)) + geom_col(position = "dodge") + facet_wrap(~variable, ncol=1) + theme_classic()+ 
  theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+x_axis_90_degrees + theme(axis.title=element_blank(), legend.position ="bottom") 
dev.off()



all_df<-left_join(corr, heri)
names(all_df)[c(2, 9)]<-c("Spearman correlation", "Broad sense heritability (H2)")

all_df%<>%select(Trait, Treatment, "Spearman correlation", "Broad sense heritability (H2)")%>%melt()
all_df<-dcast(all_df, Trait~Treatment+variable,  value.var =  "value")
all_df<-all_df%>%relocate(Trait , "control_Spearman correlation", "drought_Spearman correlation", `control_Broad sense heritability (H2)` , `drought_Broad sense heritability (H2)` )%>%arrange(desc(`control_Spearman correlation`))
fwrite(all_df, "Results/Suppl_Corr_Heri_harvestdata.txt", sep="\t")
