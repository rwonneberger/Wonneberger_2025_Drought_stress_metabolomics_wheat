# Calculate BLUEs, Spearman correlation and broad sense heritability (H2) for imaging traits

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")
library(lme4)
library(lmerTest)


exp1_c<-fread("Data/All_control_sd_NA_filt1.txt")
exp1_d<-fread("Data/All_drought_sd_NA_filt1.txt")

exp2_c<-fread("Data/All_control_sd_NA_filt2.txt")
exp2_d<-fread("Data/All_drought_sd_NA_filt2.txt")

exp1_d$exp<-"Exp1"
exp1_c$exp<-"Exp1"

exp2_d$exp<-"Exp2"
exp2_c$exp<-"Exp2"

exp1_d$trt<-"Drought"
exp2_d$trt<-"Drought"

exp1_c$trt<-"Control"
exp2_c$trt<-"control"

df_d<-rbind(exp1_d,  exp2_d)
df_c<-rbind(exp1_c,  exp2_c)


df_d$Day_Int_<-df_d$Day_Int_ + 20
df_c$Day_Int_<-df_c$Day_Int_ + 20

goi<-c("PPPW_003", "PPPW_004", "PPPW_007", "PPPW_011", "PPPW_012", "PPPW_016", "PPPW_017", "PPPW_023",
       "PPPW_025", "PPPW_033", "PPPW_034", "PPPW_046")

df_sub_d<-df_d%>%filter(genotype%in%goi)
df_sub_c<-df_c%>%filter(genotype%in%goi)


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
DataOutput$Day <- "" #fill in Day of data collection

#this empty dataframe is for variance components
DataVarComp <- data.frame()
DataVarCompOutput <- data.frame()
HeritabilityData <- data.frame()
DataOutputBLUEs<-data.frame()
DataOutputBLUPs<-data.frame()


#this empty dataframe is for dropped variance components
drops <- c("var1","var2","sdcor") 

colnames(df)
str(df[,14:ncol(df)])
#take only the columns with numerical data
colnum=c(14:(ncol(df)-2))


df%<>%filter(Day != "day 1", Day != "day 24")

correlations<-data.frame(Trait = "", Day = "", Corr = "", Pval ="")


for (i in 1:11){  #this second loop runs through each TRAIT, one at a time
  for (j in (unique(df$Day_Int_ ))){
    x=colnum[i]  #set the current [i] column as x
    trait=colnames(df)[x] #sets the current column header as the trait name
    df1 <- df #make a copy that we can use for our analysis
    df1_sub<-as.data.frame(df%>%filter(Day_Int_  == j))
    colnames(df1_sub)[x]="y"  #renames the trait variable as "y" for the model analysis below
    
    # Calculate correlations between experiments
    
    genotype_mean<-as.data.table(df1_sub)[, .(ymean=mean(y, na.rm=T)), by=list(genotype, exp)]
    genotype_mean%<>%arrange(genotype, exp)
    genotype_exp1<-genotype_mean%>%filter(exp == "Exp1")
    genotype_exp2<-genotype_mean%>%filter(exp == "Exp2")
    
    
    corr_df<-data.frame(Trait = trait , Day = j, Corr = cor(genotype_exp1$ymean, genotype_exp2$ymean, method="spearman"), Pval =cor.test(genotype_exp1$ymean, genotype_exp2$ymean, method="spearman")$p.value)
    
    correlations<-rbind(correlations, corr_df)
    
    #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
    #Knowing variability components allows us to calculate Heritability.
    #Random-effects terms are distinguished by vertical bars or pipes (|) 
    model <- lmer(y ~ (1| genotype) + (1|exp) + (1|genotype:exp),df1_sub) #our random effects mixed model
    summary(model)
    
    model_blue <- lmer(y ~  genotype + (1|exp) + (1|genotype:exp),df1_sub)
    blues<-as.data.frame(fixef(model_blue))
    #blues$Genotype<-rownames(blues)
    #blues$Genotype[1]<-"PPPW_003"
    #blues$Genotype<-gsub("genotype", "", blues$Genotype)
    blues$BLUE<-blues$`fixef(model_blue)`[1]+blues$`fixef(model_blue)`
    blues$BLUE[1]<-blues$`fixef(model_blue)`[1]
    blues$`fixef(model_blue)`<-NULL
    colnames(blues) <- "BLUE" #rename the BLUP column by the trait in our loop
    #add columns to existing dataframe   
    blues$Day<-j
    blues$Trait<-trait
    DataOutputBLUEs <- rbind(DataOutputBLUEs,blues)
    
    blup = coef(model)$genotype #coef extracts model coefficients from lmer objects returned by modeling functions
    hist(blup[,1]) #plot it out
    blup$Day<-j
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
    varComp$Day <- j
    #add columns to existing dataframe
    DataVarComp <- rbind(DataVarComp,varComp) 
    
  }}


#Reshape the Blues/Blups to something more readable
DataOutputBLUEs$Genotype<-rep(goi)
DataOutputBLUEOutput <- reshape(DataOutputBLUEs, idvar = c("Genotype", "Day"), timevar = "Trait", direction = "wide")

DataOutputBLUPs$Genotype<-rep(goi)
DataOutputBLUPOutput <- reshape(DataOutputBLUPs, idvar = c("Genotype", "Day"), timevar = "Trait", direction = "wide")

#reshape our variance components dataframe so that we can run the heritability script
DataVarCompOutput <- reshape(DataVarComp, idvar = c("Trait", "Day"), timevar = "grp", direction = "wide")

#the broad sense heritability script
#Taken from Gioia et al. 2017

HeritabilityData <- (DataVarCompOutput[,4]) / ((DataVarCompOutput[,4]) +  (DataVarCompOutput[,3]/nloc) + (((DataVarCompOutput[,6])) / (nrep*nloc)))

#summary statistics
DataColMax <- colMax(df[,colnum])
DataColMin <- colMin(df[,colnum])
DataColMean <- colMean(df[,colnum])
DataColMedian <- colMean(df[,colnum])
DataColStdev <- colStdev(df[,colnum])

#CVg <- ((DataVarCompOutput[,3])) / ((DataColMean)^2)
CVg <- (sqrt(DataVarCompOutput[,4])) / ((DataColMean))
#out <- LSD.test(model,"Entry", p.adj="bonferroni")

#bind the heritability to the variance component data frame

DataVarCompOutput <- cbind(DataVarCompOutput,HeritabilityData,CVg,DataColMin,DataColMax,DataColMean,DataColMedian,DataColStdev)

#output that beast
write.csv(correlations,paste0("Results/Corr_Exp1_Exp2_", trt, ".csv"))

write.csv(DataVarCompOutput,paste0("Results/H2_Exp1_Exp2_", trt, ".csv"))
fwrite(DataOutputBLUEOutput,paste0("Results/BLUEs_Exp1_Exp2_", trt, ".txt"), sep="\t", row.names = F)
fwrite(DataOutputBLUPOutput,paste0("Results/BLUPs_Exp1_Exp2_", trt, ".txt"), sep="\t", row.names = F)


#######################
# Plot corr between exp

corr_drought<-fread("Results/Corr_Exp1_Exp2_drought.csv")[-1,]
corr_control<-fread("Results/Corr_Exp1_Exp2_control.csv")[-1,]
corr_drought$Treatment<-"drought"
corr_control$Treatment<-"control"


corr<-rbind(corr_drought, corr_control)
#corr$Day<-as.numeric(sub("day ", "", corr$Day))

# for (i in unique(corr$Trait)){
#   ggplot(corr%>%filter(Trait == i), aes(x=Day, y=Corr, group= Treatment)) +   geom_line(aes(col = Treatment)) + geom_point() + theme_classic()+ 
#     theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ ylab("Spearman correlation") + ggtitle(i) + geom_vline(xintercept=45, linetype="dashed", color = "red")
# ggsave(paste0("Figs/Corr_", i, ".tiff"), width=5, height=3)
# }

corr$Day2<-corr$Day-21

tiff(paste0("Figs/AdditionalFile7.tiff"), width=17, height=12, compression = "lzw", res=600, units="cm")
ggplot(corr, aes(x=Day2, y=Corr, group= Treatment)) +   geom_line(aes(col = Treatment)) + geom_point(aes(col = Treatment), size=0.5) + theme_classic()+ 
  theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ facet_wrap(~Trait, scales="free_y", ncol=4) +xlab("DAD")+ylab("Spearman correlation") +  geom_vline(xintercept=45-21, linetype="dashed", color = "red") + theme(legend.position ="bottom")
dev.off()



###################################################################
# Plot h2

heri_drought<-fread("Results/H2_Exp1_Exp2_drought.csv")
heri_control<-fread("Results/H2_Exp1_Exp2_control.csv")
heri_drought$Treatment<-"drought"
heri_control$Treatment<-"control"

heri<-rbind(heri_drought, heri_control)
#heri$Day<-as.numeric(sub("day ", "", heri$Day))

# for (i in unique(heri$Trait)){
#   ggplot(heri%>%filter(Trait == i), aes(x=Day, y=HeritabilityData), group= Treatment) +   geom_line(aes(col = Treatment)) + geom_point(aes(col = Treatment)) + theme_classic()+ 
#     theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ ylab("Broad sense h2") + ggtitle(i) + geom_vline(xintercept=45, linetype="dashed", color = "red")
#   ggsave(paste0("Figs/H2_", i, ".tiff"), width=5, height=3)
# }

heri$Day2<-heri$Day-21

tiff(paste0("Figs/AdditionalFile6.tiff"), width=17, height=12, compression = "lzw", res=600, units="cm")
ggplot(heri, aes(x=Day2, y=HeritabilityData, group= Treatment)) +   geom_line(aes(col = Treatment)) + geom_point(aes(col = Treatment), size=0.5) + theme_classic()+ 
  theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ facet_wrap(~Trait, scales="free_y", ncol=4) +ylab("Broad sense heritability (H2)") +  geom_vline(xintercept=45-21, linetype="dashed", color = "red")+ theme(legend.position ="bottom")
dev.off()


