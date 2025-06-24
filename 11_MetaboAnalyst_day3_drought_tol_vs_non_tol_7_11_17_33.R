# This script reads in the outputs of the TargetSearch script and prepares different subsets for Metaboanalyst analysis. this script is not entirely automated! i.e. when you want to use a different parameter as the DAM contrast, e.g. tolerance level, you'll need to adjust automatically

# Not that the script now filters out the metabolites that have < 2 fold change between blanks and samples. This calculation is done on the dataset which is only IS-normalized (not weight-normalized), so you will need to provide both datasets


source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

# Select which file should be read in: IS-normalized, weight-normalized or both
# Note that at the moment the TargetSearch script doesn't produce an unnormalized file in the correct format
normalize.against.IS = T
normalize.against.weight = T

# Targeted or untargeted search?
targeted = T

####################################  Select which filters should be applied
time = c(3)
phenotype = c("control")
tolerance = c("tolerant", "non-tolerant")

####################################  Select contrast on which the DAM should be performed.
# Usually this is the treatment, i.e. drought vs control, but it also makes sense to compare tolerance levels or different time points. When comparing timepoints make sure that the dataset only contains two timepoints, otherwise MetaboAnalyst will run the analysis on the first two timepoints in the dataset

contrast = "Subject"

# Run this for all pairwise comparisons of interest
geno1<-"PPPW_033" # good geno
geno2<-"PPPW_011" # bad geno

####################################  Dataset prep starts here - no more user inputes required

if(normalize.against.IS) {
  if (normalize.against.weight) {
    if (targeted) {
      filename<-"PeakTable_IS_Weight_normalized"
      filterfile<-"PeakTable_IS_normalized"
    } 
    else {
      filename<-"Untargeted_PeakTable_IS_Weight_normalized"
      filterfile<-"Untargeted_PeakTable_IS_normalized"}
  }  else {
    filename<-"PeakTable_IS_normalized"
    filterfile<-"PeakTable_IS_normalized"}
} else {
  if (normalize.against.weight) {
    filename<-"PeakTable_Weight_normalized"
    filterfile<-"PeakTable_Weight_normalized"
    warning("You're going to use a weight-normalized file for filtering. Make sure that this is ok")
  } 
  else {stop("Unnormalized file does not work with this script")}  
}

# Read in peak intensity file
df<-read_csv(paste0(filename, ".csv"))
names(df)[1]<-"Metabolite"

# Read in Annotation file
pheno.table<-read_csv("Annotation_nooutliers_blanks.csv")
#pheno.table<-read_csv("../Data/Nordrome_Metabolomics/Annotation_nooutliers.csv")

# Make a lookup table to change the sample names in the peak intensity file
lookup<-pheno.table[, c(1,9)]
names(lookup)[1]<-"Sample"
lookup$Sample<-gsub(".cdf", "", lookup$Sample)


# Filter metabolites that have a > 2 fold change in the blanks, and then remove the blanks
filterfile1<-read_csv(paste0(filterfile, ".csv"))
names(filterfile1)[1]<-"Metabolite"
blanks<-lookup[grepl( "Blank", lookup$Unique_name, fixed = TRUE),]$Sample
samp<-lookup[!grepl( "Blank", lookup$Unique_name, fixed = TRUE),]$Sample

df_blanks<-filterfile1%>%dplyr::select(Metabolite, all_of(blanks))
df_samp<-filterfile1%>%dplyr::select(Metabolite, all_of(samp))

df_blanks$Mean<-rowMeans(df_blanks[,-1], na.rm=T)
df_samp$Mean<-rowMeans(df_samp[,-1], na.rm=T)

df_blanks$Miss<-rowSums(is.na(df_blanks[,-c(1, 26)]))
df_samp$Miss<-rowSums(is.na(df_samp[,-c(1, 514)]))

df_analysis<-data.frame(df_samp$Metabolite, df_blanks$Mean, df_samp$Mean, df_blanks$Miss, df_samp$Miss)
df_analysis$Fold<-df_analysis$df_samp.Mean/df_analysis$df_blanks.Mean

metab_keep<-df_analysis%>%filter(Fold > 2 | df_blanks.Miss == 24)%>%dplyr::select(df_samp.Metabolite)

df%<>%filter(Metabolite%in%metab_keep$df_samp.Metabolite)


# Prepare a meta file with samples in rows. First col = sample name

#pheno.table$rawname<-gsub(".cdf", "", pheno.table$rawname)
pheno.table$name_pot<-pheno.table$Rep<-pheno.table$`Amount [mg(FW)]  (initial weight)`<-pheno.table$type<-pheno.table$rawname<-NULL

pheno.table%<>%relocate(Unique_name, condition, day, name_geno)
pheno.table<-dplyr::filter(pheno.table, !grepl("Blank", Unique_name))

names(pheno.table)<-c("Sample", "Phenotype", "Time",  "Subject", "Tolerance", "Batch")

# change the filename sothat we know that the output we generate is for use with MetaboAnalyst
filename<-gsub("PeakTable_", "MA_", filename)
filename<-gsub("alized", "", filename)

# Prepare a file with samples in rows and metabolites in columns. First col = sample name
df1<-transpose_df(df)
df1$Sample<-rownames(df1)
df1<-left_join(df1, lookup, by="Sample")
df1%<>%relocate(Unique_name)
df1$Sample<-NULL
names(df1)[1]<-"Sample"
df1<-dplyr::filter(df1, !grepl("Blank", Sample))

# check that sample order is the same in both files, seems to be important?
df1%<>%arrange(Sample)
pheno.table%<>%arrange(Sample)

df1$Sample == pheno.table$Sample

## Subset the data
pheno.table%<>%filter(Subject%in%c(geno1, geno2), Time %in% time, Phenotype %in% phenotype, Tolerance %in% tolerance)
df1%<>%filter(Sample %in% pheno.table$Sample)





# Data wrangling to make sure the selected contrast is in the "Phenotype" column
if (contrast != "Phenotype") {
  pheno.table$Phenotype<-NULL
  names(pheno.table)[which(names(pheno.table) == contrast)]<-"Phenotype"
  pheno.table%<>%relocate(Sample, Phenotype)
}

# Throw an error if the chosen contrast does not have 2 levels (two are required for FC analysis)
if (length(table(pheno.table$Phenotype)) != 2){
  stop(paste0("The chosen contrast has ", length(table(pheno.table$Phenotype)), " levels. Two levels are required for DAM analysis."))
}


if (length(phenotype) == 2){phenotype = "all"}
if (length(tolerance) == 2){tolerance = "all"}
if (length(time) > 1){time<-paste(time, collapse="")}

if (targeted) {target = "T"} else {target = "U"}

filename<-gsub("Untargeted_", "", filename)

filestring<-paste0(target, "_", filename, "_time_", time, "_trt_", phenotype, "_tol_", tolerance, "_geno_", geno1, "_", geno2)


# Write both files
write.csv(pheno.table, paste0(filestring, "_meta.csv"),  row.names = F)
write.csv(df1, paste0(filestring, ".csv"),  row.names = F)

# Prepare a file for one factor analysis 
df2<-merge(pheno.table, df1)
write.csv(df2, paste0(filestring, "_onefactor.csv"),  row.names = F)



















































# analysis of metabolomics data using MetaboAnalystR


source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

library(MetaboAnalystR)
library(ggplot2)
library(reshape2)
library(ggforce)
library(magrittr)
library(dplyr)
library(data.table)
library(factoextra)
library(ggtree)
library(ape)
library(ggnewscale)
library(corrplot)
library(ggrepel)

###################### Select which dataset to analyse ###################### 
normalize.against.IS = T
normalize.against.weight = T

######### Select which filters should be applied
time = c(3)
phenotype = c("control")
tolerance = c("tolerant", "non-tolerant")

# Note that at the moment, the main level has to be drought/control because MA always selects the column called Phenotype for this. If we want to compare e.g. tolerance levels we will have to make this script more flexible first

# Targeted or untargeted search?
targeted = T

###################### Parameters for filtering the data ###################### 
perc_miss <- 0.2                # Remove features with more than x% missing values
imputation_method <- "min"      # Options: min, mean, median KNN, PPCA, BPCA, svdImpute
filtering_method <- "rsd"       # Options: rsd, nrsd, mean, sd, mad, iqr
perc_removed <- 25              # Percent data to remove by selected filtering method
rowNorm <- "MedianNorm"         # Row-wise normalization: Options: QuantileNorm, CompNorm, SumNorm, MedianNorm, SpecNorm
transNorm <- "LogNorm"          # Data transformation: Options: LogNorm, Crnorm
scaleNorm <- "ParetoNorm"       # scaling options: Options: MeanCenter, AutoNorm, ParetoNorm, RangeNorm

# Create the right filename
if (length(phenotype) == 2){phenotype = "all"}
if (length(tolerance) == 2){tolerance = "all"}
if (length(time) > 1){time<-paste(time, collapse="")}
if (targeted) {target="T"} else {target="U"}

if(normalize.against.IS) {
  if (normalize.against.weight) {filename<-paste0(target, "_MA_IS_Weight_norm")
  }  else {filename<-paste0(target, "_MA_IS_norm")}
} else {
  if (normalize.against.weight) {filename<-paste0(target, "_MA_Weight_norm")} 
  else {stop("Unnormalized file not found")}  
}

filestring<-paste0(filename, "_time_", time, "_trt_", phenotype, "_tol_", tolerance, "_geno_", geno1, "_", geno2)

############################################################################################################
############################## Analysis starts here ########################################################
############################################################################################################

###################### Read the normalized data (the output from the the TargetSearch script) ###################### 
mSet<-NULL # Set it to NULL at the beginning of the script, otherwise you'll have to restart R everytime you want to rerun this script
mSet<-InitDataObjects("pktable", "mf", FALSE)
mSet<-SetDesignType(mSet, "multi")
mSet<-Read.TextDataTs(mSet, paste0(filestring, ".csv"), "rowmf")
mSet<-ReadMetaData(mSet, paste0(filestring, "_meta.csv"))
mSet$dataSet$cls.type<-"disc" # This was not necessary before but suddenly it is...
mSet<-SanityCheckData(mSet)

###################### visualize proportion of missing data ###################### 
orig.data <- qs::qread("data_orig.qs")
miss<-sapply(orig.data, function(x) sum(is.na(x)))
hist(miss, breaks = 60)

###################### Remove missing data and impute ###################### 
# Remove features with more than x% missing values
mSet<-RemoveMissingPercent(mSet, percent=perc_miss)

# Impute missing data (Replace by LoDs (1/5 of the minimum positive value of each variable))
mSet<-ImputeMissingVar(mSet, method=imputation_method)
mSet<-SanityCheckMeta(mSet, 1)

###################### Filtering ######################
# Relative standard deviation (=SD/mean) filtering - remove x%
mSet<-FilterVariable(mSet, filtering_method, perc_removed, "F", 25, F)


###################### Normalization, transformation, scaling ######################
# Normalization by median, Log10 transformation, Pareto scaling
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, rowNorm, transNorm, scaleNorm, ratio=FALSE, ratioNum=20)


###################### Extract the normalized data which we're going to use for most purposes from now on ######################
norm_data<-as.data.frame(mSet$dataSet$norm)
dim(norm_data)

# write the normalized data for potential analysis later
norm_out<-norm_data
norm_out$Sample<-rownames(norm_out)
norm_out%<>%relocate(Sample)
fwrite(norm_out, paste0("Normalized_", filestring, ".txt"), sep="\t")


########################################################################################
###################### Exploratory analysis done #######################################
########################################################################################


###################### Fold-changes, t-tests and Volcano plots ######################

#For the complete list of all FC

data <- qs::qread("row_norm.qs");
fwrite(as.data.frame(data), paste0(filestring, "_rownorm.txt"),row.names=T, sep="\t")
m1 <- colMeans(data[which(mSet$dataSet$cls==levels(mSet$dataSet$cls)[1]), ]);
m2 <- colMeans(data[which(mSet$dataSet$cls==levels(mSet$dataSet$cls)[2]), ]);

ratio <- m2/m1 #a value log2(value) > 0 means that it's up in the good genotype

fc.log <- signif(log2(ratio), 5)
fc.log<-as.data.frame(fc.log)
fc.log$Metabolite<-rownames(fc.log)
fwrite(fc.log, paste0( filestring, "_fcall.txt"), row.names=T, sep="\t")

# Normal analysis 
mSet<-FC.Anal(mSet, 2, 1, FALSE) # Note it must be set to 1!
fwrite(as.data.frame(mSet$analSet$fc$sig.mat), paste0( filestring, "_fc.txt"), row.names=T, sep="\t")



# All t-test
inx1 <- which(mSet$dataSet$cls==levels(mSet$dataSet$cls)[1]);
inx2 <- which(mSet$dataSet$cls==levels(mSet$dataSet$cls)[2]);

data <- as.matrix(mSet$dataSet$norm);


univ.test <- function(x){t.test(x[inx1], x[inx2], paired = F, var.equal = T)};

my.fun <- function(x) {
  tmp <- try(univ.test(x));
  if(class(tmp) == "try-error") {
    return(c(NA, NA));
  }else{
    return(c(tmp$statistic, tmp$p.value));
  }
}
ttest_all <- apply(data, 2, my.fun);
ttest_all<-as.data.frame(t(ttest_all))
ttest_all$fdr.p <- p.adjust(ttest_all$V2, "fdr")
names(ttest_all)<-c("stat", "pval", "fdr")
ttest_all$Metabolites<-rownames(ttest_all)

mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "fdr")
ttest<-as.data.frame(mSet$analSet$tt$sig.mat)
ttest$Metabolites<-rownames(ttest)
norm_data1<-norm_data
norm_data1$Sample<-rownames(norm_data1)
meta<-as.data.frame(mSet$dataSet$meta.info)
meta$Sample<-rownames(meta)
norm_data1<-merge(meta, norm_data1)
norm_data1%<>%dplyr::select(-Sample,  -Tolerance, -Batch, -Time)
norm_data1

norm_data1_mean<-norm_data1%>%
  group_by(Phenotype) %>%
  summarise_all("mean")
norm_data1_mean<-transpose_df_Col1(norm_data1_mean, "Metabolites")
norm_data1_mean$direction<-""
norm_data1_mean<-norm_data1_mean[complete.cases(norm_data1_mean),]
names(norm_data1_mean)[2]<-"non-tolerant"
names(norm_data1_mean)[3]<-"tolerant"
norm_data1_mean[norm_data1_mean$tolerant < norm_data1_mean$`non-tolerant`,]$direction<-"down"
norm_data1_mean[norm_data1_mean$tolerant > norm_data1_mean$`non-tolerant`,]$direction<-"up"

norm_data1_mean_merged<-merge(norm_data1_mean, ttest_all, by="Metabolites")

fwrite(norm_data1_mean_merged, paste0( filestring, "_ttall.txt"), row.names=T, sep="\t")

names(fc.log)[2]<-"Metabolites"
out_all<-merge(fc.log, norm_data1_mean_merged, by="Metabolites")
names(out_all)<-c("Metabolite", "log2(FC)", "non-tolerant", "tolerant", "Direction", "Tstat", "t_pvalue", "T_FDR")

fwrite(out_all, paste0(filestring, "_ttest.txt"), row.names=T, sep="\t")

ttest1<-left_join(ttest, norm_data1_mean)


fwrite(ttest1, paste0( filestring, "_tt.txt"), row.names=T, sep="\t")
mSet<-PlotTT(mSet, paste0("../Figs/", filestring, "_Ttest"), "png", 72, width=NA)
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 1, F, 0.1, TRUE, "fdr") # note must be 1
mSet<-PlotVolcano(mSet, paste0("../Figs/", filestring, "_Volcano"),1, 0, "png", 72, width=NA)


volcano_df<-data.frame("FC" = mSet$analSet$fc$fc.log, "TT" = mSet$analSet$tt$p.log)
volcano_df$Metabolite<-rownames(volcano_df)
volcano_df$Status <- "Not sig."
volcano_df$Status[volcano_df$FC > 1 & volcano_df$TT > -log10(0.05)] <- "Up in tolerant"
volcano_df$Status[volcano_df$FC < -1 & volcano_df$TT > -log10(0.05)] <- "Down in tolerant"
volcano_df$damlabel <- NA
volcano_df$damlabel[volcano_df$Status != "Not sig."] <- volcano_df$Metabolite[volcano_df$Status != "Not sig."]

mycolors <- c("#1984c5", "#c23728",  "#e2e2e2")

names(mycolors) <- c("Down in tolerant", "Up in tolerant", "Not sig.")


p1<-ggplot(volcano_df, aes(FC, TT, color = Status, label=damlabel)) + geom_point() + geom_hline(yintercept=-log10(0.05), linetype="dashed") + geom_vline(xintercept=c(-1, 1), linetype="dashed")+ scale_colour_manual(values = mycolors)+ geom_text_repel(size = 4, show_guide=FALSE, max.overlaps = 15) + theme_minimal() + xlab("log2(fold change)") + ylab("-log10(p), FDR-adjusted") + theme(legend.position = "bottom")
p1


fwrite(as.data.frame(mSet$analSet$volcano$sig.mat), paste0(filestring, "_sig_DAM.txt"), sep="\t", row.names=T)




###################################
# Make upset plot to see which mets are in common between genotypes

library(UpSetR)
library(ComplexHeatmap)


day=3
trt="drought"

# Read in the ttest results
g1<-fread(paste0("T_MA_IS_Weight_norm_time_3_trt_", trt, "_tol_all_geno_PPPW_017_PPPW_007_ttall.txt"))
g2<-fread(paste0("T_MA_IS_Weight_norm_time_3_trt_", trt, "_tol_all_geno_PPPW_033_PPPW_007_ttall.txt"))
g3<-fread(paste0("T_MA_IS_Weight_norm_time_3_trt_", trt, "_tol_all_geno_PPPW_017_PPPW_011_ttall.txt"))
g4<-fread(paste0("T_MA_IS_Weight_norm_time_3_trt_", trt, "_tol_all_geno_PPPW_033_PPPW_011_ttall.txt"))


g1$Genotype<-"PPPW_017_PPPW_007"
g2$Genotype<-"PPPW_033_PPPW_007"
g3$Genotype<-"PPPW_017_PPPW_011"
g4$Genotype<-"PPPW_033_PPPW_011"


g1up<-g1%>%filter(direction == "up")
g2up<-g2%>%filter(direction == "up")
g3up<-g3%>%filter(direction == "up")
g4up<-g4%>%filter(direction == "up")


g1down<-g1%>%filter(direction == "down")
g2down<-g2%>%filter(direction == "down")
g3down<-g3%>%filter(direction == "down")
g4down<-g4%>%filter(direction == "down")


tt_all<-rbind(g1, g2, g3, g4)


tt_all%>%filter(fdr < 0.05)


## Read in the VIP values from OPLS-DA (run in Metaboanalyst V5 GUI)

opls_017_007<-fread(paste0("Metab_list_T_MA_IS_Weight_norm_time_3_trt_", trt, "_tol_all_PPPW_017_PPPW_007_VIP.txt")) 
opls_033_007<-fread(paste0("Metab_list_T_MA_IS_Weight_norm_time_3_trt_", trt, "_tol_all_PPPW_033_PPPW_007_VIP.txt"))
opls_017_011<-fread(paste0("Metab_list_T_MA_IS_Weight_norm_time_3_trt_", trt, "_tol_all_PPPW_017_PPPW_011_VIP.txt"))
opls_033_011<-fread(paste0("Metab_list_T_MA_IS_Weight_norm_time_3_trt_", trt, "_tol_all_PPPW_033_PPPW_011_VIP.txt"))


opls_all<-rbind(opls_017_007, opls_033_007, opls_017_011, opls_033_011)
names(opls_all)<-c("Metabolites", "VIP", "OrthoVIP", "Genotype")


all_df<-merge(tt_all, opls_all, by=c("Metabolites", "Genotype"), all=T)
all_df%<>%filter(VIP >=1)
all_df%<>%filter(fdr < 0.05)
all_df%<>%filter(pval < 0.05)

all_df_up<-all_df%>%filter(direction == "up")
all_df_down<-all_df%>%filter(direction == "down")


df_up<-list(
  PPPW_017_PPPW_007 = all_df_up[all_df_up$Genotype == "PPPW_017_PPPW_007",]$Metabolites,
  PPPW_033_PPPW_007 = all_df_up[all_df_up$Genotype == "PPPW_033_PPPW_007",]$Metabolites, 
  
  PPPW_017_PPPW_011 = all_df_up[all_df_up$Genotype == "PPPW_017_PPPW_011",]$Metabolites, 
  
  PPPW_033_PPPW_011 = all_df_up[all_df_up$Genotype == "PPPW_033_PPPW_011",]$Metabolites
  
)


m1_up = make_comb_mat(df_up, min_set_size = 2)

p3<-UpSet(m1_up , comb_order = order(comb_size(m1_up), decreasing=T))







df_down<-list(
  PPPW_017_PPPW_007 = all_df_down[all_df_down$Genotype == "PPPW_017_PPPW_007",]$Metabolites,
  PPPW_033_PPPW_007 = all_df_down[all_df_down$Genotype == "PPPW_033_PPPW_007",]$Metabolites, 
  
  PPPW_017_PPPW_011 = all_df_down[all_df_down$Genotype == "PPPW_017_PPPW_011",]$Metabolites, 
  
  PPPW_033_PPPW_011 = all_df_down[all_df_down$Genotype == "PPPW_033_PPPW_011",]$Metabolites
  
)


m1_down = make_comb_mat(df_down, min_set_size = 2)

p4<-UpSet(m1_down , comb_order = order(comb_size(m1_down), decreasing=T))




# Run the other script ManusScript_MetaboAnalyst_day3_c_d_tol_non_tol_7_11_17_33.R first and then combine all 4 Figs into one


gb_p1 = grid.grabExpr(draw(p1))
gb_p2 = grid.grabExpr(draw(p2))
gb_p3 = grid.grabExpr(draw(p3))
gb_p4 = grid.grabExpr(draw(p4))

tiff("Figs/AdditionalFile19.tiff", height=18, width=28, units="cm", res=300)
ggarrange(gb_p1, gb_p2, gb_p3, gb_p4, ncol=2, nrow=2, labels=c("a", "b", "c", "d"))
dev.off()
