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
time = c(3, 4)
phenotype = c( "drought")
tolerance = c("tolerant", "non-tolerant")

####################################  Select contrast on which the DAM should be performed.
# Usually this is the treatment, i.e. drought vs control, but it also makes sense to compare tolerance levels or different time points. When comparing timepoints make sure that the dataset only contains two timepoints, otherwise MetaboAnalyst will run the analysis on the first two timepoints in the dataset

contrast = "Time"

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
pheno.table%<>%filter(Time %in% time, Phenotype %in% phenotype, Tolerance %in% tolerance)
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

filestring<-paste0(target, "_", filename, "_time_", time, "_trt_", phenotype, "_tol_", tolerance)



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
time = c(3,4)
phenotype = c( "drought")
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

filestring<-paste0(filename, "_time_", time, "_trt_", phenotype, "_tol_", tolerance)

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

# Visualize the changes to the data

mSet<-PlotSampleNormSummary(mSet, paste0("../Figs/", filestring, "_snorm_"), "png", 72, width=NA)


###################### Extract the normalized data which we're going to use for most purposes from now on ######################
norm_data<-as.data.frame(mSet$dataSet$norm)
dim(norm_data)

# write the normalized data for potential analysis later
norm_out<-norm_data
norm_out$Sample<-rownames(norm_out)
norm_out%<>%relocate(Sample)
fwrite(norm_out, paste0("Normalized_", filestring, ".txt"), sep="\t")


# To plot all data we need a custom plot
df<-norm_data
df$Sample<-rownames(df)
df_long<-melt(df)
df_long$new_name<-paste0(substr(df_long$Sample, 12, 12), "_", substr(df_long$Sample, 10, 10), "_", substr(df_long$Sample, 1, 8), "_", substr(df_long$Sample, 14, 14))

########################################################################################
###################### Exploratory analysis done #######################################
########################################################################################


###################### Fold-changes, t-tests and Volcano plots ######################

#For the complete list of all FC

data <- qs::qread("row_norm.qs");
fwrite(as.data.frame(data), paste0(filestring, "_rownorm.txt"),row.names=T, sep="\t")
m1 <- colMeans(data[which(mSet$dataSet$cls==levels(mSet$dataSet$cls)[1]), ]);
m2 <- colMeans(data[which(mSet$dataSet$cls==levels(mSet$dataSet$cls)[2]), ]);

ratio <- m2/m1

fc.log <- signif(log2(ratio), 5)
fc.log<-as.data.frame(fc.log)
fc.log$Metabolite<-rownames(fc.log)
fwrite(fc.log, paste0( filestring, "_fcall.txt"), row.names=T, sep="\t")

# Normal analysis 
mSet<-FC.Anal(mSet, 2, 1, FALSE)
fwrite(as.data.frame(mSet$analSet$fc$sig.mat), paste0( filestring, "_fc.txt"), row.names=T, sep="\t")
mSet<-PlotFC(mSet, paste0("../Figs/", filestring, "_fc"), "png", 72, width=NA)



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
norm_data1%<>%dplyr::select(-Sample,  -Subject, -Batch, -Tolerance)
norm_data1

norm_data1_mean<-norm_data1%>%
  group_by(Phenotype) %>%
  summarise_all("mean")
norm_data1_mean<-transpose_df_Col1(norm_data1_mean, "Metabolites")
norm_data1_mean$direction<-""
norm_data1_mean<-norm_data1_mean[complete.cases(norm_data1_mean),]
norm_data1_mean[norm_data1_mean$`3` < norm_data1_mean$`4`,]$direction<-"up"
norm_data1_mean[norm_data1_mean$`3` > norm_data1_mean$`4`,]$direction<-"down"

norm_data1_mean_merged<-merge(norm_data1_mean, ttest_all, by="Metabolites")

fwrite(norm_data1_mean_merged, paste0( filestring, "_ttall.txt"), row.names=T, sep="\t")

names(fc.log)[2]<-"Metabolites"
out_all<-merge(fc.log, norm_data1_mean_merged, by="Metabolites")
names(out_all)<-c("Metabolite", "log2(FC)", "mean_Day3", "mean_Day4", "Direction", "Tstat", "t_pvalue", "T_FDR")

fwrite(out_all, paste0(filestring, "_all_FC_ttest.txt"), row.names=T, sep="\t")

ttest1<-left_join(ttest, norm_data1_mean)


fwrite(ttest1, paste0( filestring, "_tt.txt"), row.names=T, sep="\t")
mSet<-PlotTT(mSet, paste0("../Figs/", filestring, "_Ttest"), "png", 72, width=NA)
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 1, F, 0.1, TRUE, "fdr")
#mSet<-PlotVolcano(mSet, paste0("../Figs/", filestring, "_Volcano"),1, 0, "png", 72, width=NA)


volcano_df<-data.frame("FC" = mSet$analSet$fc$fc.log, "TT" = mSet$analSet$tt$p.log)
volcano_df$Metabolite<-rownames(volcano_df)
volcano_df$Status <- "Not sig."
volcano_df$Status[volcano_df$FC > 1 & volcano_df$TT > -log10(0.05)] <- "Up after rewatering"
volcano_df$Status[volcano_df$FC < -1 & volcano_df$TT > -log10(0.05)] <- "Down after rewatering"
volcano_df$damlabel <- NA
volcano_df$damlabel[volcano_df$Status != "Not sig."] <- volcano_df$Metabolite[volcano_df$Status != "Not sig."]

#Write df to file so it can be used for plotting in another script
fwrite(volcano_df, paste0(filestring, "_forvolcanoplot.txt"), sep="\t", row.names=T)


mycolors <- c("#1984c5", "#c23728",  "#e2e2e2")

names(mycolors) <- c("Down after rewatering", "Up after rewatering", "Not sig.")

tiff(paste0("Figs/", filestring,"_volcano.tiff"), height=15, width=20, units="cm", res=300)
p1<-ggplot(volcano_df, aes(FC, TT, color = Status, label=damlabel)) + geom_point() + geom_hline(yintercept=-log10(0.05), linetype="dashed")+ theme_minimal() + geom_vline(xintercept=c(-1, 1), linetype="dashed")+ scale_colour_manual(values = mycolors)+ geom_text_repel(size = 2, show_guide=FALSE, max.overlaps = 15)  + xlab("log2(fold change)") + ylab("-log10(p), FDR-adjusted") + theme(legend.position = "right", text = element_text(size = 8))
p1
dev.off()

fwrite(as.data.frame(mSet$analSet$volcano$sig.mat), paste0(filestring, "_sig_DAM.txt"), sep="\t", row.names=T)



###################### Pathway analysis ######################
# Select the metabolites that are kept in the current analysis after filtering
# metabolite_reference<-names(mSet$dataSet$norm)
# 
met_lookup<-fread("../Data/KEGG_lookup.txt")

# Get the KEGG terms of the DAMs
DAM<-as.data.frame(mSet$analSet$volcano$sig.mat)

# Add filters if you wish
#DAM%<>%filter(log2(FC) >1)

met_list<-met_lookup%>%filter(Name%in%rownames(DAM))

# Check if the number of DAMs is equal to the number of entries in met_list. If not, there is a discrepancy between metabolite names, most likely because in TargetSearch we renamed duplicates to "XX (copy y)"
length(rownames(DAM)) == dim(met_list)[1]
setdiff(rownames(DAM), met_list$Name)

fwrite(met_list, paste0("Metab_list_", filestring, "_DAMs_volcano.txt"), sep="\t")


# Do the same but for the ttest DAMs, ie the ones selected only based on significant differences, not on fc
# Get the KEGG terms of the DAMs
DAM_ttest<-as.data.frame(mSet$analSet$tt$sig.mat)

# Add filters if you wish
#DAM%<>%filter(log2(FC) >1)

met_list<-met_lookup%>%filter(Name%in%rownames(DAM_ttest))

# Check if the number of DAMs is equal to the number of entries in met_list. If not, there is a discrepancy between metabolite names, most likely because in TargetSearch we renamed duplicates to "XX (copy y)"
length(rownames(DAM_ttest)) == dim(met_list)[1]
setdiff(rownames(DAM_ttest), met_list$Name)

fwrite(met_list, paste0("Metab_list_", filestring, "_DAMs_ttest.txt"), sep="\t")

out_all_up<-out_all%>%filter(Direction == "up", T_FDR < 0.05)%>%dplyr::select(Metabolite)
met_list_up<-met_list%>%filter(Name%in%out_all_up$Metabolite)

fwrite(met_list_up, paste0("Metab_list_", filestring, "_up_DAMs_ttest.txt"), sep="\t")

out_all_down<-out_all%>%filter(Direction == "down", T_FDR < 0.05)%>%dplyr::select(Metabolite)
met_list_down<-met_list%>%filter(Name%in%out_all_down$Metabolite)

fwrite(met_list_down, paste0("Metab_list_", filestring, "_down_DAMs_ttest.txt"), sep="\t")




################ Plot volcano and enrichment

enrich<-fread("/Enrichment_recovery_FELLA_output.txt")

p2<-ggplot(enrich, aes(x=Analysis, y=KEGG.name, size=-log10(`p-value`), color=-log10(`p-value`))) + geom_point()+theme_bw()+
  guides(colour = guide_legend())+xlab("Analysis")+x_axis_90_degrees + theme(legend.position="right")+ theme(axis.title=element_blank()) + guides(color=guide_legend(title="-log10(p)"), size=guide_legend(title="-log10(p)"))


tiff("Figs/Fig3_Volcano_PathwayEnrichment_recovery.tiff", width=17, height=22, units="cm", compression = "lzw",res=600)
ggarrange(p1, p2, ncol=1, labels=c("a", "b"), heights=c(0.6, 0.4))
dev.off()


normalized<-fread(paste0("/Results/", filestring, "_rownorm.txt"), sep="\t")

volcano_df_sub<-volcano_df%>%filter((FC > 1 | FC < -1), TT > -log10(0.05))

normalized_sub<-normalized%>%dplyr::select(V1, all_of(volcano_df_sub$Metabolite))

normalized_sub$Genotype<-substr(normalized_sub$V1, start = 1, stop = 8)
normalized_sub$Day<-substr(normalized_sub$V1, start = 12, stop = 12)

normalized_sub <- normalized_sub %>% 
  mutate(Day2 = if_else(Day == 3, "22 DAD", "28 DAD"))
normalized_sub$V1<-normalized_sub$Day<-NULL

normalized_sub %<>%
  group_by(Genotype, Day2)%>%
  summarise_if(is.numeric, mean, na.rm=T)

normalized_sub_melt<-melt(normalized_sub)

p3<-ggplot(normalized_sub_melt,aes( x=Genotype, y=variable, fill=log2(value)))+geom_tile()+ theme_classic()+ 
  theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA), legend.position="bottom")+facet_grid(~Day2, scales="free_x")+
  scale_fill_gradientn(
    colours = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(255),
    values = c(1.0, (0 - min(log2(normalized_sub_melt$value))) / (max(log2(normalized_sub_melt$value)) - min(log2(normalized_sub_melt$value))), 0))+
  theme(axis.title=element_blank())+x_axis_90_degrees

tiff("Figs/suppl_recovery_rownorm_heat.tiff", width=17, height=22, units="cm", compression = "lzw",res=600)
p3

dev.off()
