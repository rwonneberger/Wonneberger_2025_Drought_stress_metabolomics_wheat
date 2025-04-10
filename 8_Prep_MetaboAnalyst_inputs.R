# This script reads in the outputs of the TargetSearch script and prepares different subsets for Metaboanalyst analysis. This script is not entirely automated! i.e. when you want to use a different parameter as the DAM contrast, e.g. tolerance level, you'll need to adjust automatically

# Note that the script now filters out the metabolites that have < 2 fold change between blanks and samples. This calculation is done on the dataset which is only IS-normalized (not weight-normalized), so you will need to provide both datasets

#setwd("/data/home/ronja.wonneberger/Projects/CresWheat/Results")

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

# Select which file should be read in: IS-normalized, weight-normalized or both
# Note that at the moment the TargetSearch script doesn't produce an unnormalized file in the correct format
normalize.against.IS = T
normalize.against.weight = T

# Targeted or untargeted search?
targeted = T

####################################  Select which filters should be applied
time = c(3)
phenotype = c("control", "drought")
tolerance = c("tolerant", "non-tolerant")

####################################  Select contrast on which the DAM should be performed.
# Usually this is the treatment, i.e. drought vs control, but it also makes sense to compare tolerance levels or different time points. When comparing timepoints make sure that the dataset only contains two timepoints, otherwise MetaboAnalyst will run the analysis on the first two timepoints in the dataset

contrast = "Phenotype"

####################################  Dataset prep starts here - no more user inputs required

if(normalize.against.IS) {
  if (normalize.against.weight) {
    if (targeted) {
      filename<-"PeakTable_IS_Weight_normalized" # Additional file 5 sheet 2
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
pheno.table<-read_csv("Annotation_nooutliers_blanks.csv")  # Additional file 5 sheet 1
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

# pheno.table[pheno.table$Subject == "PPPW_003" | pheno.table$Subject == "PPPW_007" | pheno.table$Subject == "PPPW_011" | pheno.table$Subject == "PPPW_016" | pheno.table$Subject == "PPPW_017" | pheno.table$Subject == "PPPW_023",]$Phenotype <- "Group1"
# 
# pheno.table[pheno.table$Subject == "PPPW_004" | pheno.table$Subject == "PPPW_012" | pheno.table$Subject == "PPPW_025" | pheno.table$Subject == "PPPW_033" | pheno.table$Subject == "PPPW_034" | pheno.table$Subject == "PPPW_046",]$Phenotype <- "Group2"

# Write both files
write.csv(pheno.table, paste0(filestring, "_meta.csv"),  row.names = F)
write.csv(df1, paste0(filestring, ".csv"),  row.names = F)


# Prepare a file for one factor analysis 
df2<-merge(pheno.table, df1)
write.csv(df2, paste0(filestring, "_onefactor.csv"),  row.names = F)





