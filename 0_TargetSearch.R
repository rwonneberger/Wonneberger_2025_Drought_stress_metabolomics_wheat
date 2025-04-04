# GC-MS analysis using TargetSearch (https://bioconductor.org/packages/release/bioc/html/TargetSearch.html )
# Version 2023-04-13
# written by Dr. Dariusz Kruszka and Mr. Dennis Psaroudakis, IPK Gatersleben. Modified by Dr. Ronja Wonneberger, SLU.

setwd("/data/home/ronja.wonneberger/Projects/CresWheat/")

############ 1. Load packages #################################################

library(tidyverse)
library(TargetSearch)
library(readr)
library(data.table)
library(magrittr)

############ 2. Import files ##################################################

# prepare path for calculation
cdf.path <- file.path("./Data/Nordrome_Metabolomics")
sample.file <- file.path(cdf.path,"Sample_list_outlier127.txt")
samples <- ImportSamples(sample.file, CDFpath = cdf.path, RIpath = cdf.path)
rim.file  <- file.path(cdf.path, "RIMS.txt")

# read phenodata (sample metadata and weights)
pheno.table<-read_csv(paste0(cdf.path, "/Annotation_outlier127.csv"))
pheno.table<-column_to_rownames(pheno.table,"rawname")

############ 3. Import and check FAMEs (RIMs) #################################

rimLimits <- ImportFameSettings(rim.file, mass = 87)
pdf("Results/RIMS_outliers.pdf",width=8,height=8) #for all samples
for(i in 1:length(samples)) {
  print(paste("Plotting RIMs for sample", i, "of", length(samples)))
  checkRimLim(samples[i], rimLimits)
}
dev.off()


# # Plot a few randomly selected samples:
# rimLimits <- ImportFameSettings(rim.file, mass = 87)
# pdf("Results/RIMS_subset.pdf",width=8,height=8) #for all samples
# select_plots<-sample(1:length(samples), 20)
# for(i in select_plots) {
#   print(paste("Plotting RIMs for sample", i, "of", length(samples)))
#   checkRimLim(samples[i], rimLimits)
# }
# dev.off()

############ 4. File conversion, baseline correction, Peak finding, RIcorrection ####

# This next step converts and compresses the CDF files into a format that speeds up all further analysis.
# It is a little time-consuming but saves time later on. If you're impatient, you can also skip it.
# Important! This assumes that the CDF files were already baseline corrected. If they were not,
# please skip the next line and instead uncomment and execute the one below to perform baseline correction here.
#samples <- ncdf4Convert(samples, path=".") # without baseline correction
samples <- ncdf4Convert(samples, path="Baseline_corrected", baseline=TRUE, bsline_method="quantiles") # further potentially relevant parameters: smooth=0, qntl=0.50, width=30,steps=10
# Mine are baseline corrected but it doesn't not hurt to do it again. If there are already files, remove them/make a new folder

# RI correction and peak picking
# choose good threshold. I think 10 is too small, IntThreshold = 50 
RImatrix <- RIcorrect(samples, rimLimits, IntThreshold = 50, pp.method = "ppc", Window = 15, showProgressBar=T, massRange=c(85,500))

# Outliers visualization
outliers <- FAMEoutliers(samples, RImatrix, threshold = 5, pdffile ="Results/outliers1_outliers.pdf")

############ 5. Targeted Analysis #############################################

# set library path
lib.file <- file.path(cdf.path, "1GMDB_FAMELib_TS_20110228_IS.txt")
#import library
lib <- ImportLibrary(lib.file, RI_dev = c(2000,1000,500),TopMasses = 10, SelMasses=5, ExcludeMasses = c(73, 74, 75, 147, 148, 149))

##### since the sorbitol peak was slightly outside the RI window in some samples, we manually adjust the RI window

RIdev(lib)[352,1]<-3000
RIdev(lib)[352,2]<-2000
RIdev(lib)[352,3]<-1000

# find median RI per metabolite and update the library with it
lib <- medianRILib(samples, lib, makeReport=T, showProgressBar=T,pdfFile = "Results/medianLibRep_split_outliers.pdf")
# normalize RI
cor_RI <- sampleRI(samples, lib, r_thres = 0.95, method = "dayNorm", minPairObs = 6, showProgressBar = T, 
                   makeReport = T, pdfFile = "Results/sampleRI_Rep_outliers.pdf")


# find the peaks based on lib and save RI and intensity 
peakData <- peakFind(samples, lib, cor_RI, showProgressBar =T)
met.RI <- retIndex(peakData)
met.Intensity <- Intensity(peakData)

# create profile and write it
MetabProfile <- Profile(samples, lib, peakData, r_thres = 0.95,method = "dayNorm")
finalProfile <- ProfileCleanUp(MetabProfile, timeSplit = 500, r_thres = 0.95)
Write.Results(lib,finalProfile,quantMatrix='maxint', prefix="Results/TargetSearch_outliers")

# spectra visualization
plotAllSpectra(lib, peakData, type = "ht",pdfFile="Results/spectra_outliers.pdf",width = 10, height = 8)

# write data
## correlating metabolites
met.final <- quantMatrix(lib, finalProfile)
## all metabolites in the library
met.all   <- quantMatrix(lib, MetabProfile)

met.final<-t(met.final)
met.all<-t(met.all)
dimnames(met.final)[2]<-gsub(".cdf", "", rownames(pheno.table)[1])
dimnames(met.all)[2]<-gsub(".cdf", "", rownames(pheno.table)[1])

dimnames(met.final)[1]<-attributes(met.final)[6]
dimnames(met.all)[1]<-attributes(met.all)[6]

# Rename metabolites that occur multiple times to make the metabolite names unique
# This makes it easier to work with Metaboanalyst etc. (a '(copy X)' is appended to all
# metabolites that occur multiple times)
rownames(met.all) = ave(rownames(met.all), rownames(met.all), FUN=function(x) if (length(x)>1) paste0(x[1], ' (copy ', seq_along(x), ')') else x[1])
rownames(met.final) = ave(rownames(met.final), rownames(met.final), FUN=function(x) if (length(x)>1) paste0(x[1], ' (copy ', seq_along(x), ')') else x[1])
makeDF <- function(x) cbind(quantMass=attr(x,"quantMass"),isSelMass=attr(x,"isSelMass"), x)
write.table(makeDF(met.final), file='Results/TargetSearch_outliers.final.txt', sep='\t', quote=F)
write.table(makeDF(met.all), file='Results/TargetSearch_outliers.all.txt', sep='\t', quote=F)

############ 6. Normalization ###############################

# matrix transposition 
QMAt<-t(met.final)
#dimnames(QMAt)[2]<-attributes(QMAt)[6] #Had to add this so that we can look for the standards later in the column names of the matrix

# Since we're doing this a lot we might just re-use the code later rather than repeat it
plot_histogram <- function(SDsel, data.stage) {
  #Test for %SD
  SDsel_av <- apply(SDsel, 2,mean, na.rm=T)
  SDsel_sd <- apply(SDsel, 2,sd, na.rm=T)
  SD_percent <- ( SDsel_sd / SDsel_av ) * 100
  SD_percent_av <- median(SD_percent,na.rm =T)
  SD_percent_av
  
  #histogram
  QMAt.h<-SDsel
  QMAt.h[is.na(QMAt.h)]<-0
  png(paste0("Results/histogram_",data.stage,"_outliers.png"), res = 300, units = "cm", width = 15, height = 10)
  raw.hist<-hist.default(QMAt.h, main = paste0("histogram of ",data.stage," data, RSD = ", SD_percent_av, "%"))
  dev.off()
}

# Before
plot_histogram(QMAt, "raw")

# SETTINGS FOR NORMALIZATION
# You can normalize against one or more internal standards, against sample weight, or both (or none).
normalize.against.IS = T
IS = c(
  "U-13C-sorbitol-D (internal standard)"#, "Glycerol (3TMS)"
  #"d4-alanine  (2TMS) (internal standard)" # no peak in our samples so do not use!
)
normalize.against.weight = T

############ 6.1 Internal Standard Normalization ###############################

if(normalize.against.IS) {
  # Throw an error if any of the internal standards are not found in the quantified metabolites
  is.not.found = setdiff(IS, colnames(QMAt))
  if(length(is.not.found) > 0) {
    stop(paste0("The following internal standards are not found in the quantified metabolites:\n",
                paste(is.not.found, collapse=", ")))
  }
  
  # Warn if IS are found in some of the samples
  for(i in IS) {
    is.missing = rownames(QMAt)[is.na(QMAt[,i])]
    if(length(is.missing) > 0) {
      warning(paste("Internal standard", i, "not found in the following samples. If no other IS can be found, the sample is set to NA:"))
      warning(paste(is.missing, collapse=", "))
    }
  }
  
  # Internal standard normalisation ############ 
  # Note the following lines only work if you have at least 2 standard
  # IS.medians = apply(QMAt[,IS], 2, FUN=function(x) { median(x, na.rm=T) } )
  # IS.Vec = apply(sapply(IS, FUN=function(x) { IS.medians[x] / QMAt[,x] } ), 1, FUN=function(x) { median(x, na.rm=T) })
  
  # If you only have one standard, use this instead
  IS.Vec<- c(median(QMAt[,IS], na.rm=T)/QMAt[,IS])
  
  norm_rib<-(QMAt*IS.Vec)
  
  norm_rib <- norm_rib[ ,!colnames(norm_rib) %in% c(IS)] # Eliminates IS (now identical everywhere)
  write.csv(t(norm_rib), "Results/PeakTable_IS_normalized_outliers.csv")
  
  plot_histogram(norm_rib, "IS_normalized_outliers")
} else {
  IS.Vec = 1
  norm_rib <- QMAt
}


############ 6.2 Weight Normalization ##########################################

if(normalize.against.weight) {
  # weight normalization
  weight.Vec = median(pheno.table$`Amount [mg(FW)]  (initial weight)`)/pheno.table$`Amount [mg(FW)]  (initial weight)`
  norm_weight<-(norm_rib * weight.Vec)
  if(normalize.against.IS) {
    write.csv(t(norm_weight), "Results/PeakTable_IS_Weight_normalized_outlier127.csv")
    plot_histogram(norm_weight, "IS_and_weight_normalized_outliers")
  } else {
    write.csv(t(norm_weight), "Results/PeakTable_Weight_normalized_outlier127.csv")
    plot_histogram(norm_weight, "weight_normalized_outliers")
  }
} else {
  weight.Vec = 1
  norm_weight <- norm_rib
}

