# analysis of metabolomics data using MetaboAnalystR



writetopath<-"path/"

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
phenotype = c("control",  "drought")
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
#mSet<-PlotNormSummary(mSet, paste0("../Figs/", filestring, "_norm_"), "png", 72, width=NA) # Broken??
mSet<-PlotSampleNormSummary(mSet, paste0("../Figs/", filestring, "_snorm_"), "png", 72, width=NA)


###################### Extract the normalized data which we're going to use for most purposes from now on ######################
norm_data<-as.data.frame(mSet$dataSet$norm)
dim(norm_data)

# write the normalized data for potential analysis later
norm_out<-norm_data
norm_out$Sample<-rownames(norm_out)
norm_out%<>%relocate(Sample)
fwrite(norm_out, paste0("Normalized_", filestring, ".txt"), sep="\t")


########################################################################################
###################### Preprocessing done ##############################################
########################################################################################



###################### Fold-changes, t-tests and Volcano plots ######################

#For the complete list of all FC

data <- qs::qread("row_norm.qs");
fwrite(as.data.frame(data), paste0(filestring, "_rownorm.txt"),row.names=T, sep="\t")
m1 <- colMeans(data[which(mSet$dataSet$cls==levels(mSet$dataSet$cls)[1]), ]);
m2 <- colMeans(data[which(mSet$dataSet$cls==levels(mSet$dataSet$cls)[2]), ]);

ratio <- m2/m1 #a value log2(value) > 0 means that it's up in drought

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
norm_data1%<>%dplyr::select(-Sample,  -Subject, -Batch)
norm_data1

norm_data1_mean<-norm_data1%>%
  group_by(Phenotype) %>%
  summarise_all("mean")%>%dplyr::select(-Time, -Tolerance)
norm_data1_mean<-transpose_df_Col1(norm_data1_mean, "Metabolites")
norm_data1_mean$direction<-""
norm_data1_mean<-norm_data1_mean[complete.cases(norm_data1_mean),]
norm_data1_mean[norm_data1_mean$control < norm_data1_mean$drought,]$direction<-"up"
norm_data1_mean[norm_data1_mean$control > norm_data1_mean$drought,]$direction<-"down"

norm_data1_mean_merged<-merge(norm_data1_mean, ttest_all, by="Metabolites")

fwrite(norm_data1_mean_merged, paste0( filestring, "_ttall.txt"), row.names=T, sep="\t")

names(fc.log)[2]<-"Metabolites"
out_all<-merge(fc.log, norm_data1_mean_merged, by="Metabolites")
names(out_all)<-c("Metabolite", "log2(FC)", "mean_control", "mean_drought", "Direction", "Tstat", "t_pvalue", "T_FDR")

fwrite(out_all, paste0(filestring, "_all_FC_ttest.txt"), row.names=T, sep="\t")

ttest1<-left_join(ttest, norm_data1_mean)


fwrite(ttest1, paste0( filestring, "_tt.txt"), row.names=T, sep="\t")
mSet<-PlotTT(mSet, paste0("../Figs/", filestring, "_Ttest"), "png", 72, width=NA)
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 1, F, 0.1, TRUE, "fdr")



volcano_df<-data.frame("FC" = mSet$analSet$fc$fc.log, "TT" = mSet$analSet$tt$p.log)
volcano_df$Metabolite<-rownames(volcano_df)
volcano_df$Status <- "Not sig."
volcano_df$Status[volcano_df$FC > 1 & volcano_df$TT > -log10(0.05)] <- "Up in drought"
volcano_df$Status[volcano_df$FC < -1 & volcano_df$TT > -log10(0.05)] <- "Down in drought"
volcano_df$damlabel <- NA
volcano_df$damlabel[volcano_df$Status != "Not sig."] <- volcano_df$Metabolite[volcano_df$Status != "Not sig."]

mycolors <- c("#1984c5", "#c23728",  "#e2e2e2")

names(mycolors) <- c("Down in drought", "Up in drought", "Not sig.")

p1<-ggplot(volcano_df, aes(FC, TT, color = Status, label=damlabel)) + geom_point() + geom_hline(yintercept=-log10(0.05), linetype="dashed")+ theme_minimal() + geom_vline(xintercept=c(-1, 1), linetype="dashed")+ scale_colour_manual(values = mycolors)+ geom_text_repel(size = 2, show_guide=FALSE, max.overlaps = 15)  + xlab("log2(fold change)") + ylab("-log10(p), FDR-adjusted") + theme(legend.position = "right", text = element_text(size = 8))




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


if (all(time ==c(1,2,3,4))){
## Plot selected mets over time, per geno


normalized<-fread(paste0("/Results/", filestring, "_rownorm.txt"), sep="\t")
peek(normalized)
normalized$Geno<-substr(normalized$V1, 1,8)
normalized$Treatment<-substr(normalized$V1, 10,10)
normalized$DAD<-substr(normalized$V1, 12,12)

normalized[normalized$DAD == 1,]$DAD <- "7"
normalized[normalized$DAD == 2,]$DAD <- "13"
normalized[normalized$DAD == 3,]$DAD <- "22"
normalized[normalized$DAD == 4,]$DAD <- "28"

normalized$DAD<-factor(normalized$DAD, levels =c(7, 13, 22, 28))

normalized <- normalized %>% 
  mutate(Treatment = if_else(Treatment == "c", "Control", "Drought"))

norm_data_melt<-as.data.table(melt(normalized))
norm_data_melt_avg<-norm_data_melt[,.(value_mean=mean(value, na.rm=T)),by=list(Geno, Treatment, DAD, variable)]
norm_data_melt_avg$mean_log<-log10(norm_data_melt_avg$value_mean)

norm_data_melt_avg%>%filter(variable%like%"Alanine")

mets<-c("Alanine (3TMS)", "Raffinose (11TMS)", "Proline (2TMS)", "Spermine (6TMS)", "Spermidine (4TMS)", "Histamine (3TMS)", "Putrescine (4TMS)", "Nicotinamide (2TMS)")


fontsize_text=7
knownmet<-ggplot(norm_data_melt_avg%>%filter(variable%in%mets), aes(x=DAD, y=value_mean, group=Treatment, color=Treatment)) + geom_line()+ facet_grid(variable~Geno, scales="free_y") +  theme_classic()+ 
  theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ ylab("Row-normalized peak intensity")+theme(legend.position = "bottom")+ theme(text = element_text(size = fontsize_text),axis.text.x = element_text(size = fontsize_text-1))

tiff(paste0(writetopath, "Figs/AdditionalFile14.tiff"), width=17, height=20, units="cm", compression = "lzw",res=600)

print(knownmet)
dev.off()

}




normalized<-fread(paste0("/Results/", filestring, "_rownorm.txt"), sep="\t")

volcano_df_sub<-volcano_df%>%filter((FC > 1 | FC < -1), TT > -log10(0.05))

normalized_sub<-normalized%>%dplyr::select(V1, all_of(volcano_df_sub$Metabolite))

normalized_sub$Genotype<-substr(normalized_sub$V1, start = 1, stop = 8)
normalized_sub$Trt<-substr(normalized_sub$V1, start = 10, stop = 10)

normalized_sub <- normalized_sub %>% 
  mutate(Trt2 = if_else(Trt == "c", "Control", "Drought"))
normalized_sub$V1<-normalized_sub$Trt<-NULL

normalized_sub %<>%
  group_by(Genotype, Trt2)%>%
  summarise_if(is.numeric, mean, na.rm=T)

normalized_sub_melt<-melt(normalized_sub)

p2<-ggplot(normalized_sub_melt,aes( x=Genotype, y=variable, fill=log2(value)))+geom_tile()+ theme_classic()+ 
  theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA), axis.text = element_text(size = 7))+
  facet_grid(~Trt2, scales="free_x")+
  scale_fill_gradientn(
    colours = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(255),
    values = c(1.0, (0 - min(log2(normalized_sub_melt$value))) / (max(log2(normalized_sub_melt$value)) - min(log2(normalized_sub_melt$value))), 0))+
  theme(axis.title=element_blank())+x_axis_90_degrees + theme(legend.position = "bottom")
p2




# Make a pathway enrichment plot for the last three sampling time points


enrich<-fread("Enrichment_234_FELLA_output.txt")

tiff(paste0(writetopath, "Figs/Fig9.tiff"), width=8.5, height=10, units="cm", compression = "lzw",res=600)
ggplot(enrich, aes(x=Day, y=KEGG.name, size=-log10(`p-value`), color=-log10(`p-value`))) + geom_point()+theme_bw()+
 xlab("Analysis")+x_axis_90_degrees + theme(legend.text = element_text(size=7), legend.position="bottom")+ theme(axis.title=element_blank(), axis.text=element_text(size=6), legend.title=element_text(size=7)) + guides(color=guide_legend(title="-log10(p)", title.position = "top", ncol=2), size=guide_legend(title="-log10(p)", title.position = "top", ncol=2))
dev.off()




