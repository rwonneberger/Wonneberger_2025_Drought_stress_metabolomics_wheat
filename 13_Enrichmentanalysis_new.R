# Metabolimics enrichment analysis using FELLA

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("FELLA")

library(FELLA)
library(biomaRt)

# Download wheat data from KEGG and build a graph
graph <- buildGraphFromKEGGREST(organism = "taes")
tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(
  keggdata.graph = graph,
  databaseDir = tmpdir,
  internalDir = FALSE,
  matrices = "diffusion",
  normality =  "diffusion",
  niter = 50)

fella.data <- loadKEGGdata(
  databaseDir = tmpdir,
  internalDir = FALSE,
  loadMatrix =  "diffusion"
)

cat(getInfo(fella.data))

# Read in the background set
background<-fread("../Data/KEGG_lookup_cnumbers.csv", header=F)
background<-as.vector(background)$V1
background<-background[background!=""]

outlist<-list()

files_to_run<-c("Metab_list_T_MA_IS_Weight_norm_time_2_trt_all_tol_all_DAMs_ttest.txt", "Metab_list_T_MA_IS_Weight_norm_time_3_trt_all_tol_all_DAMs_ttest.txt", "Metab_list_T_MA_IS_Weight_norm_time_4_trt_all_tol_all_DAMs_ttest.txt")

for (i in 1:3){

# Read in the metabolites of interest
dams<-fread(files_to_run[i])[,2]
dams<-as.vector(dams)$KEGG
dams<-dams[dams!=""]

# Map input metabolites to the graph
analysis_dams <- enrich(
  compounds = dams, 
  method = "diffusion", 
  approx = "normality", 
  compoundsBackground = background,
  data = fella.data)

#show(analysis_dams)

# Get metabolites that do not belong to the background. Actually quite a few...
#getExcluded(analysis_dams)

# tiff(paste0("\\\\storage-al.slu.se/home$/rawo0001/My Documents/Projects/CresWheat/Manuscripts/Metabolome_Drought_stress/Draft1/Tables_Figs/FELLA_enrichment_graph_M", i+1,".tiff"), width=50, height=50, units="cm", res=100)
# plot(
#   x = analysis_dams,
#   method = "diffusion",
#   main = i+1,
#   threshold = 0.01,
#   nlimit = 1000,
#   data = fella.data)
# dev.off()

myTable <- generateResultsTable(
  object = analysis_dams, 
  method = "diffusion", 
  threshold = 0.01, 
  nlimit = 1000,
  data = fella.data)
#myTable%>%filter(Entry.type == "pathway")%>%arrange(p.score)

outlist[[i]]<-myTable%>%filter(Entry.type == "pathway")%>%arrange(p.score)
}


day2<-outlist[[1]]
day3<-outlist[[2]]
day4<-outlist[[3]]


day2$Day<-"13 DAD"
day3$Day<-"22 DAD"
day4$Day<-"28 DAD"
alldays<-rbind(day2, day3)
alldays<-rbind(alldays, day4)
alldays$KEGG.name<-gsub(" - .*", "", alldays$KEGG.name)

names(alldays)[4]<-"p-value"

fwrite(alldays, "Enrichment_234_FELLA_output.txt", sep="\t")

ggplot(alldays, aes(x=Day, y=KEGG.name, size=`p-value`, color=`p-value`)) + geom_point()+theme_bw()+
  guides(colour = guide_legend())+xlab("Time point")+ylab("Pathway")+ x_axis_90_degrees












## Recovery stage

outlist<-list()

files_to_run<-c("Metab_list_T_MA_IS_Weight_norm_time_34_trt_drought_tol_all_DAMs_ttest.txt", "Metab_list_T_MA_IS_Weight_norm_time_34_trt_drought_tol_all_up_DAMs_ttest.txt", "Metab_list_T_MA_IS_Weight_norm_time_34_trt_drought_tol_all_down_DAMs_ttest.txt")

for (i in 1:3){
  
  # Read in the metabolites of interest
  dams<-fread(files_to_run[i])[,2]
dams<-as.vector(dams)$KEGG
dams<-dams[dams!=""]

# Map input metabolites to the graph
analysis_dams <- enrich(
  compounds = dams, 
  method = "diffusion", 
  approx = "normality", 
  compoundsBackground = background,
  data = fella.data)

#show(analysis_dams)

# Get metabolites that do not belong to the background. Actually quite a few...
#getExcluded(analysis_dams)


plot(
x = analysis_dams,
method = "diffusion",
main = "My first enrichment using the diffusion analysis in FELLA",
threshold = 0.1,
nlimit = 200,
data = fella.data)

myTable <- generateResultsTable(
  object = analysis_dams, 
  method = "diffusion", 
  threshold = 0.01, 
  nlimit = 1000,
  data = fella.data)
outlist[[i]]<-myTable%>%filter(Entry.type == "pathway")%>%arrange(p.score)
}

recovery_all<-outlist[[1]]
recovery_up<-outlist[[2]]
recovery_down<-outlist[[3]]


recovery_all$Analysis<-"All DAMs"
recovery_up$Analysis<-"DAMs up"
recovery_down$Analysis<-"DAMs down"
allrecovery<-rbind(recovery_all, recovery_up)
allrecovery<-rbind(allrecovery, recovery_down)
allrecovery$KEGG.name<-gsub(" - .*", "", allrecovery$KEGG.name)


allrecovery$Analysis<-factor(allrecovery$Analysis, levels=c("All DAMs", "DAMs up", "DAMs down"))
names(allrecovery)[4]<-"p-value"

fwrite(allrecovery, "Enrichment_recovery_FELLA_output.txt", sep="\t")

ggplot(allrecovery, aes(x=Analysis, y=KEGG.name, size=`p-value`, color=`p-value`)) + geom_point()+theme_bw()+
  guides(colour = guide_legend())+xlab("Analysis")+ylab("Pathway")+x_axis_90_degrees 




