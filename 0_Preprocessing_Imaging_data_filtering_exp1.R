# Use the statgenHTP package to filter the imaging data time series dataset and remove outliers in Experiment 1


source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

library("statgenHTP")




## time should be in YYYY-MM-DD timepoint 

NORDrought <- read.csv("imaging1.csv", header = TRUE, sep = ",") # Additional file 3 sheet 2

NORDrought%<>%filter(Genotype != "empty pot", Genotype != "dead", Genotype != "no plant", Day != "day 24") # 17/06 should be removed because there was an issue with the imaging

NORDrought$Genotype <- as.factor(NORDrought$Genotype)
NORDrought$Treatment <- as.factor(NORDrought$Treatment)

NORDrought <- NORDrought %>% group_by(Genotype, Treatment, Day) %>% mutate(Rep= row_number())

NORDrought%<>%relocate(PlantID, Genotype, Treatment, Day, Rep)

NORDrought %<>% group_by(PlantID, Date) %>% filter(duplicated(PlantID) | n()==1)

phenoTP <- createTimePoints(dat = NORDrought,
                            experimentName = "NORDrought",
                            genotype = "Genotype",
                            timePoint = "Date",
                            repId = "Rep",
                            plotId = "PlantID")

summary(phenoTP)

## Extract the time points table.
timepoint <- getTimePoints(phenoTP)

####################################RAW DATA PLOTS#############################################

trait_list=c("biomass", "height_zoom_corrected_geometry_trait_based_on_visible_light_side_view_mm_", "compactness_16_geometry_trait_based_on_visible_light_side_view_", "width_zoom_corrected_geometry_trait_based_on_visible_light_side_view_mm_", "hsv_h_mean_color_related_trait_based_on_visible_light_side_view_", "hsv_h_brown2green_color_related_trait_based_on_visible_light_side_view_", "hsv_h_yellow2green_color_related_trait_based_on_visible_light_side_view_", "hsv_h_red2green_color_related_trait_based_on_visible_light_side_view_", "lab_b_mean_color_related_trait_based_on_visible_light_side_view_", "lab_l_mean_color_related_trait_based_on_visible_light_side_view_", "lab_a_mean_color_related_trait_based_on_visible_light_side_view_")

directory="path/to/output/dir/"

subdir_list<-c("biomass", "height", "compactness", "width", "hsv_h_mean", "hsv_h_brown2green", "hsv_h_yellow2green", "hsv_h_red2green", "lab_b_mean", "lab_l_mean", "lab_a_mean")

for (dir in subdir_list){
  dir.create(file.path(directory, dir))
}


trt_list<-c("control", "drought")

## Plot the raw data, controls and treatments together
pdf(file= paste0(directory, "Rawdata_timepoints_all.pdf"), width=10, height=8)
for (i in 1:length(trait_list)){
  trait=trait_list[i]

  plot(phenoTP, 
       traits = trait,
       plotType = "raw",
       genotypes = NULL)
  
}
dev.off()


# All subsequent analysis is performed separately for control and drought
#phenoTPc
NORDroughtc<-NORDrought%>%filter(Treatment == "control")
phenoTPc <- createTimePoints(dat = NORDroughtc,
                            experimentName = "NORDrought_control",
                            genotype = "Genotype",
                            timePoint = "Date",
                            repId = "Rep",
                            plotId = "PlantID")


#phenoTPd
NORDroughtd<-NORDrought%>%filter(Treatment == "drought")
phenoTPd <- createTimePoints(dat = NORDroughtd,
                             experimentName = "NORDrought_drought",
                             genotype = "Genotype",
                             timePoint = "Date",
                             repId = "Rep",
                             plotId = "PlantID")




## Create a boxplot for "trait" and "replicates" using the default all time points.
pdf(file= paste0(directory, "Rawdata_timepoints_boxplots_reps_all_control.pdf"), width=10, height=8)
for (i in 1:length(trait_list)){
  trait=trait_list[i]
  plot(phenoTPc, 
       plotType = "box",
       traits = trait, 
       timePoints = 1:28,
       colorBy = "repId")
}
dev.off()

pdf(file= paste0(directory, "Rawdata_timepoints_boxplots_reps_all_drought.pdf"), width=10, height=8)
for (i in 1:length(trait_list)){
  trait=trait_list[i]
plot(phenoTPd, 
     plotType = "box",
     traits = trait, 
     timePoints = 1:28,
     colorBy = "repId")
}
dev.off()


## Create a correlation plot for "trait" for a selection of time points.
pdf(file= paste0(directory, "Rawdata_timepoints_corr_all_control.pdf"), width=10, height=8)
for (i in 1:length(trait_list)){
  trait=trait_list[i]
plot(phenoTPc, 
     plotType = "cor",
     traits = trait,
     timePoints = seq(from = 1, to = 28, by = 1))
}
dev.off()


pdf(file= paste0(directory, "Rawdata_timepoints_corr_all_drought.pdf"), width=10, height=8)
for (i in 1:length(trait_list)){
  trait=trait_list[i]
  plot(phenoTPd, 
       plotType = "cor",
       traits = trait,
       timePoints = seq(from = 1, to = 28, by = 1))
}
dev.off()



############################Check for single outliers###########################

# Loop through all traits and both treatments


for (i in 1:length(trait_list)){
  for (j in 1:length(trt_list)){
    trait=trait_list[i]
    subdir=subdir_list[i]
    trt = trt_list[j]

  # Make a subset based on the treatment and create the phenoTP object
    NORDrought_sub<-NORDrought%>%filter(Treatment == trt)
    phenoTP <- createTimePoints(dat = NORDrought_sub,
                                 experimentName = paste0("NORDrought_",trt),
                                 genotype = "Genotype",
                                 timePoint = "Date",
                                 repId = "Rep",
                                 plotId = "PlantID")
    

    # Since we have a drought and a rewatering phase, there will be an abrupt change in the curve, which will often be flagged as outliers. Therefore, we will need to adjust the confIntSize and nnLocfit for some traits
    # At the same time it's very difficult to decide for the other traits what is an outlier and what isn't because there are a lot of zig zags in the last days of drought.
    
    if (trait%in%c("biomass")){
      confIntSize = 7
      nnLocfit = 0.55
    } else {
      confIntSize = 5
      nnLocfit = 0.5
    }
    
    # Find all outliers from raw data
    resuVatorHTP <- detectSingleOut(TP = phenoTP,
                                trait = trait,
                                plotIds = NULL,
                                confIntSize = confIntSize,
                                nnLocfit = nnLocfit)
    
    
    # We don't want to remove outliers around end of drought/rewatering , so we set outlier=1 after 16.06. back to 0
    resuVatorHTP[(resuVatorHTP$timePoint == as.POSIXct("2021-06-16", "UTC") | resuVatorHTP$timePoint == as.POSIXct("2021-06-17", "UTC") |  resuVatorHTP$timePoint == as.POSIXct("2021-06-18", "UTC") | resuVatorHTP$timePoint == as.POSIXct("2021-06-18", "UTC") |  resuVatorHTP$timePoint == as.POSIXct("2021-06-19", "UTC") | resuVatorHTP$timePoint == as.POSIXct("2021-06-20", "UTC") |  resuVatorHTP$timePoint == as.POSIXct("2021-06-21", "UTC") |  resuVatorHTP$timePoint == as.POSIXct("2021-06-22", "UTC")) & resuVatorHTP$outlier == 1, ]$outlier<-0

    pdf(file= paste0(directory, subdir, "/Single_outliers_", subdir, trt,"_1.pdf"), width =10, height=8)
    plot(resuVatorHTP,
         outOnly = FALSE)
    dev.off()

        
    ## remove outliers detected in first check
    phenoTPOut <- try(removeSingleOut(phenoTP,
                                       resuVatorHTP))
    if (inherits(phenoTPOut, "try-error")) {
      phenoTPOut <-phenoTP
      print("No outliers detected")
    } 
    
    
    
        
    ## save outliers 1st check
    write.table(resuVatorHTP, file =paste0(directory, subdir, "/single_outliers_", subdir, trt, "_1.txt"), append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
    
    # Find remaining outliers after first round of outliers have been removed
    resuVatorHTP1 <- detectSingleOut(TP = phenoTPOut,
                                     trait = trait,
                                     plotIds = NULL,
                                     confIntSize = confIntSize,
                                     nnLocfit = nnLocfit)
    
    resuVatorHTP1[(resuVatorHTP1$timePoint == as.POSIXct("2021-06-16", "UTC") | resuVatorHTP1$timePoint == as.POSIXct("2021-06-17", "UTC") |  resuVatorHTP1$timePoint == as.POSIXct("2021-06-18", "UTC") | resuVatorHTP1$timePoint == as.POSIXct("2021-06-18", "UTC") |  resuVatorHTP1$timePoint == as.POSIXct("2021-06-19", "UTC") | resuVatorHTP1$timePoint == as.POSIXct("2021-06-20", "UTC") |  resuVatorHTP1$timePoint == as.POSIXct("2021-06-21", "UTC") |  resuVatorHTP1$timePoint == as.POSIXct("2021-06-22", "UTC")) & resuVatorHTP1$outlier == 1, ]$outlier<-0
    
    
    #save plots
    pdf(file= paste0(directory,  subdir, "/Single_outliers_", subdir, trt,"_2.pdf"), width=10, height=8)
    plot(resuVatorHTP1,
         outOnly = FALSE)
    dev.off()
    
  
    
    ## remove outliers 2nd run
   
    phenoTPOut1 <- try(removeSingleOut(phenoTPOut,
                                       resuVatorHTP1))
    if (inherits(phenoTPOut1, "try-error")) {
      phenoTPOut1 <-phenoTPOut
      print("No outliers detected")
    } 
    

    
    ## save outliers 2nd run
    write.table(resuVatorHTP1, file = paste0(directory,  subdir, "/single_outliers_", subdir, trt, "_2.txt"), append = FALSE, quote = TRUE, sep = " ",
                eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                col.names = TRUE, qmethod = c("escape", "double"),
                fileEncoding = "")
    
    
    # Make final plot after outliers have been removed twice
    resuVatorHTP2 <- detectSingleOut(TP = phenoTPOut1,
                                     trait = trait,
                                     plotIds = NULL,
                                     confIntSize = confIntSize,
                                     nnLocfit = nnLocfit)
    
    resuVatorHTP2[(resuVatorHTP2$timePoint == as.POSIXct("2021-06-16", "UTC") | resuVatorHTP2$timePoint == as.POSIXct("2021-06-17", "UTC") |  resuVatorHTP2$timePoint == as.POSIXct("2021-06-18", "UTC") | resuVatorHTP2$timePoint == as.POSIXct("2021-06-18", "UTC") |  resuVatorHTP2$timePoint == as.POSIXct("2021-06-19", "UTC") | resuVatorHTP2$timePoint == as.POSIXct("2021-06-20", "UTC") |  resuVatorHTP2$timePoint == as.POSIXct("2021-06-21", "UTC") |  resuVatorHTP2$timePoint == as.POSIXct("2021-06-22", "UTC")) & resuVatorHTP2$outlier == 1, ]$outlier<-0
    
    
    #save plots
    pdf(file= paste0(directory,  subdir, "/Single_outliers_", subdir, trt,"_final.pdf"), width=10, height=8)
    plot(resuVatorHTP2,
         outOnly = FALSE)
    dev.off()

    
    ########################Check for outliers in timeseries########################
    # Fit P-splines using data where outliers are removed twice
    phenoTP_df <- as.data.frame(phenoTPOut1)
    
    #save phenotype file after second run of outliers have been removed
    write.table(phenoTP_df, file = paste0(directory,  subdir, "/", subdir, trt, "_noout.txt"), append = FALSE, quote = TRUE, sep = " ",
                eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                col.names = TRUE, qmethod = c("escape", "double"),
                fileEncoding = "")
    
 

    rm(NORDrought_sub, phenoTP, phenoTPOut, phenoTPOut1, resuVatorHTP, resuVatorHTP1, resuVatorHTP2)

  }}    
    
    
    
#######################################################
# Read in all files with the cleaned data, select the respective trait column and merge to a new file which can be used for further analyses - manually change drought/control as desired

control_list=list()

for (i in 1:length(subdir_list)){
  control_list[[i]]<-(paste0(directory, subdir_list[i], "/", subdir_list[i], "control_noout.txt"))
}

df0<-fread(control_list[[1]])
df00<-df0%>%select(2:15)

for (i in control_list){
  df<-fread(i)
  col_select<-which(control_list == i)+15
  df1<-df%>%select(2, 5,col_select)
  df00<-left_join(df00, df1, by=c("plotId", "Day"))
 
}

names(df00)[15:25]<-subdir_list
df00$Replicate<-NULL
fwrite(df00, paste0(directory, "All_control_noout.txt"), sep="\t")

       

###################################
# Add an SD filter, and remove reps with more than x missing datapoints
df<-fread(paste0(directory, "All_control_noout.txt"))

d1<-df%>%filter(genotype == "PPPW_034")

d1_long<-melt(d1, id.vars=c("plotId", "genotype", "Treatment", "Day", "Time", "timePoint", "Day_Int_", "Day_Float_", "Weight_A_g_", "Weight_B_g_", "Water_sum_of_day_", "Water_weight.diff_", "repId"))


d1_long<-d1_long%>%filter( genotype == "PPPW_034", variable == "hsv_h_red2green")

ggplot(d1_long, aes(x=Day_Int_, y=value)) + geom_point() +  stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha = 0.2)
ggplot(d1_long, aes(x=Day_Int_, y=value, group=plotId, color=plotId )) + geom_point() + geom_line(data=d1_long[!is.na(d1_long$value),])

# Remove values with an SD > 2
# Remove plotIds if more than 20% data points missing (round(29 days * 0.15) = 5)

df_long<-melt(df, id.vars=c("plotId", "genotype", "Treatment", "Day", "repId","Time", "timePoint", "Day_Int_", "Day_Float_", "Weight_A_g_", "Weight_B_g_", "Water_sum_of_day_", "Water_weight.diff_"))
df_filtered_long <- df_long %>%  
  group_by(genotype, Day_Int_, variable) %>%                # summarize by name
  mutate(value = replace(value, abs(value - mean(value)) > 2*sd(value), NA)) %>%
  group_by(plotId, variable) %>%     
  mutate(num_na = sum(is.na(value))) %>% 
  filter(!any(num_na >=5)) %>%
  select(-num_na)%>%ungroup



df_filtered <- df_filtered_long %>% pivot_wider(names_from = variable, values_from = value)

d2<-df_filtered_long%>%filter(genotype == "PPPW_034", variable == "hsv_h_red2green")
ggplot(d2, aes(x=Day_Int_, y=value )) + geom_point() + stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha = 0.2)
ggplot(d2, aes(x=Day_Int_, y=value, group=plotId, color=plotId  )) + geom_point() + geom_line(data=d2[!is.na(d2$value),])

fwrite(df_filtered, paste0(directory, "All_control_sd_NA_filt.txt"), sep="\t") 







