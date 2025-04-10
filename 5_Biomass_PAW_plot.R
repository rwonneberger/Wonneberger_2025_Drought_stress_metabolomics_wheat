# Visualize the cleaned up imaging data- from the LemnaTec

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")
library(corrplot)
library(factoextra)



# Read in the cleaned datasets and rbind them

ct<-fread("Data/All_control_sd_NA_filt2.txt")
dt<-fread("Data/All_drought_sd_NA_filt2.txt")
ct$trt<-"control"
dt$trt<-"drought"


# Calculate PAW
all_weight<-2015.4
weight_PAW<-1151.17

ct$PAW<-(ct$Weight_A_g_-all_weight)/weight_PAW*100
dt$PAW<-(dt$Weight_A_g_-all_weight)/weight_PAW*100

#Remove lines with odd PAW values
ct%<>%filter(PAW > 0)
dt%<>%filter(PAW > 0)

# Rbind
df<-rbind(ct, dt)


ct_mean<-ct[, .(biomass=mean(biomass, na.rm=T), 
                height=mean(height, na.rm=T), 
                compactness=mean(compactness, na.rm=T),
                width=mean(width, na.rm=T),
                hsv_h_mean=mean(hsv_h_mean, na.rm=T),
                hsv_h_yellow2green=mean(hsv_h_yellow2green, na.rm=T),
                hsv_h_brown2green=mean(hsv_h_brown2green, na.rm=T),
                hsv_h_red2green=mean(hsv_h_red2green, na.rm=T),
                lab_b_mean=mean(lab_b_mean, na.rm=T),
                lab_l_mean=mean(lab_l_mean, na.rm=T),
                lab_a_mean=mean(lab_a_mean, na.rm=T),
                PAW_mean=mean(PAW, na.rm=T)
), by=list(genotype, timePoint, Day_Int_)]



dt_mean<-dt[, .(biomass=mean(biomass, na.rm=T), 
                height=mean(height, na.rm=T), 
                compactness=mean(compactness, na.rm=T),
                width=mean(width, na.rm=T),
                hsv_h_mean=mean(hsv_h_mean, na.rm=T),
                hsv_h_yellow2green=mean(hsv_h_yellow2green, na.rm=T),
                hsv_h_brown2green=mean(hsv_h_brown2green, na.rm=T),
                hsv_h_red2green=mean(hsv_h_red2green, na.rm=T),
                lab_b_mean=mean(lab_b_mean, na.rm=T),
                lab_l_mean=mean(lab_l_mean, na.rm=T),
                lab_a_mean=mean(lab_a_mean, na.rm=T),
                PAW_mean=mean(PAW, na.rm=T)
                
), by=list(genotype, timePoint, Day_Int_)]


#### Make time series plots on the genotype level
ct_mean$trt<-"control"
dt_mean$trt<-"drought"
all_mean<- rbind(dt_mean, ct_mean)
names(all_mean)

all_mean$DAS<-all_mean$Day_Int_+20

df$DAS<-df$Day_Int_+20

df%<>%filter(DAS != 21)

df_melt<-as.data.table(melt(df[,-c(1,4:13, 25:26)], id.vars=c("genotype", "DAS", "Treatment")))


df_confint<-df_melt[,.(trait_mean=mean(value, na.rm = T),
                           trait_sd=sd(value, na.rm = T),
                           trait_lower=mean(value, na.rm = T) - (qt(0.975, df=19)*sd(value, na.rm = T)/sqrt(20)),
                           trait_upper=mean(value, na.rm = T) + (qt(0.975, df=19)*sd(value, na.rm = T)/sqrt(20))
),by=list(genotype, DAS, Treatment, variable)]


PAW<-data.frame(Genotype=ct_mean$Genotype, Day=ct_mean$Day, Control = ct_mean$PAW, Drought = dt_mean$PAW)
PAW<-as.data.table(PAW)

names(PAW)<-c("genotype", "DAS", "control", "drought")

PAW_melt<-melt(PAW, id.vars=c("genotype", "DAS"), variable.name = "Treatment", value.name = "PAW")

df_sub<-df_confint%>%filter(variable == "biomass")

table(PAW_melt$genotype == df_sub$genotype)
table(PAW_melt$DAS == df_sub$DAS)


df_sub1<-merge(df_sub, PAW_melt, by=c("genotype", "DAS", "Treatment"))

df_sub1$DAD<-df_sub1$DAS-21


scaleFUN <- function(x) sprintf("%.0f", x)

tiff("Figs/AdditionalFile9.tiff", height=15, width=17, units="cm", res=600, compression = "lzw")

ggplot(df_sub1, aes(x=DAD, y=trait_mean, group=Treatment)) +     
  geom_line( aes( col = Treatment)) + 
  geom_point(aes( col = Treatment), size=0.5)+ 
  geom_ribbon( aes(ymin = trait_lower , ymax = trait_upper), alpha = 0.3)+ 
  facet_wrap(~genotype, ncol=4)  +
  geom_line( aes(y=PAW*(1/.000006), col=Treatment), linewidth=0.3,linetype="dashed")+
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Biomass",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*.000006, name="% PAW")
  )+ geom_vline(xintercept=45-21, linetype="dashed", color = "red")+
  theme_classic()+ 
  theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
  xlab("DAS") +
  theme(legend.position="bottom") + 
  geom_vline(xintercept=45-21, linetype="dashed", color = "red")
dev.off()

