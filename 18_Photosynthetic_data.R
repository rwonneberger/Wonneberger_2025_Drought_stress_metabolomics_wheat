source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

library(multcompView)
library(psych)


psi<-fread("Data/PSI.txt")
  
names(psi)
names(psi)[8]<-"Treatment"

# Remove values outside 2 SDs


psi_long<-melt(psi, id.vars=c("Measuring Date", "Measuring Time", "Timepoint", "Day", "Round Order", "Plant ID", "Genotype", "Treatment", "Size"))
psi_filtered_long <- psi_long %>%  
  group_by(Genotype, Day, `Treatment`, variable) %>%                # summarize by name
  mutate(value = replace(value, abs(value - mean(value)) > 2*sd(value), NA)) %>%
  group_by(`Plant ID`, variable) %>%     
  mutate(num_na = sum(is.na(value))) %>% 
  filter(!any(num_na >=5)) %>%
  dplyr::select(-num_na)%>%ungroup



psi_filtered <- psi_filtered_long %>% pivot_wider(names_from = variable, values_from = value)




df_out1<-data.frame()

for (i in unique(psi_filtered_long$Day)){
  for (j in unique(psi_filtered_long$Genotype)){
    
    sub_df<-psi_filtered_long%>%filter(Day == i, Genotype == j, variable == "QY_Lss1")
    
    sub_df$Treatment<-as.factor(sub_df$Treatment)
    
    
    #Runs Tukeys HSD test on all line/year combos to determine which are significantly different and assigns letters for comparison
    exp_tukey <- TukeyHSD(exp_aov <- aov(value ~ Treatment, data = sub_df))
    exp_tukey_result <- data.frame(exp_tukey$Treatment)
    
    
    #Create a matrix to hold p-values
    p_mat <- matrix(NA,length(levels(sub_df$Treatment)),length(levels(sub_df$Treatment)))
    
    
    #Extract p-values from Tukey HSD
    p_mat[lower.tri(p_mat)] <- rownames(exp_tukey_result)
    p_mat[lower.tri(p_mat)] <- round(exp_tukey_result$p.adj,3)
    p_mat[upper.tri(p_mat)] <- t(p_mat[lower.tri(p_mat)])
    rownames(p_mat) <- levels(sub_df$Treatment)
    colnames(p_mat) <- levels(sub_df$Treatment)
    
    #Assign letters to indicate if groups are significantly different (adjusted p-value <0.05)
    letters <- multcompLetters(p_mat,threshold=0.05)
    letters <- data.frame(letters$Letters)
    letters$Treatment <- row.names(letters)
    letters$Genotype<-j
    letters$Day<-i
    
    
    #Merge letters with trait data
    df_subset <- merge(sub_df,letters,by=c("Genotype", "Treatment", "Day"))
    df_out1<-rbind(df_out1, df_subset)
    
  }}


df_out1$value<-as.numeric(df_out1$value)






df_out2<-data.frame()

for (i in unique(psi_filtered_long$Day)){
  for (j in unique(psi_filtered_long$Genotype)){
    
    sub_df<-psi_filtered_long%>%filter(Day == i, Genotype == j, variable == "QY_Lss2")
    
    sub_df$Treatment<-as.factor(sub_df$Treatment)
    
    
    #Runs Tukeys HSD test on all line/year combos to determine which are significantly different and assigns letters for comparison
    exp_tukey <- TukeyHSD(exp_aov <- aov(value ~ Treatment, data = sub_df))
    exp_tukey_result <- data.frame(exp_tukey$Treatment)
    
    
    #Create a matrix to hold p-values
    p_mat <- matrix(NA,length(levels(sub_df$Treatment)),length(levels(sub_df$Treatment)))
    
    
    #Extract p-values from Tukey HSD
    p_mat[lower.tri(p_mat)] <- rownames(exp_tukey_result)
    p_mat[lower.tri(p_mat)] <- round(exp_tukey_result$p.adj,3)
    p_mat[upper.tri(p_mat)] <- t(p_mat[lower.tri(p_mat)])
    rownames(p_mat) <- levels(sub_df$Treatment)
    colnames(p_mat) <- levels(sub_df$Treatment)
    
    #Assign letters to indicate if groups are significantly different (adjusted p-value <0.05)
    letters <- multcompLetters(p_mat,threshold=0.05)
    letters <- data.frame(letters$Letters)
    letters$Treatment <- row.names(letters)
    letters$Genotype<-j
    letters$Day<-i
    
    
    #Merge letters with trait data
    df_subset <- merge(sub_df,letters,by=c("Genotype", "Treatment", "Day"))
    df_out2<-rbind(df_out2, df_subset)
    
  }}


df_out2$value<-as.numeric(df_out2$value)











df_out3<-data.frame()

for (i in unique(psi_filtered_long$Day)){
  for (j in unique(psi_filtered_long$Genotype)){
    
    sub_df<-psi_filtered_long%>%filter(Day == i, Genotype == j, variable == "Ratio_QY_LH")
    
    sub_df$Treatment<-as.factor(sub_df$Treatment)
    
    
    #Runs Tukeys HSD test on all line/year combos to determine which are significantly different and assigns letters for comparison
    exp_tukey <- TukeyHSD(exp_aov <- aov(value ~ Treatment, data = sub_df))
    exp_tukey_result <- data.frame(exp_tukey$Treatment)
    
    
    #Create a matrix to hold p-values
    p_mat <- matrix(NA,length(levels(sub_df$Treatment)),length(levels(sub_df$Treatment)))
    
    
    #Extract p-values from Tukey HSD
    p_mat[lower.tri(p_mat)] <- rownames(exp_tukey_result)
    p_mat[lower.tri(p_mat)] <- round(exp_tukey_result$p.adj,3)
    p_mat[upper.tri(p_mat)] <- t(p_mat[lower.tri(p_mat)])
    rownames(p_mat) <- levels(sub_df$Treatment)
    colnames(p_mat) <- levels(sub_df$Treatment)
    
    #Assign letters to indicate if groups are significantly different (adjusted p-value <0.05)
    letters <- multcompLetters(p_mat,threshold=0.05)
    letters <- data.frame(letters$Letters)
    letters$Treatment <- row.names(letters)
    letters$Genotype<-j
    letters$Day<-i
    
    
    #Merge letters with trait data
    df_subset <- merge(sub_df,letters,by=c("Genotype", "Treatment", "Day"))
    df_out3<-rbind(df_out3, df_subset)
    
  }}


df_out3$value<-as.numeric(df_out3$value)





pointsize=0.2
fontsize_text=7

p1<-ggplot(df_out1, aes(x=Treatment, y=value)) + geom_violin()+ geom_boxplot(width=0.3, fill="#D55E00") +geom_jitter(width=0.2, size=0.2)+ facet_grid(Day~Genotype, scales="free_y")+ theme_classic()+ 
  theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ geom_text(aes(y=max(value, na.rm=T)*1.05,label=letters.Letters), size=3, vjust="inward") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(axis.title=element_blank(), axis.text.x=element_blank()) + ggtitle("QY_Lss1")+  
  theme(text = element_text(size = fontsize_text), plot.title = element_text(size=fontsize_text)) 

p2<-ggplot(df_out2, aes(x=Treatment, y=value)) + geom_violin()+ geom_boxplot(width=0.3, fill="#D55E00") +geom_jitter(width=0.2, size=0.2)+ facet_grid(Day~Genotype, scales="free_y")+ theme_classic()+ 
  theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ geom_text(aes(y=max(value, na.rm=T)*1.05,label=letters.Letters), size=3, vjust="inward") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(axis.title=element_blank(), axis.text.x=element_blank()) + ggtitle("QY_Lss2")+  
  theme(text = element_text(size = fontsize_text), plot.title = element_text(size=fontsize_text)) 


p3<-ggplot(df_out3, aes(x=Treatment, y=value)) + geom_violin()+ geom_boxplot(width=0.3, fill="#D55E00") +geom_jitter(width=0.2, size=0.2)+ facet_grid(Day~Genotype, scales="free_y")+ theme_classic()+ 
  theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ geom_text(aes(y=max(value, na.rm=T)*1.05,label=letters.Letters), size=3, vjust="inward") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(axis.title=element_blank()) + ggtitle("Ratio_QY_LH")+  
  theme(text = element_text(size = fontsize_text), plot.title = element_text(size=fontsize_text)) 
#dev.off()




tiff(paste0("Figs/AdditionalFile13.tiff"), width=17, height=22, compression = "lzw", res=600, units="cm")

ggarrange(p1, p2, p3, ncol=1, labels=c("a", "b", "c"))
dev.off()





#

















