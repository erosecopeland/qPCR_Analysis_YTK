library(tidyverse)
library(scales)
library(multcomp)
library(RColorBrewer)
library(ggsignif)

###### Import data #####
raw <- read.csv("data/Updated_Raw_Data.csv") #exported from MicPCR
stdcurve <- read_csv("data/Std.Curve-338F.528R.csv") #exported from MicPCR

#Metadata containing:
#Sample name
#Mic qPCR Run
#Dilution factor of input template into qPCR
#Sample Info
#DNA Concentration
#Amount of starting material used (in grams or ml)
#Volume DNA was eluted (in uL)

metadata <- read.csv("data/YTK_samples_stats.csv", header = T) #Based on the data

##### Standard Curve Parameters ####
# Calculate standard curve parameters
fit <- lm(stdcurve$Cq ~ log10(stdcurve$`Given Concentration`))
std_slope <- round(fit$coefficients[[2]], 3)
std_intercept <- round(fit$coefficients[[1]], 3)
std_rsq <- round(summary(fit)$r.squared, 3)
std_efficiency <- round(-1 + 10^(-1/std_slope), 2)

##### Calculations and adding in metadata ####

#colours from previous graphs in publication
colours <- c("#FFFE03","#00FF01", "#F5941F","#F21A26","#C2C2C2")
#water, rotifer, artermia, micro, macro

#remove any samples or columns not required
raw <- raw %>%
  dplyr::filter(Sample != "NTC") %>%
  arrange(Sample) %>%
  dplyr::select(-Calculated_Concentration) %>%
  dplyr::filter(Sample != "R_T2_n2")

metadata <- metadata %>%
  dplyr::filter(Sample != "NTC" | Sample != "R_T2_n2")

# Calculate average and standard error (should be less than 0.25) of technical replicates
calc <- raw %>%
  group_by(Sample) %>%
  dplyr::summarise(avg_Cq = mean(Cq), std_dev = sd(Cq), std_err = sd(Cq)/sqrt(n())) %>%
  data.frame()

# Calculate copies per reaction
calc <- calc %>%
  dplyr::mutate(Copiesperrxn = 10^((avg_Cq-std_intercept)/(std_slope)))

# Correct for :
#1. dilution factor used when preparing qPCR
#2. extraction volume (uL eluted)
#3. DNA yield
#4. The initial amount of material extracted

Extraction_Vol <- 200

calc <- calc %>%
  left_join(metadata) %>%
  dplyr::mutate(initial_copies_peruL = Copiesperrxn*Dilution_Factor) %>%
  dplyr::mutate(Total_Yield = DNA_Conc*Extraction_Vol) %>%
  dplyr::mutate(CopiesperExtraction = initial_copies_peruL*Extraction_Vol) %>%
  dplyr::mutate(FinalCopiesperng = CopiesperExtraction/Total_Yield) %>%
  dplyr::mutate(FinalCopiesperngperg = FinalCopiesperng/Starting_amount) %>%
  dplyr::mutate(FinalCopiesperg = CopiesperExtraction/Starting_amount)

  
#selecting only columns of interest  
dat <- calc %>%
  dplyr::select(c("Sample", "Sample_Stage", "Sample_Type", "Pair", "DNA_Conc", "Total_Yield", "FinalCopiesperngperg", "FinalCopiesperg"))


#### Manipulation for graphing ####

# Coerce metadata columns to factors
cols_factor <- c("Sample", "Sample_Stage", "Sample_Type", "Pair")
dat[cols_factor] <- lapply(dat[cols_factor], factor)

ref1 <- c("L1", "L3-9", "L14-18", "L29-35", "L49-53", "Water","Rot", "Art", "Micro", "Macro" )
ref2 <- c("Yolk_Stage", "Rot_Stage", "Art_Stage", "Micro_Stage", "Macro_Stage")
ref3 <- c("Y_T1_n2", "Y_T1_n5", "Y_T1_n9", "W1", "W2", "W3", "R_T2_n4", "R_T2_n5", "R_T2_n8", "R1", "R2", "R3", 
          "A_T2_n5", "A_T2_n8", "A_T2_n9", "A1", "A2", "A3","Mi_T2_n2", "Mi_T2_n4", "Mi_T2_n10", "MI1", "MI2", "MI3",
          "Ma_T2_n3", "Ma_T2_n5", "Ma_T2_n10", "MA1", "MA2", "MA3")
ref4 <- c("Larvae", "Feed/Water")
dat$Sample_Stage <- factor(dat$Sample_Stage, levels = ref1)
dat$Pair <- factor(dat$Pair, levels = ref2)
dat$Sample <- factor(dat$Sample, levels = ref3)
dat$Sample_Type <- factor(dat$Sample_Type, levels = ref4)
levels(dat$Sample_Stage)
levels(dat$Pair)
levels(dat$Sample)

dat <- arrange(dat, Sample_Stage)

# make copy for records
write.csv(dat, "outputs/YTK_qPCR_data.csv")

#make a copy for qPCR community - Copyrighter
qPCR <- dat %>%
  slice(c(28:30, 19:21, 1,3, 13:15, 7:9)) %>%
  rename("Community Name" = Sample, "Total Abundance" = FinalCopiesperngperg) %>%
  select(c(`Community Name`, `Total Abundance`)) %>%
  as.data.frame() 
write.csv(qPCR, "outputs/qPCR_Total_Community.csv", row.names = F)


#### Graphing ####

colnames(dat)
#can try to graph first
#DNA YIELD #

#dat <- dat %>%
#  group_by(Sample_Stage) %>%
#  dplyr::mutate(DNA_Conc_sd = sd(DNA_Conc))

ggplot(dat, aes(x=Sample_Stage, y = DNA_Conc, fill = Pair)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(~Sample_Type, scales = "free", labeller = label_wrap_gen(width = 30)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
      axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=8), 
      strip.text = element_text(size = 8), legend.position = "bottom",
      panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_manual(values=colours)
#geom_errorbar(aes(ymin=DNA_Conc_sd, ymax=DNA_Conc_sd), width=.1)

ggsave("outputs/DNA_Conc_barplot.pdf")

ggplot(dat, aes(x=Sample_Stage, y = DNA_Conc, fill = Pair)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~Sample_Type, scales = "free_x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=8),
        strip.text = element_text(size = 14),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.x = element_blank(), 
        legend.position = "none") +
  ylab("DNA Concentration (ng/uL)") +
  scale_fill_manual(values=colours)

ggsave("outputs/DNA_Conc_boxplot.pdf")

#### COPY NO NORMALISED BY GRAMS MATERIAL ####
#WE ARE NOT PUBLISHING THE FOOD INFO#

temp <- expression("ANOVA, p < 0.005")

dat %>%
  filter(Sample_Type == "Larvae") %>%
ggplot(aes(x = Sample_Stage, y = FinalCopiesperg, fill = Pair)) +
  geom_boxplot() +
  scale_y_log10(#limits = c(1000000, 1000000000),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"), legend.text=element_text(size=12), 
        strip.text = element_text(size = 14),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.x = element_blank(),
        legend.position = "none") +
  #geom_hline(yintercept=272, linetype="dashed", color = "red") +
  ylab(bquote('16S rRNA Gene Copies ('* g^-1*')')) +
  scale_fill_manual(values = colours) +
  geom_signif(comparisons=list(c("L1", "L14-18")), annotation=c("**"), y_position = 9.7, vjust=0.4) +
  geom_signif(comparisons=list(c("L1", "L29-35")),annotation=c("***"), y_position = 9.9, vjust=0.4) +
  geom_signif(comparisons=list(c("L1", "L49-53")),annotation=c("**"), y_position = 10.1, vjust=0.4) +
  geom_signif(comparisons=list(c("L3-9", "L14-18")),annotation=c("*"), y_position = 8.9, vjust=0.4) +
  geom_signif(comparisons=list(c("L3-9", "L29-35")),annotation=c("**"), y_position = 9.1, vjust=0.4) +
  geom_signif(comparisons=list(c("L3-9", "L49-53")),annotation=c("*"), y_position = 9.3, vjust=0.4) +
  annotate("text", x = 4.5, y = 1600000, label = as.character(temp), size = 5)

ggsave("outputs/Copiesperg_boxplot_logaxis.pdf")

#http://www.biotechworld.it/boinf/2016/10/27/boxplots-with-ggplot2/

##### Other Graphs not Published ####

#copies per g - all samples
ggplot(dat, aes(x = Sample, y = FinalCopiesperg, fill = Pair)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Sample_Stage, scales = "free_x", labeller = label_wrap_gen(width = 40), ncol = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text=element_text(size=6),
      axis.title=element_text(size=10,face="bold"),
      legend.text=element_text(size=10), 
      strip.text = element_text(size = 8), legend.position = "none", 
      panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ylab("Final Copies (g-1)")+
  scale_fill_manual(values = colours)

ggsave("outputs/All_Samples_Copiesperg_barplot.pdf")


#COPY NO NORMALISED BY G PRODUCT ONLY #
ggplot(dat, aes(x = Sample, y = FinalCopiesperg, fill = Pair)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Sample_Stage, scales = "free_x", labeller = label_wrap_gen(width = 25), ncol = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ggtitle("Copies per g over Life Stages and Feed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=10), 
        legend.position = "bottom", strip.text = element_text(size = 8))
ggsave("Copiesperg_barplot.pdf")

ggplot(dat, aes(x = Sample_Stage, y = FinalCopiesperg, fill = Sample_Type)) +
  geom_boxplot() +
  facet_wrap(~Sample_Type, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  ggtitle("Copy Number Normalised By G Starting Material") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text=element_text(size=10),
      axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), 
      strip.text = element_text(size = 14))

ggsave("Copiesperg_boxplot_logaxis.pdf")


#### STATS #### 
library(multcomp)

####Is there an effect by life stage?
#firstly, comparisons should only be made between the life bstages, so remove feed
dat_LS <- dat %>%
  filter(Sample_Type == "Larvae")

#check normality of the data
#qqnorm(dat_LS$FinalCopiesperngperg); qqline(dat_LS$FinalCopiesperngperg)

qqnorm(dat_LS$FinalCopiesperg); qqline(dat_LS$FinalCopiesperg)
#data are right skewed. Will need to transform for normality.

#try log transformation
dat_LS$logFCng <- log(dat_LS$FinalCopiesperngperg)
qqnorm(dat_LS$logFCng); qqline(dat_LS$logFCng)

dat_LS$logFC <- log(dat_LS$FinalCopiesperg)
qqnorm(dat_LS$logFC); qqline(dat_LS$logFC)
#looks better if want to use linear regression, which assumes normality
#esp for small sample sizes

LS.lm <- lm(logFCng ~ Sample_Stage, data = dat_LS)
anova(LS.lm) # for raw data not sig and fan shape in residuals
#for transformed data sig when transformed:

#Analysis of Variance Table#

#Response: logFCng
#Df Sum Sq Mean Sq F value   Pr(>F)   
#Sample_Stage  4 38.966  9.7414  9.7342 0.001774 **
#  Residuals    10 10.007  1.0007                    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1                   

#summary(LS.lm) will compare everything to alphabetical baseline, not helpful

LS.lm <- lm(logFC ~ Sample_Stage, data = dat_LS)
anova(LS.lm) 

#Analysis of Variance Table

#Response: logFC
#Df Sum Sq Mean Sq F value    Pr(>F)    
#Sample_Stage  4 45.251 11.3127  16.354 0.0002188 ***
#  Residuals    10  6.917  0.6917                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Overall there is an effect of Life Stage both ng-1.g-1 and .g-1

#check plots to see if the assumptions of the model are met
plot(LS.lm) #qqplot looks okay, residuals look better than untransformed data
#for both normalised and unnormalised data

#post hoc testing - which comparisons are siginificant?
LScomp <- glht(LS.lm, linfct=mcp(Sample_Stage = "Tukey"))
summary(LScomp)

#Simultaneous Tests for General Linear Hypotheses

#Multiple Comparisons of Means: Tukey Contrasts


#Fit: lm(formula = logFC ~ Sample_Stage, data = dat_LS)

#Linear Hypotheses:
#  Estimate Std. Error t value Pr(>|t|)    
#L3-9 - L1 == 0         0.8530     0.6791   1.256 0.721500    
#L14-18 - L1 == 0       3.5480     0.6791   5.225 0.002730 ** 
#  L29-35 - L1 == 0       4.3829     0.6791   6.454 0.000521 ***
#  L49-53 - L1 == 0       3.7155     0.6791   5.471 0.002012 ** 
#  L14-18 - L3-9 == 0     2.6950     0.6791   3.969 0.017599 *  
#  L29-35 - L3-9 == 0     3.5299     0.6791   5.198 0.002854 ** 
#  L49-53 - L3-9 == 0     2.8625     0.6791   4.215 0.012091 *  
#  L29-35 - L14-18 == 0   0.8349     0.6791   1.229 0.735975    
#L49-53 - L14-18 == 0   0.1675     0.6791   0.247 0.999032    
#L49-53 - L29-35 == 0  -0.6674     0.6791  -0.983 0.857139    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#(Adjusted p values reported -- single-step method)

################################normalised by dna yield

#Fit: lm(formula = logFCng ~ Sample_Stage, data = dat_LS)

#Linear Hypotheses:
#  Estimate Std. Error t value Pr(>|t|)   
#L3-9 - L1 == 0         1.3484     0.8168   1.651  0.50123   
#L14-18 - L1 == 0       3.8969     0.8168   4.771  0.00521 **
#  L29-35 - L1 == 0       3.5697     0.8168   4.370  0.00962 **
#  L49-53 - L1 == 0       4.0298     0.8168   4.934  0.00418 **
#  L14-18 - L3-9 == 0     2.5485     0.8168   3.120  0.06514 . 
#L29-35 - L3-9 == 0     2.2213     0.8168   2.720  0.12028   
#L49-53 - L3-9 == 0     2.6815     0.8168   3.283  0.05060 . 
#L29-35 - L14-18 == 0  -0.3272     0.8168  -0.401  0.99369   
#L49-53 - L14-18 == 0   0.1329     0.8168   0.163  0.99981   
#L49-53 - L29-35 == 0   0.4601     0.8168   0.563  0.97760    


##is there a difference between feeds?
dat_feed <- dat %>%
  filter(Sample_Type == "Feed/Water")

#check normality of the data
#qqnorm(dat_feed$FinalCopiesperngperg); qqline(dat_feed$FinalCopiesperngperg)
qqnorm(dat_feed$FinalCopiesperg); qqline(dat_feed$FinalCopiesperg)
#data are also right skewed. Will need to transform for normality.

#try log transformation
dat_feed$logFC <- log(dat_feed$FinalCopiesperg)
qqnorm(dat_LS$logFC); qqline(dat_LS$logFC)
#looks better

feed.lm <- lm(logFC ~ Sample_Stage, data = dat_feed)
anova(feed.lm)

#Analysis of Variance Table

#Response: logFC
#Df  Sum Sq Mean Sq F value    Pr(>F)    
#Sample_Stage  4 223.662  55.916  92.152 7.543e-08 ***
#  Residuals    10   6.068   0.607                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Overall there a difference between feeds for both normalised and not norm by DNA yield

#check plots to see if the assumptions of the model are met
plot(feed.lm) #qplot looks okay, residuals look better than untransformed data

#post hoc testing - which comparisons are siginificant?
feedcomp <- glht(feed.lm, linfct=mcp(Sample_Stage = "Tukey"))
summary(feedcomp)

###### BIG NOTE HERE - WATER IS NOT IN g BUT IN mL
#I am not really sure how to deal with this besides not normalising by initial weight
#which seems a bit dodgey... 


#### Standard Curve Plot ####

# Plot standard curve

ggplot(stdcurve, aes(x = `Given Concentration`, y = Cq)) +
  geom_point(pch = 19, size = 3, alpha = 0.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_smooth(method='lm',formula = y ~ x, se = FALSE, color = "black") +
  annotate("text", x = max(stdcurve$`Given Concentration`)/10, y = max(stdcurve$Cq) - 1, hjust = 0, label = paste("Efficiency ",std_efficiency,"\nR squared ",std_rsq,"\ny =", std_slope,"x +" , std_intercept)) +
  scale_y_continuous(breaks = seq(round(min(stdcurve$Cq)), round(max(stdcurve$Cq)), 3)) +
  xlab('16S rRNA Gene Copy Number') +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

ggsave("outputs/Standard_Curve.pdf")
