library(survivalAnalysis)
library(survival)
library(RegParallel)
library(ggplot2)
library(survminer)
library(survcomp)
library(peperr)
library(SurvMetrics)
library(caret)
library(riskRegression)
library(forestplot)
library(forestploter)
library(pROC)
library(ROCit)
library(tidyverse)
library(tidytidbits)
require(RColorBrewer)
library(dplyr) 
library(readxl)
library(gridExtra)
library(circlize)
library(lattice)
library(grid)
library(waterfalls)
library(plotly)
library(knitr)
library(ComplexHeatmap)
library(cowplot)

#setwd("Your folder")
##Figure 1 The 21-gene signature in 62 primary osteosarcomas 
#data preparation
data<-read.table("Norm_count_table.txt",row.names = 1,header = TRUE,check.names = FALSE)
metadata <- read.table("metadata.txt",header = TRUE, row.names = 1)
# --check that both objects are aligned by name
rownames(metadata) %in% colnames(data)
# transform the expression data to Z score
x <- t(scale(t(data)))
# filter the Z-scores expression data to match the samples in our metadata
x <- x[,which(colnames(x) %in% rownames(metadata))]
# check that sample names match exactly between metadata and Z-scores 
all((colnames(x) == rownames(metadata)) == TRUE)
# create a merged metadata and Z-scores object
coxdata <- data.frame(metadata, t(x))
##Univariate Cox regression Model
fit_model<- RegParallel(
  data = coxdata,
  formula = 'Surv(OS_months, STATUS_OS) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE,x = T, y = T),
  FUNtype = 'coxph',
  variables = colnames(coxdata)[5:ncol(coxdata)], 
  blocksize = 770,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)
res <- fit_model[!is.na(fit_model$P),]
res
##### BUILD PROGNOSTIC RISK-SCORE genes with logRank<=0.05
###Compute for each sample PROGNOSTIC RISK-SCORE :risk score = ∑coefficient value(BETA) ∗ gene expression level.
##Got column B-coefficient and genes from res 
gene<- c('CD36','CCR3','ICAM2','IRF5','CD97','CCR6','TNFRSF18','ADA','MCAM','ETS1',
         'SPA17','ABCB1','CHUK','IL19','IL11RA','CD99','FAS','FLT3','C1QBP','ITGA1','FOXP3')
res1<-res[res$Variable %in% gene,]
res1<-res1[,c(2,3)]
coxdata1<-t(coxdata[names(coxdata) %in% gene])
##merge two dataframe
risk_score<-cbind.data.frame(res1,coxdata1)
mat_prod<-risk_score
for (i in 3:ncol(risk_score)){mat_prod[,i]=risk_score$Beta*risk_score[,i]}
mat_sum<-mat_prod
for (i in 3:ncol(mat_prod)){mat_sum[,i]=sum(mat_prod[i])}
matsum<-t(mat_sum[1,3:ncol(mat_sum)])  ##matsum=store risk-score value for each patients 
colnames(matsum)<-c("Risk_Score")
cutoff<-coxdata[,c('OS_months', 'STATUS_OS','EFS_months', 'STATUS_EFS')]
rownames(cutoff) %in% rownames(matsum)
cutoff$Risk_Score<-matsum[,c(1)]
##########################################################################################################################################
###########################Overall performance measures###########################
#Supplementary Figure S2 ROC to compute overall performance of univariate cox model in discovery cohorts
r <-roc(cutoff$STATUS_OS, cutoff$Risk_Score,ci=TRUE, boot.n=1000, ci.alpha=0.95, stratified=FALSE)
ppi <-600
png("Suppl_Fig.S2.png", width=14*ppi, height=8*ppi, res=ppi)
plot.roc(r, print.auc=TRUE, auc.polygon=FALSE, print.thres=TRUE,col=1)
dev.off() 
########################################
#############FOREST PLOT##############################
#Figure 1A Hazard ratio distribution of 21 genes across samples
model<-fit_model[fit_model$Variable %in% gene,]

model$se <- (model$HRupper - model$HRlower) / (2 * 1.96)
model$Weight <- 1/model$se^2
model$Weight <- model$Weight / sum(model$Weight)
pooled_effect <- round(sum(model$HR * model$Weight), 2)
model<-model[,c(1,6,10:12)]
model$P<-round(model$P,5)
model$HR<-round(model$HR,2)
model$HRlower<-round(model$HRlower,2)
model$HRupper<-round(model$HRupper,2)
png("Figure1A",width = 6,height =8.5,units = 'in',res = 300)
my_ticks <- c(.2,.3,.5, 1, 2,5)
attr(my_ticks, "labels") <-c(.2,.3,.5, 1, 2,5)
model |>
  mutate(mean=HR,
         lower = model$HRlower,
         upper = model$HRupper)%>%
  mutate(est = sprintf("%.2f [%.2f - %.2f]", 
                       mean, lower, upper), 
         .after = Variable)%>% 
  mutate(P= sprintf("%.3f", P),.after = est) |>
  forestplot(labeltext = c(Variable,est,P),
    xticks =my_ticks,
    xlog = TRUE,
    xlab = "HR (95% CI)",
    colgap = unit(3, "mm"),
    txt_gp=fpTxtGp(label= gpar(cex = 1)),
    shapes_gp = fpShapesGp(default = gpar(lwd = 1)),
    boxsize = 0.30,
    vertices = TRUE)|> 
  fp_add_header(Variable = c("Gene"),
                est = c("HR (95% CI)"),
                P=c("p-value")) |>
  fp_add_lines("black") |> 
  fp_set_style(box = "black",
               line = "black",
               axes = gpar(cex=1), 
               zero= gpar(col = "black", lty = "solid"))
dev.off()
###################################################################################################################################
##Figure 1B Histogram shows the distribution of risk scores, where patients classified as Low-R (Lo-R) are represented in yellow (n = 39),while patients classified as High-R (Hi-R) are represented in blue 
his<-read.table("Risk_Score_62Pts.txt",header = TRUE,check.names = FALSE)
at <- c(-15,  0, 15)
col <- ifelse(his$Risk_Score_Group== 'High', "dodgerblue", "darkgoldenrod2")
ppi <-600
png("Figure1B.png", width=10*ppi, height=6*ppi, res=ppi)
barplot(his$Risk_Score -2.4 ,yaxt = "n",col=col,border="black", space=0.0, ylim=c(-15,15))
axis(2, at = at, labels = at + 2.4)
dev.off()
##### Heatmap displays Z-score expression level of the overall survival 21-related genes
heat<-read.table("Exp_value_21genes.txt",row.names = 1,header = TRUE,check.names = FALSE)
heat=as.matrix(heat)
metadata_h<- read.table("Metadata_Heatmap.txt",header = TRUE,row.names = 1)
# --check that both objects are aligned by name
rownames(metadata_h) %in% colnames(heat)
##scale the data to Z-scores (by row)
heat <- t(scale(t(heat)))
##set colour scheme and choose breaks
hcl.pals("divergingx")
col_fun = colorRamp2(c(-2, 0, 2), hcl_palette = "Fall")
hmap <- Heatmap(heat,name = 'Gene\nZ-score',col = col_fun,
               # row (gene) parameters
                cluster_rows = TRUE,
                show_row_dend = FALSE,
                # column (sample) parameters
                cluster_columns = FALSE,
                show_column_dend = FALSE,
               show_column_names = TRUE,
               clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
               clustering_method_rows = 'ward.D2')
ppi <-600
png("Figure1B_heatmap.png", width=6*ppi, height=2*ppi, res=ppi)              
ht=draw(hmap) 
dev.off()
###################################################################################################
## Figure 1-C) Kaplan-Meier survival analysis of the risk for death of high-risk (Hi-R) patients compared to low-risk (Lo-R) patients, D) Event-free survival (EFS)
# set Risk_Score cut-offs for high and low risk group
highRisk <- 2.419
cutoff$Risk_Score <- ifelse(cutoff$Risk_Score >= highRisk , 'High','Low')
cutoff$Risk_Score<-factor(cutoff$Risk_Score)
# relevel the factors to have mid as the ref level
cutoff$Risk_Score <- factor(cutoff$Risk_Score, levels = c( 'High','Low' ))
##C-INDEX
Concordance<-coxph(Surv(OS_months, STATUS_OS) ~ Risk_Score,data = cutoff)
summary(Concordance)
##LogRank test Kaplan-Meier 
fit<-survfit(Surv(OS_months, STATUS_OS) ~ Risk_Score,data = cutoff)
#to estimate the probability of surviving to 5-year
##Percentage a 5-year overall survival
print(fit)
summary(fit,times = 60)
summary(fit)$table
##Perform log-rank test:Tests if there is a difference between two survival curves
surv_diff <- survdiff(Surv(OS_months, STATUS_OS) ~ Risk_Score,data = cutoff)
surv_diff
##Compute HR form surv_diff 
D1 <- surv_diff$obs[1]
D2 <- surv_diff$obs[2]
E1 <- surv_diff$exp[1]
E2 <- surv_diff$exp[2]
HR <- (D1/D2)/(E1/E2)
HR
#Standard errors and confidence interval 
SE_lnHR = sqrt(1/E1 + 1/E2)
SE_lnHR
L = log(HR)
lower <- exp(L - 1.96*SE_lnHR)
upper <- exp(L + 1.96*SE_lnHR)
ci95 <- c(lower=lower, upper=upper)
ci95
ggsurv<-ggsurvplot(fit, risk.table = TRUE, pval = FALSE,
                    ggtheme = theme_minimal(),
                    ylab = "Overall Survival",
                    xlab = "Follow-up Months",
                   break.x.by = 12,
                   size = 2,
                    mark.time=TRUE,
                   fontsize=4,
                    legend.labs = c("High Risk", "Low Risk"),
                    font.legend=c(12,"plain"),
                  legend.title="Risk score",
                    #surv.median.line = "hv",
                    risk.table.y.text.col = TRUE,
                    risk.table.y.text = FALSE,
                    palette = c("dodgerblue" , "darkgoldenrod2"))
 ggsurv$plot <- ggsurv$plot+ ggtitle("Training cohort")+ theme(
   plot.title = element_text(size=14, face="bold.italic"))+
   ggplot2::annotate("text", 
                     x = 15, y = 0.15, # x and y coordinates of the text
                     label = "p= 3e-06  \n HR= 10.60 (3.77-29.73) ", size = 4) 
ppi <-600
png("Figure1C.png", width=14*ppi, height=8*ppi, res=ppi)
ggsurv 
dev.off()
#####################################################
### Figure 1-D) Event-free survival (EFS)
fit<-survfit(Surv(EFS_months, STATUS_EFS) ~ Risk_Score,data = cutoff)
#to estimate the probability of surviving to 5-year
print(fit)
summary(fit,times = 60)
summary(fit)$table
## log-rank test:Tests if there is a difference between two survival curves
surv_diff <- survdiff(Surv(EFS_months, STATUS_EFS) ~ Risk_Score,data = cutoff)
surv_diff
##Compute HR form surv_diff 
D1 <- surv_diff$obs[1]
D2 <- surv_diff$obs[2]
E1 <- surv_diff$exp[1]
E2 <- surv_diff$exp[2]
HR <- (D1/D2)/(E1/E2)
HR
#Standard errors and confidence interval 
SE_lnHR = sqrt(1/E1 + 1/E2)
SE_lnHR
L = log(HR)
lower <- exp(L - 1.96*SE_lnHR)
upper <- exp(L + 1.96*SE_lnHR)
ci95 <- c(lower=lower, upper=upper)
ci95
ggsurv<-ggsurvplot(fit, risk.table = TRUE, pval = FALSE,
                    ggtheme = theme_minimal(),
                    ylab = "Event Free Survival",
                    xlab = "Follow-up Months",
                    break.x.by = 12,
                    size = 2,
                    mark.time=TRUE,
                    fontsize=4,
                    legend.labs = c("High Risk", "Low Risk"),
                    font.legend=c(12,"plain"),
                    legend.title="Risk score",
                    #surv.median.line = "hv",
                    risk.table.y.text.col = TRUE,
                    risk.table.y.text = FALSE,
                    palette = c("dodgerblue" , "darkgoldenrod2"))
ggsurv$plot <- ggsurv$plot+ ggtitle("Training cohort")+ theme(
   plot.title = element_text( size=14, face="bold.italic"))+
   ggplot2::annotate("text", 
                     x = 15, y = 0.15, # x and y coordinates of the text
                     label = "p= 3e-06  \n HR= 6.18 (2.54-15.03) ", size = 4)                  
ppi <-600
png("Figure1D.png", width=12*ppi, height=8*ppi, res=ppi)
ggsurv 
dev.off()       
###################################################################################################################################################
###############################################################################################
### Table 2 MULTIVARIATE ANALYSIS 
multv_data<-data.frame(metadata[,1:9],cutoff[5])
covariate_names <- c(SEX='SEX', Age="Age" ,Serum_alkaline_phosphatase="Serum_alkaline_phosphatase",HistologicResponse="HistologicResponse",Risk_Score="Risk_Score")

result<-multv_data %>%
  mutate(SEX=rename_factor(SEX, `1` = "M", `2` = "F"), HistologicResponse=rename_factor(HistologicResponse, `1` = "GR", `2` = "PR"),
         Serum_alkaline_phosphatase=rename_factor(Serum_alkaline_phosphatase, `1` = "Normal", `2` = "High"),
        Risk_Score=rename_factor(Risk_Score, `1` = "High", `2` = "Low")) %>%
  analyse_multivariate(vars(OS_months, STATUS_OS),
                       covariates = vars(SEX,Age,HistologicResponse,Risk_Score,Serum_alkaline_phosphatase),
                       covariate_name_dict = covariate_names)


result
#############################################################################################################################################
###############################################################################################################################################
##Figure 2: Evaluation of the 21-gene signature in the validation cohorts. 
## Figure 2A The OS and EFS in the TCGA validation cohort
TCGA_data<-read.table("TARGET_OS_Exp_value.txt",row.names = 1,header = TRUE,check.names = FALSE)
metadata_TARGET <- read.table("metadata_TARGET.txt",header = TRUE, row.names = 1)
# --check that both objects are aligned by name
rownames(metadata_TARGET) %in% colnames(TCGA_data)
# transform the expression data to Z score
x <- t(scale(t(TCGA_data)))
# filter the Z-scores expression data to match the samples in our pdata
x <- x[,which(colnames(x) %in% rownames(metadata_TARGET))]
# check that sample names match exactly between metadata and Z-scores 
all((colnames(x) == rownames(metadata_TARGET)) == TRUE)
# create a merged metadata and Z-scores object
coxdata <- data.frame(metadata_TARGET, t(x))
# tidy column names
colnames(coxdata)[1:4] <- c('EFS_Status', 'EFS_months','OS_Status','OS_months')
##Univariate Cox regression
fit_model<- RegParallel(
  data = coxdata,
  formula = 'Surv(OS_months, OS_Status) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE,x = T, y = T),
  FUNtype = 'coxph',
  variables = colnames(coxdata)[5:ncol(coxdata)], 
  blocksize = 20,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)
res <- fit_model[!is.na(fit_model$P),]
##### Compute PROGNOSTIC RISK-SCORE 
gene<- c('CD36','CCR3','ICAM2','IRF5','CD97','CCR6','TNFRSF18','ADA','MCAM','ETS1',
         'SPA17','ABCB1','CHUK','IL19','IL11RA','CD99','FAS','FLT3','C1QBP','ITGA1','FOXP3')
res1<-res[res$Variable %in% gene,]
res1<-res1[,c(2,3)]
coxdata1<-t(coxdata[names(coxdata) %in% gene])
risk_score<-cbind.data.frame(res1,coxdata1)
mat_prod<-risk_score
for (i in 3:ncol(risk_score)){mat_prod[,i]=risk_score$Beta*risk_score[,i]}
mat_sum<-mat_prod
for (i in 3:ncol(mat_prod)){mat_sum[,i]=sum(mat_prod[i])}
matsum<-t(mat_sum[1,3:ncol(mat_sum)])  
colnames(matsum)<-c("Risk_Score")
cutoff<-coxdata[,c('EFS_Status', 'EFS_months','OS_Status','OS_months')]
rownames(cutoff) %in% rownames(matsum)
cutoff$Risk_Score<-matsum[,c(1)]
#set Risk_Score cut-offs for high and low risk group
highRisk <- 2.419
cutoff$Risk_Score <- ifelse(cutoff$Risk_Score >= highRisk , 'High','Low')
cutoff$Risk_Score<-factor(cutoff$Risk_Score)
cutoff$Risk_Score <- factor(cutoff$Risk_Score, levels = c( 'High','Low' ))
##LogRank test Kaplan-Meier 
fit<-survfit(Surv(OS_months, OS_Status) ~ Risk_Score,data = cutoff)
surv_diff <- survdiff(Surv(OS_months, OS_Status) ~ Risk_Score,data = cutoff)
surv_diff
##Compute HR form surv_diff 
D1 <- surv_diff$obs[1]
D2 <- surv_diff$obs[2]
E1 <- surv_diff$exp[1]
E2 <- surv_diff$exp[2]
HR <- (D1/D2)/(E1/E2)
HR
#Standard errors and confidence interval 
SE_lnHR = sqrt(1/E1 + 1/E2)
SE_lnHR
L = log(HR)
lower <- exp(L - 1.96*SE_lnHR)
upper <- exp(L + 1.96*SE_lnHR)
ci95 <- c(lower=lower, upper=upper)
ci95
ggsurv<-ggsurvplot(fit, risk.table = TRUE, pval = FALSE,
                   ggtheme = theme_minimal(),
                   ylab = "Overall Survival",
                   xlab = "Follow-up Months",
                   break.x.by = 12,
                   size = 2,
                   mark.time=TRUE,
                   fontsize=4,
                   legend.labs = c("High Risk", "Low Risk"),
                   font.legend=c(12,"plain"),
                   legend.title="Risk score",
                   #surv.median.line = "hv",
                   risk.table.y.text.col = TRUE,
                   risk.table.y.text = FALSE,
                   palette = c("dodgerblue" , "darkgoldenrod2"))

ggsurv$plot <- ggsurv$plot+ ggtitle("Target-OS TCGA")+ theme(
  plot.title = element_text( size=14, face="bold.italic"))+
  ggplot2::annotate("text", 
                    x = 15, y = 0.15, # x and y coordinates of the text
                    label = "p= 9e-07   \n HR= 8.68 (1.85-40.59) ", size = 4)                  

ppi <-600
png("Figure2A.png", width=14*ppi, height=8*ppi, res=ppi)
ggsurv 
dev.off()      
#####################################################
###Figure 2B COMPUTE EFS
fit<-survfit(Surv(EFS_months, EFS_Status) ~ Risk_Score,data = cutoff)
##Perform log-rank test
surv_diff <- survdiff(Surv(EFS_months, EFS_Status) ~ Risk_Score,data = cutoff)
surv_diff
##Compute HR form surv_diff 
D1 <- surv_diff$obs[1]
D2 <- surv_diff$obs[2]
E1 <- surv_diff$exp[1]
E2 <- surv_diff$exp[2]
HR <- (D1/D2)/(E1/E2)
HR
#Standard errors and confidence interval 
SE_lnHR = sqrt(1/E1 + 1/E2)
SE_lnHR
L = log(HR)
lower <- exp(L - 1.96*SE_lnHR)
upper <- exp(L + 1.96*SE_lnHR)
ci95 <- c(lower=lower, upper=upper)
ci95
ggsurv<-ggsurvplot(fit, risk.table = TRUE, pval = FALSE,
                   ggtheme = theme_minimal(),
                   ylab = "Event Free Survival",
                   xlab = "Follow-up Months",
                   break.x.by = 12,
                   size = 2,
                   mark.time=TRUE,
                   fontsize=4,
                   legend.labs = c("High Risk", "Low Risk"),
                   font.legend=c(12,"plain"),
                   legend.title="Risk score",
                   #surv.median.line = "hv",
                   risk.table.y.text.col = TRUE,
                   risk.table.y.text = FALSE,
                   palette = c("dodgerblue" , "darkgoldenrod2"))

ggsurv$plot <- ggsurv$plot+ ggtitle("Target-OS TCGA")+ theme(
  plot.title = element_text( size=14, face="bold.italic"))+
  ggplot2::annotate("text", 
                    x = 8, y = 0.15, # x and y coordinates of the text
                    label = "p= 7e-04 \n HR= 3.59(1.14-11.31) ", size = 4)                  

ppi <-600
png("Figure2B.png", width=12*ppi, height=8*ppi, res=ppi)
ggsurv 
dev.off()  
#######################################################################################################     
##Figure 2: C-D) The OS and EFS in the GSE33382 validation cohort
GSE33382_data<-read.table("GSE33382_Exp_value.txt",row.names = 1,header = TRUE,check.names = FALSE)
GSE33382_metadata<- read.table("GSE33382_metadata.txt",header = TRUE, row.names = 1)
# --check that both objects are aligned by name
rownames(GSE33382_metadata) %in% colnames(GSE33382_data)
# transform the expression data to Z score
x <- t(scale(t(GSE33382_data)))
# filter the Z-scores expression data to match the samples in our pdata
x <- x[,which(colnames(x) %in% rownames(GSE33382_metadata))]
# check that sample names match exactly between metadata and Z-scores 
all((colnames(x) == rownames(GSE33382_metadata)) == TRUE)
# create a merged metadata and Z-scores object
coxdata <- data.frame(GSE33382_metadata, t(x))
# tidy column names
colnames(coxdata)[1:4] <- c('EFS_Status', 'EFS_months','OS_Status','OS_months')
##Univariate Cox regression
fit_model<- RegParallel(
  data = coxdata,
  formula = 'Surv(OS_months, OS_Status) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE,x = T, y = T),
  FUNtype = 'coxph',
  variables = colnames(coxdata)[5:ncol(coxdata)], 
  blocksize = 20,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)
res <- fit_model[!is.na(fit_model$P),]
res
##### Compute PROGNOSTIC RISK-SCORE 
gene<- c('CD36','CCR3','ICAM2','IRF5','CD97','CCR6','TNFRSF18','ADA','MCAM','ETS1',
         'SPA17','ABCB1','CHUK','IL19','IL11RA','CD99','FAS','FLT3','C1QBP','ITGA1','FOXP3')
res1<-res[res$Variable %in% gene,]
res1<-res1[,c(2,3)]
coxdata1<-t(coxdata[names(coxdata) %in% gene])
risk_score<-cbind.data.frame(res1,coxdata1)
mat_prod<-risk_score
for (i in 3:ncol(risk_score)){mat_prod[,i]=risk_score$Beta*risk_score[,i]}
mat_sum<-mat_prod
for (i in 3:ncol(mat_prod)){mat_sum[,i]=sum(mat_prod[i])}
matsum<-t(mat_sum[1,3:ncol(mat_sum)])  
colnames(matsum)<-c("Risk_Score")
cutoff<-coxdata[,c('EFS_Status', 'EFS_months','OS_Status','OS_months')]
rownames(cutoff) %in% rownames(matsum)
cutoff$Risk_Score<-matsum[,c(1)]
# set Risk_Score cut-offs for high and low risk group
highRisk <- 2.4
cutoff$Risk_Score <- ifelse(cutoff$Risk_Score >= highRisk , 'High','Low')
cutoff$Risk_Score<-factor(cutoff$Risk_Score)
# relevel the factors to have mid as the ref level
cutoff$Risk_Score <- factor(cutoff$Risk_Score, levels = c( 'High','Low' ))
##LogRank test Kaplan-Meier 
fit<-survfit(Surv(OS_months, OS_Status) ~ Risk_Score,data = cutoff)
##Perform log-rank test
surv_diff <- survdiff(Surv(OS_months, OS_Status) ~ Risk_Score,data = cutoff)
surv_diff
##Compute HR form surv_diff 
D1 <- surv_diff$obs[1]
D2 <- surv_diff$obs[2]
E1 <- surv_diff$exp[1]
E2 <- surv_diff$exp[2]
HR <- (D1/D2)/(E1/E2)
HR
#Standard errors and confidence interval 
SE_lnHR = sqrt(1/E1 + 1/E2)
SE_lnHR
L = log(HR)
lower <- exp(L - 1.96*SE_lnHR)
upper <- exp(L + 1.96*SE_lnHR)
ci95 <- c(lower=lower, upper=upper)
ci95
ggsurv<-ggsurvplot(fit, risk.table = TRUE, pval = FALSE,
                   ggtheme = theme_minimal(),
                   ylab = "Overall Survival",
                   xlab = "Follow-up Months",
                   break.x.by = 12,
                   size = 2,
                   mark.time=TRUE,
                   fontsize=4,
                   legend.labs = c("High Risk", "Low Risk"),
                   font.legend=c(12,"plain"),
                   legend.title="Risk score",
                   #surv.median.line = "hv",
                   risk.table.y.text.col = TRUE,
                   risk.table.y.text = FALSE,
                   palette = c("dodgerblue" , "darkgoldenrod2"))

ggsurv$plot <- ggsurv$plot+ ggtitle("GSE33382")+ theme(
  plot.title = element_text( size=14, face="bold.italic"))+
  ggplot2::annotate("text", 
                    x = 15, y = 0.15, # x and y coordinates of the text
                    label = "p= 3e-04 \n HR= 18.63 (1.88-193.70) ", size = 4)                  

ppi <-600
png("Figure2C.png", width=14*ppi, height=8*ppi, res=ppi)
ggsurv 
dev.off()     
#####################################################
###Figure 2D) EFS
##LogRank test Kaplan-Meier 
fit<-survfit(Surv(EFS_months, EFS_Status) ~ Risk_Score,data = cutoff)
print(fit)
summary(fit)$table
##Perform log-rank test
surv_diff <- survdiff(Surv(EFS_months, EFS_Status) ~ Risk_Score,data = cutoff)
surv_diff
##Compute HR form surv_diff 
D1 <- surv_diff$obs[1]
D2 <- surv_diff$obs[2]
E1 <- surv_diff$exp[1]
E2 <- surv_diff$exp[2]
HR <- (D1/D2)/(E1/E2)
HR
#Standard errors and confidence interval 
SE_lnHR = sqrt(1/E1 + 1/E2)
SE_lnHR
L = log(HR)
lower <- exp(L - 1.96*SE_lnHR)
upper <- exp(L + 1.96*SE_lnHR)
ci95 <- c(lower=lower, upper=upper)
ci95
ggsurv<-ggsurvplot(fit, risk.table = TRUE, pval = FALSE,
                   ggtheme = theme_minimal(),
                   ylab = "Event Free Survival",
                   xlab = "Follow-up Months",
                   break.x.by = 12,
                   size = 2,
                   mark.time=TRUE,
                   fontsize=4,
                   legend.labs = c("High Risk", "Low Risk"),
                   font.legend=c(12,"plain"),
                   legend.title="Risk score",
                   #surv.median.line = "hv",
                   risk.table.y.text.col = TRUE,
                   risk.table.y.text = FALSE,
                   palette = c("dodgerblue" , "darkgoldenrod2"))

ggsurv$plot <- ggsurv$plot+ ggtitle("GSE33382")+ theme(
  plot.title = element_text( size=14, face="bold.italic"))+
  ggplot2::annotate("text", 
                    x = 8, y = 0.15, # x and y coordinates of the text
                    label = "p= 0.004   \n HR= 6.64 (0.91-48.33) ", size = 4)                  

ppi <-600
png("Figure2D.png", width=14*ppi, height=8*ppi, res=ppi)
ggsurv 
dev.off()     

