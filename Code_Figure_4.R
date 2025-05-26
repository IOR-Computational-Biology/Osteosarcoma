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
### Figure 4 The 31-gene mifamurtide-related signature in a subgroup of 33 Pgp-positive primary 
#Osteosarcomas undergoing adjuvant mifamurtide treatment combined with chemotherapy
Pgp_positive<-read.table("33_Pgp+_Norm_count.txt",row.names = 1,header = TRUE,check.names = FALSE)
metadata_Pgp<- read.table("metadata_33_PgpPos.txt",header = TRUE, row.names = 1)
# --check that both objects are aligned by name
rownames(metadata_Pgp) %in% colnames(Pgp_positive)
# transform the expression data to Z score
x <- t(scale(t(Pgp_positive)))
# filter the Z-scores expression data to match the samples in our pdata
x <- x[,which(colnames(x) %in% rownames(metadata_Pgp))]
# check that sample names match exactly between metadata and Z-scores 
all((colnames(x) == rownames(metadata_Pgp)) == TRUE)
# create a merged metadata and Z-scores object
coxdata <- data.frame(metadata_Pgp, t(x))
##Univariate Cox regression model
fit_model<- RegParallel(
  data = coxdata,
  formula = 'Surv(OS_months, STATUS_OS) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE,x = T, y = T),
  FUNtype = 'coxph',
  variables = colnames(coxdata)[6:ncol(coxdata)], 
  blocksize = 770,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)

res <- fit_model[!is.na(fit_model$P),]
res
#####################################################################################
##### PROGNOSTIC RISK-SCORE logRank<0.05 31_gene_signature
gene<- c('CXCR2','CD36','CCR3','TXNIP','IRF5','MICB','S100A8','FPR2','CMA1','BST1','SAA1',
         'S100A12','RRAD','CD97','NRP1','IL1R2','HLA.DQB1','RORC','C7','TNFRSF10C','VCAM1',
         'CXCL10','CXCR1','STAT6','CXCR6','C6',
         'LILRA1','CCL8','RUNX1','A2M','ABCB1')
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
cutoff<-coxdata[,c('OS_months', 'STATUS_OS','EFS_months','STATUS_EFS')]
rownames(cutoff) %in% rownames(matsum)
cutoff$Risk_Score<-matsum[,c(1)]
################################################################################
###########################Overall performance measures###########################
#Supplementary Figure-S5 ROC Curve to compute overall performance of univariate cox model in 33 out of 60 samples undergoing adjuvant mifamurtide treatment
r <-roc(cutoff$STATUS_OS, cutoff$Risk_Score,ci=TRUE, boot.n=1000, ci.alpha=0.95, stratified=FALSE)
ppi <-600
png("Suppl_Fig.S4.png", width=14*ppi, height=8*ppi, res=ppi)
plot.roc(r, print.auc=TRUE, auc.polygon=FALSE, print.thres=TRUE,col=1)
dev.off() 
###########################################
##Figure 4A Forest plot showing genes significantly associated with prognosis
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
png("Figure 4A.png",width = 7,height =12,units = 'in',res = 300)
my_ticks <- c(.1,.3,.5, 1, 2,7)
attr(my_ticks, "labels") <-c(.1,.3,.5, 1, 2,7)
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
######################################################################
###Figure 4B
histog<-read.table("risk_Score_33Pgp+.txt",header = TRUE,check.names = FALSE)
at <- c(-18,  0, 40)
col <- ifelse(histog$Risk_Score_Group== 'High', "dodgerblue", "darkgoldenrod2")
ppi <-600
png("Figure4B.png", width=10*ppi, height=6*ppi, res=ppi)
barplot(histog$Risk_Score -8.69 ,yaxt = "n",col=col,border="black", space=0.0, ylim=c(-18,40))
axis(2, at = at, labels = at + 8.69)
dev.off()
########HEATMAP displays Z-score expression level of the overall survival 31-related genes
exp_value_31genes<-read.table("Exp_value_31genes.txt",row.names = 1,header = TRUE,check.names = FALSE)
exp_value_31genes=as.matrix(exp_value_31genes)
metadata_exp_value_31genes<- read.table("Metadata_Heat_33Pgp+.txt",header = TRUE,row.names = 1)
# --check that both objects are aligned by name
rownames(metadata_exp_value_31genes) %in% colnames(exp_value_31genes)
##scale the data to Z-scores (by row)
heat <- t(scale(t(exp_value_31genes)))
##set colour scheme and choose breaks
hcl.pals("divergingx")
col_fun = colorRamp2(c(-2, 0, 2), hcl_palette = "Fall")

hmap <- Heatmap(heat,name = 'Gene\nZ-score',col = col_fun,
               # row (gene) parameters
                cluster_rows = TRUE,
                show_row_dend = FALSE,
               show_row_names = TRUE,
                # column (sample) parameters
                cluster_columns = FALSE,
                show_column_dend = FALSE,
               show_column_names = FALSE,
                clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                clustering_method_rows = 'ward.D2')
ppi <-600
png("Figure_4B.png", width=6*ppi, height=2*ppi, res=ppi)              
ht=draw(hmap) 
dev.off()
  ############################################################################################################################
### Figure 4-C) Kaplan-Meier survival analysis of the risk for death of high-risk (Hi-R) patients compared to low-risk (Lo-R) patients in 33Pgp+ patients, D) Event-free survival (EFS)
# set Risk_Score cut-offs for high and low risk group
highRisk <-8.69
cutoff$Risk_Score <- ifelse(cutoff$Risk_Score >= highRisk , 'High','Low')
cutoff$Risk_Score<-factor(cutoff$Risk_Score)
# relevel the factors to have mid as the ref level
cutoff$Risk_Score <- factor(cutoff$Risk_Score, levels = c( 'High','Low' ))
##LogRank test Kaplan-Meier 
fit<-survfit(Surv(OS_months, STATUS_OS) ~ Risk_Score,data = cutoff)
#To estimate the probability of surviving to 5-year
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

ggsurv$plot <- ggsurv$plot+ ggtitle("PgP+ cohort")+ theme(
  plot.title = element_text(size=14, face="bold.italic"))+ 
  ggplot2::annotate("text", 
                    x = 8, y = 0.15, # x and y coordinates of the text
                    label = "p=1e-09\n HR= 24.51 (3.39-176.84) ", size = 4)                  

ppi <-600
png("Fig4C.png", width=14*ppi, height=8*ppi, res=ppi)
ggsurv 
dev.off()
#####################################################
#### Figure 4D Event-free survival (EFS)
##LogRank test Kaplan-Meier 
fit<-survfit(Surv(EFS_months, STATUS_EFS) ~ Risk_Score,data = cutoff)
#to estimate the probability of surviving to 5year
print(fit)
summary(fit,times = 60)
summary(fit)$table
##Perform log-rank test
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

ggsurv$plot <- ggsurv$plot+ ggtitle("PgP+ cohort")+ theme(
  plot.title = element_text(size=14, face="bold.italic"))+ 
  ggplot2::annotate("text", 
                    x = 8, y = 0.15, # x and y coordinates of the text
                    label = "p= 3e-06  \n HR= 7.80 (1.49-40.68) ", size = 4)                  

ppi <-600
png("Fig4D.png", width=14*ppi, height=8*ppi, res=ppi)
ggsurv 
dev.off()

###############################################################################################
###MULTIVARIATE ANALYSIS#######################################################################################
multv_data_Pgp<-data.frame(metadata_Pgp[,1:9],cutoff[5])
covariate_names <- c(SEX='SEX', Age="Age" ,Serum_alkaline_phosphatase="Serum_alkaline_phosphatase",HistologicResponse="HistologicResponse",Risk_Score="Risk_Score")

result<-multv_data_Pgp %>%
  mutate(SEX=rename_factor(SEX, `1` = "M", `2` = "F"), HistologicResponse=rename_factor(HistologicResponse, `1` = "GR", `2` = "PR"),
         Serum_alkaline_phosphatase=rename_factor(Serum_alkaline_phosphatase, `1` = "Normal", `2` = "High"),
        Risk_Score=rename_factor(Risk_Score, `1` = "High", `2` = "Low")) %>%
  analyse_multivariate(vars(OS_months, STATUS_OS),
                       covariates = vars(SEX,Age,HistologicResponse,Risk_Score,Serum_alkaline_phosphatase),
                       covariate_name_dict = covariate_names)
                  
 

result
write.table(result[["summaryAsFrame"]] , "Table3.txt", quote = FALSE,sep="\t")
##########################################################################################################################################























