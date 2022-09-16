suppressMessages({library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.mRNA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)})
clin <- survivalTCGA(BRCA.clinical,extract.cols="admin.disease_code")
expr_with_vital <- BRCA.mRNA %>% 
  as_tibble() %>% 
  select(bcr_patient_barcode, WDR48, PTEN) %>% 
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>% 
  inner_join(clin, by="bcr_patient_barcode") 
expr <- expressionsTCGA(BRCA.mRNA, extract.cols = c("WDR48", "PTEN"))
WDR48_expr <- ggplot(expr, aes(dataset, WDR48, fill=dataset)) + geom_boxplot()
PTEN_expr <- ggplot(expr, aes(dataset, PTEN, fill=dataset)) + geom_boxplot()
WDR48_only <-BRCA.mRNA %>% #590 observations
  as_tibble() %>% 
  select(bcr_patient_barcode, WDR48) %>% 
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>% 
  inner_join(clin, by="bcr_patient_barcode") 
PTEN_only <-BRCA.mRNA %>% #590 observations
  as_tibble() %>% 
  select(bcr_patient_barcode, PTEN) %>% 
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>% 
  inner_join(clin, by="bcr_patient_barcode") 
AKT_only <-BRCA.mRNA %>% #590 observations
  as_tibble() %>% 
  select(bcr_patient_barcode, AKT1) %>% 
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>% 
  inner_join(clin, by="bcr_patient_barcode")
GAPDH_only <-BRCA.mRNA %>% #590 observations
  as_tibble() %>% 
  select(bcr_patient_barcode, GAPDH) %>% 
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>% 
  inner_join(clin, by="bcr_patient_barcode")
TP53_only <- BRCA.mRNA %>% #590 observations
  as_tibble() %>% 
  select(bcr_patient_barcode, TP53) %>% 
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>% 
  inner_join(clin, by="bcr_patient_barcode")
##both expressions
both_genes <-BRCA.mRNA %>% #590 observations
  as_tibble() %>% 
  select(bcr_patient_barcode, PTEN, WDR48) %>% 
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>% 
  inner_join(clin, by="bcr_patient_barcode") 

akt_pten <-BRCA.mRNA %>% #590 observations
  as_tibble() %>% 
  select(bcr_patient_barcode, AKT1, PTEN) %>% 
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>% 
  inner_join(clin, by="bcr_patient_barcode") 
akt_wdr <-BRCA.mRNA %>% #590 observations
  as_tibble() %>% 
  select(bcr_patient_barcode, AKT1, WDR48) %>% 
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>% 
  inner_join(clin, by="bcr_patient_barcode") 

WDR48_only$WDR48 <- WDR48_only$WDR48/GAPDH_only$GAPDH
PTEN_only$PTEN <- PTEN_only$PTEN/GAPDH_only$GAPDH
AKT_only$AKT1 <- AKT_only$AKT1/GAPDH_only$GAPDH
akt_pten$PTEN <- akt_pten$PTEN/GAPDH_only$GAPDH
akt_pten$AKT1 <- akt_pten$AKT1/GAPDH_only$GAPDH
akt_wdr$AKT1 <- akt_wdr$AKT1/GAPDH_only$GAPDH
akt_wdr$WDR48 <- akt_wdr$WDR48/GAPDH_only$GAPDH
both_genes$PTEN <- both_genes$PTEN/GAPDH_only$GAPDH
both_genes$WDR48 <- both_genes$WDR48/GAPDH_only$GAPDH

l <- length(PTEN_only$PTEN) #length = 590 for both
bounds <- c(round(l/4), round(l-l/4))
#top5bounds <- c(round(590/20), round(590-590/20))
pctlsWDR48 <- c(sort(WDR48_only$WDR48)[bounds[1]], sort(WDR48_only$WDR48)[bounds[2]])
pctlsPTEN <- c(sort(PTEN_only$PTEN)[bounds[1]], sort(PTEN_only$PTEN)[bounds[2]])
pctlsAKT <- c(sort(AKT_only$AKT1)[bounds[1]], sort(AKT_only$AKT1)[bounds[2]])
pctlsTP <- c(sort(TP53_only$TP53)[bounds[1]], sort(TP53_only$TP53)[bounds[2]])


#25 vs 75 implementation
WDR48_only <- WDR48_only[with(WDR48_only, order(-WDR48)),][-c(bounds[1]: bounds[2]),]
PTEN_only <- PTEN_only[with(PTEN_only, order(-PTEN)),][-c(bounds[1]: bounds[2]),]
AKT_only <- AKT_only[with(AKT_only, order(-AKT1)),][-c(bounds[1]: bounds[2]),]
WDR48_only$WDR48 <- ifelse(WDR48_only$WDR48 < pctlsWDR48[1],0, 2)
PTEN_only$PTEN <- ifelse(PTEN_only$PTEN < pctlsPTEN[1],0, 2)
AKT_only$AKT1 <- ifelse(AKT_only$AKT1 < pctlsAKT[1], 0, 2)


WDR48_only$WDR48 <- ifelse(WDR48_only$WDR48 < median(WDR48_only$WDR48),0, 2)
PTEN_only$PTEN <- ifelse(PTEN_only$PTEN < median(PTEN_only$PTEN),0, 2)
AKT_only$AKT1 <- ifelse(AKT_only$AKT1 < median(AKT_only$AKT1), 0, 2)
WDR48_only$exprs <- ifelse(WDR48_only$WDR48 == 2 , "WDR48 High", "WDR48 Low")
PTEN_only$exprs <- ifelse(PTEN_only$PTEN == 2 , "PTEN High", "PTEN Low")
AKT_only$exprs <- ifelse(AKT_only$AKT1 == 2 , "AKT High", "AKT Low")

#median
both_genes$PTEN <- ifelse(both_genes$PTEN < median(both_genes$PTEN), 0, 1)
both_genes$WDR48 <- ifelse(both_genes$WDR48 < median(both_genes$WDR48), 0, 1)
both_genes <- both_genes[(both_genes$PTEN == 1 & both_genes$WDR48 == 1) | (both_genes$PTEN == 0 & both_genes$WDR48 == 0),]
both_genes$exprs<- ifelse(both_genes$PTEN == 1 & both_genes$WDR48 == 1, "Both High", "Both Low")


##25-75 
WDR48_only$WDR48 <- ifelse(WDR48_only$WDR48  < pctlsWDR48[1],0,ifelse(WDR48_only$WDR48  > pctlsWDR48[2], 3,
                                                                     ifelse((WDR48_only$WDR48 < median(WDR48_only$WDR48 ) & WDR48_only$WDR48 > pctlsWDR48[1]), 1,2)))
PTEN_only$PTEN <- ifelse(PTEN_only$PTEN < pctlsPTEN[1],0,ifelse(PTEN_only$PTEN> pctlsPTEN[2], 3,
                                                                  ifelse((PTEN_only$PTEN< median(PTEN_only$PTEN) & PTEN_only$PTEN> pctlsPTEN[1]), 1,2)))
AKT_only$AKT1 <- ifelse(AKT_only$AKT1 < pctlsAKT[1],0,ifelse(AKT_only$AKT1> pctlsAKT[2], 3,
                                                                ifelse((AKT_only$AKT1< median(AKT_only$AKT1) & AKT_only$AKT1> pctlsAKT[1]), 1,2)))
PTEN_only <- PTEN_only[PTEN_only$PTEN == 3 | PTEN_only$PTEN == 0,]
WDR48_only <- WDR48_only[WDR48_only$WDR48 == 3 | WDR48_only$WDR48 == 0,]
WDR48_only$exprs<- ifelse(WDR48_only$WDR48 ==3, "WDR48 High", "WDR48 Low")
PTEN_only$exprs<- ifelse(PTEN_only$PTEN ==3, "PTEN High", "PTEN Low")
AKT_only$exprs<- ifelse(AKT_only$AKT1 ==3, "AKT1 High", "AKT1 Low")
TP53_only$TP53 <- ifelse(TP53_only$TP53 < pctlsTP[1],0,ifelse(TP53_only$TP53> pctlsTP[2], 3,
                                                             ifelse((TP53_only$TP53< median(TP53_only$TP53) & TP53_only$TP53> pctlsTP[1]), 1,2)))
TP53_only <- TP53_only[TP53_only$TP53 == 3| TP53_only$TP53 == 0,]
TP53_only$exprs<- ifelse(TP53_only$TP53 ==3, "TP53 High", "TP53 Low")
##name all quartiles
#quartile 0,1,2,3
both_genes$WDR48 <- ifelse(both_genes$WDR48 < pctlsWDR48[1],0,ifelse(both_genes$WDR48 > pctlsWDR48[2], 3,
                      ifelse((both_genes$WDR48< median(both_genes$WDR48) & both_genes$WDR48> pctlsWDR48[1]), 1,2)))
both_genes$PTEN <- ifelse(both_genes$PTEN < pctlsPTEN[1],0,ifelse(both_genes$PTEN > pctlsPTEN[2], 3,
                      ifelse((both_genes$PTEN< median(both_genes$PTEN) & both_genes$PTEN> pctlsPTEN[1]), 1,2)))
akt_pten$AKT1 <- ifelse(akt_pten$AKT1  < pctlsAKT[1],0,ifelse(akt_pten$AKT1  > pctlsAKT[2], 3,
                                                                     ifelse((akt_pten$AKT1< median(akt_pten$AKT1) & akt_pten$AKT1> pctlsAKT[1]), 1,2)))
akt_pten$PTEN <- ifelse(akt_pten$PTEN < pctlsPTEN[1],0,ifelse(akt_pten$PTEN > pctlsPTEN[2], 3,
                                                                  ifelse((akt_pten$PTEN< median(akt_pten$PTEN) & akt_pten$PTEN> pctlsPTEN[1]), 1,2)))
akt_wdr$AKT1 <- ifelse(akt_wdr$AKT1  < pctlsAKT[1],0,ifelse(akt_wdr$AKT1  > pctlsAKT[2], 3,
                                                              ifelse((akt_wdr$AKT1< median(akt_wdr$AKT1) & akt_wdr$AKT1> pctlsAKT[1]), 1,2)))
akt_wdr$WDR48 <- ifelse(akt_wdr$WDR48 < pctlsWDR48[1],0,ifelse(akt_wdr$WDR48 > pctlsWDR48[2], 3,
                                                              ifelse((akt_wdr$WDR48< median(akt_wdr$WDR48) & akt_wdr$WDR48> pctlsWDR48[1]), 1,2)))

##fit for both
##label 1 for both high, 2 for PTEN high, WDR low, 3 for PTEN low, WDR high, 4 for both low, strictest for "both high/low"
both_genes <- both_genes[both_genes$PTEN == 3 | both_genes$PTEN == 0 | both_genes$WDR48 == 3 | both_genes$WDR48 == 0  ,]
both_genes$exprs<- ifelse(both_genes$PTEN >1 & both_genes$WDR48 >1, "Both High", 
                               ifelse(both_genes$PTEN <2 & both_genes$WDR48 <2, "Both Low", 
                               ifelse(both_genes$PTEN <2 & both_genes$WDR48 >1, "PTEN Low, WDR High"
                                      ,"PTEN High, WDR Low")))

both_genes <- both_genes[both_genes$PTEN == 3 | both_genes$PTEN == 0 | both_genes$WDR48 == 3 | both_genes$WDR48 == 0  ,]
both_genes$exprs<- ifelse(both_genes$PTEN >1 & both_genes$WDR48 >1, "Both High", 
                          ifelse(both_genes$PTEN <2 & both_genes$WDR48 <2, "Both Low", 
                                 ifelse(both_genes$PTEN <2 & both_genes$WDR48 >1, "PTEN Low, WDR High"
                                        ,"PTEN High, WDR Low")))
both_genes<- both_genes[both_genes$exprs == "Both High" | both_genes$exprs == "Both Low",]

akt_pten <- akt_pten[akt_pten$PTEN == 3 | akt_pten$PTEN == 0 | akt_pten$AKT1 == 3 | akt_pten$AKT1 == 0  ,]
akt_pten$exprs<- ifelse(akt_pten$PTEN >1 & akt_pten$AKT1 >1, "Both High", 
                          ifelse(akt_pten$PTEN <2 & akt_pten$AKT1 <2, "Both Low", 
                                 ifelse(akt_pten$PTEN <2 & akt_pten$AKT1 >1, "PTEN Low, AKT1 High"
                                        ,"PTEN High, AKT1 Low")))
akt_wdr <- akt_wdr[akt_wdr$WDR48 == 3 | akt_wdr$WDR48 == 0 | akt_wdr$AKT1 == 3 | akt_wdr$AKT1 == 0  ,]
akt_wdr$exprs<- ifelse(akt_wdr$WDR48 >1 & akt_wdr$AKT1 >1, "Both High", 
                        ifelse(akt_wdr$WDR48 <2 & akt_wdr$AKT1 <2, "Both Low", 
                               ifelse(akt_wdr$WDR48 <2 & akt_wdr$AKT1 >1, "WDR48 Low, AKT1 High"
                                      ,"WDR48 High, AKT1 Low")))

fitboth <- coxph(Surv(times, patient.vital_status)~exprs, data=both_genes)
fitWDR48 <- coxph(Surv(times, patient.vital_status)~WDR48, data=WDR48_only)
fitPTEN <- coxph(Surv(times, patient.vital_status)~PTEN, data=PTEN_only)
fitAKT <- coxph(Surv(times, patient.vital_status)~AKT1, data=AKT_only)
sfitAKW <-  survfit(Surv(times, patient.vital_status)~exprs, data=akt_wdr)
sfitWDR48 <- survfit(Surv(times, patient.vital_status)~exprs, data=WDR48_only)
sfitAKT <- survfit(Surv(times, patient.vital_status)~exprs, data=AKT_only)
sfitAKPT <- survfit(Surv(times, patient.vital_status)~exprs, data=akt_pten)
sfitTP <-  survfit(Surv(times, patient.vital_status)~exprs, data=TP53_only)
#summary(sfitWDR48, times=seq(0,365*5,365))
sfitPTEN <- survfit(Surv(times, patient.vital_status)~exprs, data=PTEN_only)
#summary(sfitPTEN, times=seq(0,365*5,365))
sfitboth <- survfit(Surv(times, patient.vital_status)~exprs, data=both_genes)
survPTEN <- ggsurvplot(sfitPTEN, conf.int=FALSE, pval=FALSE, xlab = "Time in Days", legend.labs = c("PTEN High, n = 295", "PTEN Low, n = 295"))
survWDR48 <- ggsurvplot(sfitWDR48, conf.int=FALSE, pval=FALSE, xlab = "Time in Days", legend.labs = c("WDR48 High, n = 295", "WDR48 Low, n = 295"))
survboth <- ggsurvplot(sfitboth, conf.int = FALSE, pval = FALSE, xlab = "Time in Days")
survAKT <- ggsurvplot(sfitAKT, conf.int=FALSE, pval=FALSE, xlab = "Time in Days",legend.labs = c("AKT1 High, n = 295", "AKT1 Low, n = 295"))
survAKPT <- ggsurvplot(sfitAKPT, conf.int=FALSE, pval=FALSE, xlab = "Time in Days")
survAKW <-ggsurvplot(sfitAKW, conf.int=FALSE, pval=FALSE, xlab = "Time in Days")
survTP <- ggsurvplot(sfitTP, conf.int=FALSE, pval=FALSE, xlab = "Time in Days")


