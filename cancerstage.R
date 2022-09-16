clinical <- read.delim("BRCA.clin.merged.txt", header = TRUE, sep = "\t", dec = ".")
mrna <- read.delim("BRCA.transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.data.txt", header = TRUE, sep = "\t", dec = ".")
merged <- merge(clinical, mrna, by.x="patient.samples.sample.portions.portion.analytes.analyte.2.aliquots.aliquot.bcr_aliquot_barcode", by.y="Hybridization.REF")
merged <- merged[!is.na(merged$patient.stage_event.pathologic_stage), ]
#merged <- merged[merged$patient.stage_event.pathologic_stage != "stage x",]
merged <- rename(merged, Stage = patient.stage_event.pathologic_stage)
merged <- rename(merged, Barcode = X)
merged <- rename(merged, Barcode1 = patient.samples.sample.portions.portion.analytes.analyte.2.aliquots.aliquot.bcr_aliquot_barcode)
merged$Stage[merged$Stage=="stage ia"] <- "stage i"
merged$Stage[merged$Stage=="stage ib"] <- "stage i"
merged$Stage[merged$Stage=="stage iia"] <- "stage ii"
merged$Stage[merged$Stage=="stage iib"] <- "stage ii"
merged$Stage[merged$Stage=="stage iiia"] <- "stage iii"
merged$Stage[merged$Stage=="stage iiib"] <- "stage iii"
merged$Stage[merged$Stage=="stage iiic"] <- "stage iii"

PTEN_avg <- tapply(merged$PTEN, merged$Stage, mean)
WDR48_avg <- tapply(merged$WDR48, merged$Stage, mean)
AKT_avg<- tapply(merged$AKT1, merged$Stage, mean)

fitmerged <- coxph(Surv(aggregate(d[, 3:4], list(d$Name), mean))~ Stage, data=merged)



l <- length(merged$PTEN) #length = 590 for both
bounds <- c(round(l/4), round(l-l/4))
#top5bounds <- c(round(590/20), round(590-590/20))
pctlsWDR48 <- c(sort(merged$WDR48)[bounds[1]], sort(merged$WDR48)[bounds[2]])
pctlsPTEN <- c(sort(merged$PTEN)[bounds[1]], sort(merged$PTEN)[bounds[2]])
pctlsAKT <- c(sort(merged$AKT1)[bounds[1]], sort(merged$AKT1)[bounds[2]])
merged$WDR48 <- ifelse(merged$WDR48 < pctlsWDR48[1],0,ifelse(merged$WDR48 > pctlsWDR48[2], 3,
                ifelse((merged$WDR48< median(merged$WDR48) & merged$WDR48> pctlsWDR48[1]), 1,2)))

merged$PTEN <- ifelse(merged$PTEN < pctlsPTEN[1],0,ifelse(merged$PTEN > pctlsPTEN[2], 3,
                    ifelse((merged$PTEN< median(merged$PTEN) & merged$PTEN> pctlsPTEN[1]), 1,2)))
merged$AKT1 <- ifelse(merged$AKT1 < pctlsAKT[1],0,ifelse(merged$AKT1 > pctlsAKT[2], 3,
                    ifelse((merged$AKT1< median(merged$AKT1) & merged$AKT1> pctlsAKT[1]), 1,2)))

merged <- merged[merged$PTEN == 3 | merged$PTEN == 0 | merged$WDR48 == 3 | merged$WDR48 == 0  ,]
merged$exprs<- ifelse(merged$PTEN >1 & merged$WDR48 >1, "Both High", 
                          ifelse(merged$PTEN <2 & merged$WDR48 <2, "Both Low", 
                                 ifelse(merged$PTEN <2 & merged$WDR48 >1, "PTEN Low, WDR High"
                                        ,"PTEN High, WDR Low")))
