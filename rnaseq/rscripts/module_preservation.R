setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/separate_normalization")
# module preservation
#mods = list();
# Sof thresholding powers for network definition
power = c(3, 2);
#collectGarbage();
#labels = list();
#nn = if (interactive()) nSets else 1;

exp_wtdata <- read.csv("wtfm.csv")
exp_cfdata <- read.csv("cf.csv")
exp_cmdata <- read.csv("cm.csv")
exp_rmdata <- read.csv("rm.csv")
exp_rfdata <- read.csv("rw-rf_normalized.csv")

multiExpr <- list(wtdata = list(data = as.data.frame(t(exp_wtdata[,-c(1)]))),
            cfdata = list(data = as.data.frame(t(exp_cfdata[,-c(1)]))),
            cmdata = list(data = as.data.frame(t(exp_cmdata[,-c(1)]))),
            rmdata = list(data = as.data.frame(t(exp_rmdata[,-c(1)]))),
            rfdata = list(data = as.data.frame(t(exp_rfdata[,-c(1)]))))
names(multiExpr$wtdata$data) = exp_wtdata$Gene_ID
rownames(multiExpr$wtdata$data) = names(exp_wtdata)[-c(1)]
names(multiExpr$cfdata$data) = exp_cfdata$Gene_ID
rownames(multiExpr$cfdata$data) = names(exp_cfdata)[-c(1)]
names(multiExpr$cmdata$data) = exp_cmdata$Gene_ID
rownames(multiExpr$cmdata$data) = names(exp_cmdata)[-c(1)]
names(multiExpr$rmdata$data) = exp_rmdata$Gene_ID
rownames(multiExpr$rmdata$data) = names(exp_rmdata)[-c(1)]
names(multiExpr$rfdata$data) = exp_rfdata$Gene_ID
rownames(multiExpr$rfdata$data) = names(exp_rfdata)[-c(1)]
checkSets(multiExpr)


wtmods <- blockwiseModules(multiExpr$wtdata$data, checkMissingData = FALSE, maxBlockSize = 14000, corType = "bicor",
                               maxPOutliers = 0.1, power = 12, networkType="signed",
                               checkPower = FALSE, minModuleSize = 100, TOMType = "signed",
                               networkCalibration = "full quantile", saveConsensusTOMs = TRUE, numericLabels = FALSE,
                               deepSplit=4, mergeCutHeight = 0.1, verbose=4); 
table(wtmods$colors) %>% sort(decreasing = TRUE)
module.dist <- as.data.frame(table(wtmods$colors) %>% sort(decreasing = TRUE))
colnames(module.dist) <- c("Module", "Genes")
write.csv(module.dist,"module.distribution_wt_sft121.csv")

#clamsrwmods <- blockwiseModules(exp$femdata$data, checkMissingData = FALSE, maxBlockSize = 14000, corType = "bicor",
                           # maxPOutliers = 0.1, power = 9, networkType="signed",
                           # checkPower = FALSE, minModuleSize = 100, TOMType = "signed",
                           # networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                           # deepSplit=4, mergeCutHeight = 0.1, verbose=4); 
#table(clamsrwmods$colors) %>% sort(decreasing = TRUE)
#module.dist <- as.data.frame(table(clamsrwmods$colors) %>% sort(decreasing = TRUE))
#colnames(module.dist) <- c("Module", "Genes")
#write.csv(module.dist,"module.distribution_clams-rw0.8_sft9.csv")

#clamsmods <- blockwiseModules(exp2$femdata$data, checkMissingData = FALSE, maxBlockSize = 14000, corType = "bicor",
                              #  maxPOutliers = 0.1, power = 8, networkType="signed",
                              #  checkPower = FALSE, minModuleSize = 100, TOMType = "signed",
                              #  networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                              #  deepSplit=4, mergeCutHeight = 0.1, verbose=4); 
#table(clamsmods$colors) %>% sort(decreasing = TRUE)
#module.dist <- as.data.frame(table(clamsmods$colors) %>% sort(decreasing = TRUE))
#colnames(module.dist) <- c("Module", "Genes")
#write.csv(module.dist,"module.distribution_clams-cfcm0.8_sft8.csv")

rwmods <- blockwiseModules(multiExpr$rmdata$data, checkMissingData = FALSE, maxBlockSize = 14000, corType = "bicor",
                                maxPOutliers = 0.1, power = 12, networkType="signed",
                                checkPower = FALSE, minModuleSize = 100, TOMType = "signed",
                                networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                                deepSplit=4, mergeCutHeight = 0.1, verbose=4); 
table(rwmods$colors) %>% sort(decreasing = TRUE)
module.dist <- as.data.frame(table(rwmods$colors) %>% sort(decreasing = TRUE))
colnames(module.dist) <- c("Module", "Genes")
write.csv(module.dist,"module.distribution_rm0.9_sft12.csv")

rfmods <- blockwiseModules(multiExpr$rfdata$data, checkMissingData = FALSE, maxBlockSize = 14000, corType = "bicor",
                           maxPOutliers = 0.1, power = 8, networkType="signed",
                           checkPower = FALSE, minModuleSize = 100, TOMType = "signed",
                           networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                           deepSplit=4, mergeCutHeight = 0.1, verbose=4); 
table(rfmods$colors) %>% sort(decreasing = TRUE)
module.dist <- as.data.frame(table(rfmods$colors) %>% sort(decreasing = TRUE))
colnames(module.dist) <- c("Module", "Genes")
write.csv(module.dist,"module.distribution_rf0.9_sft8.csv")

cfmods <- blockwiseModules(multiExpr$cfdata$data, checkMissingData = FALSE, maxBlockSize = 14000, corType = "bicor",
                           maxPOutliers = 0.1, power = 8, networkType="signed",
                           checkPower = FALSE, minModuleSize = 100, TOMType = "signed",
                           networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                           deepSplit=4, mergeCutHeight = 0.1, verbose=4); 
table(cfmods$colors) %>% sort(decreasing = TRUE)
module.dist <- as.data.frame(table(cfmods$colors) %>% sort(decreasing = TRUE))
colnames(module.dist) <- c("Module", "Genes")
write.csv(module.dist,"module.distribution_cf0.9_sft8.csv")

cmmods <- blockwiseModules(multiExpr$cmdata$data, checkMissingData = FALSE, maxBlockSize = 14000, corType = "pearson",
                           maxPOutliers = 0.1, power = 16, networkType="signed",
                           checkPower = FALSE, minModuleSize = 100, TOMType = "signed",
                           networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                           deepSplit=4, mergeCutHeight = 0.1, verbose=4); 
table(cmmods$colors) %>% sort(decreasing = TRUE)
module.dist <- as.data.frame(table(cmmods$colors) %>% sort(decreasing = TRUE))
colnames(module.dist) <- c("Module", "Genes")
write.csv(module.dist,"module.distribution_cm0.9_sft16.csv")

# Generate MEs
wtMEs <- orderMEs(wtmods$MEs)
#clamsrwMEs <- orderMEs(clamsrwmods$MEs)
#clamsMEs <- orderMEs(clamsmods$MEs)
rwMEs <- orderMEs(rwmods$MEs)
rfMEs <- orderMEs(rfmods$MEs)
cfMEs <- orderMEs(cfmods$MEs)
cmMEs <- orderMEs(cmmods$MEs)

# Calculate Module Membership ####

moduleMembershipwt <- bicorAndPvalue(multiExpr$wtdata$data, wtMEs, 
                              alternative = "two.sided", use = "pairwise.complete.obs", 
                                              maxPOutliers = 0.1)

#moduleMembershipcrw <- bicorAndPvalue(exp$femdata$data, clamsrwMEs, 
                                     #alternative = "two.sided", use = "pairwise.complete.obs", 
                                     #maxPOutliers = 0.1)

#moduleMembershipcfcm <- bicorAndPvalue(exp2$femdata$data, clamsMEs, 
                                      #alternative = "two.sided", use = "pairwise.complete.obs", 
                                      #maxPOutliers = 0.1)

moduleMembershiprw <- bicorAndPvalue(multiExpr$rmdata$data, rwMEs, 
                                      alternative = "two.sided", use = "pairwise.complete.obs", 
                                      maxPOutliers = 0.1)

moduleMembershiprf <- bicorAndPvalue(multiExpr$rfdata$data, rfMEs, 
                                     alternative = "two.sided", use = "pairwise.complete.obs", 
                                     maxPOutliers = 0.1)

moduleMembershipcf <- bicorAndPvalue(multiExpr$cfdata$data, cfMEs, 
                                      alternative = "two.sided", use = "pairwise.complete.obs", 
                                      maxPOutliers = 0.1)

moduleMembershipcm <- bicorAndPvalue(multiExpr$cmdata$data, cmMEs, 
                                      alternative = "two.sided", use = "pairwise.complete.obs", 
                                      maxPOutliers = 0.1)

ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host="https://useast.ensembl.org")

# Generate module membership files

MM_wt <- as.data.frame(moduleMembershipwt$bicor)
colnames(MM_wt) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_wt), fixed = TRUE)
MM_wt$Probe <- rownames(MM_wt)
MM_wt$Module <- wtmods$colors

write.table(MM_wt, "WT Module Membershipsft121.txt", sep = "\t", quote = FALSE, row.names = FALSE)

MM_rw <- as.data.frame(moduleMembershiprw$bicor)
colnames(MM_rw) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_rw), fixed = TRUE)
MM_rw$Probe <- rownames(MM_rw)
MM_rw$Module <- rwmods$colors

write.table(MM_rw, "RW Module Membershipsft12.txt", sep = "\t", quote = FALSE, row.names = FALSE)

MM_rf <- as.data.frame(moduleMembershiprf$bicor)
colnames(MM_rf) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_rf), fixed = TRUE)
MM_rf$Probe <- rownames(MM_rf)
MM_rf$Module <- rfmods$colors

write.table(MM_rf, "RF Module Membershipsft8.txt", sep = "\t", quote = FALSE, row.names = FALSE)

MM_cf <- as.data.frame(moduleMembershipcf$bicor)
colnames(MM_cf) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_cf), fixed = TRUE)
MM_cf$Probe <- rownames(MM_cf)
MM_cf$Module <- cfmods$colors

write.table(MM_cf, "CF Module Membershipsft8.txt", sep = "\t", quote = FALSE, row.names = FALSE)

MM_cm <- as.data.frame(moduleMembershipcm$bicor)
colnames(MM_cm) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_cm), fixed = TRUE)
MM_cm$Probe <- rownames(MM_cm)
MM_cm$Module <- cmmods$colors

write.table(MM_cm, "CM Module Membershipsft16.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Get Module Hub Probes and Genes ####
hubProbes_wt <- sapply(colnames(MM_wt)[!colnames(MM_wt) %in% c("Probe", "Module")], function(x){
  temp <- MM_wt[MM_wt$Module == x,]
  temp$Probe[temp[, x] == max(temp[, x])] %>% as.character %>% unique %>% sort})

hubProbes_rw <- sapply(colnames(MM_rw)[!colnames(MM_rw) %in% c("Probe", "Module")], function(x){
  temp <- MM_rw[MM_rw$Module == x,]
  temp$Probe[temp[, x] == max(temp[, x])] %>% as.character %>% unique %>% sort})

hubProbes_rf <- sapply(colnames(MM_rf)[!colnames(MM_rf) %in% c("Probe", "Module")], function(x){
  temp <- MM_rf[MM_rf$Module == x,]
  temp$Probe[temp[, x] == max(temp[, x])] %>% as.character %>% unique %>% sort})

hubProbes_cf <- sapply(colnames(MM_cf)[!colnames(MM_cf) %in% c("Probe", "Module")], function(x){
  temp <- MM_cf[MM_cf$Module == x,]
  temp$Probe[temp[, x] == max(temp[, x])] %>% as.character %>% unique %>% sort})

hubProbes_cm <- sapply(colnames(MM_cm)[!colnames(MM_cm) %in% c("Probe", "Module")], function(x){
  temp <- MM_cm[MM_cm$Module == x,]
  temp$Probe[temp[, x] == max(temp[, x])] %>% as.character %>% unique %>% sort})

ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host="https://useast.ensembl.org")

hubGenes_wt <- lapply(hubProbes_wt, function(x){
  getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", values = x, mart = ensembl, 
        verbose = TRUE) %>% unlist %>% as.character %>% unique %>% sort %>% paste(collapse = ", ")}) %>% unlist
hubGenes_wt
write.csv(hubGenes_wt, "hubGenes_wt12.csv")

hubGenes_rw <- lapply(hubProbes_rw, function(x){
  getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", values = x, mart = ensembl, 
        verbose = TRUE) %>% unlist %>% as.character %>% unique %>% sort %>% paste(collapse = ", ")}) %>% unlist
hubGenes_rw
write.csv(hubGenes_rw, "hubGenes_rw12.csv")

hubGenes_rf <- lapply(hubProbes_rf, function(x){
  getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", values = x, mart = ensembl, 
        verbose = TRUE) %>% unlist %>% as.character %>% unique %>% sort %>% paste(collapse = ", ")}) %>% unlist
hubGenes_rf
write.csv(hubGenes_rf, "hubGenes_rf8.csv")

hubGenes_cf <- lapply(hubProbes_cf, function(x){
  getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", values = x, mart = ensembl, 
        verbose = TRUE) %>% unlist %>% as.character %>% unique %>% sort %>% paste(collapse = ", ")}) %>% unlist
hubGenes_cf
write.csv(hubGenes_cf, "hubGenes_cf8.csv")

hubGenes_cm <- lapply(hubProbes_cm, function(x){
  getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", values = x, mart = ensembl, 
        verbose = TRUE) %>% unlist %>% as.character %>% unique %>% sort %>% paste(collapse = ", ")}) %>% unlist
hubGenes_cm
write.csv(hubGenes_cm, "hubGenes_cm16.csv")

### Load Traits ###

### RW ###

cov_wt <- read.csv("wtfm_traits.csv")
cov_rw <- read.csv("rw-rm_traits.csv")
cov_rw <- read.csv("rw_rm.me.traits.csv")
cov_rf <- read.csv("rf_traits.csv")
cov_cf <- read.csv("cf_traits.csv")
cov_cf <- read.csv("clams_cf.me.traits.csv")
cov_cm <- read.csv("cm_traits4.csv")
cov_cm <- read.csv("clams_cm.me2.traits.csv")
cov_cm <- read.csv("clams_cm.me3.traits.csv")
cov_all <- read.csv("clams-rw_traits.csv")

##########################################################################
##                   Setting Up Module Trait Correlations               ##
##########################################################################

##Notes## 
## Module trait correlations require the rownames of input data and trait data to match
## If rownames are not equal the output data is incorrect 
## Make sure rownames(input) == rownames(traits) by using the table function 

rownames(cov_wt) <- cov_wt$Mice
cov_wt <- cov_wt[,c("Timepoint", "SexScore")]
cov_wt <- as.matrix(cov_wt)
table(rownames(multiExpr$wtdata$data) == rownames(cov_wt)) # All TRUE
table(rownames(multiExpr$wtdata$data) == rownames(cov_wt)) # All TRUE

rownames(cov_rw) <- cov_rw$SampleID
cov_rw <- cov_rw[,c("Entrainment", "del", "ovexp", "comp", "PeriodBasal", "AmplitudeBasal", "PeriodFirst4", "AmplitudeFirst4", "PeriodLast4", "AmplitudeLast4")]
cov_rw <- as.matrix(cov_rw)
table(rownames(multiExpr$rmdata$data) == rownames(cov_rw)) # All TRUE
library(dplyr)
rownames(cov_rf) <- cov_rf$ID
cov_rf <- cov_rf[,c("Entrainment", "del", "ovexp", "comp", "PeriodBasal", "AmplitudeBasal", "PeriodFirst4", "AmplitudeFirst4", "PeriodLast4", "AmplitudeLast4")]
cov_rf %>% arrange(cov_rf$ID, decreasing = FALSE)
cov_rf <- as.matrix(cov_rf)

table(rownames(multiExpr$rfdata$data) == rownames(cov_rf)) # All TRUE

rownames(cov_rw) <- cov_rw$SampleID
cov_rw <- cov_rw[,c("Entrainment", "wt", "del", "ovexp", "comp", "delvwt", "ovexpvwt", "compvwt", "ovexpvdel", "compvdel", "compvovexp", "PeriodBasal", "AmplitudeBasal", "PeriodFirst4", "AmplitudeFirst4", "PeriodLast4", "AmplitudeLast4")]
cov_rw <- as.matrix(cov_rw)
table(rownames(multiExpr$rmdata$data) == rownames(cov_rw)) # All TRUE

rownames(cov_cf) <- cov_cf$SampleID
cov_cf <- cov_cf[,c("Entrainment", "del", "ovexp", "comp", "AvgRER", "AvgLightRER", "AvgDarkRER", "AvgHeat", "AvgLightHeat", "AvgDarkHeat", "AvgKcalIntake", "AvgLightKcalIntake", "AvgDarkKcalIntake")]
cov_cf <- as.matrix(cov_cf)
table(rownames(multiExpr$cfdata$data) == rownames(cov_cf)) # All TRUE

rownames(cov_cf) <- cov_cf$SampleID
cov_cf <- cov_cf[,c("Entrainment", "wt", "del", "ovexp", "comp", "delvwt", "ovexpvwt", "compvwt", "ovexpvdel", "compvdel", "compvovexp", "AvgRER", "AvgLightRER", "AvgDarkRER", "AvgHeat", "AvgLightHeat", "AvgDarkHeat", "AvgKcalIntake", "AvgLightKcalIntake", "AvgDarkKcalIntake")]
cov_cf <- as.matrix(cov_cf)
table(rownames(multiExpr$cfdata$data) == rownames(cov_cf)) # All TRUE

rownames(cov_cm) <- cov_cm$SampleID
cov_cm <- cov_cm[,c("Entrainment", "EntrainmentScore", "del", "ovexp", "comp", "AvgRER", "AvgLightRER", "AvgDarkRER", "AvgHeat", "AvgLightHeat", "AvgDarkHeat", "AvgKcalIntake", "AvgLightKcalIntake", "AvgDarkKcalIntake")]
cov_cm <- as.matrix(cov_cm)
table(rownames(multiExpr$cmdata$data) == rownames(cov_cm)) # All TRUE


rownames(cov_cm) <- cov_cm$sampleID
cov_cm <- cov_cm[,c("Entrainment", "wt", "del", "ovexp", "comp", "delvwt", "ovexpvwt", "compvwt", "ovexpvdel", "compvdel", "compvovexp", "AvgRER", "AvgLightRER", "AvgDarkRER", "AvgHeat", "AvgLightHeat", "AvgDarkHeat", "AvgKcalIntake", "AvgLightKcalIntake", "AvgDarkKcalIntake")]
cov_cm <- as.matrix(cov_cm)
table(rownames(multiExpr$cmdata$data) == rownames(cov_cm)) # All TRUE

rownames(cov_all) <- cov_all$SampleID
cov_all <- cov_all[,c("SampleID", "Entrainment", "del", "ovexp", "comp", "SexScore", "AvgRER", "AvgLightRER", "AvgDarkRER", "AvgHeat", "AvgLightHeat", "AvgDarkHeat", "AvgKcalIntake", "AvgLightKcalIntake", "AvgDarkKcalIntake", "PeriodBasal", "AmplitudeBasal", "PeriodFirst4", "AmplitudeFirst4", "PeriodLast4", "AmplitudeLast4")]
cov_all <- as.matrix(cov_all)
table(rownames(multiExpr$cmdata$data) == rownames(cov_all)) # All TRUE

pheno <- list(wt = list(data = cov_wt),
              rw = list(data = cov_rw),
              cf = list(data = cov_cf),
              cm = list(data = cov_cm))
pheno$cm

pheno2 <- list(wt = list(data = cov_wt),
               rw = list(data = cov_rw),
               cf = list(data = cov_cf),
               cm = list(data = cov_cm),
               all = list(data = cov_all))

# Get Meta-Analysis Correlations ####
moduleMembership <- MM_wt
Mods <- moduleMembership$Module
MEs_wt <- moduleEigengenes(t(exp_wtdata[,-c(1)]), colors = Mods)$eigengenes
rownames(MEs_wt) <- rownames(t(exp_wtdata[,-c(1)]))
MEs_wt
MM_wt
write.csv(MEs_wt, "MEs_wt12.csv")

moduleMembership <- MM_rw
Mods <- moduleMembership$Module
MEs_rw <- moduleEigengenes(t(exp_rmdata[,-c(1)]), colors = Mods)$eigengenes
rownames(MEs_rw) <- rownames(t(exp_rmdata[,-c(1)]))
MEs_rw
MM_rw
write.csv(MEs_rw, "MEs_rw12.csv")

moduleMembership <- MM_rf
Mods <- moduleMembership$Module
MEs_rf <- moduleEigengenes(t(exp_rfdata[,-c(1)]), colors = Mods)$eigengenes
rownames(MEs_rf) <- rownames(t(exp_rfdata[,-c(1)]))
MEs_rf
MM_rf
write.csv(MEs_rf, "MEs_rf8.csv")

moduleMembership <- MM_cf
Mods <- moduleMembership$Module
MEs_cf <- moduleEigengenes(t(exp_cfdata[,-c(1)]), colors = Mods)$eigengenes
rownames(MEs_cf) <- rownames(t(exp_cfdata[,-c(1)]))
MEs_cf
MM_cf
write.csv(MEs_cf, "MEs_cf8.csv")

moduleMembership <- MM_cm
Mods <- moduleMembership$Module
MEs_cm <- moduleEigengenes(t(exp_cmdata[,-c(1)]), colors = Mods)$eigengenes
rownames(MEs_cm) <- rownames(t(exp_cmdata[,-c(1)]))
MEs_cm
MM_cm
write.csv(MEs_cm, "MEs_cm16.csv")

MEs <- list(wt = list(data = MEs_wt), 
            rw = list(data = MEs_rw),
            cf = list(data = MEs_cf),
            cm = list(data = MEs_cm))

#################################################################
## Function to calculate module eigengenes and write CSV files ##
#################################################################
calculateMEsAndWriteCSV <- function(moduleMembership, expData, outputFileName) {
  Mods <- moduleMembership$Module
  MEs <- moduleEigengenes(t(expData[,-c(1)]), colors = Mods)$eigengenes
  rownames(MEs) <- rownames(t(expData[,-c(1)]))
  MM <- moduleMembership
  write.csv(MEs, outputFileName)
  list(MEs = MEs, MM = MM)
}

##################################################################
## Calculate and write module eigengenes for different datasets ##
##################################################################
MEs_wt <- calculateMEsAndWriteCSV(MM_wt, exp_wtdata, "MEs3_wt.csv")
MEs_rw <- calculateMEsAndWriteCSV(MM_rw, exp_rmdata, "MEs_rw.csv")
MEs_cf <- calculateMEsAndWriteCSV(MM_cf, exp_cfdata, "MEs_cf.csv")
MEs_cm <- calculateMEsAndWriteCSV(MM_cm, exp_cmdata, "MEs_cm.csv")

#Module Trait Correlations
lengths(MEs$cf$data)

checkSets(multiExpr)$nSamples
moduleTraitCorwt = cor(MEs$wt$data, pheno$wt$data, method = "spearman", use="p"); 
moduleTraitPvaluewt = corPvalueStudent(moduleTraitCorwt, 92)
moduleTraitCorwt = cor(WtME_reordered, pheno$wt$data, method = "spearman", use="p"); 
moduleTraitPvaluewt = corPvalueStudent(moduleTraitCorwt, 92)

moduleTraitCorrw = cor(MEs$rw$data, pheno$rw$data, method = "spearman", use="p"); 
moduleTraitPvaluerw = corPvalueStudent(moduleTraitCorrw, 30)

moduleTraitCorcf = cor(MEs$cf$data, pheno$cf$data, method = "spearman", use="p"); 
moduleTraitPvaluecf = corPvalueStudent(moduleTraitCorcf, 32)

moduleTraitCorcm = cor(MEs$cm$data, pheno$cm$data, method = "spearman", use="p"); 
moduleTraitPvaluecm = corPvalueStudent(moduleTraitCorcm, 31)

sizeGrWindow(width = 22, height = 30)

star <- apply(moduleTraitPvaluerw, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})

textMatrix=paste(signif(moduleTraitCorrw,2), 
                 star, sep=""); 
#dim(textMatrix) = dim(moduleTraitCor) 
par(mar = c(6, 7.5, 1, 0.5))
#Displaythecorrelationvalueswithinaheatmapplot 

labeledHeatmap(Matrix=moduleTraitCorcm, xLabels=colnames(pheno$cm$data), yLabels=names(MEs$cm$data), 
               ySymbols=gsub("ME", "", names(MEs$cm$data)), colorLabels=FALSE, colors=blueWhiteRed(50), 
               textMatrix=star, setStdMargins=FALSE, cex.text=0.7, textAdj = c(0.5, 0.8), 
               zlim=c(-1,1), main = "CM Module Trait Spearman Correlation", cex.lab.y = 1, plotLegend = TRUE, legendLabel = "spearman correlation")
dev.off()

labeledHeatmap(Matrix=moduleTraitCorcf, xLabels=colnames(pheno$cf$data), yLabels=names(MEs$cf$data), 
               ySymbols=gsub("ME", "", names(MEs$cf$data)), colorLabels=FALSE, colors=blueWhiteRed(50), 
               textMatrix=star, setStdMargins=FALSE, cex.text=0.7, textAdj = c(0.5, 0.8), 
               zlim=c(-1,1), main = "CF Module Trait Spearman Correlation", cex.lab.y = 1, plotLegend = TRUE, legendLabel = "spearman correlation")
dev.off()

labeledHeatmap(Matrix=moduleTraitCorrw, xLabels=colnames(pheno$rw$data), yLabels=names(MEs$rw$data), 
               ySymbols=gsub("ME", "", names(MEs$rw$data)), colorLabels=FALSE, colors=blueWhiteRed(50), 
               textMatrix=star, setStdMargins=FALSE, cex.text=0.7, textAdj = c(0.5, 0.8), 
               zlim=c(-1,1), main = "RW Module Trait Spearman Correlation", cex.lab.y = 1, plotLegend = TRUE, legendLabel = "spearman correlation")
dev.off()

labeledHeatmap(Matrix=moduleTraitCorwt, xLabels=colnames(pheno$wt$data), yLabels=names(MEs$wt$data), 
               ySymbols=gsub("ME", "", names(MEs$wt$data)), colorLabels=FALSE, colors=blueWhiteRed(50), 
               textMatrix=star, setStdMargins=FALSE, cex.text=0.7, textAdj = c(0.5, 0.8), 
               zlim=c(-1,1), main = "WT Module Trait Spearman Correlation", cex.lab.y = 1, plotLegend = TRUE, legendLabel = "spearman correlation")
dev.off()
#Module Preservation Statistics#
multiExpr$wtdata$data
colorswtdata = wtmods$colors
multiColor=list(wtdata=colorswtdata)
checkSets(multiExpr)

system.time( { 
  mp=modulePreservation(multiExpr, multiColor, 
                        referenceNetworks=1, 
                        nPermutations=20, 
                        randomSeed=1, 
                        quickCor=0, 
                        verbose=3) 
} );

ref=1 
test=4 
statsObs=cbind(mp$quality$observed[[ref]][[test]][,-1],mp$preservation$observed[[ref]][[test]][,-1]) 
statsZ=cbind(mp$quality$Z[[ref]][[test]][,-1],mp$preservation$Z[[ref]][[test]][,-1]);

#Compare preservation to quality: 
print(cbind(statsObs[,c("medianRank.pres","medianRank.qual")], 
            signif(statsZ[,c("Zsummary.pres","Zsummary.qual")],2)))

write.csv(print(cbind(statsObs[,c("medianRank.pres","medianRank.qual")], 
                      signif(statsZ[,c("Zsummary.pres","Zsummary.qual")],2))), "mp_wtvrm2.csv")


#Module labels andmodule sizes are also contained in the results 
modColors=rownames(mp$preservation$observed[[ref]][[test]]) 
moduleSizes=mp$preservation$Z[[ref]][[test]][,1];
mpModswtcf <- data.frame(modColors, moduleSizes)
length(modColors)
length(moduleSizes)
#leave grey and gold modulesout 
plotMods=!(modColors%in%c("gold"));

#Text labels forpoints 
text=modColors[plotMods];

#Auxiliary convenience variable  

plotData=cbind(mp$preservation$observed[[ref]][[test]][,2], mp$preservation$Z[[ref]][[test]][,2]) 
#Maintitlesfortheplot 
mains=c("Preservation Median rank WT-RM","Preservation Zsummary");

#Start the plot 
sizeGrWindow(10,5);

#pdf(fi="Plots/BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf",wi=10,h=5) 
par(mfrow=c(1,2)) 
par(mar=c(4.5,4.5,2.5,1)) 
for (p in 1:2) 
{ 
  min = min(plotData[, p] ,na.rm = TRUE); 
  max = max(plotData[, p] ,na.rm=TRUE);
  #Adjust ploting ranges appropriately 
  if (p==2) 
  { 
    if ( min > -max/10) min = -max/10 
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)) 
  } else 
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min)) 
  plot(moduleSizes[plotMods], plotData[plotMods, p], col=1, bg = modColors[plotMods], pch=21, 
       main = mains[p], 
       cex = 2.4, 
       ylab = mains[p], xlab="Module size", log = "x", 
       ylim = ylim, 
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.4) 
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex=1, offs=0.08);
  #ForZsummary,addthresholdlines 
  if(p==2) 
  { 
    abline(h=0) 
    abline(h=2,col="blue",lty=2) 
    abline(h=10,col="darkgreen",lty=2) 
  } 
}
dev.off()

#moduleColors is testset
#wtColors is the ref set
wtColors = MM_cm$Module
moduleColors = MM_cf$Module
#wtMEs = orderMEs(MEs)
#moduleColors = wtColors

#
#wtModuleLabels = substring(names(MEs_cm), 3)
#rwModuleLabels = substring(names(MEs_rw), 3)
#cfModuleLabels = substring(names(MEs_cf), 3)
#cmModuleLabels = substring(names(MEs_cm), 3)

wtModules = substring(names(MEs_cm), 3)
rwModules = substring(names(MEs_cf), 3)
#cfModules = cfModuleLabels
#cmModules = cmModuleLabels
#cfModules = substring(names(MEs_cf), 3)
#cmModules = substring(names(MEs_cm), 3)
# Number of modules
nwtMods = length(wtModules)
nrwMods = length(rwModules)
#ncfMods = length(cfModules)
#ncmMods = length(cmModules)

# Initialize tables of p-values and of the corresponding coutns
pTable = matrix(0, nrow = nwtMods, ncol = nrwMods)
CountTbl = matrix(0, nrow = nwtMods, ncol = nrwMods)

# Pairwise comparisons
for (wtmod in 1:nwtMods)
  for (rwmod in 1:nrwMods)
  {
    wtMembers = (wtColors == wtModules[wtmod]);
    rwMembers = (moduleColors == rwModules[rwmod]);
    pTable[wtmod, rwmod] = -log10(fisher.test(wtMembers, rwMembers, alternative = "greater")$p.value); 
    CountTbl[wtmod, rwmod]=sum(wtColors==wtModules[wtmod] & moduleColors == 
                              rwModules[rwmod])
  }

#Truncate pvalues smaller than 10^{-50} to 10^{-50} 
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]); 
pTable[pTable>50] = 50; 
#Marginalcounts(reallymodulesizes) 
wtModTotals=apply(CountTbl,1,sum) 
rwModTotals=apply(CountTbl,2,sum) 
#Actualplotting 
sizeGrWindow(10,7); 
#pdf(file="Plots/ConsensusVsFemaleModules.pdf",wi=10,he=7); 
par(mfrow=c(1,1)); 
par(cex=1.0); 
par(mar=c(8,10.4,2.7,1)+0.3);
#UsefunctionlabeledHeatmaptoproducethecolor-codedtablewithallthetrimmings
labeledHeatmap(Matrix=pTable, 
  xLabels=paste(" ", rwModules), 
  yLabels=paste(" ", wtModules), 
  colorLabels = TRUE, xSymbols=paste("CF", rwModules, ":", rwModTotals, sep=""), 
  ySymbols=paste("CM", wtModules, ":", wtModTotals, sep=""), 
  textMatrix=CountTbl, colors=greenWhiteRed(100)[50:100], 
  main="CM vs CF Modules", 
  cex.text=0.5, cex.lab=0.8, setStdMargins=FALSE); 
dev.off();

#Adjacency Analysis
adjWT = adjacency(multiExpr$wtdata$data, power = 5, type = "signed")
adjRW = adjacency(multiExpr$rmdata$data, power = 5, type = "signed")
adjCF = adjacency(multiExpr$cfdata$data, power = 3, type = "signed")
adjCM = adjacency(multiExpr$cmdata$data, power = 7, type = "signed")

#Connectivity
ConnectivityWT = ConnectivityWT/max(ConnectivityWT);
ConnectivityRW = ConnectivityRW/max(ConnectivityRW)
ConnectivityCF = ConnectivityCF/max(ConnectivityCF)
ConnectivityCM = ConnectivityCM/max(ConnectivityCM)
