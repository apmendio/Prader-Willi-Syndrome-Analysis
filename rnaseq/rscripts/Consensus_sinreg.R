#### Using sine regression to assess associations between WGCNA gene modules and time ####

library(WGCNA)
BiocManager::install("rain")
BiocManager::install("ShellChron")
library(rain)
library(ShellChron)
library(ggplot2)
BiocManager::install("openxslx")
library(openxlsx)
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Mm.eg.db")
library(AnnotationDbi)
library(org.Mm.eg.db)
library(RColorBrewer)
library(data.table)
getwd()
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/individualWGCNA")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/individualWGCNA/optimization_folder/test")

MEnames = colnames(MEs_cm)
MEnames = colnames(MEs_wt)
MEnames = colnames(wtMEs)
MEnames = colnames(MEs_female)
MEnames = colnames(consensusMEs$wtfdata$data)
MEnames = colnames(consensusMEs$wtmdata$data)
no.MEs = length(MEnames)

traits = as.data.frame(pheno$wt$data)
traits = as.data.frame(pheno$cm$data)
traits = as.data.frame(pheno2$wtf$data)
traits = as.data.frame(pheno2$wtm$data)
traits = as.data.frame(cov_wt)
#row.names(traits) <- traits$sampleID
row.names(traits)
row.names(wtMEs)
#row.names(MEs_cm)
row.names(consensusMEs$wtfdata$data)
alldata = merge(traits, MEs_cm, by="row.names")
alldata = merge(traits, MEs_wt, by="row.names")
alldata = merge(traits, consensusMEs$wtfdata$data, by="row.names")
alldata = merge(traits, consensusMEs$wtmdata$data, by="row.names")
alldata = merge(traits, wtMEs, by="row.names")

alldata = alldata[
  with(alldata, order(alldata$Timepoint)),
]

head(alldata)
write.csv(alldata, "wtm_eigenvalue_update.csv")

fitlist = as.list(1:no.MEs)
names(fitlist) <- MEnames

p.values = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(p.values) = MEnames
colnames(p.values) = "p.value"

R2 = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(R2) = MEnames
colnames(R2) = "R2"

period = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(period) = MEnames
colnames(period) = "period"

peak = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(peak) = MEnames
colnames(peak) = "peak"

amplitude = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(amplitude) = MEnames
colnames(amplitude) = "amplitude"

plotpoints = data.frame(matrix(ncol=no.MEs, nrow=nrow(alldata)))
rownames(plotpoints) = alldata$Row.names
colnames(plotpoints) = MEnames

#setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/consensus_analysis/version2/wt/universalWT")

for(i in MEnames){
  
  # print status
  print(paste("Running entity:", i, "which is", which(MEnames==i), "out of", no.MEs))
  
  #create temporary data matrix and model formula
  
  x = alldata$Timepoint
  y = alldata[,i]
  fitlist[[i]] <- sinreg(x, y, plot=FALSE)
  p.values[i,1] <- fitlist[[i]][[1]][6]
  R2[i,1] <- fitlist[[i]][[1]][5]
  period[i,1] <- fitlist[[i]][[1]][3]
  peak[i,1] <- fitlist[[i]][[1]][4]
  amplitude[i,1] <- fitlist[[i]][[1]][2]
  plotpoints[,i] <- fitlist[[i]][[2]]
}

results = cbind(p.values, R2, period, peak, amplitude)

results$FDR = p.adjust(results$p.value, method="fdr")
plotpoints$Timepoint = alldata$Timepoint

for(i in MEnames) {
  y = alldata[,i]
  y1 = plotpoints[,i]
  ggplot(data=alldata, aes(x=Timepoint, y=y)) +
    geom_point() +
    geom_line(data=plotpoints, aes(x=Timepoint, y=y1, color="red")) +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
    ggtitle(paste(i)) +
    ylab("Module EigenValue") +
    xlab("Zeitgeber Time (ZT)") +
    theme_classic() +
    theme(legend.position="none")
  ggsave(paste("Sine_regression_plots_", i, "_WTM_update.pdf"))
}

# order results by best fit (highest R2 value) and save results #
results.ordered = results[
  with(results, order(results$R2, decreasing=TRUE)),
]

results.ordered$Module = rownames(results.ordered)

openxlsx::write.xlsx(results.ordered, file="Sine_regresion_results_WTM_update.xlsx")

# Checking if hub genes within modules cycle in a similar manner to EigenValues #

exp_wtdata.1 = as.data.frame(t(wtm_norm))
#colnames(exp_wtdata.1) = wtf$Gene_ID
rownames(exp_wtdata.1)
rownames(traits)
exp_wtdata.2 = merge(traits, exp_wtdata.1, by="row.names")

exp_wtdata.2 = exp_wtdata.2[
  with(exp_wtdata.2, order(exp_wtdata.2$Timepoint)),
]

#geneSymbolswt = AnnotationDbi::mapIds(org.Mm.eg.db,
                                   keys = exp_wtdata$Gene_ID,
                                   column = "SYMBOL",
                                   keytype = 'ENSEMBL') 
#wt_genes <- exp_wtdata
#wt_genes$Gene_Name <- geneSymbolswt
#wt_genes2 <- wt_genes
#wt_genes2$Gene_Name <- toupper(wt_genes2$Gene_Name)
#wt_genes2$Gene_Name
#exp_baboon <- read.csv("baboon_allgenes_counts.csv")
#test <- match(wt_genes$Gene_Name, exp_baboon$Gene_name)
#test <- merge(wt_genes2, exp_baboon, by="Gene_Name")
#write.csv(test, "baboon_merge.csv")
#test <- WtME[match(rain.order$Module, rownames(WtME)),]

hubSymbols = AnnotationDbi::mapIds(org.Mm.eg.db,
                                   keys = test,
                                   column = "SYMBOL",
                                   keytype = 'ENSEMBL') 

for(i in hubProbes_female) {
  y = exp_wtdata.2[,i]
  sinreg_hub = sinreg(exp_wtdata.2$Timepoint, y, plot=FALSE)
  
  y1 = exp_wtdata.2[,i]
  y2 = sinreg_hub[[2]]
  
  ggplot(data=exp_wtdata.2, aes(x=Timepoint, y=y1)) +
    geom_point() +
    geom_line(data=plotpoints, aes(x=Timepoint, y=y2, color="red")) +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
    ggtitle(paste(names(hubProbes_female[which(hubProbes_female==i)]), "module hub gene: ", hubProbes_female[which(hubProbes_female == i)])) +
    ylab("Expression") +
    xlab("Zeitgeber Time (ZT)") +
    theme_classic() +
    theme(legend.position="none")
  ggsave(paste("Sine_regression_plot_", hubProbes_female[which(hubProbes_female == i)], "_ME", names(hubProbes_female[which(hubProbes_female==i)]), "_WTF_min30test.pdf", sep=""))
  
}

for(i in hubProbes_wtfemale) {
  y = exp_wtdata.2[,i]
  sinreg_hub = sinreg(exp_wtdata.2$Timepoint, y, plot=FALSE)
  
  y1 = exp_wtdata.2[,i]
  y2 = sinreg_hub[[2]]
  
  ggplot(data=exp_wtdata.2, aes(x=Timepoint, y=y1)) +
    geom_point() +
    geom_line(data=plotpoints, aes(x=Timepoint, y=y2, color="red")) +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
    ggtitle(paste(names(hubProbes_wtfemale[which(hubProbes_wtfemale==i)]), "module hub gene: ", hubSymbols[which(hubProbes_wtfemale == i)])) +
    ylab("Expression") +
    xlab("Zeitgeber Time (ZT)") +
    theme_classic() +
    theme(legend.position="none")
  ggsave(paste("Sine_regression_plot_", hubSymbols[which(hubProbes_wtfemale == i)], "_ME", names(hubProbes_female[which(hubProbes_wtfemale==i)]), "_WTF_test.pdf", sep=""))
  
}

hubSymbols = AnnotationDbi::mapIds(org.Mm.eg.db,
                                   keys = hubProbes_male,
                                   column = "SYMBOL",
                                   keytype = 'ENSEMBL') 

for(i in hubProbes_male) {
  y = exp_wtdata.2[,i]
  sinreg_hub = sinreg(exp_wtdata.2$Timepoint, y, plot=FALSE)
  
  y1 = exp_wtdata.2[,i]
  y2 = sinreg_hub[[2]]
  
  ggplot(data=exp_wtdata.2, aes(x=Timepoint, y=y1)) +
    geom_point() +
    geom_line(data=plotpoints, aes(x=Timepoint, y=y2, color="red")) +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
    ggtitle(paste(names(hubProbes_male[which(hubProbes_male==i)]), "module hub gene: ", hubSymbols[which(hubProbes_male == i)])) +
    ylab("Expression") +
    xlab("Zeitgeber Time (ZT)") +
    theme_classic() +
    theme(legend.position="none")
  ggsave(paste("Sine_regression_plot_", hubSymbols[which(hubProbes_male == i)], "_ME", names(hubProbes_male[which(hubProbes_male==i)]), "_WTM_min30.pdf", sep=""))
  
}

for(i in hubProbes_wtfemale) {
  y = exp_wtdata.2[,i]
  sinreg_hub = sinreg(exp_wtdata.2$Timepoint, y, plot=FALSE)
  
  y1 = exp_wtdata.2[,i]
  y2 = sinreg_hub[[2]]
  
  ggplot(data=exp_wtdata.2, aes(x=Timepoint, y=y1)) +
    geom_point() +
    geom_line(data=plotpoints, aes(x=Timepoint, y=y2, color="red")) +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
    ggtitle(paste(names(hubProbes_wtfemale[which(hubProbes_wtfemale==i)]), "module hub gene: ", hubSymbols[which(hubProbes_wtfemale == i)])) +
    ylab("Expression") +
    xlab("Zeitgeber Time (ZT)") +
    theme_classic() +
    theme(legend.position="none")
  ggsave(paste("Sine_regression_plot_", hubSymbols[which(hubProbes_wtfemale == i)], "_ME", names(hubProbes_female[which(hubProbes_wtfemale==i)]), "_WTF_test.pdf", sep=""))
  
}
#### Heatmaps of R2 with FDRs for associations between Modules and Time using sinreg() ####

FDRs = as.matrix(results.ordered$FDR)
pVal = as.matrix(results.ordered$p.value)
rownames(pVal) = rownames(results.ordered)
colnames(pVal) = "pVal"

R2 = as.matrix(results.ordered$R2)
rownames(R2) = rownames(results.ordered)
colnames(R2) = "R2"

textMatrix = paste(ifelse((signif(FDRs, 1))<0.05, "*", ""), sep="")

tiff("Heatmap_WTM_universal_update.tiff", res=400, height=7, width=2.5, units="in")
map1 = labeledHeatmap(Matrix = pVal,
                      xLabels = colnames(pVal),
                      yLabels = gsub("ME", "", rownames(pVal)),
                      ySymbols = rownames(pVal),
                      colorLabels = FALSE,
                      colors = brewer.pal(9, "RdBu"),
                      textMatrix = textMatrix,
                      setStdMargins = FALSE,
                      invertColors = FALSE,
                      cex.text = 1,
                      cex.lab.x = 0.75,
                      cex.lab.y = 0.65,
                      main=paste("SinReg Rhythmic WTM Modules"))
dev.off()


save(results.ordered, FDRs, pVal, plotpoints, file=glue::glue("sinreg_association_with_gene_modules_WTM_universal.RData"))
rownames(multiExpr$wtdata$data)
rownames(MEs_wt)
#rain
rain_output <- rain(consensusMEs$wtfdata$data, deltat=3, period=24, measure.sequence = c(6, 6, 6, 6, 6, 6, 5, 5), peak.border=c(0.3, 0.7), verbose=FALSE)
rain_output <- rain(wtMEs, deltat=12, period=24, measure.sequence = c(8, 8, 8, 8), peak.border=c(0.3, 0.7), verbose=FALSE)
#rain_output <- rain(multiExpr$wtdata$data, deltat=3, period=24, measure.sequence = c(12, 12, 12, 12, 12, 12, 10, 10), peak.border=c(0.3, 0.7), verbose=FALSE)

rain_output$FDR = p.adjust(rain_output$pVal, method="fdr")

frain.ordered = rain_output[
  with(rain_output, order(rain_output$pVal, decreasing=FALSE)),
]

rownames(frain.ordered) <- gsub(pattern = "ME", replacement = "", x = rownames(frain.ordered), fixed = TRUE)
frain.ordered$Module = rownames(frain.ordered)
frain.ordered$ensemblid = rownames(frain.ordered)
frain.ordered <- frain.ordered[order(frain.ordered$FDR, decreasing = FALSE),]
t <- frain.ordered[order(frain.ordered$phase, decreasing = FALSE),]
t <- t[order(t$Significance, decreasing = TRUE),]
ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
wt_f_genes_rain <- lapply(rownames(mrain.ordered), function(x){
  getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", values = x, mart = ensembl, 
        verbose = TRUE) %>% unlist %>% as.character %>% unique %>% sort %>% paste(collapse = ", ")}) %>% unlist
head(wt_fm_genes_rain)
length(wt_fm_genes_rain)
Wt_counts <- multiExpr$wtdata$data
length(Wt_counts)
Wt_counts <- as.data.frame(t(Wt_counts))
colnames(Wt_counts)
rownames(Wt_counts)
Wt_counts$ensemblid <- rownames(Wt_counts)
Wt_counts$gene_name <- wt_fm_genes_rain

###organize the significance and phase

mrain.ordered <- mrain.ordered[order(mrain.ordered$phase, decreasing = FALSE),]
mrain.ordered <- mrain.ordered[order(mrain.ordered$Significance, decreasing = TRUE),]
mrain.ordered
frain.ordered <- frain.ordered[order(frain.ordered$phase, decreasing = FALSE),]
frain.ordered <- frain.ordered[order(frain.ordered$Significance, decreasing = TRUE),]
frain.ordered
write.csv(frain.ordered, "rain_regresion_results_WTF_universal_genes_update.csv")
openxlsx::write.xlsx(mrain.ordered, file="rain_regresion_results_WTF_genes.xlsx")
openxlsx::write.xlsx(Wt_counts, file="Wtfm_genenames_counts.xlsx")

#rain heatmap
FDRs = as.matrix(mrain.ordered$FDR)
rownames(FDRs) = rownames(mrain.ordered)
colnames(FDRs) = "FDR"

pVal = as.matrix(mrain.ordered$pVal)
rownames(pVal) = rownames(mrain.ordered)
colnames(pVal) = "pVal"

textMatrix = paste(ifelse((signif(FDRs, 1))<0.05, "*", ""), sep="")

tiff("Heatmap_WTM_universal_Rain.tiff", res=400, height=7, width=2.5, units="in")
map1 = labeledHeatmap(Matrix = pVal,
                      xLabels = colnames(pVal),
                      yLabels = gsub("ME", "", rownames(pVal)),
                      ySymbols = rownames(pVal),
                      colorLabels = FALSE,
                      colors = brewer.pal(9, "RdBu"),
                      textMatrix = textMatrix,
                      setStdMargins = FALSE,
                      invertColors = FALSE,
                      cex.text = 1,
                      cex.lab.x = 0.75,
                      cex.lab.y = 0.65,
                      main=paste("RAIN Rhythmic WTM
Modules"))
dev.off()

#Female samples

MEnames = colnames(MEs_cm)
no.MEs = length(MEnames)

traits = pheno$cm$data

alldata = merge(traits, MEs_cm, by="row.names")

alldata = alldata[
  with(alldata, order(alldata$Entrainment)),
]

head(alldata)
write.csv(alldata, "cm_ev.traits.csv")
fitlist = as.list(1:no.MEs)
names(fitlist) <- MEnames

p.values = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(p.values) = MEnames
colnames(p.values) = "p.value"

R2 = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(R2) = MEnames
colnames(R2) = "R2"

period = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(period) = MEnames
colnames(period) = "period"

peak = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(peak) = MEnames
colnames(peak) = "peak"

amplitude = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(amplitude) = MEnames
colnames(amplitude) = "amplitude"

plotpoints = data.frame(matrix(ncol=no.MEs, nrow=nrow(alldata)))
rownames(plotpoints) = alldata$Row.names
colnames(plotpoints) = MEnames

for(i in MEnames){
  
  # print status
  print(paste("Running entity:", i, "which is", which(MEnames==i), "out of", no.MEs))
  
  #create temporary data matrix and model formula
  
  x = alldata$Timepoint
  y = alldata[,i]
  fitlist[[i]] <- sinreg(x, y, plot=FALSE)
  p.values[i,1] <- fitlist[[i]][[1]][6]
  R2[i,1] <- fitlist[[i]][[1]][5]
  period[i,1] <- fitlist[[i]][[1]][3]
  peak[i,1] <- fitlist[[i]][[1]][4]
  amplitude[i,1] <- fitlist[[i]][[1]][2]
  plotpoints[,i] <- fitlist[[i]][[2]]
}

results = cbind(p.values, R2, period, peak, amplitude)

results$FDR = p.adjust(results$p.value, method="fdr")
plotpoints$Timepoint = alldata$Timepoint

for(i in MEnames) {
  y = alldata[,i]
  y1 = plotpoints[,i]
  ggplot(data=alldata, aes(x=Timepoint, y=y)) +
    geom_point() +
    geom_line(data=plotpoints, aes(x=Timepoint, y=y1, color="red")) +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
    ggtitle(paste(i)) +
    ylab("Module EigenValue") +
    xlab("Zeitgeber Time (ZT)") +
    theme_classic() +
    theme(legend.position="none")
  ggsave(paste("Sine_regression_plot_", i, "_females.pdf"))
}

# order results by best fit (highest R2 value) and save results #
results.ordered = results[
  with(results, order(results$R2, decreasing=TRUE)),
]

results.ordered$Module = rownames(results.ordered)

openxlsx::write.xlsx(results.ordered, file="Sine_regresion_results_females.xlsx")

# Checking if hub genes within modules cycle in a similar manner to EigenValues #

exp_femdata.1 = as.data.frame(t(exp_femdata2[,-1]))
colnames(exp_femdata.1) = exp_femdata2$Gene_ID

exp_femdata.2 = merge(traits, exp_femdata.1, by="row.names")

exp_femdata.2 = exp_femdata.2[
  with(exp_femdata.2, order(exp_femdata.2$GenotypeScores)),
]  
#exp_femdata.2 = exp_femdata.2[
#  with(exp_femdata.2, order(exp_femdata.2$Timepoint)),
#]

hubSymbols = AnnotationDbi::mapIds(org.Mm.eg.db,
                                   keys = hubProbes_female,
                                   column = "SYMBOL",
                                   keytype = 'ENSEMBL') 

for(i in hubProbes_female) {
  y = exp_femdata.2[,i]
  sinreg_hub = sinreg(exp_femdata.2$Timepoint, y, plot=FALSE)
  
  y1 = exp_femdata.2[,i]
  y2 = sinreg_hub[[2]]
  
  ggplot(data=exp_femdata.2, aes(x=Timepoint, y=y1)) +
    geom_point() +
    geom_line(data=plotpoints, aes(x=Timepoint, y=y2, color="red")) +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
    ggtitle(paste(names(hubProbes_female[which(hubProbes_female==i)]), "module hub gene: ", hubSymbols[which(hubProbes_female == i)])) +
    ylab("Expression") +
    xlab("Zeitgeber Time (ZT)") +
    theme_classic() +
    theme(legend.position="none")
  ggsave(paste("Sine_regression_plot_", hubSymbols[which(hubProbes_female == i)], "_ME", names(hubProbes_female[which(hubProbes_female==i)]), "_females.pdf", sep=""))
  
}

#### Heatmaps of R2 with FDRs for associations between Modules and Time using sinreg() ####

FDRs = as.matrix(results.ordered$FDR)
pVal = as.matrix(results.ordered$p.value)
rownames(pVal) = rownames(results.ordered)
colnames(pVal) = "pVal"

R2 = as.matrix(results.ordered$R2)
rownames(R2) = rownames(results.ordered)
colnames(R2) = "R2"

textMatrix = paste(ifelse((signif(FDRs, 1))<0.05, "*", ""), sep="")

tiff("Heatmap_females.tiff", res=400, height=7, width=2.5, units="in")
map1 = labeledHeatmap(Matrix = pVal,
                      xLabels = colnames(pVal),
                      yLabels = gsub("ME", "", rownames(pVal)),
                      ySymbols = rownames(pVal),
                      colorLabels = FALSE,
                      colors = brewer.pal(9, "RdBu"),
                      textMatrix = textMatrix,
                      setStdMargins = FALSE,
                      invertColors = FALSE,
                      cex.text = 1,
                      cex.lab.x = 0.75,
                      cex.lab.y = 0.65,
                      main=paste("Circadian Cycling of
Modules in females"))
dev.off()


save(results.ordered, FDRs, pVal, plotpoints, file=glue::glue("sinreg_association_with_gene_modules_females.RData"))

#rain
rain_output <- rain(MEs_female, deltat=3, period=24, measure.sequence = c(6, 6, 6, 6, 6, 6, 5, 5), peak.border=c(0.3, 0.7), verbose=FALSE)
rain_output$FDR = p.adjust(rain_output$pVal, method="fdr")

frain.ordered = rain_output[
  with(rain_output, order(rain_output$FDR, decreasing=FALSE)),
]
rownames(frain.ordered) <- gsub(pattern = "ME", replacement = "", x = rownames(frain.ordered), fixed = TRUE)
frain.ordered$Module = rownames(frain.ordered)

openxlsx::write.xlsx(frain.ordered, file="rain_regresion_results_females.xlsx")

#rain heatmap
FDRs = as.matrix(frain.ordered$FDR)
rownames(FDRs) = rownames(frain.ordered)
colnames(FDRs) = "FDR"

pVal = as.matrix(frain.ordered$pVal)
rownames(pVal) = rownames(frain.ordered)
colnames(pVal) = "pVal"

textMatrix = paste(ifelse((signif(FDRs, 1))<0.05, "*", ""), sep="")

tiff("Heatmap_females.rain.tiff", res=400, height=7, width=2.5, units="in")
map1 = labeledHeatmap(Matrix = pVal,
                      xLabels = colnames(pVal),
                      yLabels = gsub("ME", "", rownames(pVal)),
                      ySymbols = rownames(pVal),
                      colorLabels = FALSE,
                      colors = brewer.pal(9, "RdBu"),
                      textMatrix = textMatrix,
                      setStdMargins = FALSE,
                      invertColors = FALSE,
                      cex.text = 1,
                      cex.lab.x = 0.75,
                      cex.lab.y = 0.65,
                      main=paste("RAIN Rhythmic
Modules in females"))
dev.off()

#name extractor and gene symbol generator
male_modulemem <- read.delim("MALES Probe Module Membership.txt", sep = "\t",
                             header = TRUE, stringsAsFactors = FALSE)
male_ens.mods <- setDT(male_modulemem[,c(11:12)])
#test <- DT[order(Module)]
#test <- DT[Module == "black"]

female_modulemem <- read.delim("FEMALES Probe Module Membership.txt", sep = "\t",
                               header = TRUE, stringsAsFactors = FALSE)
female_ens.mods <- setDT(female_modulemem[,c(11:12)])

msamples.vector <- mrain.ordered$Module
modules_interest = msamples.vector
for (i in modules_interest) {
  
  gene_id <- colnames(exp$maledata$data)[consensusMods$colors == i]
  ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  genes <- getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", 
                 values = gene_id, mart = ensembl) %>% unlist %>% as.character %>% unique %>% sort
  write.csv(genes, (paste(i,"_male_geneslist.csv")), sep = "\t", quote = FALSE)
}

fsamples.vector <- frain.ordered$Module
modules_interest = fsamples.vector
for (i in modules_interest) {
  
  gene_id <- colnames(exp$femdata$data)[fMods$colors == i]
  ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  genes <- getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", 
                 values = gene_id, mart = ensembl) %>% unlist %>% as.character %>% unique %>% sort
  write.csv(genes, (paste(i,"_female_geneslist.csv")), sep = "\t", quote = FALSE)
}

#print full list with genes#
gene_id_male <- colnames(exp$maledata$data)[consensusMods$colors == i]
#ensembl_male <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
genes_interest_male <- getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", 
                             values = probes_interest_male, mart = ensembl_male) %>% unlist %>% as.character %>% unique %>% sort
write.csv(genes_interest_male, file = paste("module_interest_", modules_interest[which(modules_interest == i)], "_males.csv", sep=""))
}
###EnrichR Pathway Analysis###

modules_interest = c("green", "grey", "black", "red", "blue", "brown", "yellow", "pink", "magenta", "turquoise")
lapply(modules_interest, function(module) {
  data = read.csv(glue::glue("{module} _male_geneslist.csv")) 
  
  data %>%
    dplyr::select(x) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2018",
                       "GO_Cellular_Component_2018",
                       "GO_Molecular_Function_2018",
                       "KEGG_2019_Mouse",
                       "Panther_2016",
                       "Reactome_2016",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>% 
    purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %T>%
    #purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05)) %>% 
    #purrr::map(~ dplyr::filter(., stringr::str_detect(Genes, ";"))) %>% 
    openxlsx::write.xlsx(file = glue::glue("Module_{module}_males_enrichr.xlsx")) %>%
    DMRichR::slimGO(tool = "enrichR",
                    annoDb = "org.Mm.eg.db",
                    plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("Module_{module}_males_rrvgo_enrichr.xlsx")) %>%
    DMRichR::GOplot() %>%
    ggplot2::ggsave(glue::glue("Module_{module}_males_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10) 
  
})

## testing
male_modulemem <- read.delim("MALES Gene Module Membership.txt", sep = "\t",
                                header = TRUE, stringsAsFactors = FALSE)
DT <- setDT(male_modulemem[,c(21:22)])
test <- DT[order(Module)]
test <- DT[Module == "black"]

exp <- list(femdata = list(data = as.data.frame(t(exp_femdata[,-c(1)]))),
            maledata = list(data = as.data.frame(t(exp_maledata[,-c(1)]))))
names(exp$femdata$data) = exp_femdata$Gene_ID
rownames(exp$femdata$data) = names(exp_femdata)[-c(1)]
names(exp$maledata$data) = exp_maledata$Gene_ID
rownames(exp$maledata$data) = names(exp_maledata)[-c(1)]
checkSets(exp)

test_list <- list(modules = list(data = as.data.frame(t(male_modulemem[,c(21:22)]))))
names(test_list$modules$data) = male_modulemem$Module
rownames(test_list$modules$data) = names(DT)
names(exp$maledata$data) = exp_maledata$Gene_ID
rownames(exp$maledata$data) = names(exp_maledata)[-c(1)]
checkSets(test_list)

seto
test2 <- male_modulemem[with(male_modulemem, Module > pink), ]
test3 <- DT[, .N, by = Module]
test4 <- DT[Module == "pink", .N, by = Module]
vector_test <- c(mrain.ordered$Module)
names(vector_test)[names(vector_test) == "c(mrain.ordered$Module)"] <- "Module"
testset_modules <- data[match(vector_test, male_modulemem$Module)]
testset_modules <- order_by(vector_test, male_modulemem$Module)
mMods2 <- moduleMembership2$Module
MEs_male <- moduleEigengenes(t(exp_maledata[,-c(1)]), colors = mMods2)$eigengenes
rownames(MEs_male) <- rownames(t(exp_maledata[,-c(1)]))
mMEs <- list(male = list(data = MEs_male))
mMEs <- orderMEs(mMEs)
order_by()
female_modulemem <- read.delim("FEMALES Gene Module Membership.txt", sep = "\t",
                                header = TRUE, stringsAsFactors = FALSE)
fMods2 <- moduleMembership2$Module
MEs_female <- moduleEigengenes(t(exp_femdata[,-c(1)]), colors = fMods2)$eigengenes
rownames(MEs_female) <- rownames(t(exp_femdata[,-c(1)]))
fMEs <- list(male = list(data = MEs_female))
fMEs <- orderMEs(fMEs)

rm(moduleMembership2)

###TRIAL

traits = as.data.frame(pheno$wtf$data)
traits = as.data.frame(pheno$wtm$data)


exp_wtdata.1 = as.data.frame(t(wtf_norm))
#colnames(exp_wtdata.1) = wtf$Gene_ID
rownames(exp_wtdata.1)
rownames(traits)
exp_wtdata.2 = merge(traits, exp_wtdata.1, by="row.names")

exp_wtdata.2 = exp_wtdata.2[
  with(exp_wtdata.2, order(exp_wtdata.2$Timepoint)),
]

genes_interest = c("Arntl", "Clock", "Cry1", "Cry2", "Per1","Per2","Per3", "Meg3", "Snhg14", "Xist", "Snrpn", "Cdkl5")
genes_interest = c("Mecp2")
#genes_interest = "Arntl"
for(i in genes_interest) {
  y.0 = exp_wtdata.2[,i]
  sinreg_interest = sinreg(exp_wtdata.2$Timepoint, y.0, plot=FALSE)
  #print(mean(exp_maledata.2$ENSMUSG00000028957))
  y.1 = exp_wtdata.2[,i]
  #y.2 = sinreg_interest[[2]]
  
  ggplot(data=exp_wtdata.2, aes(x=Timepoint, y=y.1)) +
    geom_point() +
    #geom_line(aes(x=exp_maledata.2$Timepoint, y=y.2, color="red")) +
    geom_smooth(method="loess") +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
    ggtitle(paste("Wild-Type Female", "Gene: ", genes_interest[which(genes_interest == i)], "|", "Module: ", MM_female$Module[which(MM_female$Probe == i)])) +
    ylab("Expression (logcpm)") +
    xlab("Zeitgeber Time (ZT)") +
    theme_classic() +
    theme(legend.position="none", plot.title = element_text(size = 20))
  ggsave(paste("loess_plot_", genes_interest[which(genes_interest == i)], "_WTF_test2.pdf", sep=""))
  
}
test_wtf <- t(wtf_norm)
test_wtf <- as.data.frame(t(wtf_norm))

gpt <- ggplot(test_wtf, aes(x=rownames(test_wtf), y = Xist)) +
  geom_point() +
  #geom_line(aes(x=exp_maledata.2$Timepoint, y=y.2, color="red")) +
  geom_smooth(method="loess") 
gpt

ggplot(data=exp_wtdata.2, aes(x=Timepoint, y=Mecp2)) +
  geom_point() +
  #geom_line(aes(x=exp_maledata.2$Timepoint, y=y.2, color="red")) +
  geom_smooth(method="loess") +
  scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
  #ggtitle(paste("gene: ", genes_interest[which(genes_interest)])) +
  ylab("Expression") +
  xlab("Zeitgeber Time (ZT)") +
  theme_classic() +
  theme(legend.position="none")
