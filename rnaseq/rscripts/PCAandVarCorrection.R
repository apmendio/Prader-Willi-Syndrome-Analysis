library(ggfortify)
df <- iris[1:4]
df
pca_res <- prcomp(df, scale. = TRUE)
autoplot(pca_res)
autoplot(pca_res, data = iris, colour = 'Species')
df2 <- iris
head(df2)

head(All_merge)
rownames(All_merge) <- All_merge$Row.names
All_merge <- All_merge[,-c(1)]
df <- as.data.frame(t(All_merge))
pca_res <- prcomp(All_merge, scale. = TRUE)
autoplot(pca_res)
df <- All_merge[,c(1:46, 93:124, 156:187)]
df <- All_merge[,-c(151, 159)]
df <- All_merge3
head(df)
pca_res <- prcomp(t(df), scale. = TRUE)
autoplot(pca_res)

metadata <- read.csv("metadata_exp.csv")
metadata <- read.csv("metadata_exp2.csv")
metadata <- metadata[-c(151,159),]
pca_res <- prcomp(t(df), scale. = TRUE)
autoplot(pca_res, data = metadata, colour = 'Batch')
autoplot(pca_res, data = metadata, colour = 'Experiment')
autoplot(pca_res, data = metadata, colour = 'Experiment', label = TRUE, label.size = 5)
autoplot(pca_res, data = metadata, colour = 'Sex')
dim(All_merge)
dim(metadata)
## Blythe provided code to produce residuals ##
resids2 <- t(apply(t(df), 1, function(x)resid(lm(x ~ Experiment, data = metadata)))) # normalized counts goes into df
colnames(resids2) <- colnames(t(df))
resids2 <- resids2[,-c(151,159)]
pca_res2 <- prcomp(t(resids2), scale. = TRUE)
autoplot(pca_res2, data = metadata, colour = 'Experiment')
autoplot(pca_res2, data = metadata, colour = 'Experiment', label = TRUE, label.size = 5)
autoplot(pca_res2, data = metadata, colour = 'Genotype')
autoplot(pca_res2, data = metadata, colour = 'Batch')
resids <- t(apply(resids, 1, function(x)resid(lm(x ~ Experiment, data = metadata))))
colnames(resids) <- colnames(df)
pca_res3 <- prcomp(t(resids), scale. = TRUE)
autoplot(pca_res3, data = metadata, colour = 'Experiment')
autoplot(pca_res3, data = metadata, colour = 'Batch')
autoplot(pca_res3, data = metadata, colour = 'Genotype')
### separate out data
head(resids2)
dim(resids2)
wtf_norm <- resids2[,c(1:46)]
wtm_norm <- resids2[,c(47:92)]
cf_norm <- resids2[,c(93:124)]
cm_norm <- resids2[,c(125:154)]
rf_norm <- resids2[,c(155:185)]
rm_norm <- resids2[,c(186:217)]
colnames(rf_norm)

###
BiocManager::install("bladderbatch")
library(bladderbatch)
data(bladderdata)

dat <- bladderEset[1:50,]
head(dat)
pheno = pData(dat)
edata = exprs(dat)
head(edata)
batch = pheno$batch
head(batch)

mod = model.matrix(~as.factor(cancer), data=pheno)
pca_res2 <- prcomp(t(edata), scale. = TRUE)
autoplot(pca_res2)
pheno$batch <- as.character(pheno$batch)
autoplot(pca_res2, data = pheno, colour = 'batch')
# parametric adjustment
combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
pca_res3 <- prcomp(t(combat_edata1), scale. = TRUE)
autoplot(pca_res3, data = pheno, colour = 'batch')
# non-parametric adjustment, mean-only version
combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)
pca_res3 <- prcomp(t(combat_edata2), scale. = TRUE)
autoplot(pca_res3, data = pheno, colour = 'batch')
# reference-batch version, with covariates
combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3)
pca_res3 <- prcomp(t(combat_edata3), scale. = TRUE)
autoplot(pca_res3, data = pheno, colour = 'batch')
