library(ComplexHeatmap)
library(circlize)
###set up cell function to create stars on heatmap
set = 3

padj <- moduleTraitPvalue[[set]] #creating object with p value matrix

cell_fun = function(j, i, x, y, width, height, fill) {
  if(is.na(as.matrix(padj[,])[i, j]) ){
    grid.text(sprintf("%s", ''), x, y, gp = gpar(fontsize = 12))
  }else if(as.matrix(padj[,])[i, j] <0.001 ){
    grid.text(sprintf("%s", '***'), x, y, gp = gpar(fontsize = 12))
  }else if(as.matrix(padj[,])[i, j] <0.01) {
    grid.text(sprintf("%s", '**'), x, y, gp = gpar(fontsize = 12))
  }else if(as.matrix(padj[,])[i, j] <0.05) {
    grid.text(sprintf("%s", '*'), x, y, gp = gpar(fontsize = 12 ))
  }
}

hm <- Heatmap(moduleTraitCor[[set]], cell_fun = cell_fun) # creates heatmap object

hm
hm <- Heatmap(moduleTraitCor[[3]], cell_fun = cell_fun, cluster_columns = FALSE) # removes clustering by columns

hm <- Heatmap(moduleTraitCor[[3]], cell_fun = cell_fun, cluster_columns = FALSE, name = "Spearman Correlation") # adds name to color scale

labels = gsub("ME", "", rownames(moduleTraitCor[[set]])) # substitutes rownames 

## changes scale range
hm <- Heatmap(moduleTraitCor[[set]], cell_fun = cell_fun, cluster_columns = FALSE, column_title = "Spearman Correlation RW Male", column_title_gp = gpar(fontsize = 20, fontface = "bold"), name = "Spearman Correlation", col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

hm + rowAnnotation(labels = anno_text(labels, which = "row")) # removes ME

hm <- Heatmap(moduleTraitCor[[3]], cell_fun = cell_fun, cluster_columns = FALSE, name = "Spearman Correlation", col = test)
