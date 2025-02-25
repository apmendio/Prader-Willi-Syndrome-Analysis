library(dplyr)
library(GeneOverlap)
library(VennDiagram)

###############
## load data ##
###############
pws <- read.csv("/Users/osman/Downloads/PWS_signifKEGG.csv")
pws$X <- NULL
pws <- unique(pws)
rtt <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis/filtered_enrich_df.csv")
rtt <- unique(rtt)

merged_df <- merge(pws, rtt, by = 'Term')
# Remove duplicate rows based on 'Term'
matching_df <- merged_df[!duplicated(merged_df$Term), ]

# Print the result
print(matching_df)

#############################
## Gene overlap statistics ##
#############################
go.obj <- newGeneOverlap(pws$Term, 
                         rtt$Term, 
                         genome.size=904)
go.obj
go.obj <- testGeneOverlap(go.obj)
go.obj
print(go.obj)

##################
## Venn Diagram ##
##################
pdf(glue("/Users/osman/Downloads/PWS_RETT_overlap.pdf"))
temp <- venn.diagram(
  x = list(
    PWS = pws$Term,
    RETT = rtt$Term
  ),
  category.names = c("Sig PWS KEGG Terms", "Sig RTT KEGG Terms"),
  main = 'Overlap of significant PWS and RTT KEGG Terms',
  #filename = glue("{base_path}/broad_group_analysis/venn_glutamatergic.pdf"),
  filename = NULL,
  col = c('#F79120', '#372A82'), 
  fill = c('#F79119', '#372A81'),
  cat.cex = 1.2,
  cat.fontface = "bold",
  euler.d = TRUE,
  disable.logging = TRUE,
  hyper.test = TRUE
)
grid.draw(temp)
dev.off()

