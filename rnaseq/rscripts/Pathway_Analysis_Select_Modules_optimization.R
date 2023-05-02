### Pathway Analaysis for Select Modules ###

setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/WGCNAoptimization/CircadianModules/Males/GOterms")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/WGCNA_Final/circadianmodules/GoTerms")
packages <- c("edgeR", "tidyverse", "magrittr", "RColorBrewer", "org.Mm.eg.db", "AnnotationDbi",
              "enrichR", "openxlsx", "gt", "glue", "DMRichR")
BiocManager::install("enrichR")
BiocManager::install("DMRichR")
install.packages("DMRichR")
library(enrichR)
library(DMRichR)

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE)
BiocManager::install("ben-laufer/DMRichR")
if(R.Version()$major < 4)
  install.packages("dmrseq", repos = "https://bioconductor.org/packages/3.12/bioc")

enrichR:::.onAttach() # Needed or else "EnrichR website not responding"
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

### Loading Data ###

modules_interest = c("pink", "purple", "red", "salmon", "darkturquoise")
modules_interest = c("yellow")

modules_interest = read.csv("module.distribution.csv")
modules_interest <- modules_interest$Module
modules_interest = c("grey")
modules_interest <- c("tan", "salmon", "cyan", "midnightblue", "lightcyan", "lightyellow", "royalblue", "darkred")
modules_interest <- c("greenyellow")
modules_interest
getwd()

modules_interest = module.dist$Module
modules_interest <- c("darkgrey")
#DMRichR is not optimized for 2021GO

lapply(modules_interest, function(module) {
  data = read.csv(glue::glue("{module}_module.csv")) 
  
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
    openxlsx::write.xlsx(file = glue::glue("Module_{module}_enrichr.xlsx")) %>%
    DMRichR::slimGO(tool = "enrichR",
                    annoDb = "org.Mm.eg.db",
                    plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("Module_{module}_rrvgo_enrichr.xlsx")) %>%
    DMRichR::GOplot() %>%
    ggplot2::ggsave(glue::glue("Module_{module}_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10) 
  
})

test = read.csv("grey_module.csv")
test = read.csv("module_darkgreen.csv")

test %>% 
  dplyr::select(x) %>%
  purrr::flatten() %>%
  enrichR::enrichr(c("GO_Biological_Process_2021",
                     "GO_Cellular_Component_2021",
                     "GO_Molecular_Function_2021",
                     "KEGG_2019_Mouse",
                     "Panther_2016",
                     "Reactome_2016",
                     "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>% 
  purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = ""))
set_name