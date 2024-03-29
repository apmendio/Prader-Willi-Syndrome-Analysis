### Pathway Analaysis for Select Modules ###

setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/WGCNAoptimization/CircadianModules/Males/GOterms")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/WGCNA_Final/circadianmodules/GoTerms")
packages <- c("edgeR", "tidyverse", "magrittr", "RColorBrewer", "org.Mm.eg.db", "AnnotationDbi",
              "enrichR", "openxlsx", "gt", "glue", "DMRichR")
BiocManager::install("enrichR")
BiocManager::install("DMRichR")
BiocManager::install("clusterProfiler")
BiocManager::install("glue")
install.packages("DMRichR")
library(enrichR)
library(DMRichR)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(magrittr)
library(glue)

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
modules_interest = c("yellow", "pink", "black")
lapply(modules_interestrm, function(module) {
  data = read.csv(glue::glue("{module}_moduleRM12.csv")) 
  
  data %>%
    dplyr::select(x) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2023",
                       "GO_Cellular_Component_2023",
                       "GO_Molecular_Function_2023",
                       "KEGG_2019_Mouse",
                       "Panther_2016",
                       "Reactome_2022",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>% 
    purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %T>%
    #purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05)) %>% 
    #purrr::map(~ dplyr::filter(., stringr::str_detect(Genes, ";"))) %>% 
    openxlsx::write.xlsx(file = glue::glue("Module_{module}_RM12_enrichr.xlsx")) 
})

test_modules_rw3 <- lapply(modules_interest, function(module) {
  data = read.csv(glue::glue("{module}_moduleRM12.csv")) 

  data %>%
    dplyr::select(x) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2023",
                       "GO_Cellular_Component_2023",
                       "GO_Molecular_Function_2023",
                       "KEGG_2019_Mouse",
                       "Panther_2016",
                       "Reactome_2022",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>% 
    purrr::set_names(names(.) %>% stringr::str_trunc(40, ellipsis = "")) %T>%
    purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05))
    #new_name <- paste0("test_RW",module) %>%
    #assign(new_name, test_modules)
    #purrr::map(~ dplyr::filter(., stringr::str_detect(Genes, ";"))) %>% 
    #openxlsx::write.xlsx(file = glue::glue("Module_{module}_RM12_enrichr.xlsx")) 
})

test_rw_test <- lapply(modules_interest, function(module) {
  glue('RW_{module}_KG') <- glue('test_modules_rw${module}$KEGG_2019_Mouse')
})
glue('RW_{modules_interest}_KG') <- glue('test_modules_rw${modules_interest}$KEGG_2019_Mouse')

test_rw_test <- lapply(modules_interest, function(module) {
  assign(paste0("RW_",module, "_KG"), paste0("test_modules_rw$",module, "$KEGG_2019_Mouse"))
})

names(test_modules_rw2) <- modules_interestrm
test_modules_rw$grey$KEGG_2019_Mouse
test_kg <- lapply(modules_interest, function(module) {
  data = glue::glue("test_rw_test${module}$KEGG_2019_Mouse")
})
names(test_kg) <- modules_interest
names(test_kg)
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/separate_normalization/GOTermsUpdatedCorrected")

data = read.csv(glue::glue("yellow_moduleRM12.csv")) 
RM_yellow = read.csv(glue::glue("yellow_moduleRM12.csv")) 
RM_black = read.csv(glue::glue("black_moduleRM12.csv"))
RM_pink = read.csv(glue::glue("pink_moduleRM12.csv"))
RM_red = read.csv(glue::glue("red_moduleRM12.csv"))
Red_RW <- RM_red %>%
    dplyr::select(x) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2023",
                       "GO_Cellular_Component_2023",
                       "GO_Molecular_Function_2023",
                       "KEGG_2019_Mouse",
                       "Panther_2016",
                       "Reactome_2022",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>% 
    purrr::set_names(names(.) %>% stringr::str_trunc(40, ellipsis = "")) %T>%
    purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05))
#%>% 
    #purrr::map(~ dplyr::filter(., stringr::str_detect(Genes, ";")))
    #openxlsx::write.xlsx(file = glue::glue("Module_{module}_RM12_enrichr.xlsx")) 

RW_Yellow_BP <- Yellow_RW$GO_Biological_Process_2023
RW_Yellow_BP$Sample <- "RW_yellow"
RW_Black_BP <- Black_RW$GO_Biological_Process_2023 
RW_Black_BP$Sample <- "RW_black"
test_RW_modules <- rbind(RW_Yellow_BP, RW_Black_BP)
head(test_RW_modules)
test_RW_modules2 <- filter(test_RW_modules, Adjusted.P.value <= 0.05)

RW_Yellow_KG <- Yellow_RW$KEGG_2019_Mouse
RW_Yellow_KG$Sample <- "RW_yellow"
RW_Pink_KG <- Pink_RW$KEGG_2019_Mouse
RW_Pink_KG$Sample <- "RW_pink"
RW_Red_KG <- Red_RW$KEGG_2019_Mouse
RW_Red_KG$Sample <- "RW_red"
test_RW_modules <- rbind(RW_Yellow_KG, RW_Pink_KG, RW_Red_KG)
head(test_RW_modules)
test_RW_modules2 <- filter(test_RW_modules, Adjusted.P.value <= 0.05)
#gobp1 <- test9$GO_Biological_Process_2023[test9$GO_Biological_Process_2023$Adjusted.P.value < 0.05,]
ggplot(test_RW_modules2, aes(x = Sample, y = Term, color = Adjusted.P.value, size = Odds.Ratio)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  ylab("") +
  xlab("Modules") +
  ggtitle("GO Terms")


####
    DMRichR::slimGO(tool = "enrichR",
                    annoDb = "org.Mm.eg.db",
                    plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("Module_{module}_WT121_rrvgo_enrichr.xlsx")) %>%
    DMRichR::GOplot() %>%
    ggplot2::ggsave(glue::glue("Module_{module}_WT121_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10) 
  


test2 <- read.csv("module_blue.csv")
