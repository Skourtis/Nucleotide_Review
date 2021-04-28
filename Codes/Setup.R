#### Setup Project ####
install.packages("pacman")
pacman::p_load(here,tidyverse,openxlsx,
               pheatmap,visNetwork,BiocManager,
               KEGGgraph,naniar, piggyback)
# BiocManager::install("KEGGgraph")

###Creates all files needed for project
folders <- c("Codes","Datasets","Output","Datasets/Raw","Datasets/Processed")
purrr::walk(.x = folders,~dir.create(here::here(.x)))
#piggyback::pb_new_release()

piggyback::pb_track(c("Datasets/Raw/*.txt",
                      "Datasets/Raw/*.dat",
                      "Datasets/Raw/*.zip",
                      "Datasets/Raw/*.csv",
                      "Datasets/Raw/*.tab",
                      "Datasets/Raw/*.RData",
                      "Datasets/Raw/*.xlsx",
                      "Datasets/Processed/*.RData",
                      "Output/*.pdf",
                      "app/*.RData"))%>%
    piggyback::pb_upload(repo = "Skourtis/Nucleotide_Review",
                         overwrite = "use_timestamps")
