####Importing Pathways####
pacman::p_load(here,renv,tidyverse,openxlsx,
               pheatmap,visNetwork,BiocManager,
               KEGGgraph, matrixStats)
HUMAN_9606_idmapping <- readr::read_tsv(here::here("Datasets","Raw","HUMAN_9606_idmapping.dat"),
                                        col_names = FALSE)
HUMAN_9606_idmapping_hsa <- HUMAN_9606_idmapping[str_detect(HUMAN_9606_idmapping$X3, "hsa:"),-2] %>% 
    left_join(HUMAN_9606_idmapping[HUMAN_9606_idmapping$X2 == "Gene_Name", -2], by = "X1") %>% 
    na.omit() %>% dplyr::select(-X1) %>% set_names(c("KEGG_ID", "Gene_names"))

#### Protein Families ####

Protein_Families <- read_tsv(here::here("Datasets","Raw","uniprot-organism__Homo+sapiens+(Human)+[9606]_.tab")) %>% 
    as.data.frame() %>% na.omit() %>%
    left_join(HUMAN_9606_idmapping %>% subset(X2 =="Gene_Name") %>% dplyr::select(-X2),
              by= c("Entry" = "X1")) %>% 
    mutate(Gene_names = X3)
           #Protein_families = str_match(Protein_families, pattern = "^([[:print:]]*?),")[,2]) #%>% na.omit()

###Pertubation score####
Achilles <- read.csv(here::here("Datasets","Raw","Achilles_gene_effect.csv"))
sample_info <- read.csv(here::here("Datasets","Raw","sample_info.csv")) %>% 
    dplyr::select(DepMap_ID,lineage)

##Pertubation score per tissue and gene##
Achilles_per_tissue <- Achilles %>% 
    left_join(sample_info) %>% 
    group_split(lineage, .keep = TRUE) %>%
    `names<-`({.} %>% map(~ .x$lineage[1]) %>% unlist()) %>%
    ## If you want to discard the grouping variable, do the following step as well
    map(~ .x %>% select(-c(lineage,DepMap_ID)) %>% t())

Pertubation_per_tissue <- purrr::map(.x =Achilles_per_tissue,
                                     ~data.frame(Gene_name = trimws(str_remove(rownames(.x),
                                                                               pattern = "\\.\\.[:graph:]*\\.$")),
                                     Color = rowQuantiles(.x, na.rm = T)[,2]) %>% 
                                         mutate(width = 1/rowSds(.x, na.rm = T)))
##pertubation for all linages##
Achilles_t <- t(Achilles %>% dplyr::select(-DepMap_ID))
Pertubation_score <- data.frame(Gene_name = trimws(str_remove(rownames(Achilles_t),pattern = "\\.\\.[:graph:]*\\.$")),
                                Color = rowQuantiles(Achilles_t, na.rm = T)[,2]) %>% 
    mutate(width = 1/rowSds(Achilles_t, na.rm = T))
Pertubation_per_tissue[["all"]] <- Pertubation_score
### Retrieving and Merging KEGG Folate ####
tmp <- tempfile()

pathways <-  c("hsa00670", "hsa00240", "hsa00230","hsa00270","hsa00790")
#Downloads the genes and compounds in each reaction from KEGG#
Retrieve_genes<- function(pathway){
    retrieveKGML(pathway, organism="hsa", destfile=tmp, method="auto", quiet=TRUE)
    Pathway_genes <- KEGGgraph::parseKGML2Graph(tmp)
    for (i in 1:length(Pathway_genes@nodeData@defaults[["KEGGNode"]][["nodes"]])){
        df <- data.frame(Gene_id = Pathway_genes@nodeData@defaults[["KEGGNode"]][["nodes"]][[i]]@entryID,
                         Reaction = Pathway_genes@nodeData@defaults[["KEGGNode"]][["nodes"]][[i]]@reaction,
                         pathway = pathway,
                         stringsAsFactors = FALSE)
        genes_reactions <- rbind(genes_reactions,df)
        
    }
    genes_reactions
    
    
}
Retrieve_compounds <- function(pathway){
    retrieveKGML(pathway, organism="hsa", destfile=tmp, method="auto", quiet=TRUE)
    Pathway_metabolites <- KEGGgraph::parseKGML(tmp)

    
    compound_reactions <- data.frame(Reaction = NULL,
                                     Substrate = NULL,
                                     Product = NULL,
                                     Direction = NULL,
                                     pathway = NULL,
                                     stringsAsFactors = F)
    
    for (i in 1:length(Pathway_metabolites@reactions)){
        df <- data.frame(Reaction = Pathway_metabolites@reactions[[i]]@name,
                         Substrate = Pathway_metabolites@reactions[[i]]@substrateName,
                         Product = Pathway_metabolites@reactions[[i]]@productName,
                         Direction = Pathway_metabolites@reactions[[i]]@type,
                         pathway = pathway,
                         stringsAsFactors = FALSE)
        compound_reactions <- rbind(compound_reactions,df)
        
    }
    compound_reactions
    
    
}
genes_reactions <- purrr::map_df(pathways,Retrieve_genes)
compound_reactions <-  purrr::map_df(pathways,Retrieve_compounds)

####
Enzymes_Metabolites <- read_tsv(here::here("Datasets","Raw","all_unique_KEGG_metabolites_mapped_KEGG.tsv"))
Metabolites_mapped <- read.csv(here::here("Datasets","Raw","MetaboAnalyst_mapping.csv")) %>%  na.omit() 
compound_reactions <- separate_rows(compound_reactions, Reaction, sep =" ") %>% 
    mutate(Substrate = str_remove(Substrate,"cpd:"),
           Product = str_remove(Product,"cpd:"))

#Adds Reactions which were missing from KEGGPAth
genes_reactions <- rbind(genes_reactions, data.frame(Gene_id = c("hsa:1719", "hsa:200895","hsa:4522", "hsa:441024", "hsa:10797","hsa:10797","hsa:4522","hsa:10841","hsa:10841"),
                                                     Reaction = c("rn:R00937 rn:R02235", "rn:R00937 rn:R02235","rn:R00943", "rn:R01655", "rn:R01655", "rn:R01218", "rn:R01218","rn:R02302","rn:R03189"),
                                                     pathway = "hsa00670")) %>% 
    separate_rows(Reaction, sep =" ")

Pathways <- full_join(genes_reactions,compound_reactions, by = c("Reaction", "pathway")) %>% 
    distinct() %>% 
    left_join(HUMAN_9606_idmapping_hsa, by= c("Gene_id" = "KEGG_ID")) %>%
    mutate(From = MetaboSignal::MS_changeNames(paste0("cpd:", Substrate),"hsa"),
           To = MetaboSignal::MS_changeNames(paste0("cpd:", Product),"hsa")) %>% 
    distinct(Gene_id, pathway, Substrate, Product, Direction, Gene_names, .keep_all = T) %>% 
    mutate(Gene.names = if_else(is.na(Gene_names),Gene_id,Gene_names),
           From = if_else(is.na(From),Substrate,From),
           To = if_else(is.na(To),Product,To))
####Network database Creating ####
nodes <- data.frame(id =c(Pathways$From,Pathways$To),
                    label =c(Pathways$From,Pathways$To),
                    KEGG= c(Pathways$Substrate,Pathways$Product)) %>% 
    distinct() %>% 
    mutate(title =  paste0("<p><b>",KEGG ,"</b><br></p>"))
mypal <- colorRampPalette( c( "red","purple", "#0080ff" ) )(nrow(edges_enzymes))
map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}
edges_enzymes_tissues <- purrr::map(.x =Pertubation_per_tissue,
                                    ~Pathways %>% 
    dplyr::select(c( "From","To", "Gene_names", "Direction"))%>% 
    dplyr::rename(label = `Gene_names`, from = From, to = To) %>% 
    left_join(.x, by = c("label" = "Gene_name")) %>% 
    mutate(font.color = "green", 
           font.size = 10,
           width = if_else(is.na(width),min(width, na.rm =  T),width),
           Color = if_else(is.na(Color),mean(Color, na.rm =  T),Color),
           label = if_else(is.na(label)," ",label),
           color = if_else(is.na(Color),"#a6a6a6", map2color(Color,mypal)),
           arrows.middle = if_else(Direction == "irreversible", T, F),
           dashes = if_else(label == " ", T, F) ) %>% 
    left_join(Protein_Families %>% dplyr::select(Gene_names,Protein_families) %>% 
                  distinct, by = c("label" = "Gene_names")) %>% 
    add_count(from,to,Protein_families) %>%  #No. of enzymes per family per reactions will be used for grouping
    mutate(Gene = label,
           title = if_else(label == " ",NA_character_, paste0(Gene, " :: Pertubation : ", Color)),
           label = if_else(n>2,Protein_families,label),
           label = str_replace(label,"/","/\n"),
           arrows.middle = if_else(label ==  Protein_families,arrows.middle,F)) %>% #Hide names of enzymes if multiple in the same family
    arrange(Color) %>% 
    subset(!(label == " " & duplicated(.[,c("from","to")]))) %>% 
    dplyr::group_by(Protein_families, from, to) %>% 
    dplyr::mutate(title = paste0(title, collapse = "<br>")) %>% 
    ungroup() %>% 
    distinct(from,to,label,.keep_all = T)) %>% set_names(names(Pertubation_per_tissue))

#nodes <- nodes %>% mutate(id = MetaboSignal::MS_changeNames(paste0("cpd:", id),"hsa"))
visNetwork(nodes,edges_enzymes_tissues[[1]], 
           height = 1000, width = "100%",
           main = paste0("Metabolic Vulnerabolities - ",names(edges_enzymes_tissues[1]))) %>% 
    visIgraphLayout(smooth = T)%>% 
    visNodes(
        shape = "box",
        shadow = list(enabled = TRUE, size = 10)
    ) %>%
    visEdges(smooth = list(type = "dynamic"), color = list(highlight = "#C62F4B",opacity = 0.35, border = "black"), 
             font = list(align = "middle"), arrows = list(middle = list(scaleFactor = 0.1))) %>%
    visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
               selectedBy = "KEGG") %>%
    visPhysics(stabilization = FALSE,solver = "forceAtlas2Based", 
               forceAtlas2Based = list(gravitationalConstant = -150))



#####Metabolite Heatmap####
cell_line_data=openxlsx::read.xlsx(here::here("Datasets","Raw","41591_2019_404_MOESM2_ESM.xlsx"),
                          sheet=3) %>% remove_rownames()
cell_line_annotation=read.delim(file=here::here("Datasets","Raw","mod_tissue_origins.txt"),
                                header=F)
tissue_annotation=str_replace(cell_line_annotation[,1],"HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","BLOOD")

pat_others=c("AUTONOMIC_GANGLIA|BILIARY_TRACT|BONE|KIDNEY|PLEURA|PROSTATE|SALIVARY_GLAND|SOFT_TISSUE|THYROID")
tissue_annotation=str_replace_all(tissue_annotation,pat_others,"OTHERS")

mod_cell_line_data=data.frame(column_to_rownames(cell_line_data, var = "X1"))
colnames(mod_cell_line_data)=str_replace(colnames(mod_cell_line_data),"^X","")
colnames(mod_cell_line_data)=str_replace(colnames(mod_cell_line_data),"\\.","-")

mod_tissue_annotation=data.frame(tissue_annotation)
row.names(mod_tissue_annotation) <- rownames(mod_cell_line_data)

#Colour palette for the tissue annotation#
n <- 15
qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(RColorBrewer ::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#ann_colors=sample(col_vector, n)
ann_colors=col_vector[50:(50+n-1)]
names(ann_colors) <- unique(mod_tissue_annotation$tissue_annotation)
ann_colors=list(tissue_annotation=ann_colors)

folate_metabolites=c("glycine", "serine", "adenosine", "glutamate", "NAD", "NADP", "methionine", "betaine", "cystathionine", "choline", "dimethylglycine", "homocysteine", "5-adenosylhomocysteine", "putrescine", "sarcosine", "alpha-ketoglutarate" )
folate_cell_line_data=subset(mod_cell_line_data,select= colnames(mod_cell_line_data) %in% folate_metabolites)


p <- pheatmap(folate_cell_line_data
              ,color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name ="RdYlBu")))(50),
              scale="column",clustering_distance_rows="correlation",clustering_distance_cols = "correlation",
              show_rownames =F,show_colnames =T,annotation_row=mod_tissue_annotation,
              annotation_colors=ann_colors,main="Metabolites and Tissues in the Folate Pathway")

#### Essentiality Figure #####
###Producing Essentiality for Figure 
achilles_pivoted <- pivot_longer(Achilles, -1, names_to = "Gene", values_to = "KO_Sensitivity") %>% 
    mutate(Gene = trimws(str_remove(.$Gene,pattern = "\\.\\.[:graph:]*\\.$"))) %>%
    na.omit() %>% 
    mutate(KO_Sensitivity = tanh(KO_Sensitivity))
stephan_proteins <- c("MTHFD1","MTHFD1L"
                      ,"MTHFD2" ,"MTHFD2L", 'MAT1A',
                      "GNMT", "MTHFR", "DHFR",
                      "MTR" ,"SHMT2" , "BRD4" ,"AHCY"
                      ,"NMT1","BHMT", "BHMT2", "SHMT1")

for(i in stephan_proteins){

    dataset <- achilles_pivoted[achilles_pivoted$Gene == i,]
    file_name <- here::here("Output",paste0(i,".pdf", sep = ""))
    pdf(file_name)
    plot <- ggplot(dataset,aes(x = DepMap_ID, y = Gene, fill = KO_Sensitivity))+
        geom_tile()+  scale_fill_gradientn(colours = c("#ba4747", "white", "#60A8E2"), values = c(0,abs(min(achilles_pivoted$KO_Sensitivity))/(max(achilles_pivoted$KO_Sensitivity)+abs(min(achilles_pivoted$KO_Sensitivity))),1), limits = c(min(achilles_pivoted$KO_Sensitivity), max(achilles_pivoted$KO_Sensitivity)))+
        #scale_fill_distiller(palette = "RdBu", )+
        scale_x_discrete(limits=(dataset$DepMap_ID)[order(dataset$KO_Sensitivity)])  +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
    plot(plot)
    dev.off()
}


#### Proteins and Enzymes ####
#### with Enzymes
CCLE_proteins <- openxlsx::read.xlsx(here::here("Datasets","Raw","CCLE_table_2.xlsx"), sheet =2 )
Metabolic_ids <- openxlsx::read.xlsx(here::here("Datasets","Raw","uniprot_hsa01100.xlsx"))[,3]
CCLE_enzymes <- subset(CCLE_proteins, Uniprot %in% Metabolic_ids)
CCLE_enzymes <- CCLE_enzymes[,!str_detect(colnames(CCLE_enzymes),"Peptides")]

deep_map <- read.csv(here::here("Datasets","Raw","D2_combined_gene_dep_scores.csv"))
colnames(deep_map) <- str_remove_all(colnames(deep_map), "^X")
deep_map_pivoted <- pivot_longer(deep_map, -1, names_to = "Cell_line",
                                 values_to = "Sensitivity") %>%
    mutate(Gene_name = trimws(str_remove_all(.$Gene_name, "\\([:graph:]*\\)"))) %>% 
    na.omit() %>% mutate(Cell_type =  str_match(.$Cell_line, "_([:graph:]*)")[,2])

CCLE_Mut_profile <- read_tsv(here::here("Datasets","Raw","gene_set_library_crisp.gmt"), col_names  = F)
#https://cancer.sanger.ac.uk/census
Cosmic_Tiers_mutations <- read.csv(here::here("Datasets","Raw","Census_allTue_Jul_7_15_49_21_2020.csv"))%>%
    .[.$Tier == "1",]


`%not_in%` <- purrr::negate(`%in%`)

CCLE_Mut_profile[,-c(1:2)] <- CCLE_Mut_profile[,-c(1:2)] %>% replace_with_na_all(condition = ~.x %not_in% Cosmic_Tiers_mutations$Gene.Symbol)
CCLE_Mut_profile <- data.frame(Cell_line = paste(CCLE_Mut_profile$X1, toupper(CCLE_Mut_profile$X2), sep = "_"), 
                               Mutations_Cosmic = apply(CCLE_Mut_profile[,-c(1:2)], 1, function(x) paste(x[!is.na(x)], collapse = ", ")))

colnames(CCLE_enzymes) <- str_remove_all(colnames(CCLE_enzymes),"_TenPx..")
CCLE_enzymes_pivoted <- pivot_longer(CCLE_enzymes, -c(1:6), names_to = "Cell_line",
                                     values_to = "Abundance")
CCLE_enzymes_pivoted$Cell_type <- str_match(CCLE_enzymes_pivoted$Cell_line, "_([:graph:]*)")[,2]
#CCLE_enzymes_pivoted <- CCLE_enzymes_pivoted[,-c(4,8)]%>%na.omit()

save(nodes, #done
     edges_enzymes_tissues, #done
     folate_cell_line_data, #done
     mod_tissue_annotation,#done
     Metabolites_mapped,#done
     Enzymes_Metabolites, #done
     p, #done
     deep_map_pivoted, #done
     CCLE_enzymes_pivoted, #done
     CCLE_Mut_profile, #done
     file = here::here("Datasets","Processed","Input_Data.RData"))

