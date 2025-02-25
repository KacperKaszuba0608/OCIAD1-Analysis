# INPUT: 4 files
#   1. a csv file with genes from the peroxisomeDB
#   2. a csv file with genes from the MitoCarta3.0 database
#   3. 2 files from the output of the fold change calculation scripts
# OUTPUT: a csv file with organelle assignment comparisons

# Prepare data frame with organelle assignment comparison
library("org.Hs.eg.db")
library(dplyr)

peroxisomeDB <- readr::read_csv('./data/peroxisomeDB.csv', show_col_types = F) |> # an output from peroxisomeDB.R
    select(Gene.names, Organelle)
colnames(peroxisomeDB) <- c('Gene.names', 'PeroxisomeDB')

mitocarta <- readr::read_csv('./data/Human.MitoCarta3.0.csv', show_col_types = F) |> # https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls
    select(Symbol, MitoCarta3.0_List)
colnames(mitocarta) <- c('Symbol', 'MitoCarta3.0')
    
fc_M <- readr::read_csv('./data/cleaned/OCIAD1_proteomics_mitos_process.csv', show_col_types = F)
fc_T <- readr::read_csv('./data/cleaned/OCIAD1_proteomics_totals_process.csv', show_col_types = F)

my.data <- merge(fc_M, fc_T, by = c('Protein IDs', 'Gene names'), all = TRUE)

#### ONE GO TERM SEARCH ####
this.GO.perox = "GO:0005777" #peroxisome (cellular component)
this.GO.er = "GO:0005783" #ER
this.GO.mito = "GO:0005739" #mitochondrium

GO_terms <- c(this.GO.perox, this.GO.er, this.GO.mito)
organelles <- c('Peroxisome', 'Endoplasmatic\nReticulum', 'Mitochondrion')

organ_df <- my.data |>
    select(`Protein IDs`, `Gene names`) |>
    mutate(Organelle = NA)

for (n in 1:length(GO_terms)) {
    
    this.GO <- GO_terms[n]
    organ <- organelles[n]
    
    retrieved <- AnnotationDbi::select(org.Hs.eg.db, keytype = "GOALL", keys = this.GO, columns = c("ENSEMBL", "UNIPROT"))
    
    my.data_organ <- my.data |> 
        mutate("UNIPROT" = gsub(";.*", "", `Protein IDs`)) |>
        mutate("UNIPROT" = gsub("-.*", "", UNIPROT))
    my.data_organ <- my.data_organ |> left_join(retrieved, by= "UNIPROT")
    
    my.list <- my.data_organ |> filter(!is.na(GOALL)) |> select(`Protein IDs`) |> unique()
    
    organ_df$Organelle <- ifelse(my.data$`Protein IDs` %in% my.list$`Protein IDs` & !is.na(organ_df$Organelle), paste0(organ, ";", organ_df$Organelle), ifelse(my.data$`Protein IDs` %in% my.list$`Protein IDs`, organ, organ_df$Organelle))
}

# Preparing the correct format of the organelle column
valid <- c('Peroxisome', 'Mitochondrion', 'Mitochondrion;Peroxisome', 'Endoplasmatic\nReticulum')

#Replacing all organelle combinations that we dont want to plot with Other
organ_df <- organ_df |>
    mutate(Organelle = case_when(
        Organelle %in% valid ~ Organelle, 
        TRUE ~ "Other"
    ))

organ_df <- merge(organ_df, mitocarta, by.x='Gene names', by.y='Symbol', all.x=TRUE)
organ_df <- merge(organ_df, peroxisomeDB, by.x='Gene names', by.y='Gene.names', all.x=TRUE) 
organ_df <- organ_df |>
    mutate(MitoCarta3.0 = ifelse(!is.na(MitoCarta3.0), TRUE, FALSE),
           PeroxisomeDB = ifelse(!is.na(PeroxisomeDB), TRUE, FALSE),
           "One GO Mitochondrion" = grepl('Mitochondrion', Organelle),
           "One GO Peroxisome" = grepl('Peroxisome', Organelle),
           "One GO ER" = grepl('Endoplasmatic\nReticulum', Organelle)) |>
    select(-Organelle)

write.csv(organ_df, file = './data/organelles.csv', row.names = FALSE)
