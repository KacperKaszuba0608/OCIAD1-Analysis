# INPUT: 6 files used for plotting
#   - 2 CSV files from the output of the fold change calculation files
#   - 4 RDS files from the output of the GO enrichment file

library(tidyverse)
library("org.Hs.eg.db")
library(biomaRt)
library(clusterProfiler)
library(enrichplot)
library(patchwork)
library(ggrepel)

final_table_MITOS <- read_csv('./data/cleaned/OCIAD1_proteomics_mitos_process.csv', show_col_types = FALSE)
final_table_TOTALS <- read_csv('./data/cleaned/OCIAD1_proteomics_totals_process.csv', show_col_types = FALSE)

# Cut offs for labels
label.FC.cutoff <- 2.5
label.p.cutoff <- 0.05

# downloading all annotation
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
uniprots <- Rkeys(org.Hs.egUNIPROT)
ENSEMBL_ids <- AnnotationDbi::select(org.Hs.eg.db, uniprots, "ENSEMBL", "UNIPROT")
Entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, uniprots, "ENTREZID", "UNIPROT")

#########
# https://www.bioconductor.org/help/course-materials/2014/useR2014/Integration.html
# get chr info
#listFilters(mart)  
myFilter <- c("chromosome_name")
#listFilterOptions(mart, "chromosome_name")
myValues <- listFilterOptions(mart, "chromosome_name")[c(1:22,503:505)]
myAttributes <- c("ensembl_gene_id", "external_gene_name", listAttributes(mart)[9:13,1])
## assemble and query the mart
res <- getBM(attributes =  myAttributes, filters =  myFilter,
             values =  myValues, 
             mart = mart)
##########

my.data_M <- final_table_MITOS %>% 
    mutate("UNIPROT" = gsub(";.*", "", `Protein IDs`)) %>%
    mutate("UNIPROT" = gsub("-.*", "", UNIPROT))
my.data_M <- my.data_M %>% left_join(ENSEMBL_ids, by= "UNIPROT")
# 21 still un-annotated
my.data_M %>% filter(is.na(ENSEMBL)) %>% distinct(`Protein IDs`) %>% dim()

my.data_M <- my.data_M %>% left_join(res, by=c("ENSEMBL" = "ensembl_gene_id"))
my.data_M <- my.data_M %>% distinct(`Gene names`, .keep_all = TRUE)

#### volcano plot MITOS ####
my.data_temp <- my.data_M

my.data_temp$description <- case_when(
    grepl("Fatty acid oxidation", my.data_M$MitoCarta3.0_MitoPathways) ~ 2,
    !is.na(my.data_M$MitoCarta3.0_MitoPathways) ~ 1,
    TRUE ~ 0
)

my.data_temp$description <- factor(my.data_temp$description)

my.data_temp$label <- case_when(
    abs(my.data_temp$FC_MITOS_OCIAD1) > 2 & my.data_temp$sig_MITOS_OCIAD1 == TRUE & my.data_temp$description == 1 ~ my.data_temp$`Gene names`,
    my.data_temp$sig_MITOS_OCIAD1 == TRUE & my.data_temp$description == 2 ~ my.data_temp$`Gene names`,
    abs(my.data_temp$FC_MITOS_OCIAD1) > 4 & my.data_temp$sig_MITOS_OCIAD1 == TRUE & my.data_temp$description == 0 ~ my.data_temp$`Gene names`,
    my.data_temp$`Gene names` == 'OCIAD1' ~ my.data_temp$`Gene names`,
    TRUE ~ ''
)

my.data_temp$color <- case_when(
    grepl("Fatty acid oxidation", my.data_M$MitoCarta3.0_MitoPathways) ~ 'blue',
    !is.na(my.data_M$MitoCarta3.0_MitoPathways) ~ 'darkorange',
    TRUE ~ 'darkgrey'
)

v1 <- ggplot(my.data_temp, aes(x = FC_MITOS_OCIAD1, y = -log10(p_OCIAD1_MITOS), shape = sig_MITOS_OCIAD1, color = description)) +
    geom_hline(yintercept = -log10(0.05), color = 'lightgrey', alpha = 0.6, linetype = 2) +
    geom_vline(xintercept = 1, color = 'lightgrey', alpha = 0.6, linetype = 2) +
    geom_vline(xintercept = -1, color = 'lightgrey', alpha = 0.6, linetype = 2) +
    geom_point(data = filter(my.data_temp, description == 0), stroke = 1.2, alpha = 0.4) +
    geom_point(data = filter(my.data_temp, description == 1), stroke = 1.2, alpha = 0.55) +
    geom_point(data = filter(my.data_temp, description == 2), stroke = 1.2, alpha = 0.7) +
    xlab("Fold Change KO/WT") +
    ylab("-log10(p-value)") +
    theme_bw() +
    theme(legend.position = c(0.01, 0.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=9), legend.background = element_blank(), legend.key = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    guides(shape=FALSE, alpha = FALSE) +
    scale_shape_manual(values = c(1, 16)) +
    xlim(-max(abs(my.data_temp$FC_MITOS_OCIAD1)), max(abs(my.data_temp$FC_MITOS_OCIAD1))) +
    scale_color_manual(labels = c('Non-mitochondrial', 'Other mitochondrial', 'Fatty acid oxidation'), values = c('darkgrey','darkorange', 'blue')) +
    geom_text_repel(aes(label = label), max.overlaps = Inf, verbose = TRUE, min.segment.length = 0, color = my.data_temp$color, box.padding = 1)

v1

#### volcano plot TOTALS ####
my.data_T <- final_table_TOTALS %>% 
    mutate("UNIPROT" = gsub(";.*", "", `Protein IDs`)) %>%
    mutate("UNIPROT" = gsub("-.*", "", UNIPROT))
my.data_T <- my.data_T %>% left_join(ENSEMBL_ids, by= "UNIPROT")
# 21 still un-annotated
my.data_T %>% filter(is.na(ENSEMBL)) %>% distinct(`Protein IDs`) %>% dim()

my.data_T <- my.data_T %>% left_join(res, by=c("ENSEMBL" = "ensembl_gene_id"))
my.data_T <- my.data_T %>% distinct(`Gene names`, .keep_all = TRUE)

my.data_temp <- my.data_T

my.data_temp$description <- case_when(
    grepl("Fatty acid oxidation", my.data_T$MitoCarta3.0_MitoPathways) ~ 2,
    !is.na(my.data_T$MitoCarta3.0_MitoPathways) ~ 1,
    TRUE ~ 0
)

my.data_temp$description <- factor(my.data_temp$description)

my.data_temp$label <- case_when(
    abs(my.data_temp$FC_TOTALS_OCIAD1) > 2.5 & my.data_temp$sig_TOTALS_OCIAD1 == TRUE & my.data_temp$description == 1 ~ my.data_temp$`Gene names`,
    my.data_temp$sig_TOTALS_OCIAD1 == TRUE & my.data_temp$description == 2 ~ my.data_temp$`Gene names`,
    my.data_temp$`Gene names` == 'OCIAD1' ~ my.data_temp$`Gene names`,
    my.data_temp$`Gene names` == 'TIMM17A' ~ my.data_temp$`Gene names`,
    my.data_temp$`Gene names` == 'TIMM17B' ~ my.data_temp$`Gene names`,
    my.data_temp$`Gene names` == 'TIMM23' ~ my.data_temp$`Gene names`,
    TRUE ~ ''
)

my.data_temp$color <- case_when(
    grepl("Fatty acid oxidation", my.data_T$MitoCarta3.0_MitoPathways) ~ 'blue',
    !is.na(my.data_T$MitoCarta3.0_MitoPathways) ~ '#76EE35',
    TRUE ~ 'darkgrey'
)

v1 <- ggplot(my.data_temp, aes(x = FC_TOTALS_OCIAD1, y = -log10(p_OCIAD1_TOTALS), shape = sig_TOTALS_OCIAD1, color = description)) +
    geom_hline(yintercept = -log10(0.05), color = 'lightgrey', alpha = 0.6, linetype = 2) +
    geom_vline(xintercept = 1, color = 'lightgrey', alpha = 0.6, linetype = 2) +
    geom_vline(xintercept = -1, color = 'lightgrey', alpha = 0.6, linetype = 2) +
    geom_point(data = filter(my.data_temp, description == 0), stroke = 1.2, alpha = 0.4) +
    geom_point(data = filter(my.data_temp, description == 1), stroke = 1.2, alpha = 0.55) +
    geom_point(data = filter(my.data_temp, description == 2), stroke = 1.2, alpha = 0.7) +
    xlab("Fold Change KO/WT") +
    ylab("-log10(p-value)") +
    theme_bw() +
    theme(legend.position = c(0.01, 0.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=9), legend.background = element_blank(), legend.key = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    guides(shape=FALSE, alpha = FALSE) +
    scale_shape_manual(values = c(1, 16)) +
    xlim(-max(abs(my.data_temp$FC_TOTALS_OCIAD1)), max(abs(my.data_temp$FC_TOTALS_OCIAD1))) +
    scale_color_manual(labels = c('Non-mitochondrial', 'Other mitochondrial', 'Fatty acid oxidation', 'Protein import'), values = c('darkgrey', '#76EE35', 'blue', 'darkorange')) +
    geom_text_repel(aes(label = label), max.overlaps = Inf, verbose = TRUE, min.segment.length = 0, color = my.data_temp$color, box.padding = 1)

v1

#### GO TERM ENRICHMENT ####

my.data_M <- final_table_MITOS %>% 
    mutate("UNIPROT" = gsub(";.*", "", `Protein IDs`)) %>%
    mutate("UNIPROT" = gsub("-.*", "", UNIPROT))
my.data_M <- my.data_M %>% left_join(Entrez_ids, by= "UNIPROT")

my.data_T <- final_table_TOTALS %>% 
    mutate("UNIPROT" = gsub(";.*", "", `Protein IDs`)) %>%
    mutate("UNIPROT" = gsub("-.*", "", UNIPROT))
my.data_T <- my.data_T %>% left_join(Entrez_ids, by= "UNIPROT")

# GO plot heatmap up MITOS
go_enrichx2.up <- readRDS("./data/GO_enrichment/OCIAD1_proteomics_mitos_GOterms_up.rds")

geneList <- setNames(my.data_M$FC_MITOS_OCIAD1, my.data_M$`Gene names`)
color_scale <- scale_fill_gradient2(low = "red", mid = "grey", high = "green", midpoint = 0, limits = c(-5, 5))

p_go_heat_up <- heatplot(go_enrichx2.up, foldChange = geneList, showCategory=5) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_gradient(low="grey",high="green") +
    scale_y_discrete(position = 'right', labels = function(x) str_wrap(x, width = 25)) +
    coord_flip()

p_go_heat_up

# GO plot heatmap down MITOS
go_enrichx2.down <- readRDS("./data/GO_enrichment/OCIAD1_proteomics_mitos_GOterms_down.rds")

geneList <- setNames(my.data_M$FC_MITOS_OCIAD1, my.data_M$`Gene names`)
color_scale <- scale_fill_gradient2(low = "red", mid = "grey", high = "green", midpoint = 0, limits = c(-5, 5))

p_go_heat_down <- heatplot(go_enrichx2.down, foldChange = geneList, showCategory=5) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_gradient(low="red",high="grey") +
    scale_y_discrete(position = 'right', labels = function(x) str_wrap(x, width = 25)) +
    coord_flip()

p_go_heat_down

# GO plot heatmap up TOTALS
go_enrichx2.down <- readRDS('./data/GO_enrichment/OCIAD1_proteomics_totals_GOterms_up.rds')

geneList <- setNames(my.data_T$FC_TOTALS_OCIAD1, my.data_T$`Gene names`)
color_scale <- scale_fill_gradient2(low = "red", mid = "grey", high = "green", midpoint = 0, limits = c(-5, 5))

p_go_heat_up <- heatplot(go_enrichx2.up, foldChange = geneList, showCategory=5) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_gradient(low="grey",high="green") +
    scale_y_discrete(position = 'right', labels = function(x) str_wrap(x, width = 25)) +
    coord_flip()

p_go_heat_up

# GO plot heatmap down TOTALS
go_enrichx2.down <- readRDS('./data/GO_enrichment/OCIAD1_proteomics_totals_GOterms_down.rds')

geneList <- setNames(my.data_T$FC_TOTALS_OCIAD1, my.data_T$`Gene names`)
color_scale <- scale_fill_gradient2(low = "red", mid = "grey", high = "green", midpoint = 0, limits = c(-5, 5))

p_go_heat_down <- heatplot(go_enrichx2.down, foldChange = geneList, showCategory=5) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_gradient(low="red",high="grey") +
    scale_y_discrete(position = 'right', labels = function(x) str_wrap(x, width = 25)) +
    coord_flip()

p_go_heat_down

#### one_GO_term ####

my.data <- merge(final_table_MITOS, final_table_TOTALS, by = c('Protein IDs', 'Gene names'), all.y = TRUE)

this.ont = "CC"

this.GO.perox = "GO:0005777" #peroxisome (cellular component)
this.GO.er = "GO:0005783" #ER
this.GO.mito = "GO:0005739" #mitochondrium

GO_terms <- c(this.GO.perox, this.GO.er, this.GO.mito)
organelles <- c('Peroxisome', 'Endoplasmatic\nReticulum', 'Mitochondrium')

FC_organ <- data.frame(`Gene names`=my.data$`Gene names`, `FC_MITOS_OCIAD1`=my.data$`FC_MITOS_OCIAD1`, 
                       `FC_TOTALS_OCIAD1`=my.data$`FC_TOTALS_OCIAD1`, Organelle=NA)

for (n in 1:length(GO_terms)) {
    
    this.GO <- GO_terms[n]
    organ <- organelles[n]
    
    retrieved <- AnnotationDbi::select(org.Hs.eg.db, keytype = "GOALL", keys = this.GO, columns = c("ENSEMBL", "UNIPROT"))
    
    my.data_organ <- my.data |> 
        mutate("UNIPROT" = gsub(";.*", "", `Protein IDs`)) |>
        mutate("UNIPROT" = gsub("-.*", "", UNIPROT))
    my.data_organ <- my.data_organ |> left_join(retrieved, by= "UNIPROT")
    
    my.list <- my.data_organ |> filter(!is.na(GOALL)) |> select(`Protein IDs`) |> unique()
    
    FC_organ$Organelle <- ifelse(my.data$`Protein IDs` %in% my.list$`Protein IDs` & !is.na(FC_organ$Organelle),
                                 paste0(organ, ";", FC_organ$Organelle), 
                                 ifelse(my.data$`Protein IDs` %in% my.list$`Protein IDs`, organ, FC_organ$Organelle))
}

FC_organ$Organelle[is.na(FC_organ$Organelle)] <- 'Unspecified'

#### PLOTS ####

FC_organ_copy <- data.frame(FC_organ)

FC_organ_copy <- filter(FC_organ_copy, 
                        Organelle == 'Peroxisome' | 
                            Organelle == 'Mitochondrium' | 
                            Organelle == 'Mitochondrium;Peroxisome' |
                            Organelle == 'Endoplasmatic\nReticulum'
)

# organelle volcano plot MITOS
FC_organ$pval <- my.data$p_OCIAD1_MITOS
FC_organ$sig <- my.data$sig_MITOS_OCIAD1
FC_organ$mito <- my.data$MitoCarta3.0_SubMitoLocalization.x

FC_organ$case <- case_when(
    grepl('Peroxisome', FC_organ$Organelle) & !grepl('Endoplasmatic\nReticulum', FC_organ$Organelle) & !grepl('Mitochondrium', FC_organ$Organelle) ~ 1,
    !grepl('Peroxisome', FC_organ$Organelle) & !grepl('Endoplasmatic\nReticulum', FC_organ$Organelle) & grepl('Mitochondrium', FC_organ$Organelle) ~ 4,
    grepl('Peroxisome', FC_organ$Organelle) & !grepl('Endoplasmatic\nReticulum', FC_organ$Organelle) & grepl('Mitochondrium', FC_organ$Organelle) ~ 5,
    TRUE ~ 0
)

FC_organ$case <- factor(FC_organ$case)

FC_organ$label <- case_when(
    FC_organ$sig == TRUE & grepl('Peroxisome', FC_organ$Organelle) ~ FC_organ$`Gene.names`,
    FC_organ$sig == TRUE & abs(FC_organ$FC_MITOS_OCIAD1) >= 2 & FC_organ$case == 4 ~ FC_organ$`Gene.names`,
    FC_organ$sig == TRUE & abs(FC_organ$FC_MITOS_OCIAD1) >= 4 & FC_organ$case == 0 ~ FC_organ$`Gene.names`,
    FC_organ$`Gene.names` == 'ABCD3' ~ FC_organ$`Gene.names`,
    TRUE ~ ''
)

FC_organ$color <- case_when(
    grepl('Peroxisome', FC_organ$Organelle) & !grepl('Endoplasmatic\nReticulum', FC_organ$Organelle) & !grepl('Mitochondrium', FC_organ$Organelle) ~ 'darkgreen',
    !grepl('Peroxisome', FC_organ$Organelle) & !grepl('Endoplasmatic\nReticulum', FC_organ$Organelle) & grepl('Mitochondrium', FC_organ$Organelle) ~ 'darkorange',
    grepl('Peroxisome', FC_organ$Organelle) & !grepl('Endoplasmatic\nReticulum', FC_organ$Organelle) & grepl('Mitochondrium', FC_organ$Organelle) ~ 'purple',
    TRUE ~ 'darkgrey'
)

v2 <- ggplot(FC_organ, aes(x = FC_MITOS_OCIAD1, y = -log10(pval), shape = sig, color = case)) +
    geom_hline(yintercept = -log10(0.05), color = 'lightgrey', alpha = 0.6, linetype = 2) +
    geom_vline(xintercept = 1, color = 'lightgrey', alpha = 0.6, linetype = 2) +
    geom_vline(xintercept = -1, color = 'lightgrey', alpha = 0.6, linetype = 2) +
    geom_point(data = filter(FC_organ, case == 0), stroke = 1.2, alpha = 0.4) +
    geom_point(data = filter(FC_organ, case == 4), stroke = 1.2, alpha = 0.7) +
    geom_point(data = filter(FC_organ, case == 1), stroke = 1.2, alpha = 0.7) +
    geom_point(data = filter(FC_organ, case == 5), stroke = 1.2, alpha = 0.7) +
    xlab("Fold Change KO/WT") +
    ylab("-log10(p-value)") +
    theme_bw() +
    theme(legend.position = c(0.01, 1), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=9), legend.background = element_blank(), legend.key = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    guides(shape=FALSE, alpha = FALSE) +
    scale_shape_manual(values = c(1, 16)) +
    xlim(-max(abs(FC_organ$FC_MITOS_OCIAD1)), max(abs(FC_organ$FC_MITOS_OCIAD1))) +
    scale_color_manual(labels = c('Unspecified', 'Mitochondrium', 'Peroxisome', 'Mitochondrium and Peroxisome'), values = c('darkgrey', 'darkorange', 'darkgreen', 'purple')) +
    geom_text_repel(aes(label = label), max.overlaps = Inf, verbose = TRUE, min.segment.length = 0, color = FC_organ$color, box.padding = 1)

v2

# organelle volcano plot TOTALS
FC_organ$pval <- my.data$p_OCIAD1_TOTALS
FC_organ$sig <- my.data$sig_TOTALS_OCIAD1
FC_organ$mito <- my.data$MitoCarta3.0_SubMitoLocalization.y

FC_organ$case <- case_when(
    grepl('Peroxisome', FC_organ$Organelle) & !grepl('Endoplasmatic\nReticulum', FC_organ$Organelle) & !grepl('Mitochondrium', FC_organ$Organelle) ~ 1,
    !grepl('Peroxisome', FC_organ$Organelle) & !grepl('Endoplasmatic\nReticulum', FC_organ$Organelle) & grepl('Mitochondrium', FC_organ$Organelle) ~ 4,
    grepl('Peroxisome', FC_organ$Organelle) & !grepl('Endoplasmatic\nReticulum', FC_organ$Organelle) & grepl('Mitochondrium', FC_organ$Organelle) ~ 5,
    TRUE ~ 0
)

FC_organ$case <- factor(FC_organ$case)

FC_organ$label <- case_when(
    FC_organ$sig == TRUE & grepl('Peroxisome', FC_organ$Organelle) ~ FC_organ$`Gene.names`,
    FC_organ$sig == TRUE & abs(FC_organ$FC_TOTALS_OCIAD1) >= 2 & FC_organ$case == 4 ~ FC_organ$`Gene.names`,
    FC_organ$sig == TRUE & abs(FC_organ$FC_TOTALS_OCIAD1) >= 4 & FC_organ$case == 0 ~ FC_organ$`Gene.names`,
    FC_organ$`Gene.names` == 'ABCD3' ~ FC_organ$`Gene.names`,
    TRUE ~ ''
)

FC_organ$color <- case_when(
    grepl('Peroxisome', FC_organ$Organelle) & !grepl('Endoplasmatic\nReticulum', FC_organ$Organelle) & !grepl('Mitochondrium', FC_organ$Organelle) ~ 'darkgreen',
    !grepl('Peroxisome', FC_organ$Organelle) & !grepl('Endoplasmatic\nReticulum', FC_organ$Organelle) & grepl('Mitochondrium', FC_organ$Organelle) ~ 'darkorange',
    grepl('Peroxisome', FC_organ$Organelle) & !grepl('Endoplasmatic\nReticulum', FC_organ$Organelle) & grepl('Mitochondrium', FC_organ$Organelle) ~ 'purple',
    TRUE ~ 'darkgrey'
)

v2 <- ggplot(FC_organ, aes(x = FC_TOTALS_OCIAD1, y = -log10(pval), shape = sig, color = case)) +
    geom_hline(yintercept = -log10(0.05), color = 'lightgrey', alpha = 0.6, linetype = 2) +
    geom_vline(xintercept = 1, color = 'lightgrey', alpha = 0.6, linetype = 2) +
    geom_vline(xintercept = -1, color = 'lightgrey', alpha = 0.6, linetype = 2) +
    geom_point(data = filter(FC_organ, case == 0), stroke = 1.2, alpha = 0.4) +
    geom_point(data = filter(FC_organ, case == 4), stroke = 1.2, alpha = 0.7) +
    geom_point(data = filter(FC_organ, case == 1), stroke = 1.2, alpha = 0.7) +
    geom_point(data = filter(FC_organ, case == 5), stroke = 1.2, alpha = 0.7) +
    xlab("Fold Change KO/WT") +
    ylab("-log10(p-value)") +
    theme_bw() +
    theme(legend.position = c(0.01, 1), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=9), legend.background = element_blank(), legend.key = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    guides(shape=FALSE, alpha = FALSE) +
    scale_shape_manual(values = c(1, 16)) +
    xlim(-max(abs(FC_organ$FC_TOTALS_OCIAD1)), max(abs(FC_organ$FC_TOTALS_OCIAD1))) +
    scale_color_manual(labels = c('Unspecified', 'Mitochondrium', 'Peroxisome', 'Mitochondrium and Peroxisome'), values = c('darkgrey', 'darkorange', 'darkgreen', 'purple')) +
    geom_text_repel(aes(label = label), max.overlaps = Inf, verbose = TRUE, min.segment.length = 0, color = FC_organ$color, box.padding = 1)

v2