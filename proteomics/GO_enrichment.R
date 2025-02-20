# INPUT: 2 files with fold change calculation for MITOS and TOTALS
# OUTPUT: 8 files with GO terms analysis:
#   - 4 csv files with data frames used to create supplemental excel file
#   - 4 Rdata files with result of GO enrichment used to create heatmaps in `plotting.R` file

library(tidyverse)
library("org.Hs.eg.db")
library(biomaRt)
library(clusterProfiler)
library(enrichplot)

final_table_MITOS <- read_csv('./data/cleaned/OCIAD1_proteomics_mitos_process.csv', show_col_types = FALSE)
final_table_TOTALS <- read_csv('./data/cleaned/OCIAD1_proteomics_totals_process.csv', show_col_types = FALSE)

# downloading all annotation
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
uniprots <- Rkeys(org.Hs.egUNIPROT)
ENSEMBL_ids <- AnnotationDbi::select(org.Hs.eg.db, uniprots, "ENSEMBL", "UNIPROT")

#### GO TERM ENRICHMENT ####
Entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, uniprots, "ENTREZID", "UNIPROT")

my.data_M <- final_table_MITOS %>% 
    mutate("UNIPROT" = gsub(";.*", "", `Protein IDs`)) %>%
    mutate("UNIPROT" = gsub("-.*", "", UNIPROT))
my.data_M <- my.data_M %>% left_join(Entrez_ids, by= "UNIPROT")
# 21 still un-annotated
my.data_M %>% filter(is.na(ENTREZID)) %>% distinct(`Protein IDs`) %>% dim()

my.data_T <- final_table_TOTALS %>% 
    mutate("UNIPROT" = gsub(";.*", "", `Protein IDs`)) %>%
    mutate("UNIPROT" = gsub("-.*", "", UNIPROT))
my.data_T <- my.data_T %>% left_join(Entrez_ids, by= "UNIPROT")
# 21 still un-annotated
my.data_T %>% filter(is.na(ENTREZID)) %>% distinct(`Protein IDs`) %>% dim()

################################### GO MITOS ###################################
# GO enrichment up

this.group.up <- my.data_M %>% filter(sig_MITOS_OCIAD1 & FC_MITOS_OCIAD1 > (1)) %>% select(ENTREZID) %>% pull() #no_sig > 1

this.ont = "ALL"

go_enrich.up <- enrichGO(gene = this.group.up, OrgDb="org.Hs.eg.db", 
                         #universe = unique(sort(na.omit(my.data_M$ENTREZID))),
                         pvalueCutoff = 0.05, pAdjustMethod="fdr",
                         ont=this.ont)
go_enrichx.up <- setReadable(go_enrich.up, 'org.Hs.eg.db', 'ENTREZID')
go_enrichx2.up_M <- pairwise_termsim(go_enrichx.up)

# GO enrichment down
this.group.down <- my.data_M %>% filter(sig_MITOS_OCIAD1 & FC_MITOS_OCIAD1 < (-1)) %>% select(ENTREZID) %>% pull() #no_sig > 1

this.ont = "ALL"

go_enrich.down <- enrichGO(gene = this.group.down, OrgDb="org.Hs.eg.db", 
                           #universe = unique(sort(na.omit(my.data_M$ENTREZID))),
                           pvalueCutoff = 0.05, pAdjustMethod="fdr",
                           ont=this.ont)
go_enrichx.down <- setReadable(go_enrich.down, 'org.Hs.eg.db', 'ENTREZID')
go_enrichx2.down_M <- pairwise_termsim(go_enrichx.down)

# DATA TO EXPORT
df_go_mitos_up <- go_enrichx2.up_M@result
df_go_mitos_down <- go_enrichx2.down_M@result

################################## GO TOTALS ##################################
# GO enrichment up
this.group.up <- my.data_T %>% filter(sig_TOTALS_OCIAD1 & FC_TOTALS_OCIAD1 > (1)) %>% select(ENTREZID) %>% pull() #no_sig > 1

this.ont = "ALL"

go_enrich.up <- enrichGO(gene = this.group.up, OrgDb="org.Hs.eg.db", 
                         #universe = unique(sort(na.omit(my.data_T$ENTREZID))),
                         pvalueCutoff = 0.05, pAdjustMethod="fdr",
                         ont=this.ont)
go_enrichx.up <- setReadable(go_enrich.up, 'org.Hs.eg.db', 'ENTREZID')
go_enrichx2.up_T <- pairwise_termsim(go_enrichx.up)

# GO enrichment down
this.group.down <- my.data_T %>% filter(sig_TOTALS_OCIAD1 & FC_TOTALS_OCIAD1 < (-1)) %>% select(ENTREZID) %>% pull() #no_sig > 1

this.ont = "ALL"

go_enrich.down <- enrichGO(gene = this.group.down, OrgDb="org.Hs.eg.db", 
                           #universe = unique(sort(na.omit(my.data_T$ENTREZID))),
                           pvalueCutoff = 0.05, pAdjustMethod="fdr",
                           ont=this.ont)
go_enrichx.down <- setReadable(go_enrich.down, 'org.Hs.eg.db', 'ENTREZID')
go_enrichx2.down_T <- pairwise_termsim(go_enrichx.down)

# DATA TO EXPORT
df_go_totals_up <- go_enrichx2.up_T@result
df_go_totals_down <- go_enrichx2.down_T@result

#### EXOPORTING DATA ####
write.csv(df_go_mitos_up, "./data/GO_enrichment/OCIAD1_proteomics_mitos_GOterms_up.csv")
saveRDS(go_enrichx2.up_M, "./data/GO_enrichment/OCIAD1_proteomics_mitos_GOterms_up.rds")

write.csv(df_go_mitos_down, "./data/GO_enrichment/OCIAD1_proteomics_mitos_GOterms_down.csv")
saveRDS(go_enrichx2.down_M, "./data/GO_enrichment/OCIAD1_proteomics_mitos_GOterms_down.rds")

write.csv(df_go_totals_up, "./data/GO_enrichment/OCIAD1_proteomics_totals_GOterms_up.csv")
saveRDS(go_enrichx2.up_T, "./data/GO_enrichment/OCIAD1_proteomics_totals_GOterms_up.rds")

write.csv(df_go_totals_down, "./data/GO_enrichment/OCIAD1_proteomics_totals_GOterms_down.csv")
saveRDS(go_enrichx2.down_T, "./data/GO_enrichment/OCIAD1_proteomics_totals_GOterms_down.rds")