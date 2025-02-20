# INPUT: a raw data from proteinGroups.txt
# OUTPUT: a data frame with filtered and imputed data for TOTALS

# IMPUTATION OF THE MISSING DATA for TOTALS 
# At first we assign the missingness to decide which rows should be included into
# further analysis.

# import libraries and fxns
library(tidyverse)
library(protti)
source('.//proteomics/fxn.R')

# import data
protein.groups <- read_tsv('./data/proteinGroups.txt', show_col_types = F) # 5561 proteins

# filtering the data
protein.groups <- protein.groups |> filter(is.na(`Only identified by site`), # 5432 proteins left
                                           is.na(Reverse),
                                           is.na(`Potential contaminant`))
# extracting lfq values
lfq <- protein.groups |>
    select(starts_with("LFQ Intensity") 
           & (ends_with("22") | ends_with("23") | ends_with("24")) 
           & contains("TOTALS"))

colnames(lfq) <- gsub("_TOTALS_", ".", gsub("LFQ intensity ", "", colnames(lfq)))
lfq[lfq == 0] <- NA

# Calculate number of missing values in each condition
tempKO <- apply(lfq[grep("KO", colnames(lfq))], 1, function(row) sum(is.na(row)))
tempWT <- apply(lfq[grep("WT", colnames(lfq))], 1, function(row) sum(is.na(row)))
temp <- cbind.data.frame(KO=tempKO, WT=tempWT) # concatenated via columns

rows.to.keep <- which((temp[,1] == 3 & temp[,2] == 0) | (temp[,1] == 0 & temp[,2] == 3) |
                          (temp[,1] == 1 & temp[,2] == 0) | (temp[,1] == 0 & temp[,2] == 1) |
                          (temp[,1] == 0 & temp[,2] == 0) | (temp[,1] == 1 & temp[,2] == 1) |
                          (temp[,1] == 2 & temp[,2] == 0) | (temp[,1] == 0 & temp[,2] == 2))

protein.groups <- protein.groups[rows.to.keep, ] # 3909 proteins left

################################## IMPUTATION ##################################

# extracting columns to imputation
df_miss_TOTALS <- protein.groups |>
    select(starts_with("LFQ Intensity") 
           & (ends_with("22") | ends_with("23") | ends_with("24")) 
           & contains("TOTALS"))

# changing the colnames
colnames(df_miss_TOTALS) <- gsub("LFQ intensity ", "", colnames(df_miss_TOTALS))

# Preparing data to imputation #
# making longer data frame with one column containing all LFQ values
df_miss_TOTALS <- df_miss_TOTALS |>
    mutate(prot.id = paste("prot",1:nrow(df_miss_TOTALS),sep="_")) |>
    pivot_longer(1:6, names_to = "Sample", values_to = "LFQvalue") |>
    separate(col=Sample, into=c("celltype","sampletype","rep"), sep = "_", remove = FALSE) |>
    mutate(celltype = as.factor(celltype), sampletype = as.factor(sampletype),
           rep = as.factor(rep))

# changing values with 0 to NA
df_miss_TOTALS$LFQvalue[df_miss_TOTALS$LFQvalue==0] <- NA

missingness <- assign_missing(df_miss_TOTALS$prot.id, df_miss_TOTALS$celltype, df_miss_TOTALS$LFQvalue)
df_miss_TOTALS$missingness_per_cond <- missingness$missingness_per_cond
df_miss_TOTALS$missingness_per_prot <- missingness$missingness_per_prot
df_miss_TOTALS$missingness <- missingness$missingness

df_miss_TOTALS$missingness[is.na(df_miss_TOTALS$missingness)] <- "all_NA"

# Remove unnecessary columns and reordering columns
df_miss_TOTALS <- df_miss_TOTALS |>
    select(-prot.id, -Sample) |>
    mutate(missingness_per_prot = as.factor(missingness_per_prot),
           missingness_per_cond = as.factor(missingness_per_cond),
           missingness = as.factor(missingness))

# Preparing data for imputation
df_to_protti <- protein.groups |>
    select(`Protein IDs`, `Peptide sequences`, contains("LFQ intensity") & contains("TOTALS") 
           & (ends_with("22") | ends_with("23") | ends_with("24"))) |>
    mutate(`Protein IDs`= paste("prot_", 1:nrow(protein.groups), sep="")) |>
    pivot_longer(3:8, names_to = "Sample", values_to = "Intensity")|>
    mutate(Sample = gsub("LFQ intensity ", "", Sample)) |>
    mutate(Sample = gsub("_TOTALS_", "_", Sample)) |>
    separate(col =  Sample, into = c("celltype","rep"), sep = "_", remove = F) |>
    mutate(Condition = ifelse(celltype == "KO", "treated", "control"),
           Intensity = ifelse(Intensity == 0, NA, log2(Intensity))) |> 
    select(Sample, `Protein IDs`, `Peptide sequences`, Condition, Intensity)

# assign missingness with protti fxn
data_missing <- df_to_protti |>
    assign_missingness(sample=Sample,
                       condition = Condition,
                       grouping = `Protein IDs`,
                       intensity = Intensity,
                       ref_condition = "all")

temp <- df_miss_TOTALS$missingness
temp[which(temp == "all_NA")] <- NA
data_missing$missingness2 <- temp

# impute data with protti fxn using ludovic method
imputed_TOTALS <- impute(
    data_missing,
    sample = Sample,
    grouping = `Protein IDs`,
    intensity_log2 = Intensity,
    condition = Condition,
    comparison = comparison,
    missingness = missingness2,
    method = "ludovic",
    skip_log2_transform_error = TRUE,
    retain_columns = missingness
)

rm(data_missing, temp)

# extracting KO and WT
imputed_TOTALS <- imputed_TOTALS |>
    dplyr::select(Sample, imputed_intensity, `Protein IDs`, missingness2, imputed) |>
    tidyr::pivot_wider(id_cols = everything(), names_from = Sample, values_from = imputed_intensity)

# handling with duplicates because we want to keep information about imputation
duplicates <- imputed_TOTALS[duplicated(imputed_TOTALS$`Protein IDs`), 'Protein IDs']
# imputed_TOTALS[imputed_TOTALS$`Protein IDs` %in% duplicates$`Protein IDs`,] |> View()

correct_dup_rows <- as.data.frame(t(sapply(1:nrow(duplicates),function(i) {
    false_row = imputed_TOTALS[imputed_TOTALS$`Protein IDs` == duplicates$`Protein IDs`[i]
                              & imputed_TOTALS$imputed == FALSE,]
    true_row = imputed_TOTALS[imputed_TOTALS$`Protein IDs` == duplicates$`Protein IDs`[i]
                             & imputed_TOTALS$imputed == TRUE,]
    true_row[which(is.na(true_row))] = false_row[which(is.na(true_row))]
    true_row
})))

correct_dup_rows <- lapply(correct_dup_rows, unlist) |> as.data.frame()

rows_to_drop <- which(imputed_TOTALS$`Protein IDs` %in% duplicates$`Protein IDs` & imputed_TOTALS$imputed == FALSE)
imputed_TOTALS <- imputed_TOTALS[-rows_to_drop,]

imputed_TOTALS[imputed_TOTALS$`Protein IDs` %in% duplicates$`Protein IDs`, 4:9] <- correct_dup_rows[,4:9]

# Exporting results
protein.groups.export <- protein.groups

protein.groups.export$`LFQ intensity KO_TOTALS_22` <- imputed_TOTALS$KO_22
protein.groups.export$`LFQ intensity KO_TOTALS_23` <- imputed_TOTALS$KO_23
protein.groups.export$`LFQ intensity KO_TOTALS_24` <- imputed_TOTALS$KO_24
protein.groups.export$`LFQ intensity WT_TOTALS_22` <- imputed_TOTALS$WT_22
protein.groups.export$`LFQ intensity WT_TOTALS_23` <- imputed_TOTALS$WT_23
protein.groups.export$`LFQ intensity WT_TOTALS_24` <- imputed_TOTALS$WT_24
protein.groups.export$missingness <- imputed_TOTALS$missingness2
protein.groups.export$imputed <- imputed_TOTALS$imputed

write.table(protein.groups.export, file='./data/cleaned/proteinGroupsImputed_TOTALS.txt', 
            row.names = FALSE, sep = '\t')
