# OUTPUT: a CSV file with genes from the peroxisomeDB

# Create a CSV file from multi fasta file from http://216.92.14.62/Download/Homo_sapiens.fas
peroxisomeDB <- Biostrings::readDNAStringSet('http://216.92.14.62/Download/Homo_sapiens.fas')
fasta_headers <- peroxisomeDB@ranges@NAMES
temp <- sapply(1:length(fasta_headers), 
       function(i) {
           unlist(stringr::str_split(fasta_headers[i], ' '))[c(1,3)]
       }) |> t()

unique_genes <- unique(temp[,1])
multiple_IDs <- sapply(1:length(unique_genes),
                       function(i) {
                           stringr::str_flatten(temp[temp[,1] == unique_genes[i], 2], collapse = ";")
                       })

peroxisomeDB_my <- data.frame('Gene names' = unique_genes,
                               'Protein IDs' = multiple_IDs,
                               'Organelle' = rep('peroxisome', length(unique_genes)))

write.csv(peroxisomeDB_my, file = './data/peroxisomeDB.csv', row.names = FALSE)