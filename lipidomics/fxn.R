get_sets <- function(df, ...){
  #Extract columns and their names identified by keywords
  #df = data frame containing the data to extract
  #... = keywords to identify the column names

  args_list <- list(...)
  
  #Identifying column names containing the keywords
  col_names <- lapply(args_list, function(x) grep(x, colnames(df), value = TRUE))
  
  #Extracting the columns based on the identified column names
  cols <- lapply(args_list, function(x) df[, grepl(x, colnames(df))]) 
  
  return(list(column_names = col_names, columns = cols))
}


normalize <- function(df, n){
  #Normalise the data by dividing each value by its column sum, then multiplying by the mean of all column sums, and lastly performing the log2 transformation
  #(the column sums are calculated using only identified lipids)
  #df = data frame with the data
  #n = number of rows with identified lipids
  
  #Calculating the mean of all column sums of n identified lipids
  col_sum_mean <- mean(colSums(head(df, n)))
  
  df <- df |> 
    #Scaling each value by its column sum of n identified lipids
    scale(center = FALSE, scale = colSums(head(., n))) |>
    as.data.frame(.) |>
    
    #Multiplying each value by the mean of column sums
    mutate(across(everything(), function(x) x * col_sum_mean)) |>
    
    #Applying the log2 transformation
    log2(.) |>
    as.data.frame(.)
  
  return(df)
}


fc <- function(df, n){
  #Calculate the fold change
  #df = data frame containing the data, where the first condition is represented by the first n columns and the second condition is represented by the latter n columns
  #n = number of replicates in each condition
  
  df_copy <- data.frame(df)
  
  #Calculating the difference of means of the log2 transformed intensities for each condition
  df_copy$FC <- rowMeans(df_copy[,1:n]) - rowMeans(df_copy[,(n+1):(2*n)])
  
  #Identifying if the intensities are changing up or down between the conditions
  df_copy$`FC up/down` <- ifelse(df_copy$FC > 0, 1, 0)
  
  return(df_copy)
}


ttest <- function(df, n){
  #Perform a t-test on the conditions
  #df = data frame containing the data, where the first condition is represented by the first n columns and the second condition is represented by the latter n columns
  #n = number of replicates in each condition
  
  df_copy <- data.frame(df)
  
  #Reversing the log2 transformation
  df_copy <- 2^df_copy
  
  #Performing the t-test and obtaining the pvalues
  test_result <- apply(df_copy, 1, function (x) t.test(x[1:n], x[(n+1):(2*n)], paired = T))
  pvalues <- sapply(test_result, function(x) x$p.value)
  
  #Reapplying the log2 transformation
  df_copy <- log2(df_copy)
  
  df_copy$pvalue <- pvalues
  
  #Identifying if the change between the conditions is statistically significant with the cutoff of 0.05
  df_copy$Significant <- ifelse(df_copy$pvalue < 0.05, 1, 0)
  
  return(df_copy)
}


categorize_lipids <- function(df) {
  #Categorise the lipids based on their class
  #df = data frame containing the data (it has to include the 'Lipid Class' column)
  
  #Defining the lipid categories (most common ones are included but the lists can be modified)
  fatty_acyls <- c("AC")
  glycerolipids <- c("Alkanyl-TG", "Alkenyl-DG", "Alkenyl-TG", "DG", "TG")
  phospholipids <- c("Alkanyl-TG", "CL", "GM3-NANA", "LysoPC", "LysoPE", "LysoPG", "LysoPI", "LysoPS", "PC", "PE", "PE-NMe", "PE-NMe2", "PG", "PI", "Plasmanyl-PC", "Plasmanyl-PE", "Plasmenyl-PC", "Plasmenyl-PE", "PS")
  sphingolipids <- c("Cer[AP]" , "Cer[AS]", "Cer[NP]", "Cer[NS]", "HexCer[NS]", "SHexCer", "SM", "SP")
  cholesteryl_esters <- c("CE")
  unid <- c("Unidentified")
  
  #Identifying the lipid category based on the lipid class
  df <- df |>
    mutate(`Lipid Category` = case_when(
      `Lipid Class` %in% fatty_acyls ~ "Fatty Acyl",
      `Lipid Class` %in% glycerolipids ~ "Glycerolipid",
      `Lipid Class` %in% phospholipids ~ "Phospholipid",
      `Lipid Class` %in% sphingolipids ~ "Sphingolipid",
      `Lipid Class` %in% cholesteryl_esters ~ "Cholesteryl Ester",
      `Lipid Class` %in% unid ~ "Unidentified",
      TRUE ~ "Other"
    )) |>
    select(`Lipid Category`, everything())
  
  return(df)
}


fatty_acids <- function(df){
  #Extract the information about the fatty acid chains from the lipids' identification
  #df = data frame containing the data (it has to include the 'Identification' column)
  
  df_copy <- data.frame(df)
  
  #Extracting the fatty acid chain information from the lipid identification
  chains <- regmatches(df_copy$Identification, gregexpr("\\d\\d:\\d?\\d", df_copy$Identification))
  
  #Extracting the ether information from the lipid identification
  ether_info <- regmatches(df_copy$Identification, gregexpr("P-|O-", df_copy$Identification))
  
  #Extracting the information about the bonds in the fatty acid chains
  chain_bonds <- lapply(chains, function(x) strsplit(x, ':'))
  
  #Extracting the fatty acid chain lengths
  chain_lengths <- lapply(chain_bonds, function(x) unlist(lapply(x, function (y) as.numeric(y[[1]][1]))))
  
  #Extracting the number of double bonds in the fatty acid chain
  double_bonds <- lapply(chain_bonds, function(x) unlist(lapply(x, function (y) as.numeric(y[[2]][1]))))
  
  is_even <- unlist(lapply(chain_lengths, function(x) all(x %% 2 == 0)))
  
  #Checking the saturation of the fatty acids based on the number of double bonds
  saturation <- lapply(double_bonds, function(x) if(length(x) > 1) {
    ifelse(all(x <= 1), ifelse(all(x == 0), 0, 1), 2)
  } else {
    ifelse(x > 2, 2, ifelse(x == 0, 0, 1))
  })
  
  ether <- unlist(lapply(ether_info, function(x) ifelse(length(x), TRUE, FALSE)))
  
  #Adding the information to the data frame
  df_copy$`Chain Lengths` <- chain_lengths
  df_copy$`Length Even` <- is_even
  df_copy$`Double Bonds` <- double_bonds
  df_copy$`Saturation` <- saturation
  df_copy$`Ether` <- ether
  
  return(df_copy)
}


filter_classes <- function(df, n){
  #Filter out the lipid classes with not enough data points
  #df = data frame containing the data (it has to include the 'Lipid Class' column)
  #n = threshold for the number of data points to filter the classes by
  
  #Calculating the number of data points for each lipid class
  counts <- as.data.frame(table(df$`Lipid Class`))
  
  #Filtering out the lipid classes that contain less than n data points
  counts_filtered <- filter(counts, Freq >= n)
  filtered_classes <- as.character(counts_filtered$Var1)
  
  #Filtering the data frame based on the filtered lipid classes
  filtered_df <- filter(df, `Lipid Class` %in% filtered_classes)
  
  return(list(filtered_classes = filtered_classes, filtered_df = filtered_df))
}