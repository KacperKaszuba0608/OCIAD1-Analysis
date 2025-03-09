assign_missing <- function(protein.ids, condition, lfq_intensity) {
    # NOTE: data frame must be in long format
    # protein.ids = column with protein IDs
    # condition = column conditions of samples (e.g. KO, WT)
    # lfq_intensity = column with LFQ intensity

    # verification of the fxn assumptions
    if (!is.factor(condition)) {rlang::abort("The condition are not factor!")}
    
    # occurences of protein IDs per condition
    occur <- dplyr::tibble(protein.ids) |>
        dplyr::group_by(protein.ids) |>
        dplyr::summarise(len = length(protein.ids)) |>
        dplyr::distinct(len)
    
    # occurences of condition
    occur_per_cond <- dplyr::tibble(protein.ids, condition) |>
        dplyr::group_by(protein.ids, condition) |>
        dplyr::summarise(len = length(protein.ids))
    occur_per_cond <- unique(occur_per_cond$len)
    
    if (length(protein.ids)/occur != length(unique(protein.ids))) {rlang::abort("The protein IDs are not unique!")}
    if (!is.numeric(lfq_intensity)) {rlang::abort("The lfq intensities are not numeric!")}
    
    df <- data.frame("prot.IDs"=protein.ids, "condition"=condition, "lfq"=lfq_intensity)
    
    # Calculating number of missing values for each protein for each condition
    number_missing1 <- df |>
        dplyr::group_by(prot.IDs, condition) |>
        dplyr::summarise(no_NAs = sum(is.na(lfq))) |>
        dplyr::ungroup()
    
    # Assigning the missing type per condition
    missingness_per_cond <- purrr::map(number_missing1$no_NAs, function(no) {
        if (no == occur_per_cond) {
            rep("all_NA", occur_per_cond)
        } else if (no == occur_per_cond-1) {
            rep("MNAR", occur_per_cond)
        } else if (no < occur_per_cond-1 & no != 0) {
            rep("MAR", occur_per_cond)
        } else if (no == 0){
            rep("complete", occur_per_cond)
        }
    })
    missingness_per_cond <- data.frame(do.call(c,missingness_per_cond))
    
    # Calculating number of missing values for each protein
    number_missing2 <- df |>
        dplyr::group_by(prot.IDs) |>
        dplyr::summarise(no_NAs = sum(is.na(lfq))) |>
        dplyr::ungroup()
    
    # Assigining the missing type per protein
    missingness_per_prot <- purrr::map(number_missing2$no_NAs, function(no) {
        if (no == occur) {
            rep("all_NA", occur)
        }else if (no == occur-1) {
            rep("MNAR", occur)
        } else if (no < occur-1 & no != 0) {
            rep("MAR", occur)
        } else if (no == 0){
            rep("complete", occur)
        }
    }) 
    missingness_per_prot <- data.frame(do.call(c, missingness_per_prot))
    
    prot.id.miss <- purrr::map(unique(number_missing2$prot.IDs), function(id) rep(id, occur))
    prot.id.miss <- data.frame(do.call(c,prot.id.miss))
    
    # A data frame with received missing types
    df_miss <- data.frame(prot.id.miss, missingness_per_cond, missingness_per_prot)
    colnames(df_miss) <- c("prot.IDs", "missingness_per_cond", "missingness_per_prot")
    
    missingness_per_cond <- lapply(unique(protein.ids), function(id) df_miss[which(df_miss$prot.IDs == id),"missingness_per_cond"])
    df$missingness_per_cond <- do.call(c,missingness_per_cond)
    
    missingness_per_prot <- lapply(unique(protein.ids), function(id) df_miss[which(df_miss$prot.IDs == id),"missingness_per_prot"])
    df$missingness_per_prot <- do.call(c,missingness_per_prot)
    
    # Final missing assigning
    missingness_final <- df |>
        dplyr::select(condition, missingness_per_cond) |>
        tidyr::pivot_wider(id_cols=everything(), names_from = "condition", 
                           values_from="missingness_per_cond", values_fn = list) |> 
        tidyr::unnest(cols = everything())
    
    cond_id <- seq(1,nrow(missingness_final), by=3)
    missingness_final <- missingness_final[cond_id,]
    
    missingness_final2 <- purrr::map(1:nrow(missingness_final), function(i) {
        if (all(missingness_final[i,] == "all_NA")) {
            rep(NA, occur)
        } else if ((missingness_final[i,1] == "complete" & missingness_final[i,2] == "all_NA") | (missingness_final[i,2] == "complete" & missingness_final[i,1] == "all_NA")) {
            rep("MNAR", occur)
        } else if ((missingness_final[i,1] == "MAR" & missingness_final[i,2] == "MNAR") | (missingness_final[i,2] == "MAR" & missingness_final[i,1] == "MNAR")) {
            rep(NA, occur)
        } else if (all(missingness_final[i,] == "MAR")) {
            rep("MAR", occur)
        } else if (all(missingness_final[i,] == "MNAR")) {
            rep(NA, occur)
        } else if (all(missingness_final[i,] == "complete")) {
            rep("complete", occur)
        } else if ((missingness_final[i,1] == "complete" & missingness_final[i,2] =="MNAR") | (missingness_final[i,2] == "complete" & missingness_final[i,1] =="MNAR")) {
            rep("MNAR", occur)
        } else if ((missingness_final[i,1] == "complete" & missingness_final[i,2] =="MAR") | (missingness_final[i,2] == "complete" & missingness_final[i,1] =="MAR")) {
            rep(NA, occur)  # keep NA and calculate p-value with less data
        } else {
            rep(NA, occur)
        }
    })
    
    df$missingness = do.call(c,missingness_final2)
    
    # return
    ret_list <- list(df = df, missingness_per_cond=df$missingness_per_cond,
                     missingness_per_prot = df$missingness_per_prot,
                     missingness = df$missingness)
}

ttest <- function(df, grp1, grp2){ 
    x = df[grp1]
    y = df[grp2]
    x = as.numeric((x))
    y = as.numeric((y))
    results = t.test(x,y, 
                     alternative = "two.sided", #one-sided: "greater" is x > y
                     paired = F,
                     na.action=na.omit)
    results$p.value
}

shift_dist_impute <- function(df, width = 0.5, downshift = 1, seed = 7694) {
    # df = data frame containing filtered data
    # Assumes missing data (in df) follows a narrowed and downshifted normal distribution
    
    #Seed for reproducibility (useful for later analysis)
    set.seed(seed)
    
    # Create new column indicating whether the values are imputed 
    # df$imputed = !is.finite(df$LFQvalue)
    
    # Imputation
    temp <- df
    temp[!is.finite(temp)] <- NA #make sure all non-finite values are really NA
    temp.sd <- width * sd(temp, na.rm = TRUE) #shrink sd width
    temp.mean <- mean(temp, na.rm = TRUE) - 
        downshift * sd(temp, na.rm = TRUE) #shift mean of imputed values
    n.missing <- sum(is.na(temp))
    #replace NA with values from the new distribution
    temp[is.na(temp)] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd) 
    
    df <- temp
    
    return(df)
}

is_imputed <- function(df1, df2) {
    if (!is.data.frame(df1)) {df1 <- as.data.frame(df1)}
    if (!is.data.frame(df2)) {df2 <- as.data.frame(df2)}
    
    res <- sapply(1:nrow(df1)[1], 
                  function(row) {
                      na_df1 = sum(is.na(df1[row,]))
                      na_df2 = sum(is.na(df2[row,]))
                      return(ifelse(na_df1 != na_df2, TRUE, FALSE))
                  })
    return(res)
}