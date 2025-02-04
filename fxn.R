show_cv_table <- function(df) {
    DT::datatable(df, rownames = FALSE, 
                  options=list(searching=FALSE, paging=FALSE, info=FALSE)) |>
        DT::formatStyle(columns = c("mean_KO", "mean_WT", "median_KO", "median_WT"), 
                        background = DT::styleInterval(c(15, 25), c("white","lightgreen", "white"))) |>
        DT::formatStyle(columns = "method", fontWeight = "bold")
}

assign_missing <- function(protein.ids, condition, lfq_intensity) {
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
        dplyr::summarise(no_NAs = sum(is.na(lfq)))
    
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
        dplyr::summarise(no_NAs = sum(is.na(lfq)))
    
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
    
    df_miss <- data.frame(prot.id.miss, missingness_per_cond, missingness_per_prot)
    colnames(df_miss) <- c("prot.IDs", "missingness_per_cond", "missingness_per_prot")
    
    missingness_per_cond <- lapply(unique(protein.ids), function(id) df_miss[which(df_miss$prot.IDs == id),"missingness_per_cond"])
    df$missingness_per_cond <- do.call(c,missingness_per_cond)
    
    missingness_per_prot <- lapply(unique(protein.ids), function(id) df_miss[which(df_miss$prot.IDs == id),"missingness_per_prot"])
    df$missingness_per_prot <- do.call(c,missingness_per_prot)
    
    # Testing new missing assigning
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

plothist <- function(df, title="", plot.title.and.legend=TRUE) {
    if (!is.data.frame(df)) { df <- as.data.frame(df)}
    
    # Assuming LFQ_KO is our data frame with at least 3 columns
    # Reshape the data to a long format
    df <- df |>
        pivot_longer(cols = 1:ncol(df), names_to = "Rep", values_to = "LFQValue")
    
    # Plot all histograms on the same plot using ggplot
    if (plot.title.and.legend) {
        ret_plot <- ggplot(df, aes(x = LFQValue, fill = Rep)) +
            geom_histogram(alpha = 0.4, position = "identity", bins = 30) +
            labs(title = paste(title), x = "Values", y = "Frequency") +
            theme(legend.title = element_blank())
    } else {
        ret_plot <- ggplot(df, aes(x = LFQValue, fill = Rep)) +
            geom_histogram(alpha = 0.4, position = "identity", bins = 30) +
            labs(x = "Values", y = "Frequency") +
            theme(legend.title = element_blank(), legend.position = "none")
    }
    return(ret_plot)
}

plotoneviolin <- function(object, title="") {
    object <- as.data.frame(object) |>
        pivot_longer(everything(), names_to = "Sample", values_to = "LFQ_CV")
    
    ggplot(data=object, aes(x=Sample, y=LFQ_CV, fill=Sample))+
        geom_boxplot(width=0.2, na.rm = TRUE)+
        geom_violin(alpha=0.4, na.rm = TRUE)+
        labs(title=title,x=NULL,y="LFQ CV[%]")+
        theme(legend.position = "none",
              panel.background = element_rect(fill="white", colour = "grey"),
              panel.grid = element_line(colour = "grey"))
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

imputation <- function(rf_model, df_to_model, df_ml) {
    
    mar_prots <- which(df_ml$missingness == "MAR")
    mnar_prots_long <- which(df_ml$missingness == "MNAR")
    
    mar_impute <- predict(rf_model, df_to_model[mar_prots,-ncol(df_to_model)])$predictions
    
    lfq_wide <- df_ml |>
        dplyr::select(celltype, rep, LFQvalue) |>
        dplyr::mutate(batch = paste(celltype, rep, sep=".")) |>
        dplyr::select(batch, LFQvalue) |>
        tidyr::pivot_wider(id_cols = dplyr::everything(), names_from = batch, values_from = LFQvalue, values_fn = list) |>
        tidyr::unnest(c("KO.22", "KO.23", "KO.24", "WT.22", "WT.23", "WT.24")) |>
        dplyr::mutate(missingness = df_ml$missingness[seq(1,nrow(df_ml), by=6)])
    
    mnar_prots <- which(lfq_wide$missingness == "MNAR")
    
    lfq.KO <- lfq_wide |> dplyr::select(dplyr::contains("KO"))
    lfq.KO <- log2(lfq.KO)
    lfq.WT <- lfq_wide |> dplyr::select(dplyr::contains("WT"))
    lfq.WT <- log2(lfq.WT)
    
    mnar_impute.KO <- apply(lfq.KO, 2, shift_dist_impute) # function wrote by Mateusz
    mnar_impute.WT <- apply(lfq.WT, 2, shift_dist_impute) # function wrote by Mateusz
    
    lfq.KO[mnar_prots,] <- mnar_impute.KO[mnar_prots,]
    lfq.WT[mnar_prots,] <- mnar_impute.WT[mnar_prots,]
    
    lfq_imp <- cbind.data.frame(lfq.KO, lfq.WT) |>
        tidyr::pivot_longer(cols = dplyr::everything(), names_to = "batch", values_to = "LFQ")
    
    df_final <- df_ml
    
    df_final$LFQvalue <- lfq_imp$LFQ
    
    df_final$LFQvalue[mar_prots] <- log2(mar_impute)
    
    df_final <- df_final |>
        dplyr::select(celltype, rep, missingness, LFQvalue) |>
        dplyr::mutate(batch = paste(celltype, rep, sep=".")) |>
        dplyr::select(batch, LFQvalue) |>
        tidyr::pivot_wider(id_cols = everything(), names_from = batch, values_from = LFQvalue, values_fn = list) |>
        tidyr::unnest(c("KO.22", "KO.23", "KO.24", "WT.22", "WT.23", "WT.24")) |>
        dplyr::mutate(missingness = df_ml$missingness[seq(1,nrow(df_ml), by=6)])
    
    return(df_final)
}