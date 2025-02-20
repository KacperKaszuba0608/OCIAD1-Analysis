# INPUT: 2 files from the output of the fold change calculation files

library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(grid)

#### FXN ####
mean_lfq <- function(data, value='LFQ'){
    lfq_columns <- data |> select(contains(value))
    means <- apply(lfq_columns, 1, function(row) mean(row, na.rm=TRUE))
    return(means)
}

show_posthoc <- function(data, x, y, label1, label2, label3, label4, label5, sample=NULL) {
    organ_order <- data |>
        group_by(.data[[y]]) |>
        summarise(n = n()) |>
        arrange(n)
    
    ggplot(data=data, aes(x=.data[[x]], y=.data[[y]], fill=.data[[y]])) +
        geom_violin(alpha=0.7) +
        geom_text(mapping = aes(x=31, y=5.2, label=label1, size=6)) + # Other
        geom_text(mapping = aes(x=31, y=4.2, label=label2, size=6)) + # Mitochodnrion
        geom_text(mapping = aes(x=31, y=3.2, label=label3, size=6)) + # Endoplasmatic\nReticulum
        geom_text(mapping = aes(x=31, y=2.2, label=label4, size=6)) + # Peroxisome
        geom_text(mapping = aes(x=31, y=1.2, label=label5, size=6)) + # Mitochondrion;Peroxisome
        labs(title = paste0('Post-hoc resutls for ', 
                            ifelse(stringr::str_detect(x, 'WT'), 'WT ', 'KO '), 
                            ifelse(is.null(sample), '', paste0(sample, ' sample'))),
             subtitle = "Small letters indicate statistical differences between groups at the 0.05 significance level (Dunn's test).",
             x = 'Mean log2 LFQ Intensity', y='') +
        theme(legend.position = 'None') +
        scale_y_discrete(limits = organ_order[[y]])
}

calculate_ratios <- function(data, value, grouping='Organelle') {
    means <- data |>
        group_by(.data[[grouping]]) |>
        summarise(avg_WT = mean(.data[[value]]))
    
    org <- c('Peroxisome', 'Mitochondrion', 'Mitochondrion;Peroxisome', 'Endoplasmatic\nReticulum', 'Other')
    combn_org <- t(combn(org,2)) |> as.data.frame() |> filter(V2 == 'Other')
    combn_org$ratio <- apply(combn_org, 1, function(row) {
        means[means[,1]==row[1],2] - means[means[,1]==row[2],2]
    }) |> unlist()
    return(combn_org)
}

#### LOAD DATA ####
final_table_MITOS <- readr::read_csv('./data/cleaned/OCIAD1_proteomics_mitos_process.csv')
final_table_TOTALS <- readr::read_csv('./data/cleaned/OCIAD1_proteomics_totals_process.csv')

#### PREAPARING MAIN DATA SETS ####
my.data <- merge(final_table_MITOS, final_table_TOTALS, by = c('Protein IDs', 'Gene names'), all = TRUE)

lfq_mitos <- my.data |>
    select(`Protein IDs`, `Gene names`, contains('LFQ') & contains('MITOS'))

lfq_totals <- my.data |>
    select(`Protein IDs`, `Gene names`, contains('LFQ') & contains('TOTALS'))

organ_df <- my.data |>
    select(`Protein IDs`, `Gene names`) |>
    mutate(Organelle = NA)
    
#### ONE GO TERM SEARCH ####
this.GO.perox = "GO:0005777" #peroxisome (cellular component)
this.GO.er = "GO:0005783" #ER
this.GO.mito = "GO:0005739" #mitochondrium

GO_terms <- c(this.GO.perox, this.GO.er, this.GO.mito)
organelles <- c('Peroxisome', 'Endoplasmatic\nReticulum', 'Mitochondrion')

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

#### merging ####
lfq_mitos <- merge(lfq_mitos, organ_df)
lfq_mitos_WT <- lfq_mitos |> select(`Protein IDs`, `Gene names`, Organelle, contains('WT'))

lfq_totals <- merge(lfq_totals, organ_df)
lfq_totals_WT <- lfq_totals |> select(`Protein IDs`, `Gene names`, Organelle, contains('WT'))

# calculate lfq mean for each protein
lfq_mitos_WT$mean_lfq_WT_M <- mean_lfq(lfq_mitos_WT)
lfq_totals_WT$mean_lfq_WT_T <- mean_lfq(lfq_totals_WT)

lfq_mitos <- cbind.data.frame(lfq_mitos_WT) |>
    filter(!is.na(mean_lfq_WT_M)) |>
    mutate(Organelle = as.factor(Organelle))

lfq_totals <- lfq_totals_WT |>
    filter(!is.na(mean_lfq_WT_T)) |>
    mutate(Organelle = as.factor(Organelle))

#### STATISTICAL ANALYSIS ####

# Check the key criteria for performing the ANOVA test

# 1. Normality testing
# H0: The distribution of the data in the groups is derived from a Gaussian distribution.
# H1: The distribution of the data in the groups is derived from a non-Gaussian distribution.

lfq_mitos |>
    rstatix::group_by(Organelle) |>
    rstatix::shapiro_test(mean_lfq_WT_M) |>
    mutate(norm_dist = ifelse(p < 0.05, FALSE, TRUE))

lfq_totals |>
    rstatix::group_by(Organelle) |>
    rstatix::shapiro_test(mean_lfq_WT_T) |>
    mutate(norm_dist = ifelse(p < 0.05, FALSE, TRUE))

# 2. Homoscedasticity testing
bartlett.test(mean_lfq_WT_M ~ Organelle, data = lfq_mitos)
bartlett.test(mean_lfq_WT_T ~ Organelle, data = lfq_totals)

# COMMENT:
# Our data comes from a non-Gaussian distribution, so we couldn't use the ANOVA test. 
# Instead of ANOVA, we could use a non-parametric equivalent called the Kruskal-Wallis test.

#### Kruskal-Wallis Test for MITOS ####
# H0: There is no significant difference between organelles in a MITOS sample.
# H1: There is a significant difference between organelles in a MITOS sample.

kruskal.test(mean_lfq_WT_M ~ Organelle, data = lfq_mitos)

# COMMENT:
# For both conditions there is a significant difference (p < 0.05) between the organelles.
# Now we can use the post-hoc test (also non-parametric) to check between which.

# Non-Parametric Post-hoc test for MITOS

# WT MITOS
FSA::dunnTest(mean_lfq_WT_M ~ Organelle, data = lfq_mitos, method='bonferroni')$res |>
    mutate(sig_diff = ifelse(P.adj < 0.05, TRUE, FALSE)) |> DT::datatable()

# Visualization
show_posthoc(lfq_mitos, 'mean_lfq_WT_M', 'Organelle', 'a','b', 'b', 'b', 'b', sample = 'MITOS')

#### Kruskal-Wallis Test for TOTALS ####
# H0: There is no significant difference between organelles in a TOTALS sample.
# H1: There is a significant difference between organelles in a TOTALS sample.

kruskal.test(mean_lfq_WT_T ~ Organelle, data = lfq_totals) # WT TOTALS

# COMMENT:
# There is also a significant difference (p < 0.05) between the organelles in both conditions.
# We can also check between which groups using the non-parametric post-hoc test.

# Non-Parametric Post-hoc test for TOTALS

# WT TOTALS
FSA::dunnTest(mean_lfq_WT_T ~ Organelle, data = lfq_totals, method='bonferroni')$res |>
    mutate(sig_diff = ifelse(P.adj < 0.05, TRUE, FALSE)) |> DT::datatable()

# Visualization
show_posthoc(lfq_totals, 'mean_lfq_WT_T', 'Organelle', 'a', 'b', 'b', '', '', sample='TOTALS')

#### RATIOS ####

# ratio WT MITOS
ratio_M <- calculate_ratios(lfq_mitos, 'mean_lfq_WT_M')

# ratio WT TOTALS
ratio_T <- calculate_ratios(lfq_totals, 'mean_lfq_WT_T')

ratios <- merge(ratio_M, ratio_T, by = c('V1', 'V2'), suffixes=c('_M', '_T')) |>
    mutate(enrich = paste(round(ratio_M / ratio_T * 100 - 100,2), '%'))

#### SPLIT VIOLIN PLOT ####

#Split violin plot function
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                               # Original function by Jan Gleixner (@jan-glx)
                               # Adjustments by Wouter van der Bijl (@Axeman)
                               data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                               grp <- data[1, "group"]
                               newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                               newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                               newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                               if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                   stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                                   quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
                                   aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                   aesthetics$alpha <- rep(1, nrow(quantiles))
                                   both <- cbind(quantiles, aesthetics)
                                   quantile_grob <- GeomPath$draw_panel(both, ...)
                                   ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                               }
                               else {
                                   ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                               }
                           }
)

create_quantile_segment_frame <- function(data, draw_quantiles, split = FALSE, grp = NULL) {
    dens <- cumsum(data$density) / sum(data$density)
    ecdf <- stats::approxfun(dens, data$y)
    ys <- ecdf(draw_quantiles)
    violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
    violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
    violin.xs <- (stats::approxfun(data$y, data$x))(ys)
    if (grp %% 2 == 0) {
        data.frame(
            x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
            y = rep(ys, each = 2), group = rep(ys, each = 2)
        )
    } else {
        data.frame(
            x = ggplot2:::interleave(violin.xminvs, violin.xs),
            y = rep(ys, each = 2), group = rep(ys, each = 2)
        )
    }
}

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
    layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, 
          show.legend = show.legend, inherit.aes = inherit.aes, 
          params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
    
# Implementation of the split violin plot
FC_organ_copy <- merge(lfq_mitos, lfq_totals, by=c('Protein IDs', 'Gene names', 'Organelle'), all=TRUE) |>
    mutate(Organelle = gsub(';', '\n& ', Organelle))

organ_order <- FC_organ_copy |>
    group_by(Organelle) |>
    summarise(n = n()) |>
    arrange(n) |>
    mutate(order = 1:5)

FC_organ_copy <- FC_organ_copy |>
    mutate(color = case_when(
        Organelle == 'Mitochondrion' ~ 'orange',
        Organelle == 'Mitochondrion\n& Peroxisome' ~ '#998500',
        Organelle == 'Endoplasmatic\nReticulum' ~ 'purple',
        Organelle == 'Peroxisome' ~ 'darkgreen',
        TRUE ~ 'darkgrey'
    ))

#Pivoting the df so that the Batch and LFQs are their own columns
FC_organ_copy <- FC_organ_copy |> tidyr::pivot_longer(c(mean_lfq_WT_M, mean_lfq_WT_T), names_to = "Batch", values_to = "LFQ") |>
    mutate(transparency = ifelse(Batch == 'mean_lfq_WT_M', 0.9, 1))

# create custom legend
legend_grob <- grobTree(
    # Title
    textGrob("Intensity\ndistribution in", x = 0.1, y = 0.9, just = "left",
             gp = gpar(col = "black", fontsize = 9, fontface = "bold")),
    
    # Whole Cells color square
    rectGrob(x = 0.12, y = 0.77, width = 0.05, height = 0.1,
             gp = gpar(fill = "darkgrey", col = NA)),
    rectGrob(x = 0.17, y = 0.77, width = 0.05, height = 0.1, 
             gp = gpar(fill = "orange", col = NA)),
    rectGrob(x = 0.22, y = 0.77, width = 0.05, height = 0.1, 
             gp = gpar(fill = "purple", col = NA)),
    rectGrob(x = 0.27, y = 0.77, width = 0.05, height = 0.1, 
             gp = gpar(fill = "darkgreen", col = NA)),
    rectGrob(x = 0.32, y = 0.77, width = 0.05, height = 0.1, 
             gp = gpar(fill = "#998500", col = NA)),
    
    # Whole Cells text
    textGrob("Whole Cells", x = 0.37, y = 0.77, just = "left",
             gp = gpar(col = "black", fontsize = 9)),
    
    # Mitochondrial Fraction color square
    rectGrob(x = 0.12, y = 0.61, width = 0.05, height = 0.1,
             gp = gpar(fill = "darkgrey", col = NA, alpha=0.3)),
    rectGrob(x = 0.17, y = 0.61, width = 0.05, height = 0.1, 
             gp = gpar(fill = "orange", col = NA, alpha=0.3)),
    rectGrob(x = 0.22, y = 0.61, width = 0.05, height = 0.1, 
             gp = gpar(fill = "purple", col = NA, alpha=0.3)),
    rectGrob(x = 0.27, y = 0.61, width = 0.05, height = 0.1, 
             gp = gpar(fill = "darkgreen", col = NA, alpha=0.3)),
    rectGrob(x = 0.32, y = 0.61, width = 0.05, height = 0.1, 
             gp = gpar(fill = "#998500", col = NA, alpha=0.3)),
    textGrob("Mitochondrial\nFraction", x = 0.37, y = 0.61, just = "left",
             gp = gpar(col = "black", fontsize = 9))
)

stars_df <- data.frame(x = c(4.1, 3.8, 3.1, 2.8, 2.1, 1.8, 1.1, 0.8),
                       y = rep(22,8),
                       label = c('***', '***', '***', '***', '', '***', '', '***'))

#Implemented similarily as a normal violin plot, the split is determined by the fill parameter in aes
organelle_violin_split <- ggplot() +
    geom_split_violin(data=FC_organ_copy, 
                      mapping=aes(x = reorder(`Organelle`, LFQ, FUN = function(x) -median(x)), 
                                  y = LFQ, fill = interaction(Batch, color), alpha = transparency,
                                  color = ifelse(Batch == 'mean_lfq_WT_M', color, 'black')),
                      trim = FALSE, draw_quantiles = c(0.5), linewidth = 0.7, scale = 'width') +
    ylab('Log2 of the Mean LFQ Intensity') +
    theme_bw() +
    theme(axis.text = element_text(size = 12), axis.title.y = element_blank(), 
          legend.position = c(0.99, 1), legend.justification = c("right", "top"), 
          legend.title=element_blank(), legend.text = element_text(size=9), 
          legend.background = element_blank(), legend.key = element_blank()) +
    scale_y_continuous(position="right") +
    scale_fill_manual(values = setNames(FC_organ_copy$color, interaction(FC_organ_copy$Batch, FC_organ_copy$color))) +
    scale_color_identity() +
    coord_flip() +
    scale_x_discrete(limits = organ_order$Organelle) +
    geom_text(data=organ_order, mapping=aes(x=order+0.1, y=18, label=paste0('n = ',n))) +
    guides(fill = 'none', alpha='none') +
    annotation_custom(legend_grob, xmin = 1, xmax = 3.9, ymin = 30, ymax = 35) +
    geom_text(data=stars_df, mapping = aes(x=x, y=y, label=label), size=8)

organelle_violin_split

ggsave('./proteomics/plots/split_violin.png', plot = organelle_violin_split,
       width = 6, height = 8, units = "in")
ggsave('./proteomics/plots/split_violin.pdf', plot = organelle_violin_split,
       width = 6, height = 8, units = "in")
