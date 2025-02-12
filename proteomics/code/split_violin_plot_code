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


# Implementation of the split violin plot
FC_organ_copy <- data.frame(FC_organ)

valid <- c('Peroxisome', 'Mitochondrion', 'Mitochondrion;Peroxisome', 'Endoplasmatic\nReticulum')

#Replacing all organelle combinations that we dont want to plot with Other
FC_organ_copy <- FC_organ_copy %>%
  mutate(Organelle = case_when(
    Organelle %in% valid ~ Organelle, 
    TRUE ~ "Other"
  ))

#Pivoting the df so that the Batch and LFQs are their own columns
FC_organ_copy <- FC_organ_copy %>% pivot_longer(c(MITOS_WT, TOTALS_WT), names_to = "Batch", values_to = "LFQ")


#Implemented similarily as a normal violin plot, the split is determined by the fill parameter in aes
organelle_violin_split <- ggplot(FC_organ_copy, aes(x = reorder(`Organelle`, LFQ, FUN = function(x) -median(x)), y = LFQ, fill = `Batch`)) +
  geom_split_violin(trim = FALSE ,alpha = 0.6, col = '#4B4B4B', draw_quantiles = c(0.5), linewidth = 0.7, scale = 'width') +
  ylab('Log2 of the Mean LFQ Intensity') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.title.y = element_blank(), legend.position = c(0.99, 1), legend.justification = c("right", "top"), legend.title=element_blank(), legend.text = element_text(size=9), legend.background = element_blank(), legend.key = element_blank()) +
  scale_y_continuous(position="right") +
  scale_fill_manual(labels = c('MITOS', 'TOTALS'), values = c('lightblue', 'red')) +
  coord_flip()

organelle_violin_split
