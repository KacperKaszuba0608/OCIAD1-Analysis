library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggrepel)

ociad1_pd <- readr::read_csv('./OCIAD1-PD.csv') |>
    select('Protein ID', Gene, contains('FC')) |>
    mutate(Organelle = NA)

#### ONE GO TERM SEARCH ####
this.GO.perox = "GO:0005777" #peroxisome (cellular component)
this.GO.mito = "GO:0005739" #mitochondrium

GO_terms <- c(this.GO.perox, this.GO.mito)
organelles <- c('Peroxisome', 'Mitochondrion')

for (n in 1:length(GO_terms)) {
    
    this.GO <- GO_terms[n]
    organ <- organelles[n]
    
    retrieved <- AnnotationDbi::select(org.Hs.eg.db, keytype = "GOALL", keys = this.GO, columns = c("ENSEMBL", "UNIPROT"))
    
    my.data_organ <- ociad1_pd |> 
        mutate("UNIPROT" = gsub(";.*", "", `Protein ID`)) |>
        mutate("UNIPROT" = gsub("-.*", "", UNIPROT))
    my.data_organ <- my.data_organ |> left_join(retrieved, by= "UNIPROT")
    
    my.list <- my.data_organ |> filter(!is.na(GOALL)) |> select(`Protein ID`) |> unique()
    
    ociad1_pd$Organelle <- ifelse(ociad1_pd$`Protein ID` %in% my.list$`Protein ID` & !is.na(ociad1_pd$Organelle), paste0(organ, ";", ociad1_pd$Organelle), ifelse(ociad1_pd$`Protein ID` %in% my.list$`Protein ID`, organ, ociad1_pd$Organelle))
}

# Preparing the correct format of the organelle column
valid <- c('Peroxisome', 'Mitochondrion', 'Mitochondrion;Peroxisome')

#Replacing all organelle combinations that we dont want to plot with Other
to_plot <- ociad1_pd

to_plot <- to_plot |>
    mutate(Organelle = case_when(
        Organelle %in% valid ~ Organelle, 
        TRUE ~ "Other"
    ),
    # outliers = ifelse(FC_HAO1_CTRL %in% boxplot.stats(ociad1_pd$FC_HAO1_CTRL)$out & FC_HAO1_CTRL > 0, TRUE, FALSE),
    Organelle = gsub(';','\n& ', Organelle),
    color = case_when(
        Organelle == 'Mitochondrion' ~ 'orange',
        Organelle == 'Mitochondrion\n& Peroxisome' ~ '#998500',
        Organelle == 'Peroxisome' ~ 'darkgreen',
        TRUE ~ 'darkgrey'
    ))

ociad1_pd_order <- to_plot |>
    group_by(Organelle) |>
    summarise(n = n()) |>
    arrange(n) |>
    mutate(order = 1:4)

genes_to_highlight <- c('OCIAD1', 'PHB1', 'PHB2', 'ABCD3', 'FAR1', 'PEX14')

hao1_vs_ctrl <- ggplot(data = to_plot) +
    geom_vline(xintercept=2, alpha=0.3, linetype = 'dashed')+
    geom_violin(mapping=aes(y=Organelle, x=FC_HAO1_CTRL, fill=color, color=color), alpha = 0.5)+
    geom_point(data=to_plot[to_plot$Gene %in% genes_to_highlight,], 
               mapping = aes(y=Organelle, x=FC_HAO1_CTRL, color=color),size=2)+
    scale_fill_identity()+
    scale_color_identity()+
    scale_y_discrete(limits = ociad1_pd_order$Organelle) +
    scale_x_continuous(breaks = seq(-6,8,2))+
    geom_text(data=ociad1_pd_order, mapping=aes(y=order+0.1, x=-4, label=paste0('n = ',n))) +
    geom_text_repel(data = to_plot, mapping = aes(y=Organelle, x=FC_HAO1_CTRL, label=ifelse(Gene %in% genes_to_highlight, Gene, NA)))+
    theme_bw() +
    labs(y='', x='', title='Log2(FC) HAO1 vs CTRL') +
    theme(plot.title = element_text(hjust=0.5))

hao1_vs_ctrl

ggsave(plot = hao1_vs_ctrl, filename = './proteomics/plots/hao1_vs_ctrl_violin.png',
       width = 12, height = 8)
ggsave(plot = hao1_vs_ctrl, filename = './proteomics/plots/hao1_vs_ctrl_violin.pdf',
       width = 12, height = 8)
