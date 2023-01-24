require(msigdbr)
require(fgsea)
require(ggplot2)
require(readr)
require(tibble)
require(dplyr)
GSEA <- function(diftable,
                 padj,
                 uplfc,
                 downlfc,
                 ntop)
  {
    DEGs <- diftable[which(diftable$p_val_adj < padj & abs(diftable$avg_log2FC) > 1), ]
    up <- DEGs[which(DEGs$avg_log2FC > uplfc), ]
    down <- DEGs[which(DEGs$avg_log2FC < downlfc), ]
    up_down <- rbind(up, down)
    H <- msigdbr(species = "Homo sapiens", category = 'H') 
    H.symbol.ls <- H %>% 
      dplyr::select(gs_name, gene_symbol) %>% 
      group_by(gs_name) %>% 
      summarise(all.genes = list(unique(gene_symbol))) %>% 
      deframe()
    FC <- up_down[,c(2,7)]
    FC.vec <- FC$avg_log2FC
    names(FC.vec) <- FC$gene
    gsea.H <- fgseaSimple(pathways = H.symbol.ls,
                          stats = FC.vec,
                          scoreType = "std",
                          nperm=1000)
    
    my_pal <- c("#1B9E77","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#9A7D0A")
    gsea.H$pathway <- gsub("HALLMARK_","", gsea.H$pathway)
    gsea.H$pathway <- gsub("_"," ", gsea.H$pathway)
    gsea.H <- gsea.H %>%
      mutate(Significance = "NotSignificance") %>%
      mutate(Significance = ifelse(padj <= padj & NES > 0, "Upregulated", Significance)) %>%
      mutate(Significance = ifelse(padj <= padj & NES < 0, "Downregulated", Significance))
    gsea.H <- gsea.H[which(gsea.H$Significance != "NotSignificance"), ]
    gsea.H <- na.omit(gsea.H)
    gsea.H <- gsea.H[1:ntop, ]
    p <- gsea.H %>%
      ggplot(aes(x=reorder(pathway, NES),
                 y=NES)) +
      geom_col(width = 0.3, aes(color= Significance, fill = Significance)) +
      labs(x= 'Gene set', y= "Normalized enrichment score (NES)") +
      theme_classic() +
      scale_color_manual(values=c("#1B9E77", "#7570B3")) +
      scale_fill_manual(values=c("#1B9E77", "#7570B3")) +
      coord_flip() +
      theme(axis.text = element_text(family = "Times",size = 24 , colour = "black"),
            axis.text.x = element_text(family = "Times",colour = "black", size = 11),
            axis.text.y = element_text(family = "Times",colour = "black", size = 10),
            plot.subtitle = element_text(family = "Times",size = 16, colour = "black", hjust = 0.5),
            axis.title.y = element_text(family = "Times", size = 16, angle = 90),
            axis.title.x = element_text(family = "Times", size = 16, angle = 00),
            legend.text = element_text(size = 10, family = "Times"),
            legend.title = element_text(size = 20, family = "Times")) +
      labs(subtitle = paste0("GSEA"))
    return(p)
}
diftale <- read_csv("DEGs_up.csv")
GSEA(diftable = diftale, padj = 0.05, uplfc = 1, downlfc = -1, ntop = 10)

