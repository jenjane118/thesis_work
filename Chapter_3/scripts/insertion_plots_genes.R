#functions to create plots of insertions for tn-seq data

library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)

#read in insertions
mbov_ins <- read.delim(here("in_vitro_data/hmm_bovis_new_add_19_21.wig"), 
                       sep="\t", header=FALSE, stringsAsFactors=F, 
                       comment.char = '#', col.names = c("site", "insertions", "prob_ES", "prob_GD", "prob_NE", "prob_GA", "call", "gene"))
mtb_ins <- read.delim(here("in_vitro_data/hmm_mtb_new_add_22_23.wig"), 
                      sep="\t", header=FALSE, stringsAsFactors=F, 
                      comment.char = '#', col.names = c("site", "insertions", "prob_ES", "prob_GD", "prob_NE", "prob_GA", "call", "gene"))

find_ranges <- function(mbi, mtbi, mtb_i, mtb_ii, mb_i, mb_ii){
  require(dplyr)
  #find insertion range of interest
  mtb_range <- mtbi %>% 
    filter(site %in% mtb_i:mtb_ii)%>% 
    select(insertions, call, gene) %>%
    mutate(ta_number = c(1:nrow(.)), strain="Mtb")
  mb_range <- mbi %>% 
    filter(site %in% mb_i:mb_ii) %>% 
    select(insertions, call) %>%
    mutate(ta_number = c(1:nrow(.)), strain="MBovis", gene = mtb_range$gene)
  range_df <- rbind(mb_range, mtb_range) 
  return(range_df)
}

test <- find_ranges(mbov_ins, mtb_ins, 2404166, 2405521, 2389091, 2390446 )


#function to annotate genes
add_lines <- function(gene_i, gene_ta_df){
  require(grid)
  x_min <- min(gene_ta_df %>% 
                 filter(gene==gene_i) %>% select(ta_number))
  x_max <- max(gene_ta_df %>% 
                 filter(gene==gene_i) %>% select(ta_number))
  annotate("rect", 
           xmin = x_min, xmax = x_max,  
           ymin = -40, ymax = -25, fill = "gray") 
}

add_text <- function(gene_i, gene_ta_df){
  require(grid)
  x_min <- min(gene_ta_df %>% 
                 filter(gene==gene_i) %>% select(ta_number))
  x_max <- max(gene_ta_df %>% 
                 filter(gene==gene_i) %>% select(ta_number))
  annotation_custom(grob = textGrob(sub("_.*", "", gene_i)),
                    xmin = x_min, xmax = x_max, 
                    ymin = -60, ymax = -60)
}

insertion_plot <- function(mbovis_insertions, mtb_insertions, mtb_start, mtb_end, mbovis_start, mbovis_end){
  require(ggplot2)
  require(dplyr)
  require(grid)
  require(gridExtra)
  pl_df <- find_ranges(mbovis_insertions, mtb_insertions, 
                       mtb_start, mtb_end, mbovis_start, mbovis_end)
  genes <- pl_df %>% filter(gene != "") %>% distinct(gene) %>% pull()
  no_genes <- length(genes)
  total_tas <- nrow(pl_df %>% filter(strain=="Mtb"))
  gene_tas <- pl_df %>% filter(gene %in% genes & strain=="Mtb")
  strain_cols <- c("MBovis" = "blue", "Mtb" = "darkgreen")
  lines <- lapply(genes, add_lines, pl_df)
  gene_text <- lapply(genes, add_text, pl_df)
  pl <- ggplot(pl_df, aes(colour=strain)) +
    geom_jitter(aes(x=ta_number, y=insertions, shape = call), size=6, alpha=0.8) +
    scale_color_manual(values=strain_cols) +
    guides(fill=guide_legend(title="insertions"), size=FALSE) +
    coord_cartesian(xlim=c(0, total_tas), 
                    ylim=c(-60,max(pl_df$insertions)), clip="off") +
    lines + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 12)) +
    ylab(c("Insertions"))
  gg_pl <- pl + gene_text
  return(gg_pl)
}


g1 <- insertion_plot(mbov_ins, mtb_ins, 2404166, 2405521, 2389091, 2390446)
g1
ggsave(file=here("in_vitro_data/images/test.png"), g1, width = 10)

