#insertion_plots.R
#functions to find and visualise all insertions in a specified region in both tb and bovis 


#' Find the insertions in the range of interest
#' 
#' @param insertion_df with insertions from hmm file for both conditions (merged)
#' @param start start coord
#' @param end 1 end coord
#' @returns Dataframe of insertion coordinates (rows) and number of insertions for both conditions
#' @example
#' find_ranges(men_ins, 3834596, 3835092)
#' 
find_ranges <- function(insertion_df, start, end){
  ins_range <- insertion_df %>% 
    filter(site %in% start:end) %>% 
  return(insertion_df)
}

#' Function to create rectangles to represent genes in insertion plot
#' 
#' @param gene_i String of gene of interest (locus name)
#' @param gene_ta_df Dataframe of TA sites and insertions in treated/untreated
#' @returns LayerInstance of rectangle line positions for ggplot
#' @example add_text("Rv2221c_(glnE)", insertion_df)
#' 
add_lines <- function(gene_i, gene_ta_df){
  x_min <- min(gene_ta_df %>% 
                 filter(gene==gene_i) %>% select(site))
  x_max <- max(gene_ta_df %>% 
                 filter(gene==gene_i) %>% select(site))
  annotate("rect", 
           xmin = x_min, xmax = x_max,  
           ymin = -55, ymax = -35, fill = "gray") 
}

#' Function to add Mtb locus names to annotated rectangles in insertion plot
#' 
#' @param gene_i String of gene of interest name (locus)
#' @param gene_ta_df Dataframe of TA sites and insertions in two conditions
#' @returns LayerInstance for custom annotation for ggplot
#' @example 
#' add_text("Rv2221c_(glnE)", insertion_df)
#' 
library(grid)
add_text <- function(gene_i, gene_ta_df){
  x_min <- min(gene_ta_df %>% 
                 filter(gene==gene_i & condition=="untreated") %>% select(site))
  x_max <- max(gene_ta_df %>%
                 filter(gene==gene_i & condition=="untreated") %>% select(site))
  annotation_custom(grob = textGrob(gene_i),
                    xmin = x_min, xmax = x_max, 
                    ymin = -80, ymax = -80)
}

#' Function to create a ggplot object to display number of insertions at TA 
#' coordinates in mtb and mbovis on the same plot
#' @param merged_df
#' @param cond1_start
#' @param cond1_end
#' 
#' @return ggplot object
#' @examples  
#' insertion_plot(merged_hmm, 3834596, 3835092)
#' 
insertion_plot <- function(insertion_df, start_coord, end_coord){
  require(ggplot2)
  require(dplyr)
  require(grid)
  require(gridExtra)
  pl_df <- find_ranges(insertion_df, start_coord, end_coord)
  genes <- pl_df %>% filter(is.na(gene) != T) %>% distinct(gene) %>% pull()
  no_genes <- length(genes)
  total_tas <- nrow(pl_df %>% filter(condition=="untreated"))
  gene_tas <- pl_df %>% filter(gene %in% genes)
  strain_cols <- c("treated" = "blue", "untreated" = "darkgreen")
  lines <- lapply(genes, add_lines, pl_df)
  gene_text <- lapply(genes, add_text, pl_df)
  pl <- ggplot(pl_df, aes(colour=condition)) +
    geom_jitter(aes(x=site, y=mean_ins, shape=call), size=3, alpha=0.7) +
    scale_color_manual(values=strain_cols) +
    guides(fill=guide_legend(title="insertions"), size=FALSE) +
    coord_cartesian(xlim=c(min(pl_df$site), max(pl_df$site)), 
                    ylim=c(-80,max(test_df$mean_ins)), clip="off") +
    lines + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text.y = element_text(size = 12)) +
    ylab(c("Insertions")) +
    xlab(c("TA site position"))
  gg_pl <- pl + gene_text
  return(gg_pl)
}

######examples######
#read in insertions
# hmm_pos <- read_delim(here("menadione_tnseq/output/transit/transit_hmm_pos_ttr.txt"), delim="\t", comment = "#",  col_names = F, col_select = c(1,2,7))  
# colnames(hmm_pos) = c("site", "mean_ins", "call")
# hmm_pos <- hmm_pos %>% mutate(condition = "treated")
# hmm_neg <- read_delim(here("menadione_tnseq/output/transit/transit_hmm_neg_ttr.txt"), 
#                       delim="\t",
#                       comment="#", col_names = F, col_select = c(1,2,7))
# colnames(hmm_neg) = c("site", "mean_ins", "call")
# hmm_neg <- hmm_neg %>% mutate(condition="untreated")
# merge_hmm <- rbind(hmm_pos, hmm_neg)
#find_ranges(merge_hmm, 3126128, 3127700)
 
#insertion_plot(merge_hmm, 3834596, 3835092)
