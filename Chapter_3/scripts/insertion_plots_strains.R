#insertion_plots.R
#functions to find and visualise all insertions in a specified region in both tb and bovis 


#' Find the insertions in the range of interest
#' 
#' @param mbi mbovis insertions read from hmm output file.
#' @param mtbi mtb insertions read from hmm output file.
#' @param mtb_i mtb start coord
#' @param mtb_ii mtb end coord
#' @param mb_i mbovis start coord
#' @param mb_ii mbovis end coord
#' @returns Dataframe of insertion coordinates (rows) and number of insertions for Mbovis and Mtb (columns).
#' @example
#' find_ranges(mbov_ins, mtb_ins, 3834596, 3835092, 3791637, 3792133)
#' 
find_ranges <- function(mbi, mtbi, mtb_i, mtb_ii, mb_i, mb_ii){
  mtb_range <- mtbi %>% 
    filter(site %in% mtb_i:mtb_ii)%>% 
    #select(insertions, call, gene, site) %>%
    select(insertions, gene, site) %>%
    mutate(ta_number = c(1:nrow(.)), strain="Mtb")
  mb_range <- mbi %>% 
    filter(site %in% mb_i:mb_ii) %>% 
    #select(insertions, call, site) %>%
    select(insertions, site) %>%
    mutate(ta_number = c(1:nrow(.)), strain="MBovis", gene = mtb_range$gene)
  range_df <- rbind(mb_range, mtb_range) 
  return(range_df)
}

#' Function to create rectangles to represent genes in insertion plot
#' 
#' @param gene_i String of gene of interest (locus name of Mtb ortholog)
#' @param gene_ta_df Dataframe of TA sites and insertions in Mbovis and Mtb 
#' @returns LayerInstance of rectangle line positions for ggplot
#' @example add_text("Rv2221c_(glnE)", insertion_df)
#' 
add_lines <- function(gene_i, gene_ta_df){
  x_min <- min(gene_ta_df %>% 
                 filter(gene==gene_i) %>% select(ta_number))
  x_max <- max(gene_ta_df %>% 
                 filter(gene==gene_i) %>% select(ta_number))
  annotate("rect", 
           xmin = x_min, xmax = x_max,  
           ymin = -55, ymax = -35, fill = "gray") 
}

#' Function to add Mtb locus names to annotated rectangles in insertion plot
#' 
#' @param gene_i String of gene of interest name (locus from Mtb)
#' @param gene_ta_df Dataframe of TA sites and insertions in Mbovis and Mtb 
#' @returns LayerInstance for custom annotation for ggplot
#' @example 
#' add_text("Rv2221c_(glnE)", insertion_df)
#' 
add_text <- function(gene_i, gene_ta_df){
  x_min <- min(gene_ta_df %>% 
                 filter(gene==gene_i) %>% select(ta_number))
  x_max <- max(gene_ta_df %>% 
                 filter(gene==gene_i) %>% select(ta_number))
  annotation_custom(grob = textGrob(gene_i),
                    xmin = x_min, xmax = x_max, 
                    ymin = -80, ymax = -80)
}

#' Function to create a ggplot object to display number of insertions at TA 
#' coordinates in mtb and mbovis on the same plot
#' @param mbovis_insertions
#' @param mtb_insertions
#' @param mtb_start
#' @param mtb_end
#' @param mbovis_start
#' @param mbovis_end
#' 
#' @return ggplot object
#' @examples  
#' insertion_plot(mbov_ins, mtb_ins, 2489073, 2493754, 2472631, 2477326)
#' insertion_plot(mbov_ins, mtb_ins, 906391, 908029, 907181, 908859)
#' 
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
    geom_jitter(aes(x=ta_number, y=insertions, shape = call), size=4, alpha=0.8) +
    scale_color_manual(values=strain_cols) +
    guides(fill=guide_legend(title="insertions"), size=FALSE) +
    coord_cartesian(xlim=c(0, total_tas), 
                    ylim=c(-80,max(pl_df$insertions)), clip="off") +
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
mbov_ins <- read.delim(here("in_vitro_data/hmm_bovis_new_add_19_21.wig"), 
                       sep="\t", header=FALSE, stringsAsFactors=F, 
                       comment.char = '#', col.names = c("site", "insertions", "prob_ES", "prob_GD", "prob_NE", "prob_GA", "call", "gene"))
mtb_ins <- read.delim(here("in_vitro_data/hmm_mtb_new_add_22_23.wig"), 
                      sep="\t", header=FALSE, stringsAsFactors=F, 
                      comment.char = '#', col.names = c("site", "insertions", "prob_ES", "prob_GD", "prob_NE", "prob_GA", "call", "gene"))


gln <- insertion_plot(mbov_ins, mtb_ins, 2489073, 2493754, 2472631, 2477326)
#gln

rv0812 <- insertion_plot(mbov_ins, mtb_ins, 906391, 908029, 907181, 908859)
rv0812

