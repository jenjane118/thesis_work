# Jennifer J. Stiens
# 11 February, 2022

#compare_datasets.R

#script to compare mtb and bcg datasets for cholesterol catabolism genes
library(here)
library(tidyverse)
library(colorspace)

#import bovis data (from old spreadsheet to get cholesterol genes)
bovis_data <- read_csv(here("lung_tnseq/data", "chol_catab_s4.csv"), skip=1)
#col.from <- c("M. bovis locus", "H37Rv Ortholog", "Average log2 Fold Change...7", "Average log2 Fold Change...10", "Average log2 Fold Change (both lungs and nodes)")
#col.to   <- c("bovis_locus", "H37Rv_ortholog", "bovis_lung", "bovis_LN", "bovis_all")
col.from <- c("M. bovis locus", "H37Rv Ortholog")
col.to   <- c("bovis_locus", "H37Rv_ortholog")
bovis_df <- bovis_data %>%
   select(col.from) %>%
   rename_at(vars(all_of(col.from)), ~col.to)

#indicate if gene is 'attenuated' == log2FC < -1.5 and padj<0.05
# these genes are in list
#att_data <- read_csv(here("Output/attenuated_genes_suppl.csv"), skip=1, col_names = T)
#att_genes <- att_data %>% select(`M. bovis locus`) %>% pull()
#bovis_df <- bovis_df %>% mutate(signif = bovis_locus %in% att_genes)

#import resampling dataset
#load("~/git/Mbovis_in-vivo_Tnseq/resamp_data.RData")
#better to work from sharon's spreadsheet to make sure names of genes haven't been edited, etc

# add in pooled data
lung_all <- readRDS(here("lung_tnseq/r_data/lung_pooled.RData"))
node_all <- readRDS(here("lung_tnseq/r_data/node_pooled.RData"))

test_df <- left_join(bovis_df, lung_all, by=c("bovis_locus"="Orf"))
test_df2 <- left_join(test_df, node_all, by=c("bovis_locus"="Orf"), suffix=c(".lung", ".node"))
bovis_pooled.df <- test_df2 %>% 
  select(c(bovis_locus, H37Rv_ortholog, Name.lung, log2FC.lung, log2FC.node)) %>%
  dplyr::rename(Name = Name.lung)
  

#import mendum chol catabolism data from excel sheet 12864_2019_5791_MOESM2_ESM
# cholesterol sheet (copied to text file)
bcg_data <- read_delim(here("lung_tnseq/data", "chol_catab_mendum_supp.txt"), delim="\t")
# added missing 3 Rv numbers that are included in bovis data
col.from <- c("BCG locus", "H37Rv Ortholog", "mean_lf2c_danish", "mean_lf2c_pasteur")
col.to   <- c("BCG_locus", "H37Rv_ortholog", "BCG_danish", "BCG_pasteur")
bcg_df <- bcg_data %>% select(col.from) %>%
  rename_at(vars(all_of(col.from)), ~col.to)
View(bcg_df)

all_data <- left_join(bcg_df, bovis_pooled.df, by="H37Rv_ortholog")
View(all_data)
nrow(all_data)
#86
# import bellarose tb data from two time points (this is only attenuated mutations)
#br_data <- read_csv(here("data/Bellarose_attenuating mutations.csv"), skip=9) 
# 574 but last 12 are NA
#br_data <- br_data[1:562,]

#import bellarose tb data: all of the genes
#compare 14 days post infection and 49 days post infection

library(readxl)
br_ex <- read_excel(here("lung_tnseq/data/Bellerose_2020_msystems.00396-20-st001.xlsx"), skip=9, )

#change formatting to match our data ('RVBD_0007' to 'Rv0007')

col.from <- c("gene", "log2...10", "log2...25")
col.to   <- c("H37Rv_ortholog", "Mtb_early", "Mtb_late")
br_df <- br_ex %>% select(col.from) %>%
  rename_at(vars(all_of(col.from)), ~col.to)
br_df$H37Rv_ortholog <- sub("RVBD_", "Rv", br_df$H37Rv_ortholog)
# round logfold changes to 2 decimal places
br_df$Mtb_early <- round(as.numeric(br_df$Mtb_early), 2)
br_df$Mtb_late <- round(as.numeric(br_df$Mtb_late), 2)

# highlight actual log fold changes (but without considering relative amount 
# will be different betwen datasets?)

all_data <- left_join(all_data, br_df, by="H37Rv_ortholog")
#there are a lot of missing genes in the bellarose dataset? these are attenuating 
#mutations only
nrow(all_data)

# make data longer with differences for each gene
basic_data <- all_data %>% 
  select(bovis_locus, log2FC.lung, log2FC.node, BCG_danish, BCG_pasteur, Mtb_early, Mtb_late, H37Rv_ortholog) %>%
  arrange(bovis_locus) %>%
  mutate(gene_locus = paste(bovis_locus, H37Rv_ortholog, sep="/")) %>%
  dplyr::rename(c(bovis_lung = log2FC.lung, bovis_LN = log2FC.node))
    # paste bovis and tb genes into one string

View(basic_data)


long_data <- pivot_longer(basic_data, 
                          cols=c("bovis_lung", "bovis_LN", "BCG_danish", "BCG_pasteur",
                                 "Mtb_early", "Mtb_late"), 
                          names_to=c("dataset"), values_to=c("log2FC"))
#dictate column order for plotting
dat_order <- c("bovis_lung", "bovis_LN", "BCG_danish", "BCG_pasteur",
               "Mtb_early", "Mtb_late")
axis_names <- c("Mbovis Lung", "Mbovis Node", "BCG Danish", "BCG Pasteur", "Mtb 2-weeks", "Mtb 7-weeks")
long_data$dataset <- factor(long_data$dataset, levels = c(dat_order))

#if signif, color royal blue, else black
#sig_cols <- ifelse(basic_data$signif==TRUE, "blue", "black")

ggplot(long_data, aes(x=dataset, y=gene_locus, fill=log2FC)) +
  geom_tile() + 
  geom_text(aes(x = dataset, y = gene_locus, label = log2FC), size=2, fontface="bold") +
  scale_fill_continuous_sequential(palette = "YlOrRd", name="Log2FC", trans="reverse") +
  #theme(axis.text.y = element_text(size=6, face="bold", colour = sig_cols),
  theme(axis.text.y = element_text(size=6, face="bold"),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8, face="bold"),
        axis.title = element_text(size=8),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels=axis_names) +
  ylab("M.bovis/M.tb locus")

ggsave(here("images/bcg_bovis_compare_pooled.png"), height = 12, width=8) 

#high resolution tiff
ggsave(here("images/bcg_bovis_compare.tiff"),height = 8, units= "in", dpi=700)




# with no nums, need to make with no nums on legend
ggplot(long_data, aes(x=dataset, y=gene_locus, fill=log2FC)) +
  geom_tile() + 
  #geom_text(aes(x = dataset, y = gene_locus, label = log2FC), size=2) +
  scale_fill_continuous_sequential(palette = "YlOrRd", name="Attenuation level", 
                                   trans="reverse",
                                   breaks = c(-5, 0), 
                                   labels=c("High", "Low")) +
  theme(axis.text.y = element_text(size=6, face="bold"),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=5, face="bold"),
        axis.title = element_text(size=8),
        axis.title.x = element_blank()) +
        #legend.text = element_blank()) +
        #legend.title = element_blank()) +
  ylab("M.bovis/M.tb locus")

ggsave(here("images/bcg_bovis_compare_no_nums.png"), height = 12, width=6) 

#high res tiff

ggsave(here("images/bcg_bovis_compare.tiff"), height = 10, width=8, res=600)





# may need to break into two plots for publication
library(gridExtra)
nrow(basic_data)
#86
# break into 43 loci each (43x3=129)
a<-long_data[172:1,]; b<-long_data[nrow(long_data):173,]

# function to extract legend from a ggplot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

#make basic plot WITH legend
pl <- ggplot(a, aes(x=dataset, y=bovis_locus, fill=log2FC)) +
  geom_tile() + 
  scale_fill_continuous_sequential(palette = "YlOrRd", name="Log2FC", trans="reverse") +
  theme(legend.title = element_text(size=8))

#extract legend
leg = g_legend(pl)

#plots with no legends
plota <- ggplot(a, aes(x=dataset, y=bovis_locus, fill=log2FC)) +
  geom_tile() + 
  geom_text(aes(x = dataset, y = bovis_locus, label = log2FC), size=2) +
  scale_fill_continuous_sequential(palette = "YlOrRd", name="Log2FC", trans="reverse") +
  theme(axis.text.y = element_text(size=6, face="bold"),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=5, face="bold"),
        axis.title = element_text(size=8),
        axis.title.x = element_blank()) +
  ylab("M.bovis locus") +
  theme(legend.position = "none")
plotb <- ggplot(b, aes(x=dataset, y=bovis_locus, fill=log2FC)) +
  geom_tile() + 
  geom_text(aes(x = dataset, y = bovis_locus, label = log2FC), size=2) +
  scale_fill_continuous_sequential(palette = "YlOrRd", name="Mean Log2FC", trans="reverse") +
  theme(axis.text.y = element_text(size=6, face="bold"),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=5, face="bold"),
        axis.title = element_text(size=8),
        axis.title.x = element_blank()) +
  ylab("M.bovis locus") +
  theme(legend.position = "none")
#arrange plots and add legend (ncol=3 if want legend next to it)
new_plot <- grid.arrange(plota, plotb, ncol=3, leg)

ggsave(here("images/bcg_bovis_compare_split.png"), plot=new_plot, height = 10, width=8) 

#high resolution pdf

ggsave(here("images/bcg_bovis_compare.pdf"),plot=new_plot, height = 10, width=8, res=300)


