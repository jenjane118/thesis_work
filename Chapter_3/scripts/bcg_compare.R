#script to compare mendum bcg logfold changes for cholesterol catabolism genes
library(here)
library(tidyverse)
library(colorspace)

#import bovis data
bovis_data <- read_csv(here("lung_tnseq/data", "chol_catab_s4.csv"), skip=1)
col.from <- c("M. bovis locus", "H37Rv Ortholog", "Average log2 Fold Change...7", "Average log2 Fold Change...10", "Average log2 Fold Change (both lungs and nodes)")
col.to   <- c("bovis_locus", "H37Rv_ortholog", "bovis_lung", "bovis_LN", "bovis_all")
bovis_df <- bovis_data %>%
            select(col.from) %>%
            rename_at(vars(all_of(col.from)), ~col.to)
View(bovis_df)            

#import resampling dataset
#load("~/git/Mbovis_in-vivo_Tnseq/resamp_data.RData")
#better to work from sharon's spreadsheet to make sure names of genes haven't been edited, etc
  
#import mendum chol catabolism data from excel sheet 12864_2019_5791_MOESM2_ESM
# cholesterol sheet (copied to text file)
bcg_data <- read_delim(here("lung_tnseq/data", "chol_catab_mendum_supp.txt"), delim="\t")
# added missing 3 Rv numbers that are included in bovis data

col.from <- c("BCG locus", "H37Rv Ortholog", "mean_lf2c_danish", "mean_lf2c_pasteur")
col.to   <- c("BCG_locus", "H37Rv_ortholog", "danish", "pasteur")
bcg_df <- bcg_data %>% select(col.from) %>%
          rename_at(vars(all_of(col.from)), ~col.to)
View(bcg_df)

all_data <- left_join(bcg_df, bovis_df, by="H37Rv_ortholog")
View(all_data)

# import bellarose tb data from two time points




# create row of absolute value of difference between bovis and danish
# Decided to look at mean lfc instead
#all_data <- all_data %>% mutate(diff_danish = bovis_all-danish) %>% mutate(diff_pasteur = bovis_all-pasteur)

#ggplot heatmap with geom_tile()
ggplot(all_data, aes(x=1, H37Rv_ortholog, fill=bovis_all-danish)) +
  geom_tile(aes(width = 0.5)) +   #width setting not working
  scale_fill_viridis(name = "Log2FC difference") +
  xlab("BCG Danish") +
  theme(axis.text.y = element_text(size=3, face="bold"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())


# try heatmap with two datasets
# make data longer with differences for each gene
basic_data <- all_data %>% 
  select(H37Rv_ortholog, diff_danish, diff_pasteur)

long_data <- pivot_longer(basic_data, cols=c("diff_danish", "diff_pasteur"), 
             names_to=c("BCG_strain"), values_to=c("log2FC_diff"))

ggplot(long_data, aes(x=BCG_strain, y=H37Rv_ortholog, fill=log2FC_diff)) +
  geom_tile() + 
  #scale_fill_viridis(name = "Log2FC difference", option="inferno", trans='reverse') +
  scale_fill_continuous_diverging(palette = "Purple-Green", trans='reverse', 
                                  name = "Log2FC difference") + 
  theme(axis.text.y = element_text(size=6, face="bold"),
        axis.ticks = element_blank()) 
#ggsave(here("images/bcg_bovis_compare3.png"), height = 10, width=5)  


# highlight actual log fold changes (negative green, positive purple) in each dataset
# this might be difficult as relative amount will be different betwen datasets?
# have three columns, one for each dataset
# genes on side are mbovis genes

# make data longer with differences for each gene
basic_data <- all_data %>% 
  select(bovis_locus, bovis_lung, bovis_LN, danish, pasteur) %>%
  arrange(bovis_locus)

long_data <- pivot_longer(basic_data, cols=c("bovis_lung", "bovis_LN", "danish", "pasteur"), 
                          names_to=c("dataset"), values_to=c("log2FC"))

ggplot(long_data, aes(x=Strain, y=bovis_locus, fill=log2FC)) +
  geom_tile() + 
  geom_text(aes(x = Strain, y = bovis_locus, label = log2FC), size=2) +
  #scale_fill_viridis(name = "Log2FC", option="inferno", trans='reverse') +
  #scale_fill_continuous_diverging(palette = "Purple-Green", trans='reverse', 
  #                                name = "Log2FC") + 
  scale_fill_continuous_sequential(palette = "YlOrRd", name="Log2FC", trans="reverse") +
  theme(axis.text.y = element_text(size=6, face="bold"),
        axis.ticks = element_blank()) +
  ylab("M.bovis locus") 

ggsave(here("images/bcg_bovis_compare.png"), height = 12, width=3) 


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

ggsave(here("images/bcg_bovis_compare2.png"), plot=new_plot, height = 10, width=8) 
