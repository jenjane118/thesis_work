## GENE/TA SITE MAP

# 30_09_2021

# making matrix of TA site coordinates mapped to individual genes in prot table

# import prot table for bovis

# import list of TA sites from any tpp wig file (input library?)

#end up with Orf name: list of TA site coordinates

library(GenomicRanges)
library(tidyverse)

bovis_table <- read_tsv("Mbovis_LT708304.prot_table", col_names = F)
colnames(bovis_table) <- c("Descr", "start", "end", "strand", "len", "type", "-", "name", "Orf")

gr <- GRanges(
  seqnames = "LT708304",
  ranges = IRanges(bovis_table$start, end = bovis_table$end),
  strand = Rle(strand(bovis_table$strand)))
values(gr) <- DataFrame(orf = bovis_table$Orf)
gr_df <- as_tibble(gr)
gr$orf

# read in TA sites from wig file (with non-permissive sites removed)
wig <- read_delim("lung_tnseq/tpp_genewiz2/np_removed/perm_lung_tpp_MbA27.wig", 
                  skip=2, col_names = F, delim=" ")
ta_sites <- wig[,1]
ta_gr <- GRanges(
  seqnames = "LT708304",
  ranges = IRanges(start=ta_sites$X1, end = ta_sites$X1),
  strand = "*"
)
ta_df <- as.data.frame(ta_gr)
# find overlaps with genes in prot table

ov <- findOverlaps(gr, ta_gr, type="any", ignore.strand=T)
ta_gr[subjectHits(ov)]
subjectHits(ov)
# these are the index of gene in gr that each ta site maps to.
matches <- as.tibble(gr[queryHits(ov)])

# find indices range of each orf--these are number of insertion sites
# in each gene, but not indices
#c_ov <- countOverlaps(gr, ta_gr, type="any", ignore.strand=T)

#gene_ta <- as.tibble(matrix(0, nrow=length(gr), ncol=2))
#ta_df$start[1:c_ov[1]]
#ta_df$start[c_ov[1]+1:c_ov[2]]
#ta_df$start[c_ov[2] + 1: c_ov[3]]

gene_list <- gr_df$orf


# gene_ta <- tibble()
# gene_ta <- add_column(gene_ta, gene="", ta_list=list())
# ta <- c(ta_df$start[1:c_ov[1]])
# gene_ta <- gene_ta %>% add_row(tibble_row(gene = gene_list[1], ta_list = list(ta)))
# for (i in 2:length(gr)){
#   ta <- c(ta_df$start[c_ov[i-1]+1:c_ov[i]])
#   gene_ta <- gene_ta %>% add_row(tibble_row(gene = gene_list[i], ta_list = list(ta)))
# }
# 
# gene_ta$ta_list[1]
# gene_ta$gene[1]
# # make into a list of lists
# ta_sites <- gene_ta$ta_list
# names(ta_sites) <- gene_ta$gene
# 
# 
# 
# # unnests whole tibble (gene name, each ta site)
# ta_gene_list <- as.matrix(gene_ta %>% unnest(ta_list))
# 
# #gene and list of ta sites
# gene_ta %>% group_by(gene, ta_list) %>% distinct(gene)
# 

# better to use count overlaps directly and apply to list of ta sites in python?
# this is flawed because cant count through wigs because there are intergenic sites
#ov_df <- bind_cols(gene_ta$gene, c_ov)
#colnames(ov_df) <- c("gene", "ta_count")

write_csv(ov_df, "lung_tnseq/gene_ta_list.csv")



