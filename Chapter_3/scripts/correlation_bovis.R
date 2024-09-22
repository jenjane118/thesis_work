library(here)

sample_names<-read.table(here("data/names_wig.txt"),stringsAsFactors = F)

#for (j in 1:nrow(sample_names)){
  #insertion_table<-read.csv(paste(sample_names[j,1],".wig",sep = ""))
  #non_norm<-as.data.frame(matrix(0, nrow = nrow(insertion_table), ncol = 3))
  
# If .txt tab file, use this
#B_data <- read.delim(file.choose("B_R1.wig"))  

## non normalised

B__table<-read.csv("B_R1.wig", sep = " ", skip=2)   # positions start on 3rd line
D_table<-read.csv("D_R1.wig", sep = " ", skip=2)
G_table<-read.csv("G_R1.wig", sep = " ", skip=2)

pos_vector<-B_table[,1]
pos_vector

non_norm<-as.data.frame(matrix(0, nrow = nrow(B_table), ncol =3), row.names <-pos_vector)
colnames(non_norm)<-c('B', 'D', 'G')
non_norm[,1]<-B_table[,2]
non_norm[,2]<-D_table[,2]
non_norm[,3]<-G_table[,2]

head(non_norm, 20)

cor.test(non_norm$B, non_norm$G, method = "kendall")

res_non <- cor(non_norm)
round(res_non, 2)

## TTR normalised
library(tidyverse)
data <- read_delim(here("data/bovis_combined_TTR.wig"), delim="\t", skip=7, col_names = F)
combined_ttr <- data %>% select(c(X1, X2, X3, X4))
colnames(combined_ttr) <- c("site", "B", "C", "D")
#B_ttr_table<-read.csv(here("B_R1_ttr.wig"), sep = " ", skip=2)   # positions start on 3rd line
#D_ttr_table<-read.csv(here("D_R1_ttr.wig"), sep = " ", skip=2)
#G_ttr_table<-read.csv(here("G_R1_ttr.wig"), sep = " ", skip=2)


# pos_vector<-B_ttr_table[,1]
# 
# ttr_norm<-as.data.frame(matrix(0, nrow = nrow(B_ttr_table), ncol =3), row.names <-pos_vector)
# colnames(ttr_norm)<-c('B', 'D', 'G')
# ttr_norm[,1]<-B_ttr_table[,2]
# ttr_norm[,2]<-D_ttr_table[,2]
# ttr_norm[,3]<-G_ttr_table[,2]
# 
# head(ttr_norm, 20)

#transform with log10 because pearson's assumes normal distribution

#combined_ttr <- combined_ttr %>% mutate(logB=log10(B), logC=log10(C), logD=log10(D))
#get rid of -Inf
combined_ttr <- combined_ttr %>% mutate(across(B:D, ~ ifelse(.x != 0, log10(.x), 0), 
              .names = '{paste("log", .col, sep="")}'))


ttr_mat <- as.matrix(combined_ttr %>% select(c(logB,logC,logD)))

res_ttr <- cor(ttr_mat)
round(res_ttr, 2)

#corplot
library(corrplot)
corrplot(cor(ttr_mat, method="pearson"))

#might be more effective to use mean/sum insertions per gene instead of indiv insertion sites
gene_b <- read_delim(here("in_vitro_data/bovis_hiseq/hmm_20_05_21/hmm_bovis_19_genes.wig"), skip=3, col_names = T, delim="\t")
gene_b <- rename(gene_b, ORF= `#ORF`) %>% select(ORF, mean)
gene_c <- read_delim(here("in_vitro_data/bovis_hiseq/hmm_20_05_21/hmm_bovis_20_genes.wig"), skip=3, col_names = T, delim="\t")
gene_c <- rename(gene_c, ORF= `#ORF`) %>% select(ORF, mean)
gene_d <- read_delim(here("in_vitro_data/bovis_hiseq/hmm_20_05_21/hmm_bovis_21_genes.wig"), skip=3, col_names = T, delim="\t")
gene_d <- rename(gene_d, ORF= `#ORF`) %>% select(ORF, mean)

gene_df <- tibble(gene_b$mean, gene_c$mean, gene_d$mean)
colnames(gene_df) <- c("b", "c", "d")
gene_df <- gene_df %>% mutate(across(b:d, ~ ifelse(.x != 0, log10(.x), 0), 
                                                      .names = '{paste("log", .col, sep="")}'))
gene_mat <- as.matrix(gene_df[,4:6])

res <- cor(gene_mat, method="pearson")
round(res, 2)
corrplot(cor(gene_mat, method="pearson"))

#better correlation between b and d

library(corrplot)
corrplot(cor(gene_mat))

tn_stats_df<-read.csv("data/bovis_combined_TTR.wig", header=TRUE, sep = "\t", skip=6, row.names=NULL, col.names = c('position', 'B', 'D', 'G', 'gene'))
View(tn_stats_df[1:50,])
# df for sum of reads for each ta site
comb_names<-c("position", "sum reads", "gene")
combined_df<-as.data.frame(matrix(0, nrow = nrow(tn_stats_df), ncol = 3))
combined_df[,1]<-tn_stats_df[,1]
for (i in 1:nrow(combined_df)){
  combined_df[i,2]<-sum(tn_stats_df[i,2:4])
}
combined_df[,3]<-tn_stats_df[,5]
colnames(combined_df)<-comb_names
View(combined_df[0:200,])



# find raw number of ta sites inserted for TTR normalised (sum of reads is not 0)
tot_count <- combined_df[,2]
ins_count<-0
for (i in 1:length(tot_count)){
  if (tot_count[i] != 0.0){
    ins_count <- ins_count + 1
  }
}
ins_count
# 34178
length(tot_count)
# 73536
insertion_density <- ins_count/length(tot_count)
insertion_density
# 0.4647792

# max read count
max(tot_count)
# 100024.5  (uses normalisation factor?)

# mean of no insertions per gene
mean(tot_count)
# with zeros, not useful
# 489.7486


#non zero mean (mean of read counts at all non-zero sites)
# change all zeros to na
is.na(tot_count) <- tot_count==0.0
tot_count[1:100]
mean(tot_count, na.rm=TRUE) 
# 1053.723

# do this again with non_normal

tn_non_stats_df<-read.csv("data/bovis_combined_nonorm.wig", header=TRUE, sep = "\t", skip=6, row.names=NULL, col.names = c('position', 'B', 'D', 'G', 'gene'))
View(tn_non_stats_df[1:20,])

# df for sum of reads for each ta site
comb_names<-c("position", "sum reads", "gene")
combined_non_df<-as.data.frame(matrix(0, nrow = nrow(tn_non_stats_df), ncol = 3))
combined_non_df[,1]<-tn_non_stats_df[,1]
for (i in 1:nrow(combined_non_df)){
  combined_non_df[i,2]<-sum(tn_non_stats_df[i,2:4])
}
combined_non_df[,3]<-tn_non_stats_df[,5]
colnames(combined_non_df)<-comb_names
View(combined_non_df[0:200,])

# find raw number of ta sites inserted for TTR normalised (sum of reads is not 0)
tot_non_count <- combined_non_df[,2]
ins_non_count<-0
for (i in 1:length(tot_non_count)){
  if (tot_non_count[i] != 0.0){
    ins_non_count <- ins_non_count + 1
  }
}

length(tot_non_count)
# 73535   (why is this one less than ttr?)
ins_non_count
# 34178  (exactly same as normalised with TTR?)
# same number of sites inserted but number of counts different
max(tot_non_count)
# 5913

mean(tot_non_count)
# this is with zeros
# 29.54797

#non zero mean 
is.na(tot_non_count) <- tot_non_count==0.0
mean(tot_non_count, na.rm=TRUE) 
# 63.57335





# try with one of the dejesus libraries
dejesus_40_df<-read.csv("SRR4113440_1.wig", sep = " ", skip=2, header=TRUE, row.names=NULL, col.names=c('position', 'count'))
View(dejesus_40_df[1:50,])

# find raw number of ta sites inserted for TTR normalised (sum of reads is not 0)
tot_40 <- dejesus_40_df[,2]
ins_count_40<-0
for (i in 1:length(tot_40)){
  if (tot_40[i] != 0.0){
    ins_count_40 <- ins_count_40 + 1
  }
}

length(tot_40)
# 74603
ins_count_40
# 39247

#insertion density
id_40<-ins_count_40/length(tot_40)
id_40
# 0.526078

max(tot_40)
# 4121

mean(tot_40)
# 29.79875

#non zero mean 
is.na(tot_40) <- tot_40==0.0
mean(tot_40, na.rm=TRUE) 
# 56.64321


# use same methods but get data from bwa mapping (dong's script)
# count number of insertion sites from each sample ('Total No. IS')
# add together
# summary_stats just for last sample. need to read from csv?
tot_bovis_sites<-sum(summary_stats$`No.Insertion sites`)
tot_bovis_sites


#find total number of insertions straight from csv's
features_table<-read.csv("M_bovis_features.csv",header=T)
nrow(features_table)

insertion_results<-data.frame(matrix(0,ncol=(numbersamples+4),nrow=nrow(features_table)))
insertion_results[,1:4]<-features_table[,1:4]
colnames(insertion_results)<-c("START","END","NAME","TYPE", sample_names[,1])
for (i in 1:numbersamples){
  tradis_res<-read.csv(file=paste(sample_names[i,1],"_TraDIS_summary.csv",sep = ""),header=T)
  insertion_results[,i+4]<-tradis_res[,7]
}

tradis_res$No.Insertion.sites[3:20]
View(insertion_results[1:10,])


# remove duplicate positions
total_is <- 0
for (i in 3:nrow(features_table)){
  if (insertion_results[i,4] == 'gene'){
    total_is <- total_is + sum(insertion_results[i,5:7])
  }
}

for (i in 3:nrow(features_table)){
  total_is <- total_is + sum(insertion_results[i, 5:7])
}

insertion_results[1,5:7]

total_is
nrow(insertion_results)
nrow(features_table)
sum(insertion_results[3,5:7])
insertion_results[3,4]
View(insertion_results[3:20,])
max(insertion_results[3:nrow(insertion_results),7])

which.max(insertion_results[3:8292,7])
insertion_results[743,]          
insertion_results[979,]
View(insertion_results[,7])
max(insertion_results[3:8292,7])


# make vectors for each sample insertion positions and then combine
# maybe do this in python?

B_sample<-read.csv("B_S1_TraDIS_summary.csv", header=T)
B_ins_pos<-unlist(B_sample[1,13])  

G_sample<-read.csv("G_S3_TraDIS_summary.csv", header=T)
G_ins_pos<-unlist(G_sample[1,13])

D_sample<-read.csv("D_S2_TraDIS_summary.csv", header=T)
D_ins_pos<-unlist(D_sample[1,13])


total_positions<-append(G_sample, B_sample)

length(total_positions)
total_positions<-append(total_positions, D_sample)
total_positons<-unique(total_positions)

library(dplyr)
data<-read.csv("results/new_bovis_fastqs.stats", header=T)
bovis_stats <- select(data, 1, 7, 8)
bovis_stats
View(bovis_stats)
