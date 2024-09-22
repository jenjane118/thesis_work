## make prot tables for mtb and bovis that only include orthologous genes 
## (Malone et al, 2018)

# m.bovis prot table has two of each entry--just keep the protein entries
# remove entries with 'NA' (gene)

library(rtracklayer)

bovis_or<-read.delim("~/tn_seq/data/mbovis.prot_table", header=F)
head(bovis_or)
## the headings 
colnames(bovis_or)<-c("protein", "start", "end", "strand", "length", "type", "na", "name", "locus")
nrow(bovis_or)
# keep only rows where protein name != NA
bovis_pt<-bovis_or[which(!bovis_or$protein=="NA"),]
nrow(bovis_pt)
#4027
head(bovis_pt)
# write this as a better prot table to use (won't have to find unique if use in transit again)
write.table(bovis_pt, "~/tn_seq/data/bovis_best.prot_table", sep = "\t", col.names = F, quote = F, row.names = F)

## import list of orthologous genes from malone
orths<-read.table("~/tn_seq/data/mb_mtb_orthologs_1.txt", sep=",", stringsAsFactors = F, header=T)
head(orths)
nrow(orths)
#3813

# which of bovis genes also in mtb--make new prot table (bovis_on_tb.prot_table)
bovis_on_tb<-NULL
for (i in 1:nrow(bovis_pt)){
  if (bovis_pt$locus[i] %in% orths[,1]){
    bovis_on_tb<-rbind(bovis_on_tb, bovis_pt[i,])
  }
}
nrow(bovis_on_tb)
#3803 
head(bovis_on_tb)
# this prot table has 10 columns, while bovis only has 9--add one to bovis
bovis_on_tb$na2<-rep("-")
head(bovis_on_tb)
#write.table(bovis_on_tb, "~/tn_seq/data/bovis_on_tb.prot_table", sep="\t", col.names = F, row.names = F, quote=F)


# make tb prot table with ONLY orthologous genes? or can we just leave these out of analysis?

# how long is prot table for mtb
mtb_pt<-read.delim("~/Data/mtbH37Rv.prot_table", header=F)
nrow(mtb_pt)
#3906

# change name on bovis to reflect orthologous gene loci in mtb
# files are same length (based on orthologs)
for (i in 1:nrow(bovis_on_tb)){
  bovis_on_tb$mtb[i]<-orths[which(bovis_on_tb$locus[i] == orths[,1]),2]
}
head(bovis_on_tb)
#rearrange columns, now have tb and bovis locus names in same prot_table
library(dplyr)
new_bovis_on_tb<-select(bovis_on_tb, 1,2,3,4,5,6,7,8,11,9)
head(new_bovis_on_tb)

write.table(new_bovis_on_tb, "~/Data/bovis_on_tb.prot_table", sep="\t", col.names = F, row.names = F, quote=F)

# import results
resamp_results<-read.delim("~/Data/tb_bovis_resampling_output.txt", comment.char = "#", sep="\t", header=F, stringsAsFactors=F)
head(resamp_results)
colnames(resamp_results)<-c("Orf","Name",	"Desc",	"Sites", "Mean_Ctrl",
        "Mean_Exp",	"log2FC",	"Sum_Ctrl",	"Sum_Exp", "Delta_Mean", "p_value",	"Adj_p_value")

signif_10<-resamp_results[which(resamp_results$Adj_p_value<0.10),]
nrow(signif)
#75
signif_05<-resamp_results[which(resamp_results$Adj_p_value<0.05),]
zeros<-signif_05[which(signif_05$Adj_p_value == 0),] #what does zero mean? genes only in tb?
nrow(zeros)  #29
nrow(signif_05)
#62
# which of these have log fold change >2? (Carey et al)
lfc_signif<-signif_05[which(abs(signif_05$log2FC) > 2.0),]
nrow(lfc_signif)
#57

#
## create tb_on_bovis prot table
# read in h37rv prot table
head(mtb_pt)
colnames(mtb_pt)<-c("protein", "start", "end", "strand", "length", "type", "na", "name", "locus")
nrow(mtb_pt)
# which of bovis genes also in mtb--make new prot table (bovis_on_tb.prot_table)
tb_on_bovis<-NULL
for (i in 1:nrow(mtb_pt)){
  if (mtb_pt$locus[i] %in% orths[,2]){
    tb_on_bovis<-rbind(tb_on_bovis, mtb_pt[i,])
  }
}
nrow(tb_on_bovis)
#3643
head(tb_on_bovis)
# change name on bovis to reflect orthologous gene loci in bovis
# files are same length (based on orthologs)
for (i in 1:nrow(tb_on_bovis)){
  tb_on_bovis$bovis[i]<-orths[which(tb_on_bovis$locus[i] == orths[,2]),1]
}
new_tb_on_bovis<-select(tb_on_bovis, 1,2,3,4,5,6,7,8,11,9)
head(new_tb_on_bovis)
write.table(new_tb_on_bovis, "~/Data/tb_on_bovis.prot_table", sep="\t", col.names = F, row.names = F, quote=F)

# import results
rev_resamp_results<-read.delim("~/Data/rev_resampling_output.txt", comment.char = "#", sep="\t", header=F, stringsAsFactors=F)
head(rev_resamp_results)
colnames(rev_resamp_results)<-c("Orf","Name",	"Desc",	"Sites", "Mean_Ctrl",
                                "Mean_Exp",	"log2FC",	"Sum_Ctrl",	"Sum_Exp", "Delta_Mean", "p_value",	"Adj_p_value")

# with bovis as control/tb as exp
rev_signif_10<-rev_resamp_results[which(rev_resamp_results$Adj_p_value<0.10),]
nrow(rev_signif_10)
#82
rev_signif_05<-rev_resamp_results[which(rev_resamp_results$Adj_p_value<0.05),]
#signif_05<-signif_05[which(signif_05$Adj_p_value > 0),] what does zero mean?
nrow(rev_signif_05)
#68
zeros<-rev_resamp_results[which(rev_resamp_results$Adj_p_value==0),]
nrow(zeros)  #30
# which of these have log fold change >2? (Carey et al)
rev_lfc_signif<-rev_signif_05[which(abs(rev_signif_05$log2FC) > 2.0),]
nrow(rev_lfc_signif)
#63
lfc2_signif_10<-rev_signif_10[which(abs(rev_signif_10$log2FC) > 2.0),]
#72
nrow(lfc2_signif_10)


write.table(rev_lfc_signif, "bovis_tb_resampling.txt", quote = F, row.names = F)



## make barchart to show difference in saturation vs statistical power
# library
library(ggplot2)

# create a dataset
resample<-matrix(0, 6, 3)

library <- c(rep("dj_43", 2), rep("dj_67",2), rep("hiseq_tb",2))
condition <- c("total lfc", "pv0.05", "total lfc", "pv0.05", "total lfc", "pv0.05")
#sig_df$total_lfc<-c(1107, 581, 773)
#sig_df$pv0.05<-c(26, 146, 63)

values<-c(1107, 26, 581, 146, 773, 63)

sig_df<-NULL
sig_df<-data.frame(library, condition, values)
sig_df

# Grouped
#ggplot(sig_df, aes(fill=condition, y=value, x=library)) + 
#  geom_bar(position="fill", stat="identity")

# Basic barplot
p<-ggplot(data=sig_df, aes(x=library, y=total_lfc)) +
  geom_bar(stat="identity", width=0.5, color="blue", fill="white")
p

p<-ggplot(data=sig_df, aes(x=library, y=total_lfc)) +
  geom_bar(stat="identity", width=0.5, fill="steelblue")+
  theme_minimal()
p

p<-ggplot(data=sig_df, aes(x=library, y=total_lfc, fill=pv0.05)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=c("log fold change > 2")), vjust=1.6, 
            color="white", size=3.5)
  #scale_fill_brewer(palette="Paired")+
  #theme_minimal()
p

# Stacked barplot with multiple groups
ggplot(data=sig_df, aes(x=library, y=values, fill=condition)) +
  geom_bar(stat="identity", position=position_dodge())+
  ylab("log fold change > 2")
  
# this just stacks on top, I want to overlap. dodge puts next to it


## venn diagram to show if significant genes are same in all datasets


rev_lfc_signif
head(dj39_lfc_signif)


devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
a<-dj_lfc_signif$Orf
b<-tb_signif
c<-c(dj39_lfc_signif$Orf)
x<-list(dj_67=a, hiseq_tb=b, dj_43=c)
ggVennDiagram(x) + scale_fill_gradient(low="orange",high = "red")


