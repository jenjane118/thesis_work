#install.packages("circlize")
#setwd("~/tn_seq/invitro_data/")

library(here)
library(circlize)
library(dplyr)

f.wig <- read.table(file=here("in_vitro_data/add_19_21.wig"), header = FALSE, sep="\t", skip =2)
df <- cbind(f.wig[,1],f.wig) #adds extra column with dupl position for start and end
colnames(df) <- c("start","end", "value")
df$logvalue <- log10((df$value+0.01))
quantile(df[,3])
ymax.f <- quantile(df[,3], probs=c(0.9))
df$value <- ifelse(df$value>ymax.f,ymax.f,df$value)
quantile(df[,3])
ymax.f

#read in the genes from the bed file
beddata<- read.table(file=here("in_vitro_data/circos_test/LT708304_updated.bed"), header=FALSE, sep="\t")
#keep the gene id
beddata$geneid <- gsub(pattern=".*_","",beddata$V5)
#keep the gene name (where there isn't one, the locus tag that folllows will take the place of the gene name)
beddata$genename <- gsub(pattern="gene=|;locus_tag.*", "", beddata$V5)
beddata$genename <- gsub(pattern="locus_tag.*_", "", beddata$genename)
head(beddata)
#only one genome displayed around the whole circle
bovis = c("Mbovis_LT708304.1")
colnames(beddata)[1] <- "chr"
#region.size is the number of bases shown in the circle
#region.size <- 10000
region.size <- 4349904
#how many genes in bedplus and bedminus can be shown within the region.size?
numgenes <- dim(subset(beddata, beddata[,3] < region.size))[1]
#limit the plotting of the wig to the region size we are showing
endregion <- dim(subset(df, df[,2]< region.size))[1]

# make another track inside to label gaps in coverage (pos=0 for >100 conseq)
#start_gaps<-c(1000, 3000, 4000, 5000)
#end_gaps<-c(1200, 3200, 4200, 5200)
start_gaps<-c(343084,807982, 1454144,1479475)
end_gaps<-c(357910,801766, 1448969, 1474045)
mids<-(end_gaps-start_gaps)/2+start_gaps
name_l<-c("a","b","c","d")
gap_df<-data.frame(start_gaps, end_gaps, name_l)
head(gap_df)


circos.clear()

png(here("images/test_circos2.png"), width = 960, height=960)

#1) initialise plot
circos.initialize(sectors=beddata[1:numgenes,"chr"], xlim=c(0,region.size))

#2) plot genome
circos.trackPlotRegion(factors=beddata[1:numgenes,"chr"], 
                       ylim=c(0,500), 
                       track.height = 0.3,
                       panel.fun=function(x,y){
                         circos.genomicLines(region=df[1:endregion,1:2], 
                                             value=df[1:endregion,], 
                                             numeric.column=3, 
                                             type = "h",
                                             col="blue")
                          circos.yaxis(labels.cex=0.75)
                       })

#3) draw axis with nt size markers around outside of circle
circos.genomicAxis(
  h = "top",
  major.at = NULL,
  labels = NULL,
  major.by = NULL,
  tickLabelsStartFromZero = TRUE,
  #labels.cex = 0.5*par("cex"),
  labels.cex = par("cex"),
  sector.index = get.current.sector.index(),
  track.index = get.current.track.index())

#4) plot gaps in coverage with rectangles
circos.genomicTrackPlotRegion(beddata[1:numgenes,], ylim=c(0,1), track.height= 0.1,
                               panel.fun=function(x,y, ...){
                                 chr = get.cell.meta.data("sector.index")
                                 xlim1 = gap_df[,1]
                                 xlim2 = gap_df[,2]
                                 ylim = get.cell.meta.data("ylim")
                                 circos.rect(xlim1, 0, xlim2, 0.9, col="gray")
                                 circos.text(mids, -0.5, name_l, niceFacing=F, cex=0.75)
                             }
)


dev.off()


#***************************************************************************

# make another for tb

mtb_genome <- c("NC_000962.3")
#mtb_region.size <- 10000  #test size
mtb_region.size <- 4411532

#awk '{ print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' mycobrowser_H37Rv_gff_v3.gff > NC000962.bed

## getting weird plot--plotting two sectors. maybe it's my bedfile?
# make a new .bed file for tb
#gff2bed < MtbH37RvNC_000962.3.gff > MtbH37RvNC_000962_3.bed


tb.wig <- read.table(file="add_22_23.wig", header = FALSE, sep="\t", skip =2)
tb_df <- cbind(tb.wig[,1],tb.wig) #adds extra column with dupl position for start and end
head(tb_df)
colnames(tb_df) <- c("start","end", "value")
tb_df$logvalue <- log10((tb_df$value+0.01))
quantile(tb_df[,3])
tb_ymax.f <- quantile(tb_df[,3], probs=c(0.9))
# replace max value in insertions by 90% quartile
tb_df$value <- ifelse(tb_df$value>tb_ymax.f,tb_ymax.f,tb_df$value)
quantile(tb_df[,3])
tb_ymax.f

#read in the genes from the bed file (need to change this for mtb)
beddata_tb<- read.delim(file="~/Data/NC000962.bed", header=FALSE, sep="\t")
head(beddata_tb)
#keep the gene id
beddata_tb$geneid <- gsub(pattern="comments=","",beddata_tb$V5)
beddata_tb$geneid <- gsub(pattern="%2C", "", beddata_tb$geneid)
#beddata_tb<- subset(beddata_tb, select = -c(5))
colnames(beddata_tb)[1] <- "chr"
head(beddata_tb)

numgenes <- dim(subset(beddata_tb, beddata_tb[,3] < mtb_region.size))[1]
numgenes   # number of genes in region
endregion <- dim(subset(tb_df, tb_df[,2]< mtb_region.size))[1]
endregion  # number of ta sites in region

# make another track inside to label gaps in coverage (pos=0 for >100 conseq)
#find gaps for tb
#342082-351473 
#799998-806159	 (6161) 
#1446511-1453029	(6518) (argS, lysA, thrACB)
#1471642-1476251	(4609) (ribosomal rna)
#1552610-1561266 (8656) (pyrR,B,C,F,carA,B,Rv1382)
##1646223-1652592	(6369) (various genes)
#4238041-4245025
start_tb<-c(342082,799998,1446511,1471642,1552610,1646223,4238041)
end_tb<-c(351473,806159, 1453029, 1476251,1561266,1652592, 4245025)
mid_tb<-(end_tb-start_tb)/2+start_tb
mid_tb
name_tb<-c("a","b","c","d","e","f","g")
gap_tb<-data.frame(start_tb, end_tb, name_tb)
gap_tb


circos.clear()

png("tb_circos.png", width = 960, height=960)

circos.initialize(sectors=beddata_tb[1:numgenes,"chr"], xlim=c(0,mtb_region.size))

circos.trackPlotRegion(factors=beddata_tb[1:numgenes,"chr"], 
                       ylim=c(0,200), 
                       track.height = 0.3,
                       panel.fun=function(x,y){
                         circos.genomicLines(region=tb_df[1:endregion,1:2], 
                                             value=tb_df[1:endregion,], 
                                             numeric.column=3, 
                                             type = "h",
                                             col="green")
                         circos.yaxis(labels.cex=0.75)
                       })
                      
circos.genomicAxis(
  h = "top",
  major.at = NULL,
  labels = NULL,
  major.by = NULL,
  tickLabelsStartFromZero = TRUE,
  labels.cex = par("cex"),
  sector.index = get.current.sector.index(),
  track.index = get.current.track.index())

circos.genomicTrackPlotRegion(beddata2[1:numgenes,], ylim=c(0,1), track.height= 0.1,
                              panel.fun=function(x,y, ...){
                                chr = get.cell.meta.data("sector.index")
                                xlim1 = gap_tb[,1]
                                xlim2 = gap_tb[,2]
                                ylim = get.cell.meta.data("ylim")
                                circos.rect(xlim1, 0, xlim2, 0.9, col="gray")
                                circos.text(mid_tb, -0.5, name_tb, niceFacing=T, cex=0.75)
                              }
)
dev.off()




## these are for plotting and labelling individual genes 
# circos.genomicTrackPlotRegion(bedminus[1:numgenes.minus,], ylim=c(0,1),track.height=0.1,
#                               panel.fun=function(x,y, ...){
#                                 chr = get.cell.meta.data("sector.index")
#                                 xlim1 = bedminus[1:numgenes.minus,2]
#                                 xlim2 = bedminus[1:numgenes.minus,3]
#                                 ylim = get.cell.meta.data("ylim")
#                                 circos.rect(xlim1, 0, xlim2, 0.9, col=palette.calls[bedminus[,8]])
#                                 
#                              }
#)
#not enough space for these labels too!
#circos.genomicLabels(bedminus[1:numgenes.minus,],labels.column=7,cex=0.35,col="black",line_lwd=0.5,line_col="grey80",
#                                         side="inside",connection_height=0.05,labels_height=0.04)



# #limit the plotting of the wig to the region size we are showing
# endregion <- dim(subset(df.r, df.r[,2]< region.size))[1]
# circos.trackPlotRegion(factors=beddata[1:numgenes,"chr"], 
#                        ylim=c(0,abs(ymax.r)), 
#                        panel.fun=function(x,y){
#                          circos.genomicLines(region=df.r[1:endregion,1:2], 
#                                              value=df.r[1:endregion,], 
#                                              numeric.column=3, 
#                                              type = "h",
#                                              baseline="bottom", 
#                                              col="red")
#                         
#                        circos.yaxis(labels.cex=0.3)
#                        })








