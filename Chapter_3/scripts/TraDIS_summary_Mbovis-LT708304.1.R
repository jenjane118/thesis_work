# from Dong Xia adapted from BioTraDIS


# if (!require(Rsamtools)){
 # BiocManager::install("Rsamtools"))

# }


options(scipen=999,stringsAsFactors = F)
library(Rsamtools)
features_table<-read.csv("M_bovis_features.csv",header=T)
sample_names<-read.table("names.txt",stringsAsFactors = F)

what<-c("pos")
for (j in 1:nrow(sample_names)){
  summary_stats<-as.data.frame(matrix(0, nrow = nrow(features_table), ncol = 15))
  summary_stats[,1:3]<-features_table[,3:5]
  summary_stats[,4:5]<-features_table[,1:2]
  colnames(summary_stats)<-c("Name","Type","Product","Start","End","Gene Length","No.Insertion sites", "Insertion Index", "Av. Coverage","No.Insertion sites DP>5", "Insertion Index DP>5", "Av. Coverage DP>5","Pos. IS","No.reads per IS","Total no. IS")
  for (i in 1:nrow(features_table)){
    if (summary_stats$End[i]>summary_stats$Start[i]){
    summary_stats$`Gene Length`[i]<-(as.numeric(summary_stats$End[i])-as.numeric(summary_stats$Start[i]))+1
range1 <- IRanges(features_table$START[i],features_table$END[i])
which<-IRangesList(LT708304.1=range1)

	} else {
      summary_stats$`Gene Length`[i]<-(as.numeric(summary_stats$Start[i])-as.numeric(summary_stats$End[i]))+1
range2<- IRanges(features_table$END[i],features_table$START[i])
which<-IRangesList(LT708304.1=range2)
	}
    
    # Use ScanBamParam() to create a parameter object influencing what fields and which 
    # records are imported from a (binary) BAM file. Use of which requires that a BAM index 
    # file (<filename>.bai) exists
    #
    # what Object of class character indicating what fields are to be returned.
    # which Object of class IntegerRangesList indicating which reference sequence 
    # and coordinate reads must overlap.
    param<-ScanBamParam(what=what,which = which)
    
    # scanBam:Import binary ‘BAM’ files into a list structure, with facilities for selecting 
    # what fields and which records are imported, and other operations to manipulate 
    # BAM files.
    reads<-scanBam(paste(sample_names[j,1],".sort.bam",sep = ""),param = param)[[1]]
    # unique returns a vector with duplicate elements removed, so this is removing duplicate
    # positions
    a<-unique(reads$pos)
    # doesn't include positions with no reads (a>1)
    if (length(a>=1)){
      summary_stats$`No.Insertion sites`[i]<-length(a)
      summary_stats$`Pos. IS`[i]<-paste(a,collapse = ",")
      summary_stats$`Insertion Index`[i]<-as.numeric(summary_stats$`No.Insertion sites`[i])/as.numeric(summary_stats$`Gene Length`[i])
      insertionsite<-as.data.frame(matrix(NA, nrow = 1, ncol = length(a)))
      for (k in 1:ncol(insertionsite)){
        insertionsite[1,k]<-length(which(is.element(reads$pos,a[k])))
      }
      summary_stats$`Av. Coverage`[i]<-median(as.numeric(insertionsite[1,]))
      summary_stats$`No.reads per IS`[i]<-paste(insertionsite[1,],collapse = ",")
      summary_stats$`Total no. IS`[i]<-sum(as.numeric(insertionsite[1,]))
      b<-which(insertionsite[1,]>4)
      if (length(b>=1)){
        insertionsite_more5<-as.data.frame(insertionsite[,b])
        summary_stats$`No.Insertion sites DP>5`[i]<-ncol(insertionsite_more5)
        summary_stats$`Insertion Index DP>5`[i]<-as.numeric(summary_stats$`No.Insertion sites DP>5`[i])/as.numeric(summary_stats$`Gene Length`[i])
        summary_stats$`Av. Coverage DP>5`[i]<-median(as.numeric(insertionsite_more5[1,]))
      }
    }
  }
  # write csv of results for each of 3 samples
  write.csv(summary_stats,file=paste(sample_names[j,1],"_TraDIS_summary.csv",sep = ""),row.names = F)
}

  
# create data frame with samples and insertion index for each gene
numbersamples<-nrow(sample_names)
index_results<-data.frame(matrix(0,ncol=(numbersamples*2)+3,nrow=nrow(features_table)))
index_results[,1:3]<-features_table[,3:5]
colnames(index_results)<-c("NAME","TYPE","PRODUCT",rep(sample_names[,1],2))
for (i in 1:numbersamples){
  tradis_res<-read.csv(file=paste(sample_names[i,1],"_TraDIS_summary.csv",sep = ""),header=T)
  index_results[,i+3]<-tradis_res$Insertion.Index
  index_results[,i+(numbersamples+3)]<-tradis_res$Insertion.Index.DP.5
}

nrow(index_results)



index_results_average<-data.frame(matrix(0,ncol=3,nrow=nrow(features_table)))
index_results_average[,1]<-features_table[,3]
colnames(index_results_average)<-c("NAME","Av.Index","Av.Index.DP5")

# average is mean of av index for all 3 samples
for (i in 1:nrow(index_results)){
  index_results_average[i,2]<-mean(as.numeric(as.character(index_results[i,4:(numbersamples+3)])))
  index_results_average[i,3]<-mean(as.numeric(as.character(index_results[i,(numbersamples+4):ncol(index_results)])))
}
index_results_average[0:100,]

jpeg('Insertion.index.jpg')
insertion_graph<-hist(index_results_average[,2],xlim = c(0,0.05),ylim = c(0,50),nclass = 700,xlab = "Insertion Index",ylab = "Density",main = NA)
dev.off()

# make graph for average coverage and av coverage dp5
index_results<-data.frame(matrix(0,ncol=(numbersamples*2)+3,nrow=nrow(features_table)))
index_results[,1:3]<-features_table[,3:5]
colnames(index_results)<-c("NAME","TYPE","PRODUCT",rep(sample_names[,1],2))
for (i in 1:numbersamples){
  tradis_res<-read.csv(file=paste(sample_names[i,1],"_TraDIS_summary.csv",sep = ""),header=T)
  index_results[,i+3]<-tradis_res$Av..Coverage
  index_results[,i+(numbersamples+3)]<-tradis_res$Av..Coverage.DP.5
}

index_results_average<-data.frame(matrix(0,ncol=3,nrow=nrow(features_table)))
index_results_average[,1]<-features_table[,3]
colnames(index_results_average)<-c("NAME","Av.Index","Av.Index.DP5")

for (i in 1:nrow(index_results)){
  index_results_average[i,2]<-mean(as.numeric(as.character(index_results[i,4:(numbersamples+3)])))
  index_results_average[i,3]<-mean(as.numeric(as.character(index_results[i,(numbersamples+4):ncol(index_results)])))
}

# graph of average coverage per gene, all reads
jpeg('av_coverage.jpg')
coverage_graph<-hist(index_results_average[,2],xlim = c(0,60),nclass = 300,xlab = "Av.coverage.per.gene",ylab = "Density",main = NA)
dev.off()

# make insertion graphs for each sample
what<-c("pos")
which<-IRangesList(LT708304.1=IRanges(1,4349904))
for (j in 1:nrow(sample_names)){
  param<-ScanBamParam(what=what,which = which)
  reads<-scanBam(paste(sample_names[j,1],".sort.bam",sep = ""),param = param)
  positions<-reads[[1]]$pos
  jpeg(paste0(sample_names[j,1],'_insertion.jpg'))
  hist(positions,breaks = 4349904/100)
  dev.off()
}

