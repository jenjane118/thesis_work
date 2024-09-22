
#adapted from https://github.com/sanger-pathogens/Bio-Tradis/blob/master/bin/tradis_essentiality.R

########################################################
#Actual code below#
########################################################



options(scipen=999,stringsAsFactors = F)



library("MASS")


STM_baseline <- read.table("results/BCD_bothBatch_TraDIS_summary.csv", sep=",",header=TRUE,stringsAsFactors=F, quote="\"")



#ii <- STM_baseline$Insertion.Index
ii <- STM_baseline$Insertion.Index.DP.5
#identify second maxima
h <- hist(ii, breaks=200,plot=FALSE)
maxindex <- which.max(h$density[10:length(h$density)])
maxval <- h$mids[maxindex+3]

# print pdf of loess curve and later on, histogram
#pdf(paste(input, "QC_and_changepoint_plots", "pdf", sep = "."))

#find inter-mode minima with loess
nG <- length(STM_baseline$Total.no..IS)
r <- floor(maxval *2000)
I = ii < r / 2000
h1 = hist(ii[I],breaks=(0:r/2000))
lo <- loess(h1$density ~ c(1:length(h1$density))) #loess smothing over density
plot(h1$density, main="Density")
lines(predict(lo),col='red',lwd=2)
m = h1$mids[which.min(predict(lo))]
I1 = ((ii < m)&(ii >= 0))


h = hist(ii, breaks="FD",plot=FALSE) 
I2 = ((ii >= m)&(ii < h$mids[max(which(h$counts>5))]))
f1 = (sum(I1) + sum(ii == 0))/nG
f2 = (sum(I2))/nG

d1 = fitdistr(ii[I1], "exponential")
d2 = fitdistr(ii[I2], "gamma") #fit curves


#plots - GammaFit
#hist(ii,breaks="FD", xlim=c(0,max(ii)), freq=FALSE,xlab="Insertion index", main="Gamma fits")
hist(ii,breaks="FD", xlim=c(0,max(ii)), freq=FALSE,xlab="Insertion index DP>5", main="Gamma fits")
lines(0:200/500, f1*dgamma(0:200/500, 1, d1$estimate[1]), col="green") # was [2]
lines(0:200/500, f2*dgamma(0:200/500, d2$estimate[1], d2$estimate[2]), col="red")


#save console figure


# print changepoint

#calculate log-odds ratios to choose thresholds
lower <- max(which(log((pgamma(1:1000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:1000/10000, 1,d1$e[1], lower.tail=FALSE)))/(pgamma(1:1000/10000, 1,d1$e[1], lower.tail=FALSE)*(1-pgamma(1:1000/10000, d2$e[1],d2$e[2]))) , base=2) < -2))
upper <- min(which(log((pgamma(1:1000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:1000/10000, 1,d1$e[1], lower.tail=FALSE)))/(pgamma(1:1000/10000, 1,d1$e[1], lower.tail=FALSE)*(1-pgamma(1:1000/10000, d2$e[1],d2$e[2]))) , base=2) > 2))

essen <- lower/10000
ambig <- upper/10000

lines(c(lower/10000, lower/10000), c(0,20), col="yellow")
lines(c(upper/10000, upper/10000), c(0,20), col="yellow")

mtext(paste(essen, ":", "Essential changepoint"), side=3, adj=1, padj=2)
mtext(paste(ambig, ":", "Ambiguous changepoint"), side=3, adj=1, padj=3.75)

###

#save console figure or create a output file.

dev.off()


write.csv(STM_baseline, file=paste("BCD_bothBatch_TraDIS_summary.csv"))
