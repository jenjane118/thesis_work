## make insertion plot diagram to compare insertions at TA sites in mbovis and mtb

# import insertion data (reads per insertion site)

mtb_temp<-read.table("dejesus_combined_ttr.wig", header=F, sep="\t",
                     comment.char="#", stringsAsFactors = F)
mtb_df<- mtb_temp[,1:15]
read_sums<- rowSums(mtb_df[,2:15])
read_means<- rowMeans(mtb_df[,2:15])
mtb_df<- cbind(mtb_df, read_sums, read_means)
head(mtb_df)

bovis_temp<-read.table("bovis_combined_TTR.wig", header=F, sep="\t",
                       comment.char="#", stringsAsFactors = F)
head(bovis_temp)
bovis_df <- bovis_temp[,1:4]
bovis_df <- cbind(bovis_df, rowSums(bovis_df[,2:4]), rowMeans(bovis_df[,2:4]))
head(bovis_df)
colnames(bovis_df)<-c("pos", "B", "C", "G", "sum", "mean")

max(bovis_df$sum)
#90726.9
max(bovis_df$mean)
#30242.3

max(mtb_df$read_sums)
#946019.6
max(mtb_df$read_means)
#67572.83

## making insertion plots with ggplot

install.packages("ggplot2")
library(ggplot2)
library(dplyr)
#load up wig file insertions for bovis
b19<-read.delim("~/Data/bovis_hiseq/tpp/bovis_hiseq_tpp_19.wig", header=T, comment.char = "#", sep=" ")
b20<-read.delim("~/Data/bovis_hiseq/tpp/bovis_hiseq_tpp_20.wig", header=T, comment.char = "#", sep=" ")
b21<-read.delim("~/Data/bovis_hiseq/tpp/bovis_hiseq_tpp_21.wig", header=T, comment.char = "#", sep=" ")

b21[310:320,]
insertions<-b19[,2]
insertions_df<-cbind(insertions, b20[,2], b21[,2])
head(insertions_df)

totals<-rowSums(insertions_df[, c(1, 3)])
head(totals)

# sum insertions into one df
positions<-b19[,1]
mb_ins_df<-data.frame(positions, totals)
#mb_ins_df<-cbind(positions, totals)
colnames(mb_ins_df)<-c("positions", "insertions")
head(mb_ins_df)
nrow(mb_ins_df)


max(mb_ins_df$insertions)
mb_ins_df[which.max(mb_ins_df[,2])]


plot(mb_ins_df$positions, mb_ins_df$insertions, type='l', main="insertion plot mbovis",
     xlab="TA position", ylab="number of insertions")

#ggplot(data = mb_ins_df) +
#geom_bar(mapping = aes(x = positions, y = insertions), stat = "identity")

## this works better for getting pretty axes
ggplot(data=mb_ins_df, aes(x = positions, y = insertions)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  theme_bw() +
  theme(axis.text = element_text(size = 12)) +
  scale_x_continuous(name = positions, 
                     breaks = seq(0, max(positions), 500000), 
                     labels = waiver())

## make bigger plot for export
## to make bigger plot for export
#png(file="mbo_insertions.png", width=2048, height=1536)
b<-ggplot(data=mb_ins_df, aes(x = positions, y = insertions)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title=element_text(size=16,face="bold")) +
  ylim(0, 2500) +
  scale_x_continuous(name = "position", 
                     breaks = seq(0, max(positions), 500000), 
                     labels = waiver()) +
  ggtitle("Insertion plot M.bovis (hiseq)") + 
  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20, hjust=0.5))
#dev.off()
ggsave("bovis_insertions_line.png", b, width = 15, height = 10)

# do the same with tb
# load .wig files
b22<-read.delim("~/Data/mtb_hiseq/TPP/mtb22_hiseq_tpp.wig", header=T, comment.char = "#", sep=" ")
b23<-read.delim("~/Data/mtb_hiseq/TPP/mtb23_hiseq_tpp.wig", header=T, comment.char = "#", sep=" ")

ins_df<-data.frame(b22[,2], b23[,2])
totals<-rowSums(ins_df[, c(1,2)])
positions<-b22[,1]
mtb_ins_df<-data.frame(positions, totals)
colnames(mtb_ins_df)<-c("positions", "insertions")
head(mtb_ins_df)

ggplot(data=mtb_ins_df, aes(x = positions, y = insertions)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  theme_bw() +
  theme(axis.text = element_text(size = 12))

## to make bigger plot for export
#png(file="~/Data/mtb_insertions.png", width=2048, height=1536)
mtb<-ggplot(data=mtb_ins_df, aes(x = positions, y = insertions)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title=element_text(size=16,face="bold")) +
  ylim(0, 2500) +
  scale_x_continuous(name = "position", 
                     breaks = seq(0, max(positions), 500000), 
                     labels = waiver()) +
  ggtitle("Insertion plot M.tb (hiseq)") + 
  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20, hjust=0.5))
#dev.off()
ggsave("~/Data/mtb_insertions.png", mtb, width = 15, height = 10)


# make a circular barplot with ggplot

# Set a number of 'empty bar' to add at the end of genome (adding 1000 blank entries)
# maybe should be more like 100000?
empty_bar <- 1000
to_add <- data.frame( matrix(0, empty_bar, ncol(mtb_ins_df)) )
head(to_add)
colnames(to_add) <- colnames(mtb_ins_df)

data <- rbind(mtb_ins_df, to_add)
nrow(data)
data[75604,]
#data <- data %>% arrange(group)
#data$id <- seq(1, nrow(data))


#set base data--line that goes around base
# prepare a data frame for base lines
base_data <- data %>% 
  summarize(start=min(positions), end=max(positions)) %>% 
  mutate(title=mean(c(start, end)))
# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

library(tidyverse)
#png(file="mtb_insertions_circular.png", width=2048, height=2048)
p<- ggplot(data=mtb_insertions, aes(x = positions, y = insertions)) +
  geom_bar(stat="identity", width = 1, color = "black") +
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-10000,10000) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title=element_text(size=16,face="bold")
    ) +
  scale_x_continuous(name = "position", 
                   breaks = seq(0, max(positions), 500000), 
                   labels = waiver()
                   ) +
  # this makes coordinate polar instead of cartesian (circular?)
  coord_polar(start = 0) +
  ggtitle("Insertion plot M.tb (hiseq)")
  #theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20, hjust=0.5)) +
  
# Add base line information
  #geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  #geom_text(data=base_data, aes(x = title, y = -18, label="H37Rv"), hjust=1, colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)


#p
#dev.off()

ggsave("~/Data/mtb_insertions_circular.png", p, width = 15, height = 10)  #dimensions in inches, not pixels


## make circular insertion plot for bovis

bov<- ggplot(data=mb_ins_df, aes(x = positions, y = insertions)) +
  geom_bar(stat="identity", width = 1, color = "black") +
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-1000,1000) +
  #theme_bw() +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")
  ) +
  scale_x_continuous(name = "position", 
                     breaks = seq(0, max(positions), 500000), 
                     labels = waiver()
  ) +
  # this makes coordinate polar instead of cartesian (circular?)
  coord_polar(start = 0) +
  ggtitle("Insertion plot M.bovis (hiseq)")

ggsave("~/Data/bovis_insertions_circular.png", bov, width = 15, height = 10)  #dimensions in inches, not pixels




# Upload library

install.packages("circlize")
library(circlize)

circos.clear()
circos.par("track.height" = 0.4)

# Create data
#data = data.frame(
#  factor = sample(letters[1:8], 1000, replace = TRUE),
#  x = rnorm(1000),
#  y = runif(1000)
#)

data = mb_ins_df[1:1000,]


#basic circos graphic parameters
circos.par(gap.degree =.0001, cell.padding=c(0, 0, 0, 0), track.margin=c(0,0), start.degree = 90, unit.circle.segments = 10000)

chrom<-c("mbovis")

# Step1: Initialise the chart giving factor and x-axis.
circos.initialize( factors=data$positions, xlim=c(1,1000))

#Error: Maybe your `gap.degree` is too large so that there is no space to allocate sectors.
## changed to 0.0001

# Step 2: Build the regions.
circos.trackPlotRegion(factors = data$positions, y = data$insertions)


# Step 3: Add points
circos.trackLines(factors=data$positions, data$x[order(data$x)], data$y[order(data$x)], col = rgb(0.1,0.5,0.8,0.3), lwd=2, type="h")

circos.trackLines(sectors, x, y)
