# In this script I make alignment plot of a multiplicon using the genoPlotR library
## Load in library
library(genoPlotR)

args = commandArgs(trailingOnly=TRUE)
multiplicon = args[1]
species1 = args[2]
chromosome1 = args[3]
species2 = args[4]
chromosome2 = args[5]

## Make lists of all genes in multiplicon
df1 <- read.table(sprintf("region_%s_%s.txt", species1, chromosome1), header=FALSE)
df2 <- read.table(sprintf("region_%s_%s.txt", species2, chromosome2), header=FALSE)
df2_reversed <- read.table(sprintf("region_%s_%s_reversed.txt", species2, chromosome2), header=FALSE)
colnames(df1) <- c("name", "start", "end", "strand", "col", "fill")
colnames(df2) <- c("name", "start", "end", "strand", "col", "fill")
colnames(df2_reversed) <- c("name", "start", "end", "strand", "col", "fill")

## Make coordinates for comparisons
fh <- read.table(sprintf("multiplicon_%s.txt", multiplicon), sep = "\t", header = FALSE)

## Make unique list of all genes
genes <- union(df1$name, df2$name)
genes <- unique(genes)

## Create genoPlotR dna_segs
dna_seg1 <- dna_seg(df1)
dna_seg2 <- dna_seg(df2)
dna_segs <- list(dna_seg1, dna_seg2)

dna_seg2_reversed <- dna_seg(df2_reversed)
dna_segs_reversed <- list(dna_seg1, dna_seg2_reversed)

## Obtain coordinates from dna_seg's
start1 <- sapply(fh$V4, function(x) dna_seg1[dna_seg1$name==x,]$start)
end1 <- sapply(fh$V4, function(x) dna_seg1[dna_seg1$name == x,]$end)
start2 <- sapply(fh$V5, function(x) dna_seg2[dna_seg2$name==x,]$start)
end2 <- sapply(fh$V5, function(x) dna_seg2[dna_seg2$name==x,]$end)

col <- ifelse(grepl("Pseudogene", fh$V4) | grepl("Pseudogene", fh$V5), ifelse(grepl("Pseudogene", fh$V4) & grepl("Pseudogene", fh$V5), "#99FF99", "#0505ff8a"), "#f04b4b38")

start2_reversed <- sapply(fh$V5, function(x) dna_seg2_reversed[dna_seg2_reversed$name==x,]$start)
end2_reversed <- sapply(fh$V5, function(x) dna_seg2_reversed[dna_seg2_reversed$name==x,]$end)

## Add annotation (if short enough multiplicon region)
# annot1 <- annotation(x1=dna_seg1$start, text = dna_seg1$name, rot= 60)
# annot2 <- annotation(x1=dna_seg2$start, text = dna_seg2$name, rot= 60)
# annots <- list(annot1, annot2)

## Create genoPlotR comparisons
comparison1 <- as.comparison(data.frame(start1=start1, end1=end1, start2=start2, end2=end2, col=col))
comparisons <- list(comparison1)

comparison1_reversed <- as.comparison(data.frame(start1=start1, end1=end1, start2=start2_reversed, end2=end2_reversed, col=col))
comparisons_reversed <- list(comparison1_reversed)

## Adding names to plot
names <- c(sprintf("%s|%s", species1, chromosome1), sprintf("%s|%s", species2, chromosome2))
names(dna_segs) <- names
names(dna_segs_reversed) <- names

#make figure
jpeg(sprintf("multiplicon_%s_%s_%s.jpeg", multiplicon, species1, species2), width=4000, height=400)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons, scale = FALSE, dna_seg_scale= FALSE, dna_seg_label_cex = 1, annotation_cex = 1)
dev.off()

jpeg(sprintf("multiplicon_%s_%s_%s_reversed.jpeg", multiplicon, species1, species2), width=4000, height=400)
plot_gene_map(dna_segs=dna_segs_reversed, comparisons=comparisons_reversed, scale = FALSE, dna_seg_scale= FALSE, dna_seg_label_cex = 1, annotation_cex = 1)
dev.off()
