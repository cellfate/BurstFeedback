BiocManager::install("CAGEr")
library(CAGEr)
nrow(FANTOM5mouseSamples)
head(FANTOM5mouseSamples)
x <- read.table("fig5&S9/mart_export.txt",sep=',',header = T)
g <- GRanges( paste0("chr", x$Chromosome.Name)
              , IRanges(x$Gene.Start..bp., x$Gene.End..bp.)
              , ifelse( x$Strand + 1, "+", "-"))
g$gene_name <- Rle(x$Associated.Gene.Name)
g$transcript_type <- Rle(x$Gene.Biotype)
g$type <- "gene"
g$type <- Rle(g$type)
g <- sort(unique(g))
names(g) <- g$gene_name

Samples <- FANTOM5mouseSamples[grep("primary cell",FANTOM5mouseSamples[,"type"]),]
MEF_Samples <- Samples['21',]
CAGEset.MEF <- importPublicData(source = "FANTOM5", dataset = "mouse", sample = MEF_Samples[1,"sample"])
librarySizes(CAGEset.MEF)
plotReverseCumulatives(CAGEset.MEF, fitInRange = c(5, 1000), onePlot = TRUE)
normalizeTagCount(CAGEset.MEF, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.16, T = 10^6)
clusterCTSS( object = CAGEset.MEF, 
             threshold = 1, thresholdIsTpm = TRUE,
             nrPassThreshold = 1,
             method = "distclu", 
             maxDist = 20, 
             removeSingletons = TRUE, 
             keepSingletonsAbove = 5)
cumulativeCTSSdistribution(CAGEset.MEF, clusters = "tagClusters")
quantilePositions(CAGEset.MEF, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
results <- tagClustersGR(CAGEset.MEF, "Mouse_Embryonic_fibroblasts__donor1", 
              returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)
plotInterquantileWidth(CAGEset.MEF, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)
exportToBed(object = CAGEset.MEF, what = "tagClusters", qLow = 0.1, qUp = 0.9, oneFile = FALSE)
overlap.index <- findOverlaps(results,g)
genename <- names(g)
genename <- genename[overlap.index@to]

promoter.width <- results@elementMetadata@listData$interquantile_width
promoter.width <- promoter.width[overlap.index@from]
promoter.width <- data.frame(genename = genename, width = promoter.width)
median(promoter.width$width)
promoter.width$shape[promoter.width$width > 15] <- 'board'
promoter.width$shape[promoter.width$width <= 15] <- 'sharp'

save(promoter.width, file='fig5&S9/promoter_width.Rdata')
