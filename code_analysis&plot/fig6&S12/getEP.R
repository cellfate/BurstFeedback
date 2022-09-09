EP <- read.table("fig6&S10/enhancer_gene.csv", header = T, sep = ",")
colnames(EP) <- c("genename", "intensity")
genename <- read.table("genename.txt",sep='\t',header = FALSE)
genename <- merge(genename,EP,by.x = "V1",by.y = "genename",all.x = TRUE)
genename$EP <- as.numeric(!is.na(genename$intensity))
colnames(genename) <- c("genename", "intensity","EP")
enhancer_feature <- genename[,c("genename","EP","intensity")]
save(enhancer_feature,file='fig6&S10/enhancer_feature.Rdata')