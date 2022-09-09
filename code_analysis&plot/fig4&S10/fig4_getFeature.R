library(stringr)
promofile = list.files('fig4&S8/EPD data',pattern="mouse.*.bed")
promoter_feature <- read.table("genename.txt",sep='\t',header = FALSE)

## promoter feature from EPD database
for (i in 1:length(promofile)){
  data <- read.table(paste0("fig4&S8/EPD data/",promofile[i]),sep ='\t',header = FALSE)
  data$V4 = str_replace(data$V4, "_1","")
  data = data[c('V4',"V5")]
  promoter_feature = merge(promoter_feature,data,by.x = "V1",by.y = "V4",all.x = TRUE) 
}
colnames(promoter_feature) <- c('genename','TATA','Inr','CCAAT','GC')
save(promoter_feature,file='fig4&S8/promoter_feature.Rdata')

