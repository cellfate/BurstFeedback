library("ggpubr")
library("ggplot2")
library("reshape2")

####
data.all <- read.table("../code_Inference/Hierarchical_model/data/MEF_QC_all.csv", sep = ',')
data.all[data.all == -1] <- NA
genename <- read.table("genename.txt", sep = '\t', header = FALSE)
row.names(data.all) <- c(as.character(genename$V1))

# read result file
result.all <- read.table("../code_Inference/Hierarchical_model/results_MEF.csv", header = FALSE, sep = ',')
result.all <- cbind(genename, result.all)
colnames(result.all) <- c('genename', 'a', 'b', 'k', 'H', 'mean', 'noise', 'fval', 'goodFit')
result.all$mean_data <- apply(data.all,1,function(x) {mean(x, na.rm = T)})

# feedback type
result.all$Hsign[result.all$H > 0] <- 'negative'
result.all$Hsign[result.all$H < 0] <- 'positive'
result.all$Hsign[result.all$H == 0] <- 'non-feedback'

#### Goodness fit ####
p_A <- ggplot(result.all, aes(x = log10(mean_data), fill = as.factor(goodFit))) +
  geom_histogram(bins = 30) + 
  scale_fill_manual(values = c("#E26463", "#71A1C6")) +
  labs(x = expression(log[10](mean)), y = "Numbers of gene") +
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

p_A
#### error of mean and noise ####
data.mean <- apply(data.all, 1, function(x){
  x <- x[!is.na(x)]
  x <- sort(x)
  x <- x[1:floor(0.95*length(x))]
  mean(x)
})

data.var <- apply(data.all, 1, function(x){
  x <- x[!is.na(x)]
  x <- sort(x)
  x <- x[1:floor(0.95*length(x))]
  var(x)
})

data.noise <- data.var/(data.mean^2)
mat.error <- data.frame(mean = data.mean / result.all$mean,
                        noise = data.noise / result.all$noise,
                        goodFit = result.all$goodFit)
mat.error.melt <- melt(mat.error, id.vars = c("goodFit"))
mat.error.melt <- mat.error.melt[mat.error.melt$goodFit == 1,]

p_B <- ggplot(mat.error.melt, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = variable), width = 0.2, outlier.size = 0, size = 0.15, fatten = 1) +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed", size = 0.5) +
  labs(y = 'Data/Model') +
  ylim(0,4) +
  theme_bw()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

p_B
pic1 <- ggarrange(p_A, p_B, ncol = 2, nrow = 1, widths = c(2), heights = c(2), align = "hv")
pic1

ggsave('fig2&S7/figS7.pdf', width = 4, height = 2, useDingbats = FALSE)
