library('ggplot2')
library('ggsignif')
library("zoo")
library('reshape2')

load('fig6&S12/enhancer_feature.Rdata')
genename <- read.table("genename.txt", sep = '\t', header = FALSE)
result.all <- read.table("../code_Inference/Hierarchical_model/results_MEF.csv", header = FALSE, sep = ',')
result.all <- cbind(genename, result.all)
colnames(result.all) <- c('genename', 'a', 'b', 'k', 'H', 'mean', 'noise', 'fval', 'goodFit')
result.all$Hsign[result.all$H > 0] <- 'negative'
result.all$Hsign[result.all$H < 0] <- 'positive'
result.all$Hsign[result.all$H == 0] <- 'non-feedback'
result.all$rcv <- result.all$noise - (1/log2(result.all$mean))

analysis.mat = merge(enhancer_feature, result.all, by.x = "genename", by.y = "genename", all.x = TRUE)
analysis.mat = analysis.mat[analysis.mat$a > 0 & analysis.mat$goodFit == 1 & analysis.mat$EP == 1,]
analysis.mat$EP.rank <- rank(-analysis.mat$intensity)

## E-P interaction
analysis.mat2 <- analysis.mat[analysis.mat$EP.rank >= 1613 | analysis.mat$EP.rank <= 100,]
analysis.mat2[analysis.mat2$EP.rank <= 100, c("EP.level")] <- "high" 
analysis.mat2[analysis.mat2$EP.rank >= 1613, c("EP.level")] <- "low"

## Effects of EP intensity

# rolling functional smoothing
analy.mat.po <- analysis.mat[analysis.mat$Hsign == "positive" & !is.na(analysis.mat$intensity), c("a","b","mean","rcv","intensity")]
analy.mat.po <- analy.mat.po[order(analy.mat.po$intensity),]
roll.po <-  data.frame(rollmean(analy.mat.po,200))
roll.po$Hsign <- "positive" 

analy.mat.non <- analysis.mat[analysis.mat$Hsign == "non-feedback" & !is.na(analysis.mat$intensity), c("a","b","mean","rcv","intensity")]
analy.mat.non <- analy.mat.non[order(analy.mat.non$intensity),]
roll.non <- data.frame(rollmean(analy.mat.non,200))
roll.non$Hsign <- "non-feedback"

analy.mat.ne <- analysis.mat[analysis.mat$Hsign == "negative" & !is.na(analysis.mat$intensity), c("a","b","mean","rcv","intensity")]
analy.mat.ne <- analy.mat.ne[order(analy.mat.ne$intensity),]
roll.ne <- data.frame(rollmean(analy.mat.ne,200))
roll.ne$Hsign <- "negative"

analy.roll.mat <- data.frame(rbind(roll.po,roll.non,roll.ne))

color.feedback <- c("#3AA438", "grey50", "#E71419")
p_B <- ggplot(analy.roll.mat,aes(x = intensity, y = rcv)) +
  geom_smooth(aes(group = Hsign, colour = Hsign, fill = Hsign), size = 0.5, alpha = 0.2) +
  labs(y = expression(rcv), x = 'E-P intensity') + 
  theme_bw() +
  scale_colour_manual(values = color.feedback) +
  scale_fill_manual(values = color.feedback) +
  theme(legend.position = 'none',
        axis.title.x = element_text(colour = 'black', size = 8),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 6),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_B

ggsave('fig6&S12/fig6B.pdf', width = 1.8, height = 1.5, useDingbats = FALSE)

p_SD <- ggplot(analy.roll.mat,aes(x = intensity, y = rcv)) +
  geom_smooth(size = 0.5, alpha = 0.2, colour = '#F9AE42', fill = '#FBD9AD') +
  labs(y = expression(rcv), x = 'E-P intensity') + 
  theme_bw() +
  theme(legend.position = 'none',
        axis.title.x = element_text(colour = 'black', size = 8),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 6),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_SD

p_SA <- ggplot(analy.roll.mat,aes(x = intensity, y = mean)) +
  geom_smooth(aes(group = Hsign, colour = Hsign, fill = Hsign), size = 0.5, alpha = 0.2) +
  labs(y = expression(mean), x = 'E-P intensity') + 
  theme_bw() +
  scale_colour_manual(values=color.feedback)+
  scale_fill_manual(values=color.feedback)+
  theme(legend.position = 'none',
        axis.title.x = element_text(colour = 'black', size = 8),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 6),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_SA

p_SE <- ggplot(analy.roll.mat,aes(x = intensity, y = mean)) +
  geom_smooth(size = 0.5, alpha = 0.2, colour = '#F9AE42', fill = '#FBD9AD') +
  labs(y = expression(mean), x = 'E-P intensity') + 
  theme_bw() +
  scale_colour_manual(values=color.feedback)+
  scale_fill_manual(values=color.feedback)+
  theme(legend.position = 'none',
        axis.title.x = element_text(colour = 'black', size = 8),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 6),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_SE

p_SB <- ggplot(analy.roll.mat,aes(x = intensity, y = a)) +
  geom_smooth(aes(group = Hsign, colour = Hsign, fill = Hsign), size = 0.5, alpha = 0.2) +
  labs(y = expression(bf), x = 'E-P intensity') + 
  theme_bw() +
  scale_colour_manual(values=color.feedback)+
  scale_fill_manual(values=color.feedback)+
  theme(legend.position = 'none',
        axis.title.x = element_text(colour = 'black', size = 8),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 6),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_SB

p_SF <- ggplot(analy.roll.mat,aes(x = intensity, y = a)) +
  geom_smooth(size = 0.5, alpha = 0.2, colour = '#F9AE42', fill = '#FBD9AD') +
  labs(y = expression(bf), x = 'E-P intensity') + 
  theme_bw() +
  scale_colour_manual(values=color.feedback)+
  scale_fill_manual(values=color.feedback)+
  theme(legend.position = 'none',
        axis.title.x = element_text(colour = 'black', size = 8),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 6),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_SF

p_SC <- ggplot(analy.roll.mat,aes(x = intensity, y = b)) +
  geom_smooth(aes(group = Hsign, colour = Hsign, fill = Hsign), size = 0.5, alpha = 0.2) +
  labs(y = expression(bs), x = 'E-P intensity') + 
  theme_bw() +
  scale_colour_manual(values=color.feedback)+
  scale_fill_manual(values=color.feedback)+
  theme(legend.position = 'none',
        axis.title.x = element_text(colour = 'black', size = 8),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 6),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_SC

p_SG <- ggplot(analy.roll.mat,aes(x = intensity, y = b)) +
  geom_smooth(size = 0.5, alpha = 0.2, colour = '#F9AE42', fill = '#FBD9AD') +
  labs(y = expression(bs), x = 'E-P intensity') + 
  theme_bw() +
  scale_colour_manual(values=color.feedback)+
  scale_fill_manual(values=color.feedback)+
  theme(legend.position = 'none',
        axis.title.x = element_text(colour = 'black', size = 8),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 6),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_SG


pic1 <- ggarrange(p_SA, p_SB, p_SC, ncol = 3, nrow = 1, widths=c(1.8), heights=c(1.5), align = "v")
pic1

ggsave('fig6&S12/figS12A-C.pdf', width = 5.4, height = 1.5, useDingbats = FALSE)

pic2 <- ggarrange(p_SD, p_SE, p_SF, p_SG, ncol = 4, nrow = 1, widths = c(1.5), heights = c(1.5), align = "v")
pic2
ggsave('fig6&S12/figS12D-G.pdf', width = 6.5, height = 1.5, useDingbats = FALSE)

analy.roll.mat2 <- analy.roll.mat[analy.roll.mat$intensity < 35 & analy.roll.mat$Hsign == "positive",]
analy.roll.mat2$a <- analy.roll.mat2$a - mean(analy.roll.mat2$a)
analy.roll.mat2$b <- analy.roll.mat2$b - mean(analy.roll.mat2$b)
analy.roll.mat2$intensity <- analy.roll.mat2$intensity - mean(analy.roll.mat2$intensity)
analy.roll.melt.mat2 <- melt(analy.roll.mat2, id = c("mean","rcv","intensity","Hsign"))

p_C <- ggplot(analy.roll.melt.mat2)+
  geom_line(aes(x = intensity, y = value, group = variable, colour = variable), size = 0.4) +
  labs(y = "Burst kinetic (normalized)", x = 'E-P intensity (normalized)') + 
  theme_bw() +
  xlim(-5,10) +
  scale_colour_manual(values=c("#D26A50",'#5893CF'))+
  theme(legend.position = c(0.9,0.2),
        legend.title = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 8),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 6),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid.minor = element_blank())

p_C

ggsave('fig6&S12/fig6C.pdf', width = 1.85, height = 1.5, useDingbats = FALSE)