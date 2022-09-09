rm(list = ls())
gc()
library('ggplot2')
library('ggsignif')
library('reshape2')
library('ggpubr')
library('ggsci')
library("zoo")

# from CAGE dataset
load('fig5&S11/promoter_width.Rdata')

p_A <- ggplot(promoter.width, aes(x = width, fill = shape)) + 
  geom_histogram(boundary = 15) +
  geom_vline(aes(xintercept = 15), linetype = "dashed", size = 0.5, colour = "#545357") +
  scale_fill_manual(values = c("#00A28A", '#F9AE42')) + 
  labs(y = 'Number', x = 'Promoter width(bp)') + 
  xlim(0,150) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_A
ggsave('fig5&S11/fig5A.pdf', width = 2.8, height = 1.8, useDingbats = FALSE)

promoter.width <- promoter.width[!duplicated(promoter.width$genename),]
genename <- read.table("genename.txt", sep = '\t', header = FALSE)
result.all <- read.table("../code_Inference/Hierarchical_model/results_MEF.csv", header = FALSE, sep = ',')
result.all <- cbind(genename, result.all)
colnames(result.all) <- c('genename', 'a', 'b', 'k', 'H', 'mean', 'noise', 'fval', 'goodFit')
result.all$Hsign[result.all$H > 0] <- 'negative'
result.all$Hsign[result.all$H < 0] <- 'positive'
result.all$Hsign[result.all$H == 0] <- 'non-feedback'
result.all$rcv <- result.all$noise - (1/log2(result.all$mean))

analysis.mat = merge(result.all, promoter.width, by.x = "genename", by.y = "genename", all.x = TRUE) 
analysis.mat <- analysis.mat[!is.na(analysis.mat$width),]
analysis.mat = analysis.mat[analysis.mat$a > 0 & analysis.mat$goodFit == 1,]

# mean
colourbar <- c('#888584', "#31479B")
fillbar <- c('#E3E2E0', "#98BDD9")

p_SA <- ggplot(analysis.mat, aes(x = factor(Hsign, levels = c("positive", "non-feedback","negative")),
                                 fill = factor(shape, levels = c('sharp', 'board')),
                                 colour = factor(shape, levels = c('sharp', 'board')),
                                 y = log10(mean))) +
  geom_boxplot(outlier.size = 0, outlier.colour = 'black', size = 0.2, notch = TRUE) +
  labs(y = expression(log[10](mean))) +
  theme_bw() +
  ylim(0,2) +
  scale_colour_manual(values = colourbar) +
  scale_fill_manual(values = fillbar) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 8),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_SA

p_SB <- ggplot(analysis.mat, aes(x = factor(shape, levels = c('sharp', 'board')),
                                  fill = factor(shape, levels = c('sharp', 'board')),
                                  colour = factor(shape, levels = c('sharp', 'board')),
                                  y = log10(mean))) +
  geom_boxplot(outlier.size = 0, outlier.colour = 'black', size = 0.2, notch = TRUE)+
  geom_signif(comparisons=list(c('sharp','board')),map_signif_level=F,test=wilcox.test,textsize = 2,color="black") +
  labs(y = expression(log[10](mean))) +
  theme_bw() +
  ylim(0,2) +
  scale_colour_manual(values = colourbar) +
  scale_fill_manual(values = fillbar) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 8),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_SB

# rcv
colourbar <- c('#888584', "#159E85")
fillbar <- c('#E3E2E0', "#CAE6DE")

p_SC <- ggplot(analysis.mat, aes(x = factor(Hsign, levels = c("positive", "non-feedback","negative")),
                                 fill = factor(shape, levels = c('sharp', 'board')),
                                 colour = factor(shape, levels = c('sharp', 'board')),
                                 y = rcv)) +
  geom_boxplot(outlier.size = 0, outlier.colour = 'black', size = 0.2, notch = TRUE)+
  labs(y = expression(rcv)) +
  theme_bw() +
  scale_colour_manual(values = colourbar) +
  scale_fill_manual(values = fillbar) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 8),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_SC

p_SD <- ggplot(analysis.mat, aes(x = factor(shape, levels = c('sharp', 'board')),
                                  fill = factor(shape, levels = c('sharp', 'board')),
                                  colour = factor(shape, levels = c('sharp', 'board')),
                                  y = rcv)) +
  geom_boxplot(outlier.size = 0, outlier.colour = 'black', size = 0.2, notch = TRUE)+
  geom_signif(comparisons=list(c('sharp','board')),map_signif_level=F,test=wilcox.test,textsize = 2,color="black") +
  labs(y = expression(rcv)) +
  theme_bw() +
  scale_colour_manual(values = colourbar) +
  scale_fill_manual(values = fillbar) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 8),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_SD

pic1 <- ggarrange(p_SA, p_SB, p_SC, p_SD, ncol = 2, nrow = 2, widths=c(3,1.2), heights=c(2), align = "hv")
pic1

ggsave('fig5&S11/figS11A-D.pdf', width = 4, height = 3, useDingbats = FALSE)

# bf
colourbar <- c("#00A28A", '#F9AE42')
fillbar <- c("#A0D1CA", '#FBD9AD')

p_C1 <- ggplot(analysis.mat, aes(x = factor(Hsign, levels = c("negative", "non-feedback", "positive")),
                         fill = factor(shape, levels = c('board', 'sharp')),
                         colour = factor(shape, levels = c('board', 'sharp')),
                         y = log10(a))) +
  geom_boxplot(outlier.size = 0, outlier.colour = 'black', size = 0.3, notch = TRUE)+
  labs(y = expression(log[10](bf))) +
  ylim(-0.7,1.5) +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = colourbar) +
  scale_fill_manual(values = fillbar) +
  theme(legend.position = 'none',
        axis.line.x = element_line(size = 0.15),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 8),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.border = element_blank(),
        panel.grid = element_blank())
p_C1

colourbar <- c('#888584', "#EA5345")
fillbar <- c('#E3E2E0', "#FDEBE7")

p_SE <- ggplot(analysis.mat, aes(x = factor(shape, levels = c('sharp', 'board')),
                                 fill = factor(shape, levels = c('sharp', 'board')),
                                 colour = factor(shape, levels = c('sharp', 'board')),
                                 y = log10(a))) +
  geom_boxplot(outlier.size = 0, outlier.colour = 'black', size = 0.3, notch = TRUE)+
  geom_signif(comparisons=list(c('sharp','board')),map_signif_level=F,test=wilcox.test,textsize = 2,color="black") +
  labs(y = expression(log[10](bf))) +
  theme_bw() +
  scale_colour_manual(values = colourbar) +
  scale_fill_manual(values = fillbar) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 8),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_SE

# bs
colourbar <- c("#00A28A", '#F9AE42')
fillbar <- c("#A0D1CA", '#FBD9AD')
p_C2 <- ggplot(analysis.mat, aes(x = factor(Hsign, levels = c("negative", "non-feedback", "positive")),
                                 fill = factor(shape, levels = c('board', 'sharp')),
                                 colour = factor(shape, levels = c('board', 'sharp')),
                                 y = log10(b))) +
  geom_boxplot(outlier.size = 0, outlier.colour = 'black', size = 0.3, notch = TRUE)+
  labs(y = expression(log[10](bs))) +
  ylim(0,1.5) +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = colourbar) +
  scale_fill_manual(values = fillbar) +
  theme(legend.position = 'none',
        axis.line.x = element_line(size = 0.15),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 8),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.border = element_blank(),
        panel.grid = element_blank())
p_C2

colourbar <- c('#888584', "#F9AE42")
fillbar <- c('#E3E2E0', "#FBD9AD")

p_SF <- ggplot(analysis.mat, aes(x = factor(shape, levels = c('sharp', 'board')),
                                 fill = factor(shape, levels = c('sharp', 'board')),
                                 colour = factor(shape, levels = c('sharp', 'board')),
                                 y = log10(b))) +
  geom_boxplot(outlier.size = 0, outlier.colour = 'black', size = 0.3, notch = TRUE)+
  geom_signif(comparisons=list(c('sharp','board')),map_signif_level=F,test=wilcox.test,textsize = 2,color="black") +
  labs(y = expression(log[10](bs))) +
  theme_bw() +
  scale_colour_manual(values = colourbar) +
  scale_fill_manual(values = fillbar) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 8),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_SF

pic2 <- ggarrange(p_C1, p_C2, ncol = 2, nrow = 1, widths=c(2), heights=c(2), align = "v")
pic2

ggsave('fig5&S11/fig5C.pdf', width = 3.5, height = 1.8, useDingbats = FALSE)

pic3 <- ggarrange(p_SE, p_SF, ncol = 1, nrow = 2, widths=c(1.2), heights=c(2), align = "hv")
pic3
ggsave('fig5&S11/figS11E-F.pdf', width = 1.2, height = 3, useDingbats = FALSE)

# rolling functional smoothing
analy.mat.po <- analysis.mat[analysis.mat$Hsign == "positive", c("a", "b", "mean", "rcv", "width")]
analy.mat.po <- analy.mat.po[order(analy.mat.po$width),]
roll.po <-  data.frame(rollmean(analy.mat.po, 200))
roll.po$Hsign <- "positive" 

analy.mat.non <- analysis.mat[analysis.mat$Hsign == "non-feedback", c("a", "b", "mean", "rcv", "width")]
analy.mat.non <- analy.mat.non[order(analy.mat.non$width),]
roll.non <- data.frame(rollmean(analy.mat.non, 200))
roll.non$Hsign <- "non-feedback"

analy.mat.ne <- analysis.mat[analysis.mat$Hsign == "negative", c("a", "b", "mean", "rcv", "width")]
analy.mat.ne <- analy.mat.ne[order(analy.mat.ne$width),]
roll.ne <- data.frame(rollmean(analy.mat.ne, 200))
roll.ne$Hsign <- "negative"

analy.roll.mat <- data.frame(rbind(roll.po, roll.non, roll.ne))
analy.roll.mat$a <- log10(analy.roll.mat$a)
analy.roll.mat$b <- log10(analy.roll.mat$b)
analy.roll.mat$mean <- log10(analy.roll.mat$mean)

color.feedback <- c("#3AA438", "grey50", "#E71419")
p_B <- ggplot(analy.roll.mat, aes(x = width, y = rcv)) +
  geom_line(aes(group = Hsign, colour = Hsign), size = 0.5) +
  labs(y = expression(rcv), x = 'Promoter width (bp)') + 
  theme_bw() +
  scale_colour_manual(values = color.feedback)+
  scale_fill_manual(values = color.feedback)+
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid.minor.y = element_blank())
p_B
ggsave('fig5&S11/fig5B.pdf', width = 2.8, height = 1.8, useDingbats = FALSE)
