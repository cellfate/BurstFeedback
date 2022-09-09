library("ggplot2")
library("pROC")
library("ggstatsplot")
library('ggsignif')
library('ggpubr')
load('fig4&S10/promoter_feature.Rdata')
promoter_feature[is.na(promoter_feature)] <- 0

genename <- read.table("genename.txt", sep = '\t', header = FALSE)
result.all <- read.table("../code_Inference/Hierarchical_model/results_MEF.csv", header = FALSE, sep = ',')
result.all <- cbind(genename, result.all)
colnames(result.all) <- c('genename', 'a', 'b', 'k', 'H', 'mean', 'noise', 'fval', 'goodFit')
result.all$Hsign[result.all$H > 0] <- 'negative'
result.all$Hsign[result.all$H < 0] <- 'positive'
result.all$Hsign[result.all$H == 0] <- 'non-feedback'

result.all$rcv <- result.all$noise - (1/log2(result.all$mean))
analysis.mat = merge(promoter_feature, result.all, by.x = "genename", by.y = "genename", all.x = TRUE)
analysis.mat = analysis.mat[analysis.mat$a > 0 & analysis.mat$goodFit == 1,]

###################################################################
# Hierarchical Linear Regression
## rcv
rcvfit <- lm(rcv ~ (TATA * Inr + GC * CCAAT) : Hsign, data = analysis.mat)
rcvfit <- summary(rcvfit)$coef
name <- strsplit(rownames(rcvfit), ':Hsign')
feature <- sapply(name, function(x) {x[1]})
feedback <- sapply(name, function(x) {x[2]})
rcv.result <- data.frame(cbind(rcvfit[2:19,c(3,4)], feature = feature[2:19], feedback = feedback[2:19]))
rcvindex <- as.numeric(as.character(rcv.result[,2])) > 5e-2
rcv.result$sign[rcvindex] <- 'non-significant' 
rcv.result$sign[!rcvindex & as.numeric(as.character(rcv.result$t.value))>0] <- 'positive effect'
rcv.result$sign[!rcvindex & as.numeric(as.character(rcv.result$t.value))<0] <- 'negative effect'
colnames(rcv.result) <- c('tvalue', 'pvalue', 'feature', 'feedback', 'significant')

maxt.rcv <- max(abs(as.numeric(as.character(rcv.result$tvalue)))) + 0.5
p_B <- ggplot(rcv.result, aes(x = factor(feature, levels = c('TATA', 'Inr', 'TATA:Inr', 'GC', 'CCAAT', 'GC:CCAAT')),
                             y = as.numeric(as.character(tvalue)))) + 
  geom_hline(aes(yintercept = 1.96), linetype = "dashed", size = 0.5, colour = "#F7C3B4") +
  geom_hline(aes(yintercept = -1.96), linetype = "dashed", size = 0.5, colour = "#B8CEBA") +
  geom_hline(aes(yintercept = 0), size = 0.5, colour = "#EBEAEA") +
  geom_point(aes(shape = factor(feedback, levels = c('positive', 'negative', 'non-feedback')),
                 colour = factor(significant, levels = c('non-significant', 'negative effect', 'positive effect'))),
             position = position_dodge(0.5), size = 1.5) +
  labs(y = 'Effect size(rcv)', x = 'Promoter architecture') + 
  ylim(-maxt.rcv, maxt.rcv) + 
  theme_bw() +
  scale_shape_manual(values = c(15,16,17)) +
  scale_colour_manual(values = c("grey70","#F8766D",'#3C977B')) +
  scale_fill_manual(values = c("grey70","#F8766D",'#3C977B')) +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank())
p_B

## bf
afit <- lm(log10(a) ~ (TATA * Inr + GC * CCAAT) : Hsign, data = analysis.mat)
afit <- summary(afit)$coef
name <- strsplit(rownames(afit), ':Hsign')
feature <- sapply(name, function(x) {x[1]})
feedback <- sapply(name, function(x) {x[2]})
a.result <- data.frame(cbind(afit[2:19,c(3,4)], feature = feature[2:19], feedback = feedback[2:19]))
aindex <- as.numeric(as.character(a.result[,2])) > 5e-2
a.result$sign[aindex] <- 'non-significant' 
a.result$sign[!aindex & as.numeric(as.character(a.result$t.value))>0] <- 'positive effect'
a.result$sign[!aindex & as.numeric(as.character(a.result$t.value))<0] <- 'negative effect'
colnames(a.result) <- c('tvalue', 'pvalue', 'feature', 'feedback', 'significant')

maxt.a <- max(abs(as.numeric(as.character(a.result$tvalue)))) + 0.5
p_C <- ggplot(a.result, aes(x = factor(feature, levels = c('TATA', 'Inr', 'TATA:Inr', 'GC', 'CCAAT', 'GC:CCAAT')),
                           y = as.numeric(as.character(tvalue)))) + 
  geom_hline(aes(yintercept = 1.96), linetype = "dashed", size = 0.5, colour = "#F7C3B4") +
  geom_hline(aes(yintercept = -1.96), linetype = "dashed", size = 0.5, colour = "#B8CEBA") +
  geom_hline(aes(yintercept = 0), size = 0.5, colour = "#EBEAEA") +
  geom_point(aes(shape = factor(feedback, levels = c('positive', 'negative', 'non-feedback')),
                 colour = factor(significant, levels = c('non-significant', 'negative effect', 'positive effect'))),
             position = position_dodge(0.5), size = 1.5) +
  labs(y = 'Effect size(bf)', x = 'Promoter architecture') + 
  ylim(-maxt.a, maxt.a) + 
  theme_bw() +
  scale_shape_manual(values = c(15,16,17)) +
  scale_colour_manual(values = c("grey70",'#3C977B',"#F8766D")) +
  scale_fill_manual(values = c("grey70",'#3C977B',"#F8766D")) +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank())
p_C

## bs
bfit <- lm(log10(b) ~ (TATA * Inr + GC * CCAAT) : Hsign, data = analysis.mat)
bfit <- summary(bfit)$coef
name <- strsplit(rownames(bfit), ':Hsign')
feature <- sapply(name, function(x) {x[1]})
feedback <- sapply(name, function(x) {x[2]})
b.result <- data.frame(cbind(bfit[2:19,c(3,4)], feature = feature[2:19], feedback = feedback[2:19]))
bindex <- as.numeric(as.character(b.result[,2])) > 5e-2
b.result$sign[bindex] <- 'non-significant' 
b.result$sign[!bindex & as.numeric(as.character(b.result$t.value))>0] <- 'positive effect'
b.result$sign[!bindex & as.numeric(as.character(b.result$t.value))<0] <- 'negative effect'
colnames(b.result) <- c('tvalue', 'pvalue', 'feature', 'feedback', 'significant')

maxt.b <- max(abs(as.numeric(as.character(b.result$tvalue)))) + 0.5
p_D <- ggplot(b.result, aes(x = factor(feature, levels = c('TATA', 'Inr', 'TATA:Inr', 'GC', 'CCAAT', 'GC:CCAAT')),
                           y = as.numeric(as.character(tvalue)))) + 
  geom_hline(aes(yintercept = 1.96), linetype = "dashed", size = 0.5, colour = "#F7C3B4") +
  geom_hline(aes(yintercept = -1.96), linetype = "dashed", size = 0.5, colour = "#B8CEBA") +
  geom_hline(aes(yintercept = 0), size = 0.5, colour = "#EBEAEA") +
  geom_point(aes(shape = factor(feedback, levels = c('positive', 'negative', 'non-feedback')),
                 colour = factor(significant, levels = c('non-significant', 'negative effect', 'positive effect'))),
             position = position_dodge(0.5), size = 1.5) +
  labs(y = 'Effect size(bs)', x = 'Promoter architecture') + 
  ylim(-maxt.b, maxt.b) + 
  theme_bw() +
  scale_shape_manual(values = c(15,16,17)) +
  scale_colour_manual(values = c("grey70",'#3C977B',"#F8766D")) +
  scale_fill_manual(values = c("grey70",'#3C977B',"#F8766D")) +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank())
p_D

pic1 <- ggarrange(p_B, p_C, p_D, ncol = 1, nrow = 3, widths=c(6), heights=c(3.5), align = "v")
pic1

ggsave('fig4&S10/fig4B-D.pdf', width = 3, height = 4, useDingbats = FALSE)

###################################################################
# all gene regression
## rcv
rcvfit.all <- lm(rcv ~ TATA * Inr + GC * CCAAT, data = analysis.mat)
rcvfit.all <- summary(rcvfit.all)$coef
feature <- rownames(rcvfit.all)
rcv.result <- data.frame(cbind(rcvfit.all[2:7,c(3,4)], feature = feature[2:7]))
rcvindex <- as.numeric(as.character(rcv.result[,2])) > 5e-2
rcv.result$sign[rcvindex] <- 'non-significant' 
rcv.result$sign[!rcvindex & as.numeric(as.character(rcv.result$t.value)) > 0] <- 'positive effect'
rcv.result$sign[!rcvindex & as.numeric(as.character(rcv.result$t.value)) < 0] <- 'negative effect'
colnames(rcv.result) <- c('tvalue', 'pvalue', 'feature', 'significant')

maxt.rcv <- max(abs(as.numeric(as.character(rcv.result$tvalue)))) + 0.5
p_SC <- ggplot(rcv.result, aes(x = factor(feature, levels = c('TATA','Inr','TATA:Inr','GC','CCAAT','GC:CCAAT')),
                       y = as.numeric(as.character(tvalue)))) +
  geom_hline(aes(yintercept = 1.96), linetype = "dashed", size = 0.5, colour = "#F7C3B4") +
  geom_hline(aes(yintercept = -1.96), linetype = "dashed", size = 0.5, colour = "#B8CEBA") +
  geom_hline(aes(yintercept = 0), size = 0.5, colour = "#EBEAEA") +
  geom_point(aes(colour = factor(significant, levels = c('non-significant','negative effect','positive effect'))),
             position = position_dodge(0.5), size = 1.5) +
  labs(y = 'Effect size(rcv)', x = 'Promoter architecture') + 
  ylim(-maxt.rcv, maxt.rcv) + 
  theme_bw() +
  scale_colour_manual(values = c("grey70","#F8766D",'#3C977B')) +
  scale_fill_manual(values = c("grey70","#F8766D",'#3C977B')) +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank())
p_SC

## bf
afit.all <- lm(log10(a) ~ TATA * Inr + GC * CCAAT, data = analysis.mat)
afit.all <- summary(afit.all)$coef
feature <- rownames(afit.all)
a.result <- data.frame(cbind(afit.all[2:7,c(3,4)], feature = feature[2:7]))
aindex <- as.numeric(as.character(a.result[,2])) > 5e-2
a.result$sign[aindex] <- 'non-significant' 
a.result$sign[!aindex & as.numeric(as.character(a.result$t.value)) > 0] <- 'positive effect'
a.result$sign[!aindex & as.numeric(as.character(a.result$t.value)) < 0] <- 'negative effect'
colnames(a.result) <- c('tvalue', 'pvalue', 'feature', 'significant')

maxt.a <- max(abs(as.numeric(as.character(a.result$tvalue)))) + 0.5
p_SD <- ggplot(a.result, aes(x = factor(feature, levels = c('TATA','Inr','TATA:Inr','GC','CCAAT','GC:CCAAT')),
                       y = as.numeric(as.character(tvalue)))) +
  geom_hline(aes(yintercept = 1.96), linetype = "dashed", size = 0.5, colour = "#F7C3B4") +
  geom_hline(aes(yintercept = -1.96), linetype = "dashed", size = 0.5, colour = "#B8CEBA") +
  geom_hline(aes(yintercept = 0), size = 0.5, colour = "#EBEAEA") +
  geom_point(aes(colour = factor(significant, levels = c('non-significant','negative effect','positive effect'))),
             position = position_dodge(0.5), size = 1.5) +
  labs(y = 'Effect size(bf)', x = 'Promoter architecture') + 
  ylim(-maxt.a, maxt.a) + 
  theme_bw() +
  scale_colour_manual(values = c("grey70","#F8766D",'#3C977B')) +
  scale_fill_manual(values = c("grey70","#F8766D",'#3C977B')) +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank())
p_SD

## bs
bfit.all <- lm(log10(b) ~ TATA * Inr + GC * CCAAT, data = analysis.mat)
bfit.all <- summary(bfit.all)$coef
feature <- rownames(bfit.all)
b.result <- data.frame(cbind(bfit.all[2:7,c(3,4)], feature = feature[2:7]))
bindex <- as.numeric(as.character(b.result[,2])) > 5e-2
b.result$sign[bindex] <- 'non-significant' 
b.result$sign[!bindex & as.numeric(as.character(b.result$t.value)) > 0] <- 'positive effect'
b.result$sign[!bindex & as.numeric(as.character(b.result$t.value)) < 0] <- 'negative effect'
colnames(b.result) <- c('tvalue', 'pvalue', 'feature', 'significant')

maxt.b <- max(abs(as.numeric(as.character(b.result$tvalue)))) + 0.5
p_SE <- ggplot(b.result, aes(x = factor(feature, levels = c('TATA','Inr','TATA:Inr','GC','CCAAT','GC:CCAAT')),
                       y = as.numeric(as.character(tvalue)))) +
  geom_hline(aes(yintercept = 1.96), linetype = "dashed", size = 0.5, colour = "#F7C3B4") +
  geom_hline(aes(yintercept = -1.96), linetype = "dashed", size = 0.5, colour = "#B8CEBA") +
  geom_hline(aes(yintercept = 0), size = 0.5, colour = "#EBEAEA") +
  geom_point(aes(colour = factor(significant, levels = c('non-significant','negative effect','positive effect'))),
             position = position_dodge(0.5), size = 1.5) +
  labs(y = 'Effect size(bs)', x = 'Promoter architecture') + 
  ylim(-maxt.b, maxt.b) + 
  theme_bw() +
  scale_colour_manual(values = c("grey70","#F8766D",'#3C977B')) +
  scale_fill_manual(values = c("grey70","#F8766D",'#3C977B')) +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank())
p_SE

pic2 <- ggarrange(p_SC, p_SD, p_SE, ncol = 1, nrow = 3, widths=c(6), heights=c(3.5), align = "v")
pic2

ggsave('fig4&S10/figS10C-E.pdf', width = 3, height = 4, useDingbats = FALSE)

###################################################################
# the effect of TATA-box on rcv
analysis.mat.po = analysis.mat[analysis.mat$Hsign == "positive",]
analysis.mat.po$rank.rcv <- rank(analysis.mat.po$rcv)
logit.fit <- glm(TATA ~ rank.rcv,
                 family = binomial(link = 'logit'),
                 data = analysis.mat.po)
analysis.mat.po$predict <- predict(logit.fit, type = "response", newdata = analysis.mat.po)
roc.po <- roc(analysis.mat.po, TATA, predict)
auc(roc.po)

analysis.mat.no = analysis.mat[analysis.mat$Hsign == "non-feedback",]
analysis.mat.no$rank.rcv <- rank(analysis.mat.no$rcv)
logit.fit <- glm(TATA ~ rank.rcv,
                 family = binomial(link = 'logit'),
                 data = analysis.mat.no)
analysis.mat.no$predict <- predict(logit.fit, type = "response", newdata = analysis.mat.no)
roc.no <- roc(analysis.mat.no, TATA, predict)
auc(roc.no)

analysis.mat.ne = analysis.mat[analysis.mat$Hsign == "negative",]
analysis.mat.ne$rank.rcv <- rank(analysis.mat.ne$rcv)
logit.fit <- glm(TATA ~ rank.rcv,
                 family = binomial(link = 'logit'),
                 data = analysis.mat.ne)
analysis.mat.ne$predict <- predict(logit.fit, type = "response", newdata = analysis.mat.ne)
roc.ne <- roc(analysis.mat.ne, TATA, predict)
auc(roc.ne)

p_E <- ggroc(list(positive = roc.po, non = roc.no, negative = roc.ne), size = 0.2) + 
  theme_bw() +
  scale_colour_manual(values = c("#FF0000", "#4595D1", "#39AE36")) + 
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_E

ggsave('fig4&S10/fig4E.pdf', width = 1.7, height = 1.6, useDingbats = FALSE)

# the effect of TATA-box on bf and bs 
mean.notata <- apply(log10(analysis.mat[analysis.mat$TATA == 0, c("a", "b")]), 2, mean)
mean.tata.po <- apply(log10(analysis.mat[analysis.mat$TATA == 1 & analysis.mat$Hsign == "positive", c("a", "b")]), 2, mean)
mean.tata.no <- apply(log10(analysis.mat[analysis.mat$TATA == 1 & analysis.mat$Hsign == "non-feedback", c("a", "b")]), 2, mean)
mean.tata.ne <- apply(log10(analysis.mat[analysis.mat$TATA == 1 & analysis.mat$Hsign == "negative", c("a", "b")]), 2, mean)

se.notata <- apply(log10(analysis.mat[analysis.mat$TATA == 0, c("a", "b")]), 2, function(x) {sd(x)/sqrt(length(x))})
se.tata.po <- apply(log10(analysis.mat[analysis.mat$TATA == 1 & analysis.mat$Hsign == "positive", c("a", "b")]), 2, function(x) {sd(x)/sqrt(length(x))})
se.tata.no <- apply(log10(analysis.mat[analysis.mat$TATA == 1 & analysis.mat$Hsign == "non-feedback", c("a", "b")]), 2, function(x) {sd(x)/sqrt(length(x))})
se.tata.ne <- apply(log10(analysis.mat[analysis.mat$TATA == 1 & analysis.mat$Hsign == "negative", c("a", "b")]), 2, function(x) {sd(x)/sqrt(length(x))})

mean.mat <- t(data.frame(mean.notata, mean.tata.po, mean.tata.no, mean.tata.ne))
se.mat <- t(data.frame(se.notata, se.tata.po, se.tata.no, se.tata.ne))
meanse.mat <- data.frame(cbind(mean.mat, se.mat))
colnames(meanse.mat) <- c("a.mean", "b.mean", "a.se", "b.se")
meanse.mat$group = c("noTATA", "TATA.po", "TATA.no", "TATA.ne")

p_F <- ggplot(meanse.mat, aes(x = a.mean, y = b.mean, colour = group, fill = group)) +
  geom_linerange(aes(xmin = a.mean - a.se, xmax = a.mean + a.se, colour = group)) + 
  geom_linerange(aes(ymin = b.mean - b.se, ymax = b.mean + b.se, colour = group)) +
  geom_point(size = 1, shape = 21) +
  labs(x = expression(log[10](bf)), y = expression(log[10](bs))) + 
  scale_fill_manual(values = c("gray80", "#C2DFAE", "#B2C2CE", "#F29E96")) +
  scale_colour_manual(values = c("gray50", "#39AE36", "#4595D1", "#FF0000")) + 
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_F

ggsave('fig4&S10/fig4F.pdf', width = 1.7, height = 1.7, useDingbats = FALSE)


# the relation between promoter architecture and mean level of expression
colourbar <- c('#888584', "#F9AD42")
fillbar <- c('#E3E2E0', "#FFDB8A")
p_SB1 <- ggplot(analysis.mat,aes(x = factor(TATA), y = log10(mean), fill = factor(TATA), colour = factor(TATA))) + 
  geom_boxplot(outlier.size = 0, outlier.colour = 'black', size = 0.3, notch = TRUE) +
  geom_signif(comparisons = list(c('0','1')), map_signif_level = F, size = 0.2, test = wilcox.test, textsize = 2, color = "black") +
  labs(y = expression(log[10](mean)), x = 'TATA') +
  ylim(0,2.1) +
  scale_colour_manual(values = colourbar) +
  scale_fill_manual(values = fillbar) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_SB1

p_SB2 <- ggplot(analysis.mat,aes(x = factor(Inr), y = log10(mean), fill = factor(Inr), colour = factor(Inr))) + 
  geom_boxplot(outlier.size=0, outlier.colour = 'black', size = 0.3, notch = TRUE) +
  geom_signif(comparisons = list(c('0','1')), map_signif_level = F, size = 0.2, test = wilcox.test, textsize = 2, color = "black") +
  labs(y = expression(log[10](mean)), x = 'Inr') +
  ylim(0,2.1) +
  scale_colour_manual(values = colourbar) +
  scale_fill_manual(values = fillbar) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_SB2

p_SB3 <- ggplot(analysis.mat,aes(x = factor(GC), y = log10(mean), fill = factor(GC), colour = factor(GC))) + 
  geom_boxplot(outlier.size=0, outlier.colour = 'black', size = 0.3, notch = TRUE) +
  geom_signif(comparisons = list(c('0','1')), map_signif_level = F, size = 0.2, test = wilcox.test, textsize = 2, color = "black") +
  labs(y = expression(log[10](mean)), x = 'GC') +
  ylim(0,2.1) +
  scale_colour_manual(values = colourbar) +
  scale_fill_manual(values = fillbar) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_SB3

p_SB4 <- ggplot(analysis.mat,aes(x = factor(CCAAT), y = log10(mean), fill = factor(CCAAT), colour = factor(CCAAT))) + 
  geom_boxplot(outlier.size=0, outlier.colour = 'black', size = 0.3, notch = TRUE) +
  geom_signif(comparisons = list(c('0','1')), map_signif_level = F, size = 0.2, test = wilcox.test, textsize = 2, color = "black") +
  labs(y = expression(log[10](mean)), x = 'CCAAT') +
  ylim(0,2.1) +
  scale_colour_manual(values = colourbar) +
  scale_fill_manual(values = fillbar) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_SB4

pic4 <- ggarrange(p_SB1, p_SB2, p_SB3, p_SB4, ncol = 4, nrow = 1, widths=c(1.2), heights=c(2), align = "v")

ggsave('fig4&S10/figS10A.pdf', width = 3.9, height = 1.7, useDingbats = FALSE)

###################################################################
# the relation between promoter architecture and feedback
attach(analysis.mat)
cont.TATA <- table(TATA, Hsign)
prop.TATA <- as.data.frame(prop.table(cont.TATA, 1))
p_SA1 <- ggplot(prop.TATA, aes(TATA, Hsign)) + 
  geom_tile(aes(fill = Freq), colour = "white") + 
  scale_fill_gradient(low = "white",high = "#4595D1", limits = c(0.1,0.6))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position= 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_SA1

cont.Inr <- table(Inr, Hsign)
prop.Inr <- as.data.frame(prop.table(cont.Inr, 1))
p_SA2 <- ggplot(prop.Inr, aes(Inr, Hsign)) + 
  geom_tile(aes(fill = Freq), colour = "white") + 
  scale_fill_gradient(low = "white",high = "#EB9A18", limits = c(0.1,0.6))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position= 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_SA2

cont.GC <- table(GC, Hsign)
prop.GC <- as.data.frame(prop.table(cont.GC, 1))
p_SA3 <- ggplot(prop.GC, aes(GC, Hsign)) + 
  geom_tile(aes(fill = Freq), colour = "white") + 
  scale_fill_gradient(low = "white",high = "#269B59", limits = c(0.1,0.6))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position= 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_SA3

cont.CCAAT <- table(CCAAT, Hsign)
prop.CCAAT <- as.data.frame(prop.table(cont.CCAAT, 1))
p_SA4 <- ggplot(prop.CCAAT, aes(CCAAT, Hsign)) + 
  geom_tile(aes(fill = Freq), colour = "white") + 
  scale_fill_gradient(low = "white",high = "#DF4C4F", limits = c(0.1,0.6))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position= 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid.minor = element_blank())
p_SA4

detach(analysis.mat)
pic3 <- ggarrange(p_SA1, p_SA2, p_SA3, p_SA4, ncol = 4, nrow = 1, widths=c(2), heights=c(3), align = "v")
pic3

ggsave('fig4&S10/figS10B.pdf', width = 4, height = 2, useDingbats = FALSE)

