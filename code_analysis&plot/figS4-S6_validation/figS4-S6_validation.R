library(ggplot2)
library(reshape2)
library(ggpubr)

########################################################################
result.cellnumber <- read.table('../code_Inference/Hierarchical_model/results/validation_CellNumber_results.csv',sep = ',',header = F)
colnames(result.cellnumber)<-c('times','cell.number','a.true','b.true','k.true','h.true','a.est','b.est','k.est','h.est','mean','noise','fval')
result.cellnumber$h.type[result.cellnumber$h.est == 0] <- 'non-feedback' 
result.cellnumber$h.type[result.cellnumber$h.est < 0] <- 'positive'
result.cellnumber$h.type[result.cellnumber$h.est > 0] <- 'negative'


result.sensitive <- read.table('../code_Inference/Hierarchical_model/results/validation_Sensitive_results.csv',sep = ',',header = F)
colnames(result.sensitive)<-c('times','sensitive','a.true','b.true','k.true','h.true','a.est','b.est','k.est','h.est','mean','noise','fval')
result.sensitive$h.type[result.sensitive$h.est == 0] <- 'non-feedback' 
result.sensitive$h.type[result.sensitive$h.est < 0] <- 'positive'
result.sensitive$h.type[result.sensitive$h.est > 0] <- 'negative'

theta.true.set <- rbind(c(3,4,0,0), c(5,10,0,0), c(2,15,2,-3), c(4,10,10,-3), c(3,18,2,3), c(2,10,10,3))
colnames(theta.true.set) <- c('a.true', 'b.true', 'k.true', 'h.true')

colour.feedback <- c("#E71419", "#547CBE", "#3AA438")
fill.feedback <- c("#ED9899", "#B0C1DC", "#A6D2A6")

# theta_true = [3 4 0 0] 
temp <- 1
result.cellnumber.temp <- result.cellnumber[result.cellnumber[,3] == theta.true.set[temp,1] & result.cellnumber[,4] == theta.true.set[temp,2] & 
  result.cellnumber[,5] == theta.true.set[temp,3] & result.cellnumber[,6] == theta.true.set[temp,4],]

p1 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200','300','500','1000','5000')), y = log10(a.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,1])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Cell number', y = expression(log[10](bf))) +
  ylim(0.3,0.7) + 
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p1

p2 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200','300','500','1000','5000')), y = log10(b.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,2])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Cell number', y = expression(log[10](bs))) +
  ylim(0.4,0.8) + 
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p2

p3 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200', '300', '500', '1000', '5000')),
                        fill = factor(h.type, level = c('positive', 'non-feedback', 'negative')),
                        colour = factor(h.type, level = c('positive', 'non-feedback', 'negative')))) +
  geom_bar(position = "fill", width = 0.7, size = 0.15) + 
  labs(x = 'Cell number', y = 'Model selection') + 
  scale_fill_manual(values = fill.feedback) +
  scale_colour_manual(values = colour.feedback) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p3

result.sensitive.temp <- result.sensitive[result.sensitive[,3] == theta.true.set[temp,1] & result.sensitive[,4] == theta.true.set[temp,2] & 
                                              result.sensitive[,5] == theta.true.set[temp,3] & result.sensitive[,6] == theta.true.set[temp,4],]

p4 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')), y = log10(a.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,1])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Sensitive', y = expression(log[10](bf))) +
  ylim(0.3,0.7) + 
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p4

p5 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')), y = log10(b.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,2])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Sensitive', y = expression(log[10](bs))) +
  ylim(0.4,0.8) + 
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p5

p6 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')),
                                         fill = factor(h.type, level = c('positive', 'non-feedback', 'negative')),
                                         colour = factor(h.type, level = c('positive', 'non-feedback', 'negative')))) +
  geom_bar(position = "fill", width = 0.7, size = 0.15) + 
  labs(x = 'Sensitive', y = 'Model selection') + 
  scale_fill_manual(values = fill.feedback) +
  scale_colour_manual(values = colour.feedback) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p6


pic1 <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, widths=c(3,3,3), heights=c(3.5), align = "hv")
pic1

ggsave('figS4-S6_validation/figS_validation_3_4_0_0.pdf', width = 3, height = 3.2, useDingbats = FALSE)
  
# theta_true = [5 10 0 0] 
temp <- 2
result.cellnumber.temp <- result.cellnumber[result.cellnumber[,3] == theta.true.set[temp,1] & result.cellnumber[,4] == theta.true.set[temp,2] & 
                                              result.cellnumber[,5] == theta.true.set[temp,3] & result.cellnumber[,6] == theta.true.set[temp,4],]

p1 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200','300','500','1000','5000')), y = log10(a.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,1])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Cell number', y = expression(log[10](bf))) +
  ylim(0.5,0.9) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p1

p2 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200','300','500','1000','5000')), y = log10(b.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,2])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Cell number', y = expression(log[10](bs))) +
  ylim(0.8,1.2) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p2

p3 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200', '300', '500', '1000', '5000')),
                                         fill = factor(h.type, level = c('positive', 'non-feedback', 'negative')),
                                         colour = factor(h.type, level = c('positive', 'non-feedback', 'negative')))) +
  geom_bar(position = "fill", width = 0.7, size = 0.15) + 
  labs(x = 'Cell number', y = 'Model selection') + 
  scale_fill_manual(values = fill.feedback) +
  scale_colour_manual(values = colour.feedback) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p3

result.sensitive.temp <- result.sensitive[result.sensitive[,3] == theta.true.set[temp,1] & result.sensitive[,4] == theta.true.set[temp,2] & 
                                            result.sensitive[,5] == theta.true.set[temp,3] & result.sensitive[,6] == theta.true.set[temp,4],]

p4 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')), y = log10(a.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,1])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Sensitive', y = expression(log[10](bf))) +
  ylim(0.5,0.9) + 
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p4

p5 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')), y = log10(b.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,2])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Sensitive', y = expression(log[10](bs))) +
  ylim(0.8,1.2) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p5

p6 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')),
                                        fill = factor(h.type, level = c('positive', 'non-feedback', 'negative')),
                                        colour = factor(h.type, level = c('positive', 'non-feedback', 'negative')))) +
  geom_bar(position = "fill", width = 0.7, size = 0.15) + 
  labs(x = 'Sensitive', y = 'Model selection') + 
  scale_fill_manual(values = c( "#B0C1DC", "#A6D2A6", "#ED9899")) +
  scale_colour_manual(values = c("#547CBE", "#3AA438", "#E71419")) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p6


pic1 <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, widths=c(3,3,3), heights=c(3.5), align = "hv")
pic1

ggsave('figS4-S6_validation/figS_validation_5_10_0_0.pdf', width = 3, height = 3.2, useDingbats = FALSE)

# theta_true = [2 15 2 -3] 
temp <- 3
result.cellnumber.temp <- result.cellnumber[result.cellnumber[,3] == theta.true.set[temp,1] & result.cellnumber[,4] == theta.true.set[temp,2] & 
                                              result.cellnumber[,5] == theta.true.set[temp,3] & result.cellnumber[,6] == theta.true.set[temp,4],]

p1 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200','300','500','1000','5000')), y = log10(a.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,1])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Cell number', y = expression(log[10](bf))) +
  ylim(0.2,0.75) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p1

p2 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200','300','500','1000','5000')), y = log10(b.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,2])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Cell number', y = expression(log[10](bs))) +
  ylim(0.8,1.3) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p2

p3 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200', '300', '500', '1000', '5000')),
                                         fill = factor(h.type, level = c('positive', 'non-feedback', 'negative')),
                                         colour = factor(h.type, level = c('positive', 'non-feedback', 'negative')))) +
  geom_bar(position = "fill", width = 0.7, size = 0.15) + 
  labs(x = 'Cell number', y = 'Model selection') + 
  scale_fill_manual(values = fill.feedback) +
  scale_colour_manual(values = colour.feedback) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p3

result.sensitive.temp <- result.sensitive[result.sensitive[,3] == theta.true.set[temp,1] & result.sensitive[,4] == theta.true.set[temp,2] & 
                                            result.sensitive[,5] == theta.true.set[temp,3] & result.sensitive[,6] == theta.true.set[temp,4],]

p4 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')), y = log10(a.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,1])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Sensitive', y = expression(log[10](bf))) +
  ylim(0.2,0.75) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p4

p5 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')), y = log10(b.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,2])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Sensitive', y = expression(log[10](bs))) +
  ylim(0.8,1.3) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p5

p6 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')),
                                        fill = factor(h.type, level = c('positive', 'non-feedback', 'negative')),
                                        colour = factor(h.type, level = c('positive', 'non-feedback', 'negative')))) +
  geom_bar(position = "fill", width = 0.7, size = 0.15) + 
  labs(x = 'Sensitive', y = 'Model selection') + 
  scale_fill_manual(values = fill.feedback) +
  scale_colour_manual(values = colour.feedback) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p6


pic1 <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, widths=c(3,3,3), heights=c(3.5), align = "hv")
pic1

ggsave('figS4-S6_validation/figS_validation_2_15_2_-3.pdf', width = 3, height = 3.2, useDingbats = FALSE)

# theta_true = [4 10 10 -3] 
temp <- 4
result.cellnumber.temp <- result.cellnumber[result.cellnumber[,3] == theta.true.set[temp,1] & result.cellnumber[,4] == theta.true.set[temp,2] & 
                                              result.cellnumber[,5] == theta.true.set[temp,3] & result.cellnumber[,6] == theta.true.set[temp,4],]

p1 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200','300','500','1000','5000')), y = log10(a.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,1])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Cell number', y = expression(log[10](bf))) +
  ylim(0.2,0.75) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p1

p2 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200','300','500','1000','5000')), y = log10(b.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,2])), colour = '#D95319', linetype = 'dashed',size = 0.2) + 
  labs(x = 'Cell number', y = expression(log[10](bs))) +
  ylim(0.8,1.3) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p2

p3 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200', '300', '500', '1000', '5000')),
                                         fill = factor(h.type, level = c('positive', 'non-feedback', 'negative')),
                                         colour = factor(h.type, level = c('positive', 'non-feedback', 'negative')))) +
  geom_bar(position = "fill", width = 0.7, size = 0.15) + 
  labs(x = 'Cell number', y = 'Model selection') + 
  scale_fill_manual(values = fill.feedback) +
  scale_colour_manual(values = colour.feedback) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p3

result.sensitive.temp <- result.sensitive[result.sensitive[,3] == theta.true.set[temp,1] & result.sensitive[,4] == theta.true.set[temp,2] & 
                                            result.sensitive[,5] == theta.true.set[temp,3] & result.sensitive[,6] == theta.true.set[temp,4],]

p4 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')), y = log10(a.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,1])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Sensitive', y = expression(log[10](bf))) +
  ylim(0.2,0.75) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p4

p5 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')), y = log10(b.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,2])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Sensitive', y = expression(log[10](bs))) +
  ylim(0.8,1.3) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p5

p6 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')),
                                        fill = factor(h.type, level = c('positive', 'non-feedback', 'negative')),
                                        colour = factor(h.type, level = c('positive', 'non-feedback', 'negative')))) +
  geom_bar(position = "fill", width = 0.7, size = 0.15) + 
  labs(x = 'Sensitive', y = 'Model selection') + 
  scale_fill_manual(values = fill.feedback) +
  scale_colour_manual(values = colour.feedback) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p6


pic1 <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, widths=c(3,3,3), heights=c(3.5), align = "hv")
pic1

ggsave('figS4-S6_validation/figS_validation_4_10_10_-3.pdf', width = 3, height = 3.2, useDingbats = FALSE)

# theta_true = [3 18 2 3] 
temp <- 5
result.cellnumber.temp <- result.cellnumber[result.cellnumber[,3] == theta.true.set[temp,1] & result.cellnumber[,4] == theta.true.set[temp,2] & 
                                              result.cellnumber[,5] == theta.true.set[temp,3] & result.cellnumber[,6] == theta.true.set[temp,4],]

p1 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200','300','500','1000','5000')), y = log10(a.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,1])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Cell number', y = expression(log[10](bf))) +
  ylim(0,1.5) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p1

p2 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200','300','500','1000','5000')), y = log10(b.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,2])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Cell number', y = expression(log[10](bs))) +
  ylim(0.7,1.4) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p2

p3 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200', '300', '500', '1000', '5000')),
                                         fill = factor(h.type, level = c('positive', 'non-feedback', 'negative')),
                                         colour = factor(h.type, level = c('positive', 'non-feedback', 'negative')))) +
  geom_bar(position = "fill", width = 0.7, size = 0.15) + 
  labs(x = 'Cell number', y = 'Model selection') + 
  scale_fill_manual(values = fill.feedback) +
  scale_colour_manual(values = colour.feedback) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p3

result.sensitive.temp <- result.sensitive[result.sensitive[,3] == theta.true.set[temp,1] & result.sensitive[,4] == theta.true.set[temp,2] & 
                                            result.sensitive[,5] == theta.true.set[temp,3] & result.sensitive[,6] == theta.true.set[temp,4],]

p4 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')), y = log10(a.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,1])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Sensitive', y = expression(log[10](bf))) +
  ylim(0,1.5) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p4

p5 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')), y = log10(b.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,2])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Sensitive', y = expression(log[10](bs))) +
  ylim(0.7,1.4) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p5

p6 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')),
                                        fill = factor(h.type, level = c('positive', 'non-feedback', 'negative')),
                                        colour = factor(h.type, level = c('positive', 'non-feedback', 'negative')))) +
  geom_bar(position = "fill", width = 0.7, size = 0.15) + 
  labs(x = 'Sensitive', y = 'Model selection') + 
  scale_fill_manual(values = c("#B0C1DC", "#A6D2A6", "#ED9899")) +
  scale_colour_manual(values = c("#547CBE", "#3AA438", "#E71419")) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p6


pic1 <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, widths=c(3,3,3), heights=c(3.5), align = "hv")
pic1

ggsave('figS4-S6_validation/figS_validation_3_18_2_3.pdf', width = 3, height = 3.2, useDingbats = FALSE)


# theta_true = [5 10 10 3] 
temp <- 6
result.cellnumber.temp <- result.cellnumber[result.cellnumber[,3] == theta.true.set[temp,1] & result.cellnumber[,4] == theta.true.set[temp,2] & 
                                              result.cellnumber[,5] == theta.true.set[temp,3] & result.cellnumber[,6] == theta.true.set[temp,4],]

p1 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200','300','500','1000','5000')), y = log10(a.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,1])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Cell number', y = expression(log[10](bf))) +
  ylim(0,0.5) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p1

p2 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200','300','500','1000','5000')), y = log10(b.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,2])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Cell number', y = expression(log[10](bs))) +
  ylim(0.5,1.2) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p2

p3 <- ggplot(result.cellnumber.temp, aes(x = factor(cell.number, level = c('200', '300', '500', '1000', '5000')),
                                         fill = factor(h.type, level = c('positive', 'non-feedback', 'negative')),
                                         colour = factor(h.type, level = c('positive', 'non-feedback', 'negative')))) +
  geom_bar(position = "fill", width = 0.7, size = 0.15) + 
  labs(x = 'Cell number', y = 'Model selection') + 
  scale_fill_manual(values = c( "#B0C1DC", "#A6D2A6", "#ED9899")) +
  scale_colour_manual(values = c("#547CBE", "#3AA438", "#E71419")) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p3

result.sensitive.temp <- result.sensitive[result.sensitive[,3] == theta.true.set[temp,1] & result.sensitive[,4] == theta.true.set[temp,2] & 
                                            result.sensitive[,5] == theta.true.set[temp,3] & result.sensitive[,6] == theta.true.set[temp,4],]

p4 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')), y = log10(a.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,1])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Sensitive', y = expression(log[10](bf))) +
  ylim(0,0.5) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p4

p5 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')), y = log10(b.est))) +
  geom_boxplot(outlier.colour = '#D8D8D8', outlier.size = 0, size = 0.15, colour = '#3E3938', fill = '#EEEEEE', fatten = 1) +
  geom_hline(aes(yintercept = log10(theta.true.set[temp,2])), colour = '#D95319', linetype = 'dashed', size = 0.2) + 
  labs(x = 'Sensitive', y = expression(log[10](bs))) +
  ylim(0.5,1.2) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p5

p6 <- ggplot(result.sensitive.temp, aes(x = factor(sensitive, level = c('0.1','0.3','0.5','0.7','0.9')),
                                        fill = factor(h.type, level = c('positive', 'non-feedback', 'negative')),
                                        colour = factor(h.type, level = c('positive', 'non-feedback', 'negative')))) +
  geom_bar(position = "fill", width = 0.7, size = 0.15) + 
  labs(x = 'Sensitive', y = 'Model selection') + 
  scale_fill_manual(values = fill.feedback) +
  scale_colour_manual(values = colour.feedback) +
  theme_bw() +
  theme(legend.position= 'none',
        axis.title = element_text(colour = 'black', size = 6),
        axis.text = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(angle = 45, colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p6


pic1 <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, widths=c(3,3,3), heights=c(3.5), align = "hv")
pic1

ggsave('figS4-S6_validation/figS_validation_2_10_10_3.pdf', width = 3, height = 3.2, useDingbats = FALSE)
