library("ggplot2")
library("ggpubr")

genename <- read.table("genename.txt", sep = '\t', header = FALSE)
result.all <- read.table("../code_Inference/Hierarchical_model/results_MEF.csv", header = FALSE, sep = ',')
result.all <- cbind(genename, result.all)
colnames(result.all) <- c('genename', 'a', 'b', 'k', 'H', 'mean', 'noise', 'fval', 'goodFit')
result.all$Hsign[result.all$H > 0] <- 'negative'
result.all$Hsign[result.all$H < 0] <- 'positive'
result.all$Hsign[result.all$H == 0] <- 'non-feedback'
result.all = result.all[result.all$a > 0 & result.all$goodFit == 1,]
breaks <- unname(quantile(log10(result.all$mean),seq(0.0,1,0.2)))
breaks[1] <- breaks[1]-0.1
result.all$level <- cut(log10(result.all$mean), breaks, 1:5, ordered_result = T)
result.statis <- data.frame(Hsign = result.all$Hsign, level = result.all$level, counts = 1)

result.statis$Hsign <- as.factor(result.statis$Hsign)
result.statis$level <- as.factor(result.statis$level)
result.statis$counts <- as.numeric(result.statis$counts)

attach(result.statis)
result.statis2 <- aggregate(result.statis[,c('counts')], by = list(Hsign, level), FUN = sum, na.rm = TRUE)
result.statis3 <- aggregate(result.statis[,c('counts')], by = list(Hsign), FUN = sum, na.rm = TRUE)
detach(result.statis)

# color.feedback <- c("#FD7F32", "#547CBE", "#00B74A")
color.feedback <- c("#3AA438", "#547CBE", "#E71419")
mat <- log10(result.all[,c("noise", "a", "b")])
mat$Hsign <- result.all$Hsign
attach(mat)
median.mat <- aggregate(mat, by = list(Hsign), FUN = median)

# pdf of noise 
p_A <- ggplot(result.all, aes(x = log10(noise), colour = factor(Hsign, levels = c('negative', 'non-feedback', 'positive'))))+
  geom_line(stat = "density", adjust = 1.5, size = 0.45) +
  geom_vline(aes(xintercept = noise, colour = Group.1), data = median.mat, linetype = "dashed", size = 0.45) +
  labs(x = expression(log[10](CV2)), y = "density") +
  theme_bw() +
  scale_colour_manual(values = color.feedback)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank()) 
p_A

# pdf of bf
p_B <- ggplot(result.all, aes(x = log10(a), colour = factor(Hsign, levels = c('negative', 'non-feedback', 'positive'))))+
  geom_line(stat = "density", adjust = 1.5, size = 0.45) +
  geom_vline(aes(xintercept = a, colour = Group.1), data = median.mat, linetype = "dashed", size = 0.45) +
  labs(x = expression(log[10](bf)), y = "density") +
  theme_bw() +
  scale_colour_manual(values = color.feedback)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank()) 
p_B

# pdf of bs
p_C <- ggplot(result.all, aes(x = log10(b), colour = factor(Hsign, levels = c('negative', 'non-feedback', 'positive'))))+
  geom_line(stat = "density", adjust = 1.5, size = 0.45) +
  geom_vline(aes(xintercept = b, colour = Group.1), data = median.mat, linetype = "dashed", size = 0.45) +
  labs(x = expression(log[10](bs)), y = "density") +
  theme_bw() +
  scale_colour_manual(values = color.feedback)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank()) 
p_C
detach(mat)

# boxplot of noise as the mean changes
attach(result.all)
c_mean <- aggregate(log10(noise), by = list(level, Hsign), mean)
colnames(c_mean) <- c('level', 'Hsign', 'y')
c_mean <- cbind(c_mean, x = c(0.75, 1.75, 2.75, 3.75, 4.75, 1, 2, 3, 4, 5, 1.25, 2.25, 3.25, 4.25, 5.25))

p_D <- ggplot(result.all, aes(x = level, y = log10(noise),
                              fill = Hsign, colour = Hsign)) +
  geom_boxplot(outlier.size = 0, outlier.colour = 'gray70', fatten = 1, size =  0.3, alpha = 0.45) +
  geom_line(data = c_mean, aes(x = x, y = y, colour = Hsign), size = 0.3, linetype = "dashed") +
  labs(y = expression(log[10](CV2))) + 
  theme_bw()  + 
  ylim(-1,1) +
  scale_fill_manual(values = color.feedback)+
  scale_colour_manual(values = color.feedback)+
  theme(legend.title = element_blank(),
        legend.position = 'none',
        axis.title.x =  element_blank(),
        axis.text = element_text(size=6, color="black"),
        axis.line.x = element_line(size = 0.2, colour = "black"),
        axis.line.y = element_line(size = 0.2, colour = "black"),
        axis.title = element_text(size=8, color="black"),
        axis.ticks = element_line(color="black",size=0.25,lineend = 10),
        panel.border = element_blank(),
        panel.grid = element_blank())
p_D

# boxplot of bf as the mean changes
b_mean <- aggregate(log10(a), by = list(level, Hsign), mean)
colnames(b_mean) <- c('level', 'Hsign', 'y')
b_mean <- cbind(b_mean, x = c(0.75, 1.75, 2.75, 3.75, 4.75, 1, 2, 3, 4, 5, 1.25, 2.25, 3.25, 4.25, 5.25))
p_E <- ggplot(result.all, aes(x = level, y = log10(a), fill = Hsign, colour = Hsign)) +
  geom_boxplot(outlier.size = 0, outlier.colour = 'gray70', fatten = 1, size = 0.3, alpha = 0.45) +
  geom_line(data = b_mean, aes(x = x, y = y, colour = Hsign), size = 0.3, linetype = "dashed") +
  labs(y = expression(log[10](bf))) + 
  theme_bw()  + 
  ylim(-0.6,1.5) +
  scale_fill_manual(values = color.feedback)+
  scale_colour_manual(values = color.feedback)+
  theme(legend.title = element_blank(),
        legend.position = 'none',
        axis.title.x =  element_blank(),
        axis.text = element_text(size=6, color="black"),
        axis.line.x = element_line(size = 0.2, colour = "black"),
        axis.line.y = element_line(size = 0.2, colour = "black"),
        axis.title = element_text(size=8, color="black"),
        axis.ticks = element_line(color="black",size=0.25,lineend = 10),
        panel.border = element_blank(),
        panel.grid = element_blank())
p_E

# boxplot of bs as the mean changes
b_mean <- aggregate(log10(b), by = list(level, Hsign), mean)
colnames(b_mean) <- c('level', 'Hsign', 'y')
b_mean <- cbind(b_mean, x = c(0.75, 1.75, 2.75, 3.75, 4.75, 1, 2, 3, 4, 5, 1.25, 2.25, 3.25, 4.25, 5.25))
p_F <- ggplot(result.all, aes(x = level, y = log10(b), fill = Hsign, colour = Hsign)) +
  geom_boxplot(outlier.size = 0, outlier.colour = 'gray70', fatten = 1, size =  0.3, alpha = 0.45) +
  geom_line(data = b_mean, aes(x = x, y = y, colour = Hsign), size = 0.3, linetype = "dashed") +
  labs(y = expression(log[10](bs))) + 
  theme_bw()  + 
  ylim(0,1.2) +
  scale_fill_manual(values = color.feedback)+
  scale_colour_manual(values = color.feedback)+
  theme(legend.title = element_blank(),
        legend.position = 'none',
        axis.title.x =  element_blank(),
        axis.text = element_text(size=6, color="black"),
        axis.line.x = element_line(size = 0.2, colour = "black"),
        axis.line.y = element_line(size = 0.2, colour = "black"),
        axis.title = element_text(size=8, color="black"),
        axis.ticks = element_line(color="black",size=0.25,lineend = 10),
        panel.border = element_blank(),
        panel.grid = element_blank())
p_F

detach(result.all)

pic3 <- ggarrange(p_A, p_D, p_B, p_E, p_C, p_F, ncol = 2, nrow = 3, widths = c(2.5,4.5), heights = c(2), align = "hv")
ggsave('fig3&S9/fig3.pdf', width = 5.6, height = 4,useDingbats = FALSE)
