### supplementary figure - mean vs param
genename <- read.table("genename.txt", sep = '\t', header = FALSE)
result.all <- read.table("../code_Inference/PoissonFriedman_model/results_MEF.csv", header = FALSE, sep = ',')
result.all <- cbind(genename, result.all)
colnames(result.all) <- c('genename', 'a', 'b', 'k', 'H', 'mean', 'noise', 'fval', 'goodFit')
result.all$Hsign[result.all$H > 0] <- 'negative'
result.all$Hsign[result.all$H < 0] <- 'positive'
result.all$Hsign[result.all$H == 0] <- 'non-feedback'
result.all$rcv <- result.all$noise-(1/log2(result.all$mean))
result.all = result.all[result.all$a > 0 & result.all$goodFit == 1,]

# mean vs noise bf bs
p1 <- ggplot(result.all, aes(x = log10(mean), y = log10(noise))) + 
  geom_point(colour = 'grey70', size = 0.2)+
  geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), alpha = 0.2, size = 0.5, colour = "#008CD6", fill = "#008CD6") + 
  labs(y = expression(log[10](CV2)), x = expression(log[10](mean))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank())
p1

p2 <- ggplot(result.all, aes(x = log10(mean), y = rcv)) + 
  geom_point(colour = 'grey70', size = 0.2)+
  geom_smooth(method = "lm", alpha = 0.2, size = 0.5, colour = "#008CD6", fill = "#008CD6") + 
  labs(y = expression(rCV2), x = expression(log[10](mean))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank())
p2

p3 <- ggplot(result.all, aes(x = log10(mean), y = log10(a))) + 
  geom_point(colour = 'grey70', size = 0.2)+
  geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), alpha = 0.2, size = 0.5, colour = "#008CD6", fill = "#008CD6") + 
  labs(y = expression(log[10](bf)), x = expression(log[10](mean))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank())
p3

p4 <- ggplot(result.all, aes(x = log10(mean), y = log10(b))) + 
  geom_point(colour = 'grey70', size = 0.2)+
  geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), alpha = 0.2, size = 0.5, colour = "#008CD6", fill = "#008CD6") + 
  labs(y = expression(log[10](bs)), x = expression(log[10](mean))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank())
p4

pic1 <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, widths=c(4), heights=c(4), align = "v")

pic1
ggsave('fig3&S9/figS9.pdf', width = 4, height = 3.7, useDingbats = FALSE)
