# import R package
library("ggplot2")
library("reshape2")
library("ggsci")
library("ggpubr")

# read data file 
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



#### Gene fitting figure ####
# load gene expression data
data.all.T <- data.frame(t(data.all))
color.histogram <- "#D6D8DB"
color.PF <- "#E84622"
color.OnOff <- "#00837E"
# Mbnl2
prob.Mbnl2.PF <- read.table("../code_Inference/Hierarchical_model/results/results_example/prob_PF_Mbnl2.csv", sep = ',')
prob.Mbnl2.PF <- prob.Mbnl2.PF * length(na.omit(data.all.T$Mbnl2))

prob.Mbnl2.OnOff <- read.table("../code_Inference/OnOff_model/results/results_example/prob_OnOFF_Mbnl2.csv", sep = ',')
prob.Mbnl2.OnOff <- data.frame(prob.Mbnl2.OnOff * length(na.omit(data.all.T$Mbnl2)))

prob.Mbnl2 <- data.frame(cbind(1:length(prob.Mbnl2.PF)-1, t(prob.Mbnl2.OnOff[1,]), t(prob.Mbnl2.PF[1,])))
colnames(prob.Mbnl2) <- c("x_axis", "prob.Mbnl2.OnOff", "prob.Mbnl2.PF")

p_A1 <- ggplot(data.all.T, aes(x = Mbnl2)) + 
  geom_histogram(aes(x = Mbnl2), bins = length(prob.Mbnl2.PF), na.rm = TRUE, fill = color.histogram) +
  geom_line(data = prob.Mbnl2, aes(x = x_axis, y = prob.Mbnl2.OnOff), colour = color.OnOff, size = 0.5, linetype = 1) +
  geom_line(data = prob.Mbnl2, aes(x = x_axis, y = prob.Mbnl2.PF), colour = color.PF, size = 0.5, linetype = 1) +
  labs(x = "Mbnl2", y = "Cell number") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_A1

# Prr13
prob.Prr13.PF <- read.table("../code_Inference/Hierarchical_model/results/results_example/prob_PF_Prr13.csv", sep = ',')
prob.Prr13.PF <- prob.Prr13.PF * length(na.omit(data.all.T$Prr13))

prob.Prr13.OnOff <- read.table("../code_Inference/OnOff_model/results/results_example/prob_OnOFF_Prr13.csv", sep = ',')
prob.Prr13.OnOff <- data.frame(prob.Prr13.OnOff * length(na.omit(data.all.T$Prr13)))

prob.Prr13 <- data.frame(cbind(1:length(prob.Prr13.PF)-1, t(prob.Prr13.OnOff[1,]), t(prob.Prr13.PF[1,])))
colnames(prob.Prr13) <- c("x_axis", "prob.Prr13.OnOff", "prob.Prr13.PF")

p_A2 <- ggplot(data.all.T, aes(x = Prr13)) + 
  geom_histogram(aes(x = Prr13), bins = length(prob.Prr13.PF), na.rm = TRUE, fill = color.histogram) +
  geom_line(data = prob.Prr13, aes(x = x_axis, y = prob.Prr13.OnOff), colour = color.OnOff, size = 0.5, linetype = 1) +
  geom_line(data = prob.Prr13, aes(x = x_axis, y = prob.Prr13.PF), colour = color.PF, size = 0.5, linetype = 1) +
  labs(x = "Prr13", y = "Cell number") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_A2

# Ralb
prob.Ralb.PF <- read.table("../code_Inference/Hierarchical_model/results/results_example/prob_PF_Ralb.csv", sep = ',')
prob.Ralb.PF <- prob.Ralb.PF * length(na.omit(data.all.T$Ralb))

prob.Ralb.OnOff <- read.table("../code_Inference/OnOff_model/results/results_example/prob_OnOFF_Ralb.csv", sep = ',')
prob.Ralb.OnOff <- data.frame(prob.Ralb.OnOff * length(na.omit(data.all.T$Ralb)))

prob.Ralb <- data.frame(cbind(1:length(prob.Ralb.PF)-1, t(prob.Ralb.OnOff[1,]), t(prob.Ralb.PF[1,])))
colnames(prob.Ralb) <- c("x_axis", "prob.Ralb.OnOff", "prob.Ralb.PF")

p_A3 <- ggplot(data.all.T, aes(x = Ralb)) + 
  geom_histogram(aes(x = Ralb), bins = length(prob.Ralb.PF), na.rm = TRUE, fill = color.histogram) +
  geom_line(data = prob.Ralb, aes(x = x_axis, y = prob.Ralb.OnOff), colour = color.OnOff, size = 0.5, linetype = 1) +
  geom_line(data = prob.Ralb, aes(x = x_axis, y = prob.Ralb.PF), colour = color.PF, size = 0.5, linetype = 1) +
  labs(x = "Ralb", y = "Cell number") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_A3

# Plod1
prob.Plod1.PF <- read.table("../code_Inference/Hierarchical_model/results/results_example/prob_PF_Plod1.csv", sep = ',')
prob.Plod1.PF <- prob.Plod1.PF * length(na.omit(data.all.T$Plod1))

prob.Plod1.OnOff <- read.table("../code_Inference/OnOff_model/results/results_example/prob_OnOFF_Plod1.csv", sep = ',')
prob.Plod1.OnOff <- data.frame(prob.Plod1.OnOff * length(na.omit(data.all.T$Plod1)))

prob.Plod1 <- data.frame(cbind(1:length(prob.Plod1.PF)-1, t(prob.Plod1.OnOff[1,]), t(prob.Plod1.PF[1,])))
colnames(prob.Plod1) <- c("x_axis", "prob.Plod1.OnOff", "prob.Plod1.PF")

p_A4 <- ggplot(data.all.T, aes(x = Plod1)) + 
  geom_histogram(aes(x = Plod1), bins = length(prob.Plod1.PF), na.rm = TRUE, fill = color.histogram) +
  geom_line(data = prob.Plod1, aes(x = x_axis, y = prob.Plod1.OnOff), colour = color.OnOff, size = 0.5, linetype = 1) +
  geom_line(data = prob.Plod1, aes(x = x_axis, y = prob.Plod1.PF), colour = color.PF, size = 0.5, linetype = 1) +
  labs(x = "Plod1", y = "Cell number") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_A4

# bf vs bs vs mean
result.all.good = result.all[result.all$a > 0 & result.all$goodFit == 1,]
p_A <- ggplot(result.all.good, aes(x = log10(a), y = log10(b), colour = log10(mean))) + 
  geom_point(size = 1)+
  scale_colour_gradient(low = "#9CECFB", high = "#0052D4") +
  labs(y = expression(log[10](bs)), x = expression(log[10](bf))) +
  annotate("text", x = log10(result.all$a[result.all$genename == "Mbnl2"]),
           y = log10(result.all$b[result.all$genename == "Mbnl2"]), label = expression(A[1]), size = 3, colour = color.PF, face = "bold") + 
  annotate("text", x = log10(result.all$a[result.all$genename == "Prr13"]),
           y = log10(result.all$b[result.all$genename == "Prr13"]), label = expression(A[2]), size = 3, colour = color.PF, face = "bold") +
  annotate("text", x = log10(result.all$a[result.all$genename == "Ralb"]),
           y = log10(result.all$b[result.all$genename == "Ralb"]), label = expression(A[3]), size = 3, colour = color.PF, face = "bold") +
  annotate("text", x = log10(result.all$a[result.all$genename == "Plod1"]),
           y = log10(result.all$b[result.all$genename == "Plod1"]), label = expression(A[4]), size = 3, colour = color.PF, face = "bold") +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank())
p_A

pic1 <- ggarrange(p_A1, p_A2, p_A3, p_A4, ncol = 2, nrow = 2, widths=c(4), heights=c(4), align = "v")
pic1

pic2 <- ggarrange(p_A, pic1, ncol = 2, nrow = 1, widths=c(4,5), heights=c(4), align = "v")
pic2

ggsave('fig2&S7/fig2A.pdf', width = 6.6, height = 2.2, useDingbats = FALSE)


## PF model compare with on-off model
result.all.onoff <- read.table("../code_Inference/OnOff_model/results_onoff_MEF.csv", header = FALSE, sep = ',')
result.all.onoff <- cbind(genename, result.all.onoff)
colnames(result.all.onoff) <- c('genename', 'kon', 'koff', 'mu', 'fval', 'bf', 'bs', 'mean', 'noise')

# noise
noise.compare <- data.frame(cbind(result.all$noise, result.all.onoff$noise))
noise.compare$Hsign <- result.all$Hsign
colnames(noise.compare) <- c('PF', 'onoff', "Hsign")
noise.compare <- noise.compare[result.all$a > 0 & result.all$goodFit == 1,]

color.feedback <- c("#3AA438", "#547CBE", "#E71419")
p_B <- ggplot(noise.compare, aes(x = log10(PF),y = log10(onoff), colour = Hsign))+
  geom_point(size = 0.05) +
  geom_abline(slope = 1, color = '#717070',linetype = "dashed") +
  stat_cor(method = "pearson", size = 2) + 
  labs(x = "Hierarchical model", y = "Telegraph model", title = expression(log[10](noise))) + 
  theme_bw() +
  scale_color_manual(values = color.feedback) +
  xlim(-1,1) +
  ylim(-1,1) +
  theme(title = element_text(colour = 'black', size = 6),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()) 
p_B

# bf
bf.compare <- data.frame(cbind(result.all$a, result.all.onoff$bf))
bf.compare$Hsign <- result.all$Hsign
colnames(bf.compare) <- c('PF', 'onoff', "Hsign")
bf.compare <- bf.compare[result.all$a > 0 & result.all$goodFit == 1,]

p_C <- ggplot(bf.compare, aes(x = log10(PF),y = log10(onoff), colour = Hsign))+
  geom_point(size = 0.05) +
  geom_abline(slope = 1, color = '#717070',linetype = "dashed") +
  stat_cor(method = "pearson", size = 2) + 
  labs(x = "Hierarchical model", y = "Telegraph model", title = expression(log[10](bf))) + 
  theme_bw() +
  scale_color_manual(values = color.feedback) +
  xlim(-0.8,1.6) +
  ylim(-1.6,0.8) +
  theme(title = element_text(colour = 'black', size = 6),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()) 
p_C

# bs
bs.compare <- data.frame(cbind(result.all$b, result.all.onoff$bs))
bs.compare$Hsign <- result.all$Hsign
colnames(bs.compare) <- c('PF', 'onoff', "Hsign")
bs.compare <- bs.compare[result.all$a > 0 & result.all$goodFit == 1,]

p_D <- ggplot(bs.compare, aes(x = log10(PF),y = log10(onoff), colour = Hsign))+
  geom_point(size=0.05) +
  geom_abline(slope = 1, color = '#717070',linetype = "dashed") +
  stat_cor(method = "pearson", size = 2) + 
  labs(x = "Hierarchical model", y = "Telegraph model", title = expression(log[10](bs))) + 
  theme_bw() +
  scale_color_manual(values = color.feedback) +
  theme(title = element_text(colour = 'black', size = 6),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank())  
p_D

pic3 <- ggarrange(p_B, p_C, p_D, ncol = 3, nrow = 1, widths=c(4), heights=c(4), align = "v")
pic3

ggsave('fig2&S7/fig2B-D.pdf', width = 6.3, height = 2.2, useDingbats = FALSE)

