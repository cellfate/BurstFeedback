library(ggplot2)
library(ggpubr)

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/a_1.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(a = rep(c(1:10/10), each = 11), prob = data)
p1 <- ggplot(data, aes(x = rep(c(0:10), time = 10), y = prob, group = a, colour = a)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/a_2.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(a = rep(c(1:10/10), each = 21), prob = data)
p2 <- ggplot(data, aes(x = rep(c(0:20), time = 10), y = prob, group = a, colour = a)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/a_3.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(a = rep(c(1:10/10), each = 21), prob = data)
p3 <- ggplot(data, aes(x = rep(c(0:20), time = 10), y = prob, group = a, colour = a)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/a_4.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(a = rep(c(1:10*3), each = 301), prob = data)
p4 <- ggplot(data, aes(x = rep(c(0:300), time = 10), y = prob, group = a, colour = a)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/a_5.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(a = rep(c(1:10*3), each = 51), prob = data)
p5 <- ggplot(data, aes(x = rep(c(0:50), time = 10), y = prob, group = a, colour = a)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/a_6.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(a = rep(c(1:10*3), each = 301), prob = data)
p6 <- ggplot(data, aes(x = rep(c(0:300), time = 10), y = prob, group = a, colour = a)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

###################################################################################
data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/b_1.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(b = rep(c(1:10*2), each = 11), prob = data)
p7 <- ggplot(data, aes(x = rep(c(0:10), time = 10), y = prob, group = b, colour = b)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/b_2.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(b = rep(c(1:10*2), each = 11), prob = data)
p8 <- ggplot(data, aes(x = rep(c(0:10), time = 10), y = prob, group = b, colour = b)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/b_3.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(b = rep(c(1:10*2), each = 11), prob = data)
p9 <- ggplot(data, aes(x = rep(c(0:10), time = 10), y = prob, group = b, colour = b)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/b_4.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(b = rep(c(1:10*2), each = 101), prob = data)
p10 <- ggplot(data, aes(x = rep(c(0:100), time = 10), y = prob, group = b, colour = b)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  ylim(0,0.2) + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/b_5.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(b = rep(c(1:10*2), each = 51), prob = data)
p11 <- ggplot(data, aes(x = rep(c(0:50), time = 10), y = prob, group = b, colour = b)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/b_6.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(b = rep(c(1:10*2), each = 151), prob = data)
p12 <- ggplot(data, aes(x = rep(c(0:150), time = 10), y = prob, group = b, colour = b)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p12

pic1 <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
                  ncol=3,nrow=4,widths=c(3),heights=c(2),align = "v")
pic1

ggsave('figS2-S3_paramSensitive/figS2_paramSensitive.pdf', width = 7.3, height = 6)

#############################################################################
data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/h_1.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(h = rep(c(1:10), each = 11), prob = data)
p13 <- ggplot(data, aes(x = rep(c(0:10), time = 10), y = prob, group = h, colour = h)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/h_2.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(h = rep(c(-10:-1), each = 11), prob = data)
p14 <- ggplot(data, aes(x = rep(c(0:10), time = 10), y = prob, group = h, colour = h)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/h_3.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(h = rep(c(1:10), each = 31), prob = data)
p15 <- ggplot(data, aes(x = rep(c(0:30), time = 10), y = prob, group = h, colour = h)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/h_4.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(h = rep(c(-10:-1), each = 101), prob = data)
p16 <- ggplot(data, aes(x = rep(c(0:100), time = 10), y = prob, group = h, colour = h)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

#############################################################################
data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/k_1.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(k = rep(c(log10(1), log10(5), log10(10), log10(50), log10(100)), each = 11), prob = data)
p17 <- ggplot(data, aes(x = rep(c(0:10), time = 5), y = prob, group = k, colour = k)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/k_2.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(k = rep(c(log10(1), log10(5), log10(10), log10(50), log10(100)), each = 11), prob = data)
p18 <- ggplot(data, aes(x = rep(c(0:10), time = 5), y = prob, group = k, colour = k)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/k_3.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(k = rep(c(log10(1), log10(5), log10(10), log10(50), log10(100)), each = 81), prob = data)
p19 <- ggplot(data, aes(x = rep(c(0:80), time = 5), y = prob, group = k, colour = k)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  ylim(0,0.1) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

data <- read.table('../code_Inference/Hierarchical_model/misc/data_paramchange/k_4.csv', sep = ',')
data <- t(data)
data <- as.vector(data)
data <- data.frame(k = rep(c(log10(1), log10(5), log10(10), log10(50), log10(100)), each = 71), prob = data)
p20 <- ggplot(data, aes(x = rep(c(0:70), time = 5), y = prob, group = k, colour = k)) + 
  geom_line(size = 0.2) +
  scale_color_gradient(low = "#62C4FF",high = "#1F5879") +
  labs(x = 'mRNA expression', y = 'Prob.') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())

pic2 <- ggarrange(p14,p13,p16,p15,p17,p18,p19,p20,
                  ncol=2,nrow=4,widths=c(3),heights=c(2),align = "v")
pic2

ggsave('figS2-S3_paramSensitive/figS3_paramSensitive.pdf', width = 4.9, height = 6)

