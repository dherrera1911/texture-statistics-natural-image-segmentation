library(dplyr)
library(tidyr)
library(ggplot2)

# gaussian function
gauss <- function(x, mean, sigma) {
  y <- (1/(sqrt(2*pi*sigma^2))) * exp(-(x-mean)^2/(2*sigma^2))
  return(y)
} 

# model 1 good separation
same1 <- gauss(seq(-2.5, 2.5, 0.01), mean = -0.7, sigma = 0.54)
different1 <- gauss(seq(-2.5, 2.5, 0.01), mean = 0.7, sigma = 0.54)
df1 <- data.frame(same = same1, different = different1,
                  predictor = seq(0, 5, 0.01)) %>%
  pivot_longer(., c("different", "same"), names_to = "Stats", values_to = "value")


# model 2 better separation
same2 <- gauss(seq(-2.5, 2.5, 0.01), mean = -1.05, sigma = 0.54)
different2 <- gauss(seq(-2.5, 2.5, 0.01), mean = 1.05, sigma = 0.54)
df2 <- data.frame(same = same2, different = different2,
                  predictor = seq(0, 5, 0.01)) %>%
  pivot_longer(., c("different", "same"), names_to = "Stats", values_to = "value")

# model 3 bad separation
same3 <- gauss(seq(-2.5, 2.5, 0.01), mean = -0.1, sigma = 0.54)
different3 <- gauss(seq(-2.5, 2.5, 0.01), mean = 0.1, sigma = 0.54)
df3 <- data.frame(same = same3, different = different3,
                  predictor = seq(0, 5, 0.01)) %>%
  pivot_longer(., c("different", "same"), names_to = "Stats", values_to = "value")


goodSeparation <- ggplot(data = df1, aes(x = predictor, y = value,
                                         group = Stats, fill = Stats)) +
  geom_area(position = "identity", alpha = 0.7) +
  theme_bw() +
  scale_fill_manual(name = "Segment",
                     labels = c("different", "same"),
                       values = c("#b33018", "#14b74b"),
                       guide = "none") +
  theme(axis.text.x = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  xlab("Difference") +
  ylab("Density")

betterSeparation <- ggplot(data = df2, aes(x = predictor, y = value,
                                         group = Stats, fill = Stats)) +
  geom_area(position = "identity", alpha = 0.7) +
  theme_bw() +
  scale_fill_manual(name = "Statistic groups",
                     labels = c("different", "same"),
                       values = c("#b33018", "#14b74b"),
                     guide = "none") +
  theme(axis.text.x = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  xlab("Difference") +
  ylab("Density")

badSeparation <- ggplot(data = df3, aes(x = predictor, y = value,
                                         group = Stats, fill = Stats)) +
  geom_area(position = "identity", alpha = 0.7) +
  theme_bw() +
  scale_fill_manual(name = "Statistic groups",
                     labels = c("different", "same"),
                       values = c("#b33018", "#14b74b"),
                     guide = "none") +
  theme(axis.text.x = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  xlab("Difference") +
  ylab("Density")


ggsave("goodSeparation.png", goodSeparation,  height = 2.2, width = 4, units = "cm")
ggsave("betterSeparation.png", betterSeparation,  height = 2.2, width = 4, units = "cm")
ggsave("badSeparation.png", badSeparation,  height = 2.2, width = 4, units = "cm")


# model peformance plots

model = factor(c("FAS", "HOS", "FAS + HOS"), levels = c("FAS", "HOS", "FAS + HOS"))
performance1 <- c(0.3, 0.7, 0.3)
performance2 <- c(0.3, 0.3, 0.3)
performance3 <- c(0.3, 0.3, 0.1)

dfP1 <- data.frame(model = model, performance = performance1)
dfP2 <- data.frame(model = model, performance = performance2)
dfP3 <- data.frame(model = model, performance = performance3)

performance1 <- ggplot(dfP1, aes(x = model, y = performance, color = model)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(name = "Statistic groups",
                     labels = c("different", "same"),
                       values = c("#2691d4", "#b400e5", "#ba2229"),
                     guide = "none") +
  theme(panel.grid.major = element_blank(), axis.text.y = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 15, size = 8, hjust = 0.8, vjust=0.9,
                                   face = "bold")) +
  scale_x_discrete(labels = c("Spec", "HOS", "Spec+HOS")) +
  ylab("Error") +
  ylim(0, 0.9)

performance2 <- ggplot(dfP2, aes(x = model, y = performance, color = model)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(name = "Statistic groups",
                     labels = c("different", "same"),
                       values = c("#2691d4", "#b400e5", "#ba2229"),
                     guide = "none") +
  theme(panel.grid.major = element_blank(), axis.text.y = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 15, size = 8, hjust = 0.8, vjust=0.9,
                                   face = "bold")) +
  scale_x_discrete(labels = c("Spec", "HOS", "Spec+HOS")) +
  ylab("Error") +
  ylim(0, 0.9)

performance3 <- ggplot(dfP3, aes(x = model, y = performance, color = model)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(name = "Statistic groups",
                     labels = c("different", "same"),
                       values = c("#2691d4", "#b400e5", "#ba2229"),
                     guide = "none") +
  theme(panel.grid.major = element_blank(), axis.text.y = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 15, size = 8, hjust = 0.8, vjust=0.9,
                                   face = "bold")) +
  scale_x_discrete(labels = c("Spec", "HOS", "Spec+HOS")) +
  ylab("Error") +
  ylim(0, 0.9)


ggsave("scenario1.png", performance1,  height = 2.4, width = 4.3, units = "cm")
ggsave("scenario2.png", performance2,  height = 2.4, width = 4.3, units = "cm")
ggsave("scenario3.png", performance3,  height = 2.4, width = 4.3, units = "cm")


