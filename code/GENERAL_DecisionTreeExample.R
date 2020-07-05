"Minimalistic Example for a better explanation of the decision tree"
# Load the needed libraries
library(rpart)
library(rpart.plot)
library(ggplot2)
library(rattle)
library(grid)
require(gridExtra)

# ReCreate the data from table1 in the MS Thesis
example <- data.frame("weight" = c(65.4, 83.9, 67.4, 105.2, 71.5, 73.0),
                      "height" = c(187, 192, 167, 175, 173, 169),
                      "Y"      = c(1, 0, 1, 0, 0, 1))

# Recode the response to a factor
example$Y <- as.factor(example$Y)

# Fit a Tree and look at the results
fit <- rpart(Y ~ ., data = example, method = 'class', 
             control = rpart.control(minbucket = 1))
a1 <- fancyRpartPlot(fit, cex = 1.7)	

# Plot how the Decision Tree partionates the feature space!
a2 <- ggplot(data = example, aes(x = weight, y = height, col = Y)) + 
  geom_point(size = 8) +
  theme_bw() + 
  theme(text = element_text(size = 28),
        plot.title = element_text(size = 30),
        legend.position = "none") +
  annotate("text", x = 88, y = 180, label = "N1", size = 15, col = "seagreen") +
  scale_x_continuous(limits = c(63, 106), expand = c(0, 0)) +
  ggtitle("All obs. in the root node N1")

a3 <- ggplot(data = example, aes(x = weight, y = height, col = Y)) + 
  geom_point(size = 8) +
  theme_bw() + 
  theme(text = element_text(size = 28),
        plot.title = element_text(size = 30),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  geom_vline(xintercept = 69, col = "red", lwd = 1.7) +
  scale_x_continuous(limits = c(63, 106), expand = c(0, 0)) +
  annotate("text", x = 88, y = 180, label = "N2", size = 15, col = "seagreen") +
  annotate("text", x = 66, y = 180, label = "N3", size = 15, col = "seagreen") +
  ggtitle("Splitting N1 into child nodes N2 & N3") + 
  ylab("")

a4 <- ggplot(data = example, aes(x = weight, y = height, col = Y)) + 
  geom_point(size = 8) +
  theme_bw() + 
  theme(text = element_text(size = 28),
        plot.title = element_text(size = 30),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  geom_vline(xintercept = 69, col = "red", , lwd = 1.7)  +
  geom_segment(aes(x = 69, xend = 107, y = 171, yend = 171), col = "red", lwd = 1.7) + 
  scale_x_continuous(limits = c(62, 107), expand = c(0, 0),) +
  annotate("text", x = 88, y = 180, label = "N4", size = 15, col = "seagreen") +
  annotate("text", x = 80, y = 169, label = "N5", size = 15, col = "seagreen") +
  annotate("text", x = 65.5, y = 180, label = "N3", size = 15, col = "seagreen") +
  ggtitle("Splitting N2 into child nodes N4 & N5")


# Paste all plots into a single figure
grid.arrange(a2, a3, a4, nrow = 1)
