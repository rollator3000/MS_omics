"Minimalistic Example for a better explanation of the decision tree!"
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
a1 <- fancyRpartPlot(fit)	

# Plot how the Decision Tree partionates the feature space!
a2 <- ggplot(data = example, aes(x = weight, y = height, col = Y)) + 
  geom_point(size = 5) +
  theme_bw() + 
  theme(text = element_text(size=20),
        plot.title = element_text(size=18),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  annotate("text", x=88, y=180, label= "Node1", size = 7, col = "seagreen") +
  scale_x_continuous(limits = c(63, 109), expand = c(0, 0)) +
  ggtitle("Raw Data - all observations in Node1") +
  xlab("")

a3 <- ggplot(data = example, aes(x = weight, y = height, col = Y)) + 
  geom_point(size = 5) +
  theme_bw() + 
  theme(text = element_text(size=20),
        plot.title = element_text(size=18),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = 69, col = "red", lwd = 2) +
  scale_x_continuous(limits = c(63, 109), expand = c(0, 0)) +
  annotate("text", x=88, y=180, label= "Node2", size = 7, col = "seagreen") +
  annotate("text", x=66, y=180, label= "Node3", size = 7, col = "seagreen") +
  ggtitle("Splitting Node1 into Node2 & Node3") + 
  xlab("")

a4 <- ggplot(data = example, aes(x = weight, y = height, col = Y)) + 
  geom_point(size = 5) +
  theme_bw() + 
  theme(text = element_text(size=20),
        plot.title = element_text(size=18)) +
  geom_vline(xintercept = 69, col = "red", , lwd = 2)  +
  geom_segment(aes(x = 69, xend = 109, y = 171, yend = 171), col = "red", lwd = 2) + 
  scale_x_continuous(limits = c(63, 109), expand = c(0, 0),) +
  annotate("text", x=88, y=180, label= "Node4", size = 7, col = "seagreen") +
  annotate("text", x=80, y=169, label= "Node5", size = 7, col = "seagreen") +
  annotate("text", x=66, y=180, label= "Node3", size = 7, col = "seagreen") +
  ggtitle("Splitting Node 2 into Node4 & Node5")

grid.arrange(a2, a3, a4, nrow = 3, top = textGrob("Single segmentation steps of a decision tree", 
                                                  gp = gpar(fontsize=25)))