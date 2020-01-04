"Create a dataset w/ factors, characters & numeric features - for investigation!"

setwd("C:/Users/kuche_000/Desktop/MS-Thesis")
data("iris")
iris$factor                              <-  sample(c("good", "ok", "legend"), 
                                                    size = nrow(iris), 
                                                    replace = TRUE, 
                                                    prob = c(0.35, 0.4, 0.25))
iris$factor[iris$Species == "virginica"] <- "legend"
iris$factor                              <-  as.factor(iris$factor)
iris$char                                <- sample(c("lol", "not_lol"),
                                                   size = nrow(iris),
                                                   replace = T,
                                                   prob = c(0.6, 0.4))
iris$char                                <- as.character(iris$char)

write.csv2(iris, file = "./data/external/example_data/iris_example.csv",
           row.names = F)