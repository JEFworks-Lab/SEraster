library(ggplot2)
data(iris)
head(iris)
ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width, col=Species)) +
  coord_fixed() +
  geom_point(size = 2) +
  theme_bw() + scale_color_manual(values = rainbow(length(unique(iris$Species))))
