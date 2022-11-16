library(tidyverse)

ProbY <- function(y,C) {(C^y * exp(-C))/factorial(y)} 

coverage = c(1:30)

data <- data.frame(
  coverage = coverage,
  prob = unlist(lapply(coverage,function(x){sum(ProbY(0:5,x))}))
  
)

ggplot(data, aes(x=coverage, y=prob))+
  geom_point(colour="blue")+
  geom_smooth()