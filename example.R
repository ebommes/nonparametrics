setwd('') # set your working directory here
source('local_polynomial.R')

n <- 500
x <- runif(n, 1.5, 20)
eps <- rnorm(n, 0)
y <- x^2 + x + (x^(1.2) * eps)

data.l <- np(x, y)
model.l <- npsmoother(data.l, 1, 0.1)

plot(model.l, xlab = 'x', ylab = 'y')
