# install and load packages
pack <- c("doParallel", "foreach", "KernSmooth")

lapply(pack, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})

lapply(pack, library, character.only = TRUE)

# class
np <- function(x, y, xmin = NULL, xmax = NULL, x.unscaled = NULL) {
  # order observations
  y <- y[order(x)]
  x <- x[order(x)]

  # min and max
  if (is.null(xmin)) {
    xmin <- min(x)
  }

  if (is.null(xmax)) {
    xmax <- max(x)
  }

  # scale data
  x.scale <- npscale(x, xmin, xmax)

  # class np.object
  np <- list(x = x.scale$x, y = y, xmin = xmin, xmax = xmax,
             x.unscaled = x.scale$x.unscaled, x.grid = NULL, y.hat = NULL)
  class(np) <- "np.object"
  return(np)

}

npscale = function(x, xmin, xmax) {
  # check if data already scaled
  if (min(x) != xmin | max(x) != xmax) {
    if (min(x) == 0 & max(x) == 1) {
      # data already scaled, compute unscaled x
      if (is.null(x.unscaled)) {
        x.unscaled <- (x * (xmax - xmin)) + xmin
        }
    } else {
      stop("Error: Something seems to be wrong with the data (Scaling)")
    }
  } else {
    # data not scaled yet, scale it
    x.unscaled <- x
    x <- (x - xmin) / (xmax - xmin)
  }

  return(list(x = x, x.unscaled = x.unscaled))
}

plot.np.object <- function(x, ...) {
  plot(x$x.unscaled, x$y, 
       pch = 20, cex = 0.5, col = 'grey',
       bty  = 'l', xaxt = 'n', yaxt = 'n', col.lab = "grey",
       ...)
  box(lwd = 2, bty = 'l', col = 'grey')
  axis(side = 2, x$y, tick = FALSE, cex.lab = 1.5, pos = 1, 
       col.axis = "grey", labels = NULL, at = pretty(x$y, n = 5))

  axis(side = 1, x$x.unscaled, tick = FALSE, cex.lab = 1.5, pos = 1, 
       col.axis = "grey", labels = NULL, at = pretty(x$x.unscaled, n = 7))

  lines(x$y.hat$x.grid.scaled, x$y.hat$y.hat, col = 'blue', lwd = 2)
}

npsmoother = function(data.np, deg, h, xmin = NULL, xmax = NULL, 
                      x.grid = NULL) {

  if (class(data.np) != "np.object") {
    stop("Create np object before using npsmoother.")
  }
  
  kernelq <- function(u) {
    dnorm(u, mean = 0, sd = 1)
  }
    
  if (is.null(x.grid)) {
    x.grid <- seq(0, 1, length.out = 1000)
  }

  n.grid <- length(x.grid)
  n.val <- length(data.np$x)

  x <- data.np$x
  y <- data.np$y

  y.hat <- vector()

  # use all but one core in parallel computation
  cl <- detectCores() - 1
  cl <- makeCluster(cl)
  registerDoParallel(cl)

  y.hat <- foreach(i = 1:n.grid)%dopar%{

    # construct X matrix with polynomials
    if (deg == 0) {
      x.mat <- rep(1, n.grid)
    } else {
      x.mat <- matrix(rep(c(1, 1:deg), length(x)), ncol = deg + 1, byrow = TRUE)

      ind <- 2:(deg + 1)
      x.mat[, ind] <- (x - x.grid[i])^x.mat[, ind]
    }

    # weights
    w <- diag( kernelq((x - x.grid[i])/h)/h)

    # WLS
    beta <- solve((t(x.mat) %*% w %*% x.mat)) %*% (t(x.mat) %*% w %*% y)

    # local prediction is just constant of WLS
    beta[1]
  }

  stopCluster(cl)

  y.hat <- unlist(y.hat)

  x.grid.scaled <- (x.grid * (data.np$xmax - data.np$xmin)) + data.np$xmin
  data.np$y.hat <- list(y.hat = y.hat, x.grid = x.grid, deg = deg, h = h, 
                        x.grid.scaled = x.grid.scaled)

  return(data.np)
}
