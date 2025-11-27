set.seed(22)
rows <- 20000
x1 <- rbinom(rows, 1, 0.5)
x2 <- runif(rows, 0, 1)
x3 <- rnorm(rows, 0, 3)
X <- cbind(x1,x2,x3)
t_beta  <- c(0.5, 2, -1)
t_sigma <- 1

#' ######################
#' #                    #
#' # Extreme value case #
#' #                    #
#' ######################
#'
t_lambda <- 1
error <- rglg(rows, 0, 1, t_lambda)
y1 <- error
y1 <- X %*%t_beta + t_sigma*error
data.example <- data.frame(y1,X)
#library(microbenchmark)

#microbenchmark(.subset2(data.example,1)[1000:2000], data.example$y1[1000:2000], .subset(data.example,1)$y1[1000:2000], times=500)
#Unit: microseconds
#expr                                     min    lq     mean median     uq    max  neval
#.subset2(data.example, 1)[1000:2000]   2.001 4.800 5.387644  5.201 5.7010 47.200   500  ***
#data.example$y1[1000:2000]             2.600 5.701 6.358200  6.201 6.7005 18.301   500
#.subset(data.example, 1)$y1[1000:2000] 2.500 5.301 5.878398  5.701 6.2010 19.101   500

#microbenchmark(.subset2(data.example,2)[1000:2000], data.example$x1[1000:2000], .subset(data.example,2)$x1[1000:2000], times=500)
#Unit: microseconds
#expr                                    min     lq     mean  median     uq    max  neval
#.subset2(data.example, 2)[1000:2000]   1.901 2.3005 2.816772 2.4020 2.8000 26.101   500  ***
#data.example$x1[1000:2000]             2.501 2.8010 3.930204 3.0015 3.5005 91.901   500
#.subset(data.example, 2)$x1[1000:2000] 2.301 2.6020 3.432598 2.8020 3.2010 36.602   500

#microbenchmark(.subset2(data.example,3)[1000:2000], data.example$x2[1000:2000], .subset(data.example,3)$x2[1000:2000], times=500)
#Unit: microseconds
#expr                                   min    lq     mean median    uq    max  neval
#.subset2(data.example, 3)[1000:2000]   1.9 2.102 2.383780  2.300 2.402  6.901   500 ***
#data.example$x2[1000:2000]             2.4 2.800 3.037354  2.901 3.101 13.000   500
#.subset(data.example, 3)$x2[1000:2000] 2.2 2.500 2.734598  2.601 2.801 11.101   500

#microbenchmark(.subset2(data.example,4)[1000:2000], data.example$x3[1000:2000], .subset(data.example,4)$x3[1000:2000], times=500)
#Unit: microseconds
#expr                                     min    lq     mean median    uq    max  neval
#.subset2(data.example, 4)[1000:2000]   1.901 2.201 2.428998  2.302 2.501  6.001   500 ***
#data.example$x3[1000:2000]             2.400 2.800 3.190628  2.901 3.201 87.801   500
#.subset(data.example, 4)$x3[1000:2000] 2.201 2.501 2.751608  2.701 2.900  7.401   500





