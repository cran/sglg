# Normalizing function LL
c_l <- function(lambd){
       if (abs(lambd) < 0.15){
           if (lambd > 0) lambd <- 0.15
           else lambd <- -0.15
    }
  invlambdos <- 1/lambd^2
  c <- abs(lambd)/gamma(invlambdos)
  output <- c * (invlambdos^invlambdos)
  #lc <- log(abs(lambd)) - lgamma(invlambdos)
  #output <- lc + invlambdos*log(invlambdos)
  return(output)
}

# Some special functions
u_lambda <- function(lambd) {
  invlamb <- 1/lambd^2
  output <- (1/lambd) * (digamma(1 + invlamb) - log(invlamb))
  return(output)
}

v_lambda <- function(lambd) {
  invlamb <- 1/lambd^2
  output <- invlamb * trigamma(1 + invlamb) + u_lambda(lambd)^2
  return(output)
}

# K_1_lambda and K_2_lambda functions

K_1 <- function(lambd) {
  invlamb2 <- 1/lambd^2
  part1 <- 4 * (1 + digamma(1 + invlamb2) - digamma(invlamb2) - invlamb2 *
                  trigamma(invlamb2))
  part2 <- trigamma(1 + invlamb2) + (digamma(invlamb2 + 1) - log(invlamb2))^2
  output <- 1 - invlamb2 * (part1 - part2)
  return(output)
}

K_2 <- function(lambd) {
  invlamb <- 1/lambd
  invlamb2 <- invlamb^2
  output <- invlamb * ((digamma(1 + invlamb2) - digamma(invlamb2)) -
                         trigamma(1 + invlamb2) - (digamma(1 + invlamb2) - log(invlamb2))^2)
  return(output)
}

## Defining the components of the FIM

I_22 <- function(n,sigm,lambd) {
  output <- (n/(sigm^2)) * (1 + v_lambda(lambd))
  return(output)
}

I_23 <- function(n,sigm,lambd) {
  output <- (n/(sigm * lambd^2)) * K_2(lambd)
  return(output)
}

#

gfit <- function(resid, lambd) {
  Fs <- pglg(resid, shape = lambd)
  equantil <- qnorm(Fs)
  diff <- qqnorm(equantil, plot.it = FALSE)
  output <- mean(abs(diff$x - diff$y))
  return(output)
}

interval_median <- function(size,sample,alpha){
  output <- sort(sample)[qbinom(c(alpha/2,1-alpha/2),size,0.5)]
  return(output)
}
