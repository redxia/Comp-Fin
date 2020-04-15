# Starting the initial value or global functions 
library(qrandom) # library to use the qrandomnorm- Polar-Marsaglia
set.seed(0) # setting the seed to stabilize the randomness
n <- 1000 # this is a global n that is used in the function. n = Number of random draws from distribution


## Q1 ##
# function that only takes a (the covariance) and outputs correlation
Q1 <- function(a){
  z1 <- suppressMessages(qrandomnorm(n, method = "polar"))
  z2 <- suppressMessages(qrandomnorm(n, method = "polar"))
  # Mean zero for normal
  mu <- rep(0,2)
  # The covariance matrix
  covMat <- matrix(c(3,a,a,5),2,2)
  # Generating the (x,y) Sequences of Bivariate normal
  x <- mu[1] + covMat[1,1] * z1
  # correlation is
  y <- mu[2] + covMat[1,2] / covMat[1,1] * z1 + covMat[2,2] * 
    sqrt(1-(covMat[1,2]/(covMat[1,1] * covMat[2,2]))^2) * z2
  rho <- (1 / (n - 1) * sum((x - sum(x) / n) * (y - sum(y) / n))) /
    (sqrt(1 / (n - 1) * sum((x - mean(x))^2)) * sqrt(1 / (n - 1) * sum((y - mean(y))^2)))
  return(rho)
}
a <- -0.7
rho = Q1(a)
cat("rho is: ",rho,'\n')

####

## Q2 ##
# The function should only take rho and outputs E, which is the expected value of
# max(0,(X^3 + sin(Y) + X^2 * Y))
Q2 <- function(rho) {
  # Draws the normal N(0,1) distribution through the polar methods
  z1 <- qrandomnorm(n, method = "polar")
  z2 <- qrandomnorm(n, method = "polar")
  
  # Generates the X and Y distributions
  x <- z1
  y <- rho * z1 + sqrt(1 - rho^2) * z2
  
  # Compute the expected value by taking the mean
  E <- mean(pmax(0, x^3 + sin(y) + x^2 * y))
  return(E)
}
rho <- 0.6
E = Q2(rho)
cat("The expected value of max(0,(X^3 + sin(Y) + X^2 * Y)) = ",E,'\n')

####

## Q3 ##

# 3a
t <- c(1,3,5)
n <- 5000 # using a larger n

# The function takes a time and spits out the expected values of A(t)
# A(t) is defined as E(W_t^2 + sin(W_t)), W_t is the standard wiener process
# t should be a single input
Q3_At <- function(t){
  z1 <- qrandomnorm(n, method = "polar")
  W_t <- sqrt(t) * z1
  At <- (W_t^2 + sin(W_t))
  return(mean(At))
}

Q3_Bt <- function(t) {
  z1 <- qrandomnorm(n, method = "polar")
  W_t <- sqrt(t) * z1
  Bt <- (exp(t / 2) * cos(W_t))
  return(Bt)
}
E_A1 <- Q3_At(t[1]) # When t = 1
E_A3 <- Q3_At(t[2]) # When t = 3
E_A5 <- Q3_At(t[3]) # When t = 5
E_B1 <- Q3_Bt(t[1]) # When t = 1
E_B3 <- Q3_Bt(t[2]) # When t = 3
E_B5 <- Q3_Bt(t[3]) # When t = 5

# Expected Value Matrix, storing the means. 2 columns for A(t) & B(t)
E_Mat <- matrix(0,2,length(t))
rownames(E_Mat) <- c("A(t)","B(t)")
colnames(E_Mat) <- c('t = 1','t = 3','t = 5')
E_Mat[1,1] <- mean(E_A1)
E_Mat[1,2] <- mean(E_A3)
E_Mat[1,3] <- mean(E_A5)
E_Mat[2,1] <- mean(E_B1)
E_Mat[2,2] <- mean(E_B3)
E_Mat[2,3] <- mean(E_B5)

library(knitr)
kable(E_Mat, caption = "The expected values by simulation of the functions for each time steps")

# 3b
E_Mat <- rbind(E_Mat, c(var(E_B1),var(E_B3),var(E_B5)))
rownames(E_Mat)[3] <- "Variance of B(t)"
kable(E_Mat, caption = "Adding the variance of B(t)")

# 3c
# Input is all t
Q3c <- function(t){
  z1 <- qrandomnorm(n, method = "polar")
  Wa_t <- matrix(0,nrow = n,ncol = length(t))
  Wb_t <- matrix(0,nrow = n,ncol = length(t))
  At <- matrix(0, nrow = n, ncol = length(t))
  Bt <- matrix(0, nrow = n, ncol = length(t))
  colnames(Wa_t) <- c('t = 1','t = 3','t = 5')
  colnames(Wb_t) <- c('t = 1','t = 3','t = 5')
  colnames(At) <- c('t = 1','t = 3','t = 5')
  colnames(Bt) <- c('t = 1','t = 3','t = 5')
  for(i in 1:length(t)){
    Wa_t[,i] <- sqrt(t[i]) * z1 # uses the global vairable, z1, polar method of random normal
    Wb_t[,i] <- sqrt(t[i]) * -z1
    AW_t <- cbind(Wa_t[,i]^2 + sin(Wa_t[,i]),Wb_t[,i]^2 + sin(Wb_t[i]) )
    At[,i] <- apply(AW_t,1,mean)
    BW_t <- cbind(exp(t[i] / 2) * cos(Wa_t[,i]),exp(t[i] / 2) * cos(Wb_t[,i]))
    Bt[,i] <- apply(BW_t,1,mean)
  }
  return(cbind(At,Bt)) #Returns
}

varReduce <- Q3c(t)
E_varRed <- apply(varReduce,2,mean)
Var_varRed <- apply(varReduce,2,var)
Bt_reducedVar <- rbind(E_varRed[4:length(E_varRed)],Var_varRed[4:length(Var_varRed)])
rownames(Bt_reducedVar) <- c("mean","variance")
kable(Bt_reducedVar,caption = "Antithetic Variates variance reduction of B(t)")
cat("The expected value of B(t) is ",Bt_reducedVar[1,3],'\n')
####

## Q4 ##
# 4a
Q4 <- function(r,sigma,S0,t,X){
  z1 <- qrandomnorm(n, method = "polar")
  S_T <- S0 * exp((r - sigma^2/2) * t + sigma * sqrt(t) * z1) # stock formula
  payoff <- pmax(0,S_T - X) # Call function applied
  return(exp(-r*t) * payoff) #discounting the payoffs
}
r <- .04 # interest rate
sigma <- .2 # volatility
S0 <- 88 # stock price today
t <- 5 # time to maturity
x <- 100 # strikeprice
callPRC <- Q4(r,sigma,S0,t,x)
cat("The call price by monte carlo simulation is ", mean(callPRC),'\n')

# 4b
library(qrmtools) # To use the black scholes formula
callPRC_BS <- Black_Scholes(0,S0,r,sigma,x,t, type = 'call')
cat("Black Scholes call price: ",callPRC_BS,'\n')

# 4c
Q4c <- function(r,sigma,S0) {
  z1 <- qrandomnorm(n, method = "polar")
  S_Ta <- S0 * exp((r - sigma^2/2) * t + sigma * sqrt(t) * z1) # stock formula
  S_Tb <- S0 * exp((r - sigma^2 / 2) * t + sigma * sqrt(t) * (-z1))
  payoffa <- pmax(0,S_Ta - x) #derivative payoffs
  payoffb <- pmax(0,S_Tb - x)
  PV_payoffs <- exp(-r * t) * cbind(payoffa,payoffb)# The present value of payoffs
  Cal <- apply(PV_payoffs,1,mean)
  return(Cal)
}

CalPRC <- Q4c(r,sigma,S0)
cat("The estimated call price of red", mean(CalPRC),'\n')
cat("The reduced variance call price ", var(CalPRC),'\n')
cat("The original variance of call price ", var(callPRC),'\n')
####

## Q5 ##

# 5a
# This function builds the randome walk
Q5 <- function(r,sigma,S0,t){
  z1 <- rnorm(n)
  S_T <- S0 * exp((r - sigma^2/2) * t + sigma * sqrt(t) * z1) # stock formula
  return(mean(S_T))
}
n <- 1000
sigma <- .18
randomWalk <- numeric(10)
for(t in 1:10){
  randomWalk[t] <- Q5(r,sigma,S0,t)
}
plot(1:10, randomWalk, xlab = "Time", ylab = "Stock Price Simulated", main = "Stock Price", type = 'l')

# 5b
# We generate 6 random walk paths
randomWalkPaths <- matrix(0, ncol = 6, nrow = n)
for(i in 1:6) {
  timeStep <- 1
  for(t in seq(1,10,(10-1)/(n - 1))) {
    randomWalkPaths[timeStep,i] <- Q5(r,sigma,S0,t)
    timeStep <- timeStep + 1
  }
}
cat('\n')

# 5c
time <- matrix(rep(seq(1,10,(10 - 1) / (n - 1) ),6),n,6)
matplot(time,randomWalkPaths, type = 'l', xlab = 'Time', ylab = 'Stock Price')
lines(randomWalk, lwd = 3, col = 'black')

# 5d
sigma <- .35
# The code is a repeat from a. We just change the sigma
randomWalk <- numeric(10)
for(t in 1:10){
  randomWalk[t] <- Q5(r,sigma,S0,t)
}
# The code is a repeat from b. We just change the sigma
randomWalkPaths <- matrix(0, ncol = 6, nrow = n)
for(i in 1:6) {
  timeStep <- 1
  for(t in seq(1,10,(10-1)/(n - 1))) {
    randomWalkPaths[timeStep,i] <- Q5(r,sigma,S0,t)
    timeStep <- timeStep + 1
  }
}
time <- matrix(rep(seq(1,10,(10 - 1) / (n - 1) ),6),n,6)
matplot(time,randomWalkPaths, type = 'l', xlab = 'Time', ylab = 'Stock Price')
lines(randomWalk, lwd = 3, col = 'black')
####

## Q6 ##

# 6a
Q6a <- function(x){
  pi <- 4 * sqrt(1 - x^2)
  return(pi)
}

piVal <- 0
n <- 10000 # We use higher n to get a more accurate result
for(i in seq(0,1,1/n)){
  piVal <- piVal + Q6a(i) * (1 / n)
}
cat("The value of the integral is ",piVal, '\n')

# 6b
piValMC <- 4 * sqrt(1 - runif(n)^2)
cat("The Monte Carlo Estimated Pi: ", mean(piValMC))

# 6c
# Wei Cai model for importance sampling
x <- runif(n)
a <- 0.74  # the alpha parameter
g_x <- 4 * sqrt(1 - x^2)
t_x <- (1-0.74 * x^2) / (1 - a / 3)

ISmethod <- g_x / t_x
cat("Importance Sampling (IS) Method Estimated Value ", mean(ISmethod),'\n')
cat("Comparing the variances.\n","Importance sampling standard deviation: ",sd(ISmethod),'\n')
cat("MonteCarlo Estimated standard deviation: ",sd(piValMC),'\n')
cat("The absolute difference between actual pi and IS method ", abs(mean(ISmethod)-pi))
cat("The absolute difference between actual pi and Monte Carlo method ", abs(mean(piValMC)-pi))
####

