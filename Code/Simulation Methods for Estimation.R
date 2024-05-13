# Function is (4-x^2)e^(3-x^2)
# Theoretical value is 56.19551
# Theoretical variance is unable to be calculated

# Regular Monte Carlo

h <- function(x) (4-x^2)*exp(3-x^2)
rmc <- function(n){
  u <- runif(n)
  hu <- h(u)
  Ihat <- mean(hu)
  SE2_hu <- var(hu)/n
  return(list(Ihat=Ihat,SE2=SE2_hu))
}

# Importance Sampling 

f <- function(x) (2-x) * (2/3)
is <- function(n){
  u <- runif(n)
  x <- -sqrt(4-3*u) + 2
  hf <- (h(x)/f(x))
  Ihat <- mean(hf)
  SE2_hf <-var(hf)/n 
  return(list(Ihat=Ihat,SE2=SE2_hf))
}

# Stratified Sampling 

ss <- function(n,k){
  nn <- n/k
  hh <- matrix(0,k,nn)
  x <- matrix(0,k,nn)
  u <- matrix(0,k,nn)
  t <- 0
  for (i in 1:k) {
    a_j <- (1/k)*(i-1)
    b_j <- (1/k)*i
    u[i,] <- runif(nn,a_j,b_j)
  }
  for (i in 1:k) {
    hh[i,] <- h(u[i,])
  }
  for (i in 1:k) {
    t[i] <- mean(hh[i,])
  }
  Ihat <- sum(t)*(1/k)
  SE2 <- 0
  for (i in 1:k) {
    SE2 <- var(hh[i,])*(1/k^2)/nn + SE2
  }
  return(list(Ihat=Ihat,SE2=SE2))
}


# Stratified Importance Sampling 

sis <- function(n,k){
  nn <- n/k
  t <- rep(0,k)
  hh <- matrix(0,k,nn)
  u <- matrix(0,k,nn)
  ff <- matrix(0,k,nn)
  x <- matrix(0,k,nn)
  for (i in 1:k) {
    a_j <- (1/k)*(i-1)
    b_j <- (1/k)*i
    u[i,] <- runif(nn,0,1)
    x[i,] <- 2 - sqrt(b_j^2*u[i,] + 4*a_j*u[i,] - a_j^2*u[i,] - 4*b_j*u[i,] + a_j^2 + 4 - 4*a_j)
    
  }
  for (i in 1:k) {
    hh[i,] <- h(x[i,])
    ff[i,] <- f(x[i,])
  }
  for (i in 1:k) {
    t[i] <- mean(hh[i,]/ff[i,])
  }
  Ihat <- 0
  SE2 <- 0
  for (i in 1:k) {
    a_j <- (1/k)*(i-1)
    b_j <- (1/k)*i
    Ihat <- t[i]*(4*b_j/3 - 4*a_j/3 - (b_j^2 - a_j^2)/3) + Ihat
    SE2 <- var(hh[i,]/ff[i,])*(4*b_j/3 - 4*a_j/3 - (b_j^2 - a_j^2)/3)^2 / nn  + SE2
  }
  return(list(Ihat=Ihat,SE2=SE2))
}

# Stratified Sampling with optimal n_i

ss_sd <- function(n,k){
  t <- rep(0,k)
  sum_var <- 0
  n_n <- n/k
  nn <- rep(0,k)
  sd <- rep(0,k)
  u <- list()
  hh <- list()
  for (i in 1:k) {
    a_j <- (1/k)*(i-1)
    b_j <- (1/k)*i
    u[[i]] <- runif(n_n,a_j,b_j)
  }
  for (i in 1:k) {
    hh[[i]] <- h(u[[i]])
  }
  sd_total <- 0
  for (i in 1:k) {
    a_j <- (1/k)*(i-1)
    b_j <- (1/k)*i
    sd[i] <- sd(hh[[i]])
    sd_total <- sd[i]*(1/k)+sd_total
  }
  for (i in 1:k) {
    a_j <- (1/k)*(i-1)
    b_j <- (1/k)*i
    nn[i] <- n*sd[i]*(1/k) / sd_total
  }
  for (i in 1:k) {
    a_j <- (1/k)*(i-1)
    b_j <- (1/k)*i
    u[[i]] <- runif(nn[i],a_j,b_j)
  }
  for (i in 1:k) {
    hh[[i]] <- h(u[[i]]) 
  }
  for (i in 1:k) {
    t[i] <- mean(hh[[i]])
  }
  Ihat <- 0
  SE2 <- 0
  for (i in 1:k) {
    Ihat <- t[i]*(1/k) + Ihat
    SE2 <- var(hh[[i]])*(1/k^2)/nn[i] + SE2
  }
  return(list(Ihat=Ihat,SE2=SE2))
}

# Stratified Importance Sampling with optimal n_i

sis_sd <- function(n,k){
  t <- rep(0,k)
  n_n <- n/k
  nn <- rep(0,k)
  sd <- rep(0,k)
  u <- list()
  hh <- list()
  ff <- list()
  x <- list()
  for (i in 1:k) {
    a_j <- (1/k)*(i-1)
    b_j <- (1/k)*i
    u[[i]] <- runif(n_n,0,1)
    x[[i]] <- 2 - sqrt(b_j^2*u[[i]] + 4*a_j*u[[i]] - a_j^2*u[[i]] - 4*b_j*u[[i]] + a_j^2 + 4 - 4*a_j)
  }
  for (i in 1:k) {
    hh[[i]] <- h(x[[i]])
    ff[[i]] <- f(x[[i]])
  }
  sd_total <- 0
  for (i in 1:k) {
    a_j <- (1/k)*(i-1)
    b_j <- (1/k)*i
    sd[i] <- sqrt(var(hh[[i]]/ff[[i]])/n_n)
    sd_total <- sd[i]*(4*b_j/3 - 4*a_j/3 - (b_j^2 - a_j^2)/3) + sd_total
  }
  for (i in 1:k) {
    a_j <- (1/k)*(i-1)
    b_j <- (1/k)*i
    nn[i] <- n*sd[i]*(4*b_j/3 - 4*a_j/3 - (b_j^2 - a_j^2)/3) / sd_total
  }
  for (i in 1:k) {
    a_j <- (1/k)*(i-1)
    b_j <- (1/k)*i
    u[[i]] <- runif(nn[i])
    x[[i]] <- 2 - sqrt(b_j^2*u[[i]] + 4*a_j*u[[i]] - a_j^2*u[[i]] - 4*b_j*u[[i]] + a_j^2 + 4 - 4*a_j)
  }
  for (i in 1:k) {
    hh[[i]] <- h(x[[i]])
    ff[[i]] <- f(x[[i]])
  }
  for (i in 1:k) {
    t[i] <- mean(hh[[i]]/ff[[i]])
  }
  Ihat <- 0
  SE2 <- 0
  for (i in 1:k) {
    a_j <- (1/k)*(i-1)
    b_j <- (1/k)*i
    Ihat <- t[i]*(4*b_j/3 - 4*a_j/3 - (b_j^2 - a_j^2)/3) + Ihat
    SE2 <- var(hh[[i]]/ff[[i]])*(4*b_j/3 - 4*a_j/3 - (b_j^2 - a_j^2)/3)^2/nn[i] + SE2
  }
  return(list(Ihat=Ihat,SE2=SE2))
}

# n = 100 | k = 10

B <- 1000
res_mc <- rep(0,B)
res_is <- rep(0,B)
res_ss <- rep(0,B)
res_sis <- rep(0,B)
res_ss_sd <- rep(0,B)
res_sis_sd <- rep(0,B)

n <- 100
k <- 10

for (i in 1:B) {
  res_mc[i] <- rmc(n)$Ihat
  res_is[i] <- is(n)$Ihat
  res_ss[i]<- ss(n,k)$Ihat
  res_sis[i] <- sis(n,k)$Ihat
  res_ss_sd[i] <- ss_sd(n,k)$Ihat
  res_sis_sd[i] <- sis_sd(n,k)$Ihat
}

boxplot(res_mc,res_is,res_ss,res_ss_sd,res_sis,res_sis_sd,
        col=c("orange","yellow","blue","lightblue","green","darkgreen"),
        range = 0,main="Simulation Study with n = 100 | k = 10")
abline(h=56.19551,col="red")
text(x=1,y=quantile(res_mc,0.75),labels="Regular Monte Carlo",pos=3)
text(x=2,y=quantile(res_is,0.75),labels="Importance Sampling",pos=3)
text(x=3,y=quantile(res_ss,0.75),labels="Stratified Sampling",pos=3)
text(x=4,y=quantile(na.omit(res_ss_sd),0.75),labels="Stratified Sampling Opt",pos=3)
text(x=5,y=quantile(res_sis,0.75),labels="Stratified Importance Sampling",pos=3)
text(x=6.1,y=quantile(na.omit(res_sis_sd),0.75),labels="Stratified Importance Sampling Opt",pos=3)

# n = 2000 k = 20 

B <- 1000
res_mc <- rep(0,B)
res_is <- rep(0,B)
res_ss <- rep(0,B)
res_sis <- rep(0,B)
res_ss_sd <- rep(0,B)
res_sis_sd <- rep(0,B)

n <- 2000
k <- 20

for (i in 1:B) {
  res_mc[i] <- rmc(n)$Ihat
  res_is[i] <- is(n)$Ihat
  res_ss[i]<- ss(n,k)$Ihat
  res_sis[i] <- sis(n,k)$Ihat
  res_ss_sd[i] <- ss_sd(n,k)$Ihat
  res_sis_sd[i] <- sis_sd(n,k)$Ihat
}

boxplot(res_mc,res_is,res_ss,res_ss_sd,res_sis,res_sis_sd,
        col=c("orange","yellow","blue","lightblue","green","darkgreen"),
        range = 0,main="Simulation Study with n = 2000 | k = 20")
abline(h=56.19551,col="red")
text(x=1,y=quantile(res_mc,0.75),labels="Regular Monte Carlo",pos=3)
text(x=2,y=quantile(res_is,0.75),labels="Importance Sampling",pos=3)
text(x=3,y=quantile(res_ss,0.75),labels="Stratified Sampling",pos=3)
text(x=4,y=quantile(na.omit(res_ss_sd),0.75),labels="Stratified Sampling Opt",pos=3)
text(x=5,y=quantile(res_sis,0.75),labels="Stratified Importance Sampling",pos=3)
text(x=6.1,y=quantile(na.omit(res_sis_sd),0.75),labels="Stratified Importance Sampling Opt",pos=3)