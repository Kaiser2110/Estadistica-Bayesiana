library(VGAM)
library(ggplot2)
library(gridExtra)
library(mvtnorm)
library(EfficientMaxEigenpair)

load("diabetes.RData")
rownames(X.train) <- NULL
rownames(X.test) <- NULL

X <- X.train[, 2:11]
y <- X.train[, 1]
X.test <- X.test[, 2:11]
y.test <- X.test[, 1]
rm(X.train)

XtX <- (t(X)%*%X)
Xty <- t(X)%*%y

n <- 342
p <- 10

SSR <- function(beta, sigma = diag(rep(1, n))){
  return(t(y - X%*%beta)%*%sigma%*%(y - X%*%beta))
}
mod <- lm(y ~ -1 + X)
#=============================================================================================

#                                              MODELO 1

#=============================================================================================
B <- 11000

LL1 <- 0
par_M1 <- 0
M1 <- NULL

beta_b <- mod$coefficients
lambda_b <- 1
sigma2_b <- summary(mod)$sigma^2

sigma2_0 <- sigma2_b
set.seed(2110)
tictoc::tic()
for(i in 1:B){
  sig <- solve(XtX + diag(rep(lambda_b, 10)))
  beta_b <- c(rmvnorm(1, mean = sig%*%Xty, sigma = sigma2_b*sig))
  sigma2_b <- 1/rgamma(1, shape = 353/2, rate = (1/2)*(sum((y - X%*%beta_b)^2) + lambda_b*sum(beta_b^2) + sigma2_0))
  lambda_b <- rgamma(1, shape =  6, rate = (1/(2*sigma2_b))*(sum(beta_b^2)) + 2)
  
  if(i > 1000){
    LL1 <- LL1 +  (1/(B - 1000))*sum(dnorm(y, mean = X%*%beta_b, sd = sqrt(sigma2_b), log = T))
    par_M1 <- par_M1 + (1/(B - 1000))*c(beta_b, sigma2_b)
    y.est <- rnorm(n, mean = X%*%beta_b, sd = sqrt(sigma2_b))
    M1 <- rbind(M1, c(mean(y.est), sqrt(var(y.est))))
  }
}
tictoc::toc()


rm(beta_b, lambda_b, sigma2_b, sig)

#=============================================================================================

#                                              MODELO 2

#=============================================================================================
B <- 11000

LL2 <- 0
par_M2 <- 0
M2 <- NULL

beta_b <- mod$coefficients
phi2_b <- beta_b^2
sigma2_b <- summary(mod)$sigma^2

sigma2_0 <- sigma2_b

set.seed(2110)
tictoc::tic()
for(i in 1:B){
  sig <- solve(XtX + diag(sigma2_b/phi2_b))
  beta_b <- c(rmvnorm(1, mean = sig%*%Xty, sigma = sigma2_b*sig))
  phi2_b <- 1/rinv.gaussian(10, mu = 1/5*sqrt(beta_b^2), lambda = 1/25)
  sigma2_b <- 1/rgamma(1, shape = 343/2, rate = (1/2)*(sum((y - X%*%beta_b)^2) + sigma2_0))
  
  
  if(i > 1000){
    par_M2 <- par_M2 + (1/(B - 1000))*c(beta_b, sigma2_b)
    LL2 <- LL2 + (1/(B - 1000))*sum(dnorm(y, mean = X%*%beta_b, sd = sqrt(sigma2_b), log = T))
    y.est <- rnorm(n, mean = X%*%beta_b, sd = sqrt(sigma2_b))
    M2 <- rbind(M2, c(mean(y.est), sqrt(var(y.est))))
  }
}
tictoc::toc()
rm(beta_b, phi2_b, sigma2_b, sig, y.est)


#=============================================================================================

#                                              MODELO 3

#=============================================================================================
B <- 11000
L <- 20

LL3 <- 0
par_M3 <- 0
M3 <- NULL

beta_b <- mod$coefficients
ro_b <- 0.5
sigma2_b <- summary(mod)$sigma^2
sigma2_0 <- sigma2_b
t2_0 <- 50
ac <- 0
C_ro <- (abs(outer((1:n),(1:n) ,"-")))

set.seed(2110)
tictoc::tic()
for(i in 1:B){
  IC_ro <- tridiag(rep(-ro_b, n - 1), rep(-ro_b, n - 1), c(1, rep(1 + ro_b^2, n-2), 1))
  sig <- solve((1/(1 - ro_b^2))*t(X)%*%IC_ro%*%X + diag(rep(sigma2_b/t2_0, 10)))
  beta_b <- c(rmvnorm(1, mean = sig%*%Xty, sigma = sigma2_b*sig))
  res <- y - X%*%beta_b
  
  sigma2_b <- 1/rgamma(1, shape = (1 + n)/2, rate = (1/(2*(1 - ro_b^2)))*(t(res)%*%IC_ro%*%(res) + sigma2_0))
  
  # Actualizar ro
  ro.ham <- ro_b
  
  phi.p <- rnorm(1, mean = 0, sd = 1/3)
  phi <- phi.p  + (1/(2*L))*(((n -1)*ro.ham/(ro.ham^2 - 1)) - (1/(2*sigma2_b*(1 - ro.ham^2)^2))*(t(res)%*%tridiag(rep(-(1 + ro.ham^2), n - 1), rep(-(1 + ro.ham^2), n - 1), c(2*ro.ham, rep(4*ro.ham, n-2), 2*ro.ham))%*%(res)))
  
  for(j in 1:L){
    ro.ham <- ro.ham + (1/(3*L))*drop(phi)
    phi <- phi  + (1/(2*L))*(((n -1)*ro.ham/(ro.ham^2 - 1)) - (1/(2*sigma2_b*(1 - ro.ham^2)^2))*(t(res)%*%tridiag(rep(-(1 + ro.ham^2), n - 1), rep(-(1 + ro.ham^2), n - 1), c(2*ro.ham, rep(4*ro.ham, n-2), 2*ro.ham))%*%(res)))
  }
  ro.ham <- ifelse(ro.ham > 1, 2 - ro.ham, ifelse(ro.ham < 0, -ro.ham, ro.ham))
  logr <- -((n - 1)/2)*(log(1 - ro.ham^2) - log(1 - ro_b^2)) - (1/(2*sigma2_b))*(t(res)%*%((1/(1 - ro.ham^2))*tridiag(rep(-ro_b, n - 1), rep(-ro_b, n - 1), c(1, rep(1 + ro_b^2, n-2), 1)) - (1/(1 - ro_b^2))*IC_ro)%*%(res)) + (dnorm(phi, mean = 0, sd = 1/3) - dnorm(phi.p, mean = 0, sd = 1/3))
                                                        
  if(runif(1) < exp(logr)){
    ro_b <- ro.ham
    ac <- ac + 1
  }
  
  if(i > 1000){
    par_M3 <- par_M3 + (1/(B - 1000))*c(beta_b, sigma2_b, ro_b)
    LL3 <- LL3 + (1/(B - 1000))*dmvnorm(y, mean = y - res, sigma = sigma2_b*ro_b^C_ro, log = T)
    y.est <- c(rmvnorm(1, mean = y - res, sigma = sigma2_b*ro_b^C_ro))
    M3 <- rbind(M3, c(mean(y.est), sqrt(var(y.est))))
  }
}
tictoc::toc()

rm(beta_b, sigma2_b, sig, ro.ham, ro_b, logr, phi, phi.p, C_ro, y.est)


#=============================================================================================

#                                              MODELO 4

#=============================================================================================
B <- 11000

LL4 <- 0
par_M4 <- 0
M4 <- NULL

b_zb <- mod$coefficients
sigma2_b <- summary(mod)$sigma^2
sigma2_0 <- sigma2_b
z_b <- rep(1, p)


logp <- function(z){
  X_z <- X%*%diag(z_b)[, z_b == 1]
  s20 <- try(summary(lm(y ~ -1 + X_z))$sigma^2, silent = TRUE)
  p_z <- sum(z)
  
  if(p_z == 0){
    Hg <- 0
    s20 <- mean(y^2)
  }else if(p_z > 0){
    Hg <- (n/(n + 1))*X_z%*%solve(t(X_z)%*%X_z)%*%t(X_z)
  }
  SSRg <- t(y)%*%(diag(1, nrow = n) - Hg)%*%y 
  return(-0.5*(p_z*log(1 + n) + (1 + n)*log(s20 + SSRg) - log(s20)))
}

SSR_n <- function(z){
  X_z <- X%*%diag(z_b)[, z_b == 1]
  s20 <- try(summary(lm(y ~ -1 + X_z))$sigma^2, silent = TRUE)
  Hg <- (n/(n + 1))*X_z%*%solve(t(X_z)%*%X_z)%*%t(X_z)
  return((t(y)%*%(diag(1, nrow = n) - Hg)%*%y) + s20)
}

logp_b <- logp(z_b)

set.seed(2110)
tictoc::tic()
for(i in 1:B){
  for(j in sample(1:p)){
    zp <- z_b
    zp[j] <- 1 - zp[j]
    logp.p <- logp(zp)
    r <- (logp.p - logp_b)*(-1)^(zp[j] == 0)
    z_b[j] <- rbinom(1, 1, (1/(1 + exp(-r))))
    
    if (z_b[j] == zp[j]){
      logp_b <- logp.p
    }
  }
  
  C_zb <- diag(z_b)[, z_b == 1]
  sigma2_b <- 1/rgamma(1, shape = (1 + n)/2, rate = (1/2)*(SSR_n(z_b)))
  
  sig <- ((n)/(n + 1))*solve(t(C_zb)%*%XtX%*%C_zb)
  b_zb <- c(mvtnorm::rmvnorm(1, mean = sig%*%t(C_zb)%*%Xty, sigma = sigma2_b*sig))
  
  beta_b <- C_zb%*%b_zb
  if(i > 1000){
    par_M4 <- par_M4 + (1/(B - 1000))*c(beta_b, sigma2_b)
    LL4 <- LL4 + (1/(B - 1000))*sum(dnorm(y, mean = X%*%beta_b, sd = sqrt(sigma2_b), log = T))
    y.est <- rnorm(n, mean = X%*%beta_b, sd = sqrt(sigma2_b))
    M4 <- rbind(M4, c(mean(y.est), sqrt(var(y.est))))
  }
}
tictoc::toc()
rm(b_zb, sig, z_b, sigma2_b, C_zb, zp, logp.p, logp_b, r, beta_b, IC_ro, res)




#=============================================================================================

#                                              PUNTO 1

#=============================================================================================
beta_est <- matrix(c(par_M1[1:10], par_M2[1:10], par_M3[1:10], par_M4[1:10]), nrow = 10, ncol = 4)

y.test_est <- X.test%*%beta_est
EAM <- colMeans(abs(y.test_est - y.test))


y_res <- cbind(y.test, y.test_est)
colnames(y_res) <- c('Obs', 'y_1', 'y_2', 'y_3', 'y_4')
y_res <- data.frame(y_res)

p_1 <- ggplot(data = y_res) +
  geom_point(mapping = aes(x = y_1, y = Obs)) +
  geom_abline(intercept = 0, slope = 1, col = 'red') +
  labs(title = paste('Modelo 1 EAM =', round(EAM[1], 3)), x = expression(hat(y)[list(test, 1)]), y = expression(y[test]))
p_2 <- ggplot(data = y_res) +
  geom_point(mapping = aes(x = y_2, y = Obs)) +
  geom_abline(intercept = 0, slope = 1, col = 'blue') +
  labs(title = paste('Modelo 2 EAM =', round(EAM[2], 3)), x = expression(hat(y)[list(test, 2)]), y = expression(y[test]))
p_3 <- ggplot(data = y_res) +
  geom_point(mapping = aes(x = y_3, y = Obs)) +
  geom_abline(intercept = 0, slope = 1, col = 'purple') +
  labs(title = paste('Modelo 3 EAM =', round(EAM[3], 3)), x = expression(hat(y)[list(test, 3)]), y = expression(y[test]))
p_4 <- ggplot(data = y_res) +
  geom_point(mapping = aes(x = y_4, y = Obs)) +
  geom_abline(intercept = 0, slope = 1, col = 'green') +
  labs(title = paste('Modelo 4 EAM =', round(EAM[4], 3)), x = expression(hat(y)[list(test, 4)]), y = expression(y[test]))

grid.arrange(p_1, p_2, p_3, p_4, nrow = 2)
rm(p_1, p_2, p_3, p_4, beta_est, y.test_est, EAM)

#=============================================================================================

#                                              PUNTO 2

#=============================================================================================
med_obs <- mean(y)
sd_obs <- sqrt(var(y))

y_prb <- cbind(M1, M2, M3, M4)
colnames(y_prb) <- c('med1', 'sd1', 'med2', 'sd2', 'med3', 'sd3', 'med4', 'sd4')
y_prb <- data.frame(y_prb)

x_range <- range(M1[, 1], M2[, 1], M3[, 1], M4[, 1])
y_range <- range(M1[, 2], M2[, 2], M3[, 2], M4[, 2])

p_1 <- ggplot(data = y_prb) +
  geom_point(mapping = aes(x = med1, y = sd1)) +
  geom_point(mapping = aes(x = med_obs, y = sd_obs), col = 'red', size = 2, shape = 15) +
  geom_vline(xintercept = med_obs, col = 'red') +
  geom_hline(yintercept = sd_obs, col = 'red') +
  labs(title = paste('Modelo 1'), x = expression(Media(bar(y))), y = expression(Desv.Est. (sigma))) +
  coord_cartesian(xlim = x_range, ylim = y_range)


p_2 <- ggplot(data = y_prb) +
  geom_point(mapping = aes(x = med2, y = sd2)) +
  geom_point(mapping = aes(x = med_obs, y = sd_obs), col = 'blue', size = 2, shape = 15) +
  geom_vline(xintercept = med_obs, col = 'blue') +
  geom_hline(yintercept = sd_obs, col = 'blue') +
  labs(title = paste('Modelo 2'), x = expression(Media (bar(y))), y = expression(Desv.Est. (sigma))) +
  coord_cartesian(xlim = x_range, ylim = y_range)


p_3 <- ggplot(data = y_prb) +
  geom_point(mapping = aes(x = med3, y = sd3)) +
  geom_point(mapping = aes(x = med_obs, y = sd_obs), col = 'purple', size = 2, shape = 15) +
  geom_vline(xintercept = med_obs, col = 'purple') +
  geom_hline(yintercept = sd_obs, col = 'purple') +
  labs(title = paste('Modelo 3'), x = expression(Media (bar(y))), y = expression(Desv.Est. (sigma))) +
  coord_cartesian(xlim = x_range, ylim = y_range)


p_4 <- ggplot(data = y_prb) +
  geom_point(mapping = aes(x = med4, y = sd4)) +
  geom_point(mapping = aes(x = med_obs, y = sd_obs), col = 'green', size = 2, shape = 15) +
  geom_vline(xintercept = med_obs, col = 'green') +
  geom_hline(yintercept = sd_obs, col = 'green') +
  labs(title = paste('Modelo 4'), x = expression(Media (bar(y))), y = expression(Desv.Est. (sigma))) +
  coord_cartesian(xlim = x_range, ylim = y_range)

grid.arrange(p_1, p_2, p_3, p_4, nrow = 2)
p_4

#=============================================================================================

#                                              PUNTO 3

#=============================================================================================

ll_ver <- c(sum(dnorm(y, mean = X%*%par_M1[1:10], sd = sqrt(par_M1[11]), log = T)), sum(dnorm(y, mean = X%*%par_M2[1:10], sd = sqrt(par_M2[11]), log = T)), dmvnorm(y, mean = X%*%par_M3[1:10], sigma = par_M3[11]*par_M3[12]^(abs(outer((1:n),(1:n) ,"-"))), log = T), sum(dnorm(y, mean = X%*%par_M4[1:10], sd = sqrt(par_M4[11]), log = T)))


p_DIC <- 2*(ll_ver - c(LL1, LL2, LL3, LL4))
DIC <- -2*ll_ver + 2*p_DIC

knitr::kable(x = rbind(p_DIC, DIC), digits = 3, align = "c", format = 'latex')


# Optimizar sobretodo el 3 (1 hora)
# Prueba matrix C_ro (?)
