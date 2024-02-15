##--------------------------------------- Librerias Necesarias ----------------------------------------
library(tidyverse)
library(metRology)
library(tictoc)
library(ggdendro)
library(factoextra)

##-------------------------------- Cargar datos y crear el data frame ---------------------------------
load("D:/Users/USUARIO/Downloads/Personas.RData")
df <- data.frame(dominio = dat$dominio, ingtot = log(dat$ingtot))
rm(dat)

##-------------------------------- Resultados de la muestra necesarios --------------------------------
m <- length(unique(df$dominio))
n <- length(df$ingtot)
y <- mean(df$ingtot)
s2 <- var(df$ingtot)
n_j <- tapply(df$ingtot, df$dominio, length)
y_j <- tapply(df$ingtot, df$dominio, mean)
s2_j <- tapply(df$ingtot, df$dominio, var)

##----------------------------- Definición de constantes e hiperparámetros ----------------------------
B <- 11000
mu0 <- 13.495
g20 <- 11.382
nu0 <- 1
s20 <- 1.182
eta0 <- 1
t20 <- 1.182
alpha0 <- 1
beta0 <- 0.846
k <- 3
a <- 1
b <- 1.182

#===============================================================================================
#                            Muestreador de Gibbs para cada modelo
#===============================================================================================
set.seed(2110)
##----------------------------------------- M1 Normal ------------------------------------------
theta_b.1 <- rnorm(1, mean = mu0, sd = sqrt(g20))
sigma2_b.1 <- 1/rgamma(1, shape = nu0/2, rate = nu0*s20/2)

##----------------------------- M2 Normal con medias especificas -------------------------------
mu_b.2 <- rnorm(1, mean = mu0, sd = sqrt(g20))
t2_b.2 <- 1/rgamma(1, shape = eta0/2, rate = eta0*t20/2)
theta_b.2 <- rnorm(m, mean = mu_b.2, sd = sqrt(t2_b.2))
sigma2_b.2 <- 1/rgamma(1, shape = nu0/2, rate = nu0*s20/2)

##---------------------- M3 Normal con medias y varianzas especificas ---------------------------
mu_b.3 <- rnorm(1, mean = mu0, sd = sqrt(g20))
t2_b.3 <- 1/rgamma(1, shape = eta0/2, rate = eta0*t20/2)
theta_b.3 <- rnorm(m, mean = mu_b.3, sd = sqrt(t2_b.3))
sigma20_b.3 <- rgamma(1, shape = alpha0/2, rate = beta0/2)
sigma2_b.3 <- 1/rgamma(m, shape = nu0/2, rate = nu0*sigma20_b.3/2)

##-------------------------------------- M4 t ---------------------------------------------------
theta_b.4 <- rnorm(1, mean = mu0, sd = sqrt(g20))
sigma2_b.4 <- rgamma(1, shape = alpha0/2, rate = beta0/2)
zeta2_ijb.4 <- 1/rgamma(n, shape = k/2, rate = k*sigma2_b.4/2)

##--------------------------- M5 t con medias especificas ---------------------------------------
mu_b.5 <- rnorm(1, mean = mu0, sd = sqrt(g20))
t2_b.5 <- 1/rgamma(1, shape = eta0/2, rate = eta0*t20/2)
theta_b.5 <- rnorm(m, mean = mu_b.5, sd = sqrt(t2_b.5))
sigma2_b.5 <- rgamma(1, shape = alpha0/2, rate = beta0/2)
zeta2_ijb.5 <- 1/rgamma(n, shape = k/2, rate = k*sigma2_b.5/2)

##----------------------- M6 t con medias y varianzas especificas -------------------------------
mu_b.6 <- rnorm(1, mean = mu0, sd = sqrt(g20))
t2_b.6 <- 1/rgamma(1, shape = eta0/2, rate = eta0*t20/2)
theta_b.6 <- rnorm(m, mean = mu_b.6, sd = sqrt(t2_b.6))
beta_b.6 <- rgamma(1, shape = a/2, rate = b/2)
sigma2_b.6 <- rgamma(m, shape = alpha0/2, rate = beta_b.6/2)
zeta2_ijb.6 <- 1/rgamma(n, shape = k/2, rate = k*rep(sigma2_b.6, n_j)/2)

#================================================================================================
#                                    Ajuste de los modelos
#================================================================================================
for(i in 1:6){
  assign(paste0('theta_bayes', i), 0)
}
M_6 <- NULL
LL_two <- rep(NA, 60000)
tictoc::tic()
set.seed(2110)
for(i in 1:11000){
  
  theta_b.1 <- rnorm(1, mean = (mu0/g20 + (n*y)/sigma2_b.1)*(1/(1/g20 + n/sigma2_b.1)), sd = sqrt(1/(1/g20 + n/sigma2_b.1)))
  theta_b.2 <- rnorm(m, mean = (mu_b.2/t2_b.2 + (n_j*y_j)/sigma2_b.2)*(1/(1/t2_b.2 + n_j/sigma2_b.2)), sd = sqrt(1/(1/t2_b.2 + n_j/sigma2_b.2)))
  theta_b.3 <- rnorm(m, mean = (mu_b.3/t2_b.3 + (n_j*y_j)/sigma2_b.3)*(1/(1/t2_b.3 + n_j/sigma2_b.3)), sd = sqrt(1/(1/t2_b.3 + n_j/sigma2_b.3)))
  theta_b.4 <- rnorm(1, mean = (mu0/g20 + sum(df$ingtot/zeta2_ijb.4))*(1/(1/g20 + sum(1/zeta2_ijb.4))), sd = sqrt(1/(1/g20 + sum(1/zeta2_ijb.4))))
  theta_b.5 <- rnorm(m, mean = (mu_b.5/t2_b.5 + tapply(df$ingtot/zeta2_ijb.5, df$dominio, sum))*(1/(1/t2_b.5 + tapply(1/zeta2_ijb.5, df$dominio, sum))), sd = sqrt(1/(1/t2_b.5 + tapply(1/zeta2_ijb.5, df$dominio, sum))))
  theta_b.6 <- rnorm(m, mean = (mu_b.6/t2_b.6 + tapply(df$ingtot/zeta2_ijb.6, df$dominio, sum))*(1/(1/t2_b.6 + tapply(1/zeta2_ijb.6, df$dominio, sum))), sd = sqrt(1/(1/t2_b.6 + tapply(1/zeta2_ijb.6, df$dominio, sum))))
  
  
  sigma20_b.3 <- rgamma(1, shape = (alpha0 + m*nu0)/2, rate = beta0/2 + (nu0/2)*sum(1/sigma2_b.3))
  
  
  sigma2_b.1 <- 1/rgamma(1, shape = (nu0 + n)/2, rate = (nu0*s20 + (n-1)*s2 + n*(y - theta_b.1)^2)/2)
  sigma2_b.2 <- 1/rgamma(1, shape = (nu0 + n)/2, rate = (nu0*s20 + sum((n_j-1)*s2_j + n_j*(y_j - theta_b.2)^2))/2)
  sigma2_b.3 <- 1/rgamma(m, shape = (nu0 + n_j)/2, rate = (nu0*sigma20_b.3 + (n_j-1)*s2_j + n_j*(y_j - theta_b.3)^2)/2)
  sigma2_b.4 <- rgamma(1, shape = (alpha0 + n*k)/2, rate = (beta0/2) + (k/2)*sum(1/zeta2_ijb.4))
  sigma2_b.5 <- rgamma(1, shape = (alpha0 + n*k)/2, rate = (beta0/2) + (k/2)*sum(1/zeta2_ijb.5))
  sigma2_b.6 <- rgamma(m, shape = (alpha0 + (n_j*k))/2, rate = (beta_b.6/2) + (k/2)*tapply(1/zeta2_ijb.6, df$dominio, sum))
  
  
  beta_b.6 <- rgamma(1, shape = (a + (m*alpha0))/2, rate = (b/2) + (1/2)*sum(sigma2_b.6))
  
  
  zeta2_ijb.4 <- 1/rgamma(n, shape = (k + 1)/2, rate = (k*sigma2_b.4 + (df$ingtot - theta_b.4)^2)/2)
  zeta2_ijb.5 <- 1/rgamma(n, shape = (k + 1)/2, rate = (k*sigma2_b.5 + (df$ingtot - rep(theta_b.5, n_j))^2)/2)
  zeta2_ijb.6 <- 1/rgamma(n, shape = (k + 1)/2, rate = (k*rep(sigma2_b.6, n_j) + (df$ingtot - rep(theta_b.6, n_j))^2)/2)
  
  
  mu_b.2 <- rnorm(1, mean = (mu0/g20 + (m*mean(theta_b.2))/t2_b.2)*(1/(1/g20 + m/t2_b.2)), sd = sqrt(1/(1/g20 + m/t2_b.2)))
  mu_b.3 <- rnorm(1, mean = (mu0/g20 + (m*mean(theta_b.3))/t2_b.3)*(1/(1/g20 + m/t2_b.3)), sd = sqrt(1/(1/g20 + m/t2_b.3)))
  mu_b.5 <- rnorm(1, mean = (mu0/g20 + (m*mean(theta_b.5))/t2_b.5)*(1/(1/g20 + m/t2_b.5)), sd = sqrt(1/(1/g20 + m/t2_b.5)))
  mu_b.6 <- rnorm(1, mean = (mu0/g20 + (m*mean(theta_b.6))/t2_b.6)*(1/(1/g20 + m/t2_b.6)), sd = sqrt(1/(1/g20 + m/t2_b.6)))
  
  
  t2_b.2 <- 1/rgamma(1, shape = (eta0 + m)/2, rate = ((eta0*t20) + (m-1)*var(theta_b.2) + m*(mean(theta_b.2) - mu_b.2)^2)/2)
  t2_b.3 <- 1/rgamma(1, shape = (eta0 + m)/2, rate = ((eta0*t20) + (m-1)*var(theta_b.3) + m*(mean(theta_b.3) - mu_b.3)^2)/2)
  t2_b.5 <- 1/rgamma(1, shape = (eta0 + m)/2, rate = ((eta0*t20) + (m-1)*var(theta_b.5) + m*(mean(theta_b.5) - mu_b.5)^2)/2)
  t2_b.6 <- 1/rgamma(1, shape = (eta0 + m)/2, rate = ((eta0*t20) + (m-1)*var(theta_b.6) + m*(mean(theta_b.6) - mu_b.6)^2)/2)
  
  if(i > 1000){
    LL_two[(6*(i - 1000)) - 5] <- sum(dnorm(df$ingtot, mean = theta_b.1, sd = sqrt(sigma2_b.1), log = T))
    LL_two[(6*(i - 1000)) - 4 ] <- sum(dnorm(df$ingtot, mean = rep(theta_b.2, n_j), sd = sqrt(sigma2_b.2), log = T))
    LL_two[(6*(i - 1000)) - 3] <- sum(dnorm(df$ingtot, mean = rep(theta_b.3, n_j), sd = sqrt(rep(sigma2_b.3, n_j)), log = T))
    LL_two[(6*(i - 1000)) - 2] <- sum(dt.scaled(df$ingtot,df = k, mean = theta_b.4, sd = sqrt(sigma2_b.4), log = T))
    LL_two[(6*(i - 1000)) - 1] <- sum(dt.scaled(df$ingtot,df = k, mean = rep(theta_b.5, n_j), sd = sqrt(sigma2_b.5), log = T))
    LL_two[(6*(i - 1000))] <- sum(dt.scaled(df$ingtot,df = k, mean = rep(theta_b.6, n_j), sd = sqrt(rep(sigma2_b.6, n_j)), log = T))
    
    for(j in 1:6){
      assign(paste0('theta_bayes', j), get(paste0('theta_bayes', j)) + c(get(paste0('theta_b.', j)), get(paste0('sigma2_b.', j))))
    }
    
    M_6 <- rbind(M_6, c(theta_b.6, sigma2_b.6))
  }
  
}
tictoc::toc()

rm(theta_b.1, sigma2_b.1)
rm(mu_b.2, t2_b.2, theta_b.2, sigma2_b.2)
rm(mu_b.3, t2_b.3, theta_b.3,sigma20_b.3, sigma2_b.3)
rm(theta_b.4, sigma2_b.4, zeta2_ijb.4)
rm(mu_b.5, t2_b.5, theta_b.5, sigma2_b.5, zeta2_ijb.5)
rm(mu_b.6, t2_b.6, theta_b.6, beta_b.6, sigma2_b.6, zeta2_ijb.6)


##------------------------------------------------- Punto 2 -------------------------------------------------------
LL <- data.frame('LL_n' = LL_two, 'mod' = rep(c('M1', 'M2', 'M3', 'M4', 'M5', 'M6'), 10000))


ggplot(data = LL) + 
  geom_point(mapping = aes(x = seq(1:60000), y = LL_n, color = mod)) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange", "black")) +
  theme(axis.text.x = element_blank()) +
  xlab("Iteración")
  ylab("Log-verosimilitud")
  
rm(LL_two)
##------------------------------------------------- Punto 3 -------------------------------------------------------
for(i in 1:6){
  assign(paste0('theta_bayes', i), get(paste0('theta_bayes', i))/(B-1000))
}


ll_k <- c(sum(dnorm(df$ingtot, mean = theta_bayes1[1], sd = sqrt(theta_bayes1[2]), log = T)),
          sum(dnorm(df$ingtot, mean = rep(theta_bayes2[1:m],n_j), sd = sqrt(theta_bayes2[m + 1]), log = T)),
          sum(dnorm(df$ingtot, mean = rep(theta_bayes3[1:m], n_j), sd = sqrt(rep(theta_bayes3[(m+1):(2*m)], n_j)), log = T)),
          sum(dt.scaled(df$ingtot,df = k, mean = theta_bayes4[1], sd = sqrt(theta_bayes4[2]), log = T)),
          sum(dt.scaled(df$ingtot,df = k, mean = rep(theta_bayes5[1:m],n_j), sd = sqrt(theta_bayes5[m + 1]), log = T)),
          sum(dt.scaled(df$ingtot,df = k, mean = rep(theta_bayes6[1:m], n_j), sd = sqrt(rep(theta_bayes6[(m+1):(2*m)], n_j)), log = T)))

p_DIC <- 2*(ll_k - tapply(LL$LL_n, LL$mod, mean))
DIC <- -2*ll_k + 2*p_DIC
knitr::kable(x = cbind(p_DIC, DIC), digits = 3, align = "c")#, format = 'latex')
knitr::kable(x = DIC, digits = 3, align = "c", format = 'latex')
max(DIC)
rm(theta_bayes1, theta_bayes2, theta_bayes3, theta_bayes4, theta_bayes5, ll_k, p_DIC, DIC)

##------------------------------------------------- Punto 4 -------------------------------------------------------

med_bog <- median(df[df$dominio == 'BOGOTA',]$ingtot)
sd_bog <- sd(df[df$dominio == 'BOGOTA',]$ingtot)
coef_bog <- sd_bog/y_j[3]
rango_bog <- max(df[df$dominio == 'BOGOTA',]$ingtot) - min(df[df$dominio == 'BOGOTA',]$ingtot)
ran.int_bog <- quantile(df[df$dominio == 'BOGOTA',]$ingtot, 0.75) - quantile(df[df$dominio == 'BOGOTA',]$ingtot, 0.25)

est <- c(0, 0, 0, 0, 0, 0)

tic()
for(i in 1:10000){
  X <- rt.scaled(n_j[3],df = k, mean = M_6[i, 3], sd = sqrt(M_6[i, 28]))
  
  est <- est + (1/B)*c((mean(X) > y_j[3]), (median(X) > med_bog) , ((sqrt(3)*sd(X)) > sd_bog), ((sqrt(3)*sd(X))/mean(X) > coef_bog), (max(X) - min(X) > rango_bog), (quantile(X, 0.75) - quantile(X, 0.25) > ran.int_bog))
}
toc()


names(est) <- c('Media', 'Mediana','Desv. Est.', 'Coef. Var', 'Rango', 'Rango Intercuart.')
knitr::kable(x = cbind(est, c(y_j[3], med_bog, sd_bog, coef_bog, rango_bog, ran.int_bog)), digits = 3, align = "c")#, format = 'latex')

rm(med_bog, sd_bog, coef_bog, rango_bog, ran.int_bog, est)

##------------------------------------------------- Punto 5 -------------------------------------------------------

nom <- unique(df$dominio)
that  <- colMeans(M_6[,1:m])
ic1   <- apply(X = M_6[,1:m], MARGIN = 2, FUN = function(x) quantile(x, c(0.025,0.975)))
ranking <- order(that)
nom <- nom[ranking]
that <- that[ranking]
ic1  <- ic1 [, ranking]
colo <- rep(2,m)
colo[which(ic1[1,]>13.83)] <- 1
colo[which(ic1[2,]<13.83)] <- 3
colo <- c("green","black","red")[colo]

df_plot <- data.frame(x = that, y = nom)

ggplot(df_plot, aes(x = x, y = y)) +
  geom_segment(aes(x = ic1[1,], xend = ic1[2,], y = y, yend = y), color = colo) +
  geom_point(size = 2, color = colo) +
  geom_vline(xintercept = 13.83, color = "gray", size = 1.5) +
  xlim(12.85, 14.15) +
  scale_y_discrete(limits = nom) +
  labs(x = "Promedio", y = "", title = "Ranking Bayesisano: Modelo 6") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, hjust = 1))

rm(that, ic1, colo, df_plot)

##------------------------------------------------- Punto 6 -------------------------------------------------------

media_top <- colMeans((M_6[,1:25])[,ranking[21:25]])
var_top <- colMeans((M_6[,25:50])[,ranking[21:25]])

knitr::kable(x = cbind(nom[21:25], round(exp(media_top), 3), round(exp(sqrt(var_top)), 3), round(sqrt(var_top)/media_top, 3) * 100), digits = 3, align = "c", format = 'latex')

##------------------------------------------------- Punto 7 -------------------------------------------------------

M_6 <- M_6[ ,c(ranking, ranking + 25)]
minmax <- function(X){
  return((X - min(X))/(max(X) - min(X)))
}

seg <- NULL
tic()
for(i in 1:10000){
  dat <- matrix(c(minmax(M_6[i, 1:25]), minmax(M_6[i, 26:50])), nrow = 25, ncol = 2)
  
  group_k <- cutree(hclust(dist(dat)), k = 4)
  
  seg <- rbind(seg, group_k)
}
toc()

incid <- NULL
for(i in 1:25){
  for(j in 1:25){
    incid <- c(incid, mean(seg[, i] == seg[, j]))
  }
}


grid <- expand.grid(nom[ranking], nom[ranking])
grid$incid <- incid
ggplot(grid) + 
  geom_tile(mapping = aes(x = Var1, y = Var2, fill = incid)) +
  scale_fill_gradientn(colors = hcl.colors(25, "Spectral")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Matriz de Incidencia") +
  labs(y= "", x = "")


M_6_mm <- colMeans(M_6)
  

dat <- matrix(M_6_mm, nrow = 2, ncol = 25, byrow = T)
dat1 <- matrix(c(minmax(dat[1,]), minmax(dat[2,])), nrow = 25, ncol = 2)


rownames(dat1) <- nom[ranking]

ggplot() + 
  geom_point(mapping = aes(x = dat1[,1], y = dat1[,2], color = group_k)) + 
  geom_text(mapping = aes(x = dat1[,1], y = dat1[,2], label = nom[ranking], color = group_k)) +
  scale_color_gradientn(colors = hcl.colors(25, "Spectral")) +
  ggtitle("") +
  labs(y= expression(sigma), x = expression(theta))


clust <-  hclust(dist(dat1))
group_k <- cutree(clust, k = 4)

plot(clust, main = 'Dendrograma', xlab = 'Dominios')
rect.hclust(clust, k = 4)
