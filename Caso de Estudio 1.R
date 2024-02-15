library(readr)
library(foreach)
setwd("D:/Users/USUARIO/Downloads")
## Pre-procesamiento de la base de datos
df <- read_csv("victimas.csv", col_types = cols(RUPTURA = col_skip(), 
                                                CONEXO = col_skip(), ETAPA = col_skip(), 
                                                LEY = col_skip(), MUNICIPIO = col_skip(), 
                                                SECCIONAL = col_skip(), DELITO = col_skip(), 
                                                IMPUTACION = col_skip(), CONDENA = col_skip(), 
                                                ATIPICIDAD_INEXISTENCIA = col_skip(), 
                                                ACUSACION = col_skip(), CAPTURA = col_skip(), 
                                                HOMICIDIO_DOLOSO_CONSUMADO = col_skip()))

df <- df[!is.na(df$SEXO_VICTIMA),]
df <- df[!is.na(df$TOTAL_VICTIMAS),]
df <- df[df$HECHO == 'SI',]
df <- df[!(names(df) %in% c('HECHO'))]
df <- df[df$ESTADO_NOTICIA == 'ACTIVO',]
df <- df[!(names(df) %in% c('ESTADO_NOTICIA'))]
df <- df[df$ANIO_DENUNCIA == 2022,]
df <- df[!(names(df) %in% c('ANIO_DENUNCIA'))]
df <- df[df$ANIO_ENTRADA == 2022,]
df <- df[!(names(df) %in% c('ANIO_ENTRADA'))]
df <- df[df$ANIO_HECHO == 2022,]
df <- df[!(names(df) %in% c('ANIO_HECHO'))]
df <- df[df$PAIS == 'Colombia',]
df <- df[!(names(df) %in% c('PAIS'))]
df <- df[df$DEPARTAMENTO == 'BOGOTÁ, D. C.',]
df <- df[!(names(df) %in% c('DEPARTAMENTO'))]
df <- df[df$GRUPO_DELITO == 'DELITOS SEXUALES',]
df <- df[!(names(df) %in% c('GRUPO_DELITO'))]
df <- df[df$PAIS_NACIMIENTO == 'Colombia',]
df <- df[!(names(df) %in% c('PAIS_NACIMIENTO'))]
df <- df[df$GRUPO_EDAD_VICTIMA %in% c('PRIMERA INFANCIA 0 - 5', 'INFANCIA 6 - 11', 'PRE-ADOLESCENTE 12 - 13', 'ADOLESCENTE 14 - 17'),]
df <- df[!(names(df) %in% c('GRUPO_EDAD_VICTIMA'))]
df <- df[df$TOTAL_VICTIMAS <= (5/2)*quantile(df$TOTAL_VICTIMAS,0.75) + (3/2)*quantile(df$TOTAL_VICTIMAS,0.25) & 
          df$TOTAL_VICTIMAS >= (5/2)*quantile(df$TOTAL_VICTIMAS,0.25) - (3/2)*quantile(df$TOTAL_VICTIMAS,0.75)
         ,]

gc()

table(df$SEXO_VICTIMA)
# Se obtienen conteos de 115 hombres y 237 mujeres


y_h <- as.numeric(df[df$SEXO_VICTIMA == 'MASCULINO',]$TOTAL_VICTIMAS) # Total de victimas hombres
y_m <- as.numeric(df[df$SEXO_VICTIMA == 'FEMENINO',]$TOTAL_VICTIMAS) # Total de victimas mujeres


#                                 Primera visualización de los datos
#-----------------------------------------------------------------------------------------------
par(mfrow = c(1,1), mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
y <- 0:12
plot(y - .07, table(factor(x = y_h, levels = y))/length(y_h), col = 2, type = "h", ylim = c(0,.7), lwd = 3, ylab = "F. Relativa", xlab = "# de Victimas", main = "Proporción de # de Victimas por género", yaxt = "n")
points(y + .07, table(factor(x = y_m, levels = y))/length(y_m), col = 4, type = "h", lwd = 3)
axis(side = 2)
legend("topright", legend = c("Masculino", "Femenino"), bty = "n", lwd = 2, col = c(2,4))
#------------------------------------------------------------------------------------------------

#                           Tamaños de muestra y estaísticos suficientes
#------------------------------------------------------------------------------------------------
n_m <- length(y_m)
n_h <- length(y_h)
s_m <- sum(y_m)
s_h <- sum(y_h)

sd_m <- sd(y_m)
sd_h <- sd(y_h)

rm(df,y_m, y_h)
gc()
#------------------------------------------------------------------------------------------------


#===========================================================================
#                          Punto 1 A
#===========================================================================

# Hiperparámetros
a <- 0.01
b <- 0.01

# Parámetrso distribución posterior
a_post <- a + c(s_h, s_m)
b_post <- b + c(n_h, n_m)

#                                                   Gráfico de las posteriores y la previa
#------------------------------------------------------------------------------------------
par(mfrow = c(1,2), mar = c(3,3,1.4,1.4), mgp = c(1.75,.75,0))
theta <- seq(0, 5, length = 1000)
plot(NA, NA, xlim = c(0,4), ylim = c(0,5.5), xlab = expression(theta), ylab = expression(paste("p","(",theta," | ",y,")",sep="")), main = "Posterior y previa")
lines(theta, dgamma(theta, shape = a_post[1], rate = b_post[1]), col = 2, lwd = 2)
lines(theta, dgamma(theta, shape = a_post[2], rate = b_post[2]), col = 4, lwd = 2)
lines(theta, dgamma(theta, shape = a  , rate = b  ), col = 1, lwd = 1)
legend("topright", legend = c("Masculino", "Femenino", "Previa"), bty = "n", lwd = 2, col = c(2, 4, 1))
y <- 0:12
plot(y - .07, dnbinom(y, size = a_post[1], mu = a_post[1]/b_post[1]), col = 2, type = "h", ylim = c(0,.3), lwd = 3, ylab = "p(y* | y )", xlab = "y*", main = "Predictiva posterior")
points(y + .07, dnbinom(y, size = a_post[2], mu = a_post[2]/b_post[2]), col = 4, type = "h", lwd = 3)
#----------------------------------------------------------------------------------------------



#===========================================================================
#                          Punto 1 B
#===========================================================================



set.seed(2110)
B <- 10000 # # Número de muestras elegido

#    Simulación por Monte Carlo de la distribución empirica de eta
eta <- (rgamma(B, shape = a_post[2], rate = b_post[2])/rgamma(B, shape = a_post[1], rate = b_post[1]))-1

#                     Gráfico de la distribución empirica de eta por Monte Carlo
#----------------------------------------------------------------------------------------------------------
par(mfrow = c(1,1), mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
hist(eta, freq = F, col = "gray90", border = "gray90", ylab = expression(eta),xlab = "", main = paste0("Distribución de ", expression(eta),' por MonteCarlo'))
lines(density(eta), col = 4, lwd = 2)
abline(v = mean(eta), col = 1, lwd = 2, lty = 1)
abline(v = quantile(x = eta, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 2)
legend("topright", legend = c(expression(eta), "IC 95%", "Media"), col = c(4, 2, 1), lty = 1, lwd = 2, bty = "n")
#----------------------------------------------------------------------------------------------------------

#             Tabla con la media, coeficiente de variación y un intervalo de credibilidad al 95% para eta
#-----------------------------------------------------------------------------------------------------------
out <- c(mean(eta), sd(eta)/mean(eta),quantile(eta,c(0.025, 0.975)))
names(out) <- c('Media', 'Coeficiente de Variación', 'LI 95%', 'LS 95%')
knitr::kable(x = out, digits = 3, align = "c")
#------------------------------------------------------------------------------------------------------------
rm(eta)


#===========================================================================
#                          Punto 1 C
#===========================================================================
set.seed(2110)

# Estados de información externos dados
a_k <- c(0.01, 0.1, 1, 1, 1, 1)
b_k <- c(0.01, 0.1, 1, 1/2, 1/3, 1/4)

# Parámetros de la distribución posterior de cada distribución previa dada
a_post_kh <- a_k + s_h
a_post_km <- a_k + s_m
b_post_kh <- b_k + n_h
b_post_km <- b_k + n_m

#                Gráficos de la distribución empirica por Monte Carlo por cada distribución previa
#------------------------------------------------------------------------------------------------------
out <- NULL
par(mfrow = c(3,2), mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
for(i in 1:6){
  # Cálculo de eta por cada distribución previa
  eta_k <- (rgamma(B, shape = a_post_km[i], rate = b_post_km[i])/rgamma(B, shape = a_post_kh[i], rate = b_post_kh[i]))-1
  
  # Gráfico de cada distribución empirica por cada distribución previa
  hist(eta_k, freq = F, col = "gray90", border = "gray90", ylab = expression(eta),xlab = "", main = paste0("Distribución de ", expression(eta),' para la distribución previa ', i))
  lines(density(eta_k), col = 4, lwd = 2)
  abline(v = mean(eta_k), col = 1, lwd = 2, lty = 1)
  abline(v = quantile(x = eta_k, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 2)
  
  # Registro de la media, coeficiente de variación y coefciente de credibilidad por cada distribución previa
  out <- rbind(out, c(mean(eta_k), sd(eta_k)/mean(eta_k), quantile(eta_k,c(0.025, 0.975))))

  rm(eta_k)
  gc()
}
legend("topright", legend = c(expression(eta), "IC 95%", "Media"), col = c(4, 2, 1), lty = 1, lwd = 2, bty = "n")
#---------------------------------------------------------------------------------------------------------------------

#Tabla con la media apriori, cv apriori, media, cv y un intervalo de confianza de la distribución empirica por Monte Carlo para cada distribución
out <- cbind(out, a_k/b_k, 1/sqrt(a_k))
colnames(out) <- c('Media', 'Coeficiente de Variación', 'LI 95%', 'LS 95%', 'Media apriori', 'CV apriori')
rownames(out) <- c(paste0('Distr. Previa', 1:6))
knitr::kable(x = out, digits = 3, align = "c")


#===========================================================================
#                          Punto 1 D
#===========================================================================

set.seed(2110)
out <- NULL
# Matriz para los estadísticos de prueba, primera fila para hombres y segunda para mujeres
est_prueba <- matrix(0,nrow = 4, ncol = B)
rownames(est_prueba) <- c('Media Hombres', 'Media Mujeres', 'Desv. Est. Hombres', 'Desv. Est. Mujeres')

# Estadísticos observados de la muestra
t_obs <- c(s_h/n_h, s_m/n_m, sd_h, sd_m)

#                         Bondad de Ajuste
#--------------------------------------------------------------------------------------
for(i in 1:6){
  # Simulación de las muestras y cálculo de los estadísticos de prueba
  foreach(j = 1:B) %do% {
    prueba_h <- rpois(n = n_h, lambda = rgamma(1, shape = a_post_kh[i], rate = b_post_kh[i]))
    prueba_m <- rpois(n = n_m, lambda = rgamma(1, shape = a_post_km[i], rate = b_post_km[i]))
    est_prueba[,j] <- c(mean(prueba_h), mean(prueba_m), sd(prueba_h), sd(prueba_m))
    rm(prueba_h)
    rm(prueba_m)
  }
  par(mfrow = c(2,2), mar = c(3,3,1.4,1.4), mgp = c(1.75,.75,0))
  # Visualización de las distribuciones empíricas de los estadísticos de prueba por cada distribución previa
  for(j in 1:4){
    hist(x = est_prueba[j,], freq = F, col = "gray90", border = "gray90", xlab = "t", ylab = "p(t | y)", main = paste(rownames(est_prueba)[j], i))
    lines(density(est_prueba[j,]), col = 4, lwd = 2)
    abline(v = t_obs[j], col = 1, lwd = 2, lty = 1)
    abline(v = quantile(x = est_prueba[j,], probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 2)
  }
  legend("topright", legend = c("Posterior", "IC 95%", "t obs"), col = c(4, 2, 1), lty = 1, lwd = 2, bty = "n")
  
  # Cálculo de los valores p predictivos posteriores
  out <- rbind(out, rowMeans(est_prueba > t_obs))
  
  
  est_prueba <- matrix(0,nrow = 4, ncol = B)
  rownames(est_prueba) <- c('Media Hombres', 'Media Mujeres', 'Desv. Est. Hombres', 'Desv. Est. Mujeres')
}
#--------------------------------------------------------------------------------------------
rm(est_prueba)
gc()

#                     Tabla de valores p predicitivos posteriores para los estadísticos de prueba
#---------------------------------------------------------------------------------------------------
rownames(out) <- c(paste0('Distr. Previa', 1:6))
knitr::kable(x = out, digits = 3, align = "c")
#---------------------------------------------------------------------------------------------------


#===========================================================================
#                          Punto 2 A
#===========================================================================

set.seed(2110)
B <- 50000 # Número de remuestras elegido

# Cálculo de la distribución empirica de eta usando Bootstrap paramétrico
eta <- ((rpois(B, lambda = s_m)*(1/n_m))/(rpois(B, lambda = s_h)*(1/n_h)))-1


#                     Gráfico de la distribución empirica de eta por Bootstrap Paramétrico
#----------------------------------------------------------------------------------------------------------
par(mfrow = c(1,1), mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
hist(eta, freq = F, col = "gray90", border = "gray90", ylab = expression(eta),xlab = "", main = paste0("Distribución de ", expression(eta),' por Bootstrap Paramétrico'))
lines(density(eta), col = 4, lwd = 2)
abline(v = mean(eta), col = 1, lwd = 2, lty = 1)
abline(v = quantile(x = eta, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 2)
legend("topright", legend = c(expression(eta), "IC 95%", "Media"), col = c(4, 2, 1), lty = 1, lwd = 2, bty = "n")
#----------------------------------------------------------------------------------------------------------


#             Tabla con la media, coeficiente de variación y un intervalo de credibilidad al 95% para eta
#-----------------------------------------------------------------------------------------------------------
out <- c(mean(eta), sd(eta)/mean(eta),quantile(eta,c(0.025, 0.975)))
names(out) <- c('Media', 'Coeficiente de Variación', 'LI 95%', 'LS 95%')
knitr::kable(x = out, digits = 3, align = "c")
#------------------------------------------------------------------------------------------------------------
rm(eta)
gc()

#===========================================================================
#                          Punto 2 B
#===========================================================================

# Función indicadora de un intervalo
isin <- function(x, ic){
  if(ic[1] <= x & x <= ic[2]){return(T)}else{return(F)}
}

out <- NULL
B <- 2000 # Número de remuestras para el Bootstrap y de simulaciones para Monte Carlo elegido
n_k <- c(10, 20, 50, 100) # Tamaños de muestra de cada escenario

# Valor verdadero de eta
eta_real <- ((s_m/n_m)/(s_h/n_h)) - 1

for(j in n_k){
  # Cobertura de ambos métodos
  cob_bayes <- 0
  cob_frec <- 0
  
  # Simulación de las cien mil muestras aleatorias
  foreach(i = 1:100000) %do% {
    # Estadístico suficiente simulado de la distribución muestral
    Ssim_h <- rpois(n = 1, lambda = j * (s_h/n_h))
    Ssim_m <- rpois(n = 1, lambda = j * (s_m/n_m))
    
    # Intervalo de credibilidad al 95% para eta 
    ic_bayes <- quantile((rgamma(B, shape = a + Ssim_m, rate = b + j)/rgamma(B, shape = a + Ssim_h, rate = b + j))-1, c(0.025, 0.975))
    
    # Intervalo de confianza al 95% para eta
    ic_frec <- quantile(((rpois(B, lambda = Ssim_m)*(1/j))/(rpois(B, lambda = Ssim_h)*(1/j)))-1, c(0.025, 0.975))
    
    # Evaluación de la cobertura para los intervalos de credibilidad y confianza 
    if(isin(eta_real, ic_bayes)){cob_bayes <- cob_bayes + 1}
    if(isin(eta_real, ic_frec)){cob_frec <- cob_frec + 1}
  }
  
  # Cálculo de los porcentajes de cobertura de los intervalos de credibilidad y confianza
  out <- rbind(out, c(j, cob_bayes/100000, cob_frec/100000))
}
rm(ic_bayes, ic_frec, cob_bayes, cob_frec)

#                  Tabla de los porcentajes de cobertura de los intervalos de credibilidad y confianza
#------------------------------------------------------------------------------------------------------
colnames(out) <- c('Tamaño de muestra', 'Cobertura Int. Credibilidad', 'Cobertura Int. Confianza')
rownames(out) <- c(paste0('Escenario', 1:4))
knitr::kable(x = out, digits = 3, align = "c")
#-------------------------------------------------------------------------------------------------------


#===========================================================================
#                          Punto 3
#===========================================================================

# Pre-procesamiento de la base de datos
df <- read_csv("victimas.csv", col_types = cols(RUPTURA = col_skip(), 
                                                CONEXO = col_skip(), ETAPA = col_skip(), 
                                                LEY = col_skip(), MUNICIPIO = col_skip(), 
                                                SECCIONAL = col_skip(), DELITO = col_skip(), 
                                                IMPUTACION = col_skip(), CONDENA = col_skip(), 
                                                ATIPICIDAD_INEXISTENCIA = col_skip(), 
                                                ACUSACION = col_skip(), CAPTURA = col_skip(), 
                                                HOMICIDIO_DOLOSO_CONSUMADO = col_skip()))

# Aplicación de filtros
df <- df[!is.na(df$SEXO_VICTIMA),]
df <- df[df$HECHO == 'SI',]
df <- df[!(names(df) %in% c('HECHO'))]
df <- df[df$ESTADO_NOTICIA == 'ACTIVO',]
df <- df[!(names(df) %in% c('ESTADO_NOTICIA'))]
df <- df[df$ANIO_ENTRADA == df$ANIO_DENUNCIA & df$ANIO_DENUNCIA == df$ANIO_HECHO,]
df <- df[!(names(df) %in% c('ANIO_DENUNCIA'))]
df <- df[!(names(df) %in% c('ANIO_ENTRADA'))]
df <- df[df$PAIS == 'Colombia',]
df <- df[!(names(df) %in% c('PAIS'))]
df <- df[df$DEPARTAMENTO == 'BOGOTÁ, D. C.',]
df <- df[!(names(df) %in% c('DEPARTAMENTO'))]
df <- df[df$GRUPO_DELITO == 'DELITOS SEXUALES',]
df <- df[!(names(df) %in% c('GRUPO_DELITO'))]
df <- df[df$PAIS_NACIMIENTO == 'Colombia',]
df <- df[!(names(df) %in% c('PAIS_NACIMIENTO'))]
df <- df[df$GRUPO_EDAD_VICTIMA %in% c('PRIMERA INFANCIA 0 - 5', 'INFANCIA 6 - 11', 'PRE-ADOLESCENTE 12 - 13', 'ADOLESCENTE 14 - 17'),]
df <- df[!(names(df) %in% c('GRUPO_EDAD_VICTIMA'))]
df <- df[df$ANIO_HECHO <= 2022 & df$ANIO_HECHO >= 2012,]
df <- df[df$TOTAL_VICTIMAS <= (5/2)*quantile(df$TOTAL_VICTIMAS,0.75) + (3/2)*quantile(df$TOTAL_VICTIMAS,0.25) & 
           df$TOTAL_VICTIMAS >= (5/2)*quantile(df$TOTAL_VICTIMAS,0.25) - (3/2)*quantile(df$TOTAL_VICTIMAS,0.75)
         ,]

gc()

# Tabla con los tamaños de muestra por sexo y año
table(df$SEXO_VICTIMA, df$ANIO_HECHO)


B <- 10000 # Número de remuestras para el Bootstrap y de simulaciones para Monte Carlo elegido

# Hiperparámetros de la distribución previa
a <- 1
b <- 1
out <- NULL
set.seed(2110)

# Ajuste de los modelos por cada año desde el 2012 hasta el 2022
for (k in 2012:2022) {
  # Datos observados en cada año para hombres y mujeres
  y_h <- as.numeric(df[df$SEXO_VICTIMA == 'MASCULINO' & (df$ANIO_HECHO == k),]$TOTAL_VICTIMAS)
  y_m <- as.numeric(df[df$SEXO_VICTIMA == 'FEMENINO' & (df$ANIO_HECHO == k),]$TOTAL_VICTIMAS)
  
  # Estadisticos suficientes y tamaño de muestra para los datos observados 
  n_h <- length(y_h)
  n_m <- length(y_m)
  s_h <- sum(y_h)
  s_m <- sum(y_m)

  rm(y_h, y_m)
  
  # Ajuste de los modelos y cálculo de eta por el método de Monte Carlo y por Bootstrap paramétrico
  eta_bayes <- (rgamma(B, shape = a + s_m, rate = b + n_m)/rgamma(B, shape = a + s_h, rate = b + n_h))-1
  eta_frec <- ((rpois(B, lambda = s_m)*(1/n_m))/(rpois(B, lambda = s_h)*(1/n_h)))-1
  
  # Cálculo de las estimaciones y intervalos de confianza/ credibilidad al 95% y al 99% 
  out <- rbind(out, c(mean(eta_bayes),
                      mean(eta_frec),
                      quantile(eta_bayes,probs = c(.025, .975)),
                      quantile(eta_frec, probs = c(.025, .975)),
                      quantile(eta_bayes, probs = c(.005, .995)), 
                      quantile(eta_frec, probs = c(.005, .995))
                      )
               )
  rm(eta_bayes, eta_frec)
}
rm(df, est, ic, n_h, n_m, s_h, s_m)
gc()

#                           Tabla con las estimaciones y intervalos de confianza/ credibilidad al 95% y al 99%
#---------------------------------------------------------------------------------------------------------------
colnames(out) <- c("Estimación Bayes", "Estimación Frec", "LI B95", "LS B95", "LI F95", "LS F95", "LI B99", "LS B99", "LI F99", "LS F99")
rownames(out) <- paste0("Año ", 2012:2022)
knitr::kable(x = out, digits = 3, align = "c")
#---------------------------------------------------------------------------------------------------------------

#                           Gráfico con las estimaciones y intervalos de confianza/ credibilidad al 95% y al 99%
#-----------------------------------------------------------------------------------------------------------------
plot(x = 1:nrow(out) - 0.05, y = out[,1],ylim = c(-0.5, 2.3), pch = 19, cex = 1.1, xaxt = "n", xlab = "# de Vicitmas por año", ylab = expression(eta), main = "Análisis Bayesiano y Frecuentista 2012 - 2022")

# Grafico con las estimaiones bayesianas y los intervalos de credibilidad al 95% y al 99% por cada año 
lines(x = 1:nrow(out) - 0.05, y = out[,1], type = "l")
abline(h = 0, col = 1, lty = 2)
segments(x0 = 1:nrow(out) - 0.05, y0 = out[,7], x1 = 1:nrow(out) - 0.05, y1 = out[,8])
segments(x0 = 1:nrow(out) - 0.05, y0 = out[,3], x1 = 1:nrow(out) - 0.05, y1 = out[,4], lwd = 4)

# Grafico con las estimaiones frecuentistas y los intervalos de confianza al 95% y al 99% por cada año 
lines(x = 1:nrow(out) + 0.05, y = out[,2], type = 'p', pch = 19, cex = 1.1, xaxt = "n", col = 2)
lines(x = 1:nrow(out) + 0.05, y = out[,2], type = "l", col = 2)
segments(x0 = 1:nrow(out) + 0.05, y0 = out[,9], x1 = 1:nrow(out) + 0.05, y1 = out[,10], col = 2)
segments(x0 = 1:nrow(out) + 0.05, y0 = out[,5], x1 = 1:nrow(out) + 0.05, y1 = out[,6], col = 2, lwd = 4)


axis(side = 1, at = 1:nrow(out), labels = rownames(out), las = 1)
legend("topright", legend = c("Estimación Bayesiana", "Estimación Frecuentista"), col = c(1, 2), lty = 1, lwd = 2, bty = "n")
#------------------------------------------------------------------------------------------------------------------