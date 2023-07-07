library(ggplot2)
library(FuzzyStatTra)
#Gráfica del fsbp
n <- 1:20
f1 <- 1/n
f2 <- floor((n+1)/2)/n
f3 <- (floor(0.3*n)+1)/n
f4 <- (floor(0.5*n)+1)/n
datos <- data.frame(n=n, f1=f1, f2=f2, f3=f3, f4=f4)

ggplot(data=datos, aes(x=n)) +
  geom_line(aes(y=f1, color='1')) +
  geom_line(aes(y=f2, color='2')) +
  geom_line(aes(y=f3, color='3')) +
  geom_line(aes(y=f4, color='4')) +
  theme_bw() +
  ylab('fsbp') +
  scale_x_continuous(n.breaks=10) +
  ylim(0,1) +
  scale_color_manual('', values = c('red', 'blue', 'darkgreen', 'orange'), labels=c('Aumann', 'Me/M-est', expression(''*beta*'=0.30'), expression(''*beta*'=0.50')))

#Función de densidad en el Caso 1
ggplot() +
  stat_function(fun = dbeta, args = c(1,1), n=1000, size=1) +
  xlab('') +
  ylab('') +
  theme_classic() +
  annotate(geom = 'text', parse=T, label=as.character(expression(beta(1,1))), x=0.2, y=1.25) +
  xlim(0,1) +
  ylim(0,1.5)

#Ejemplo de muestra del Caso 1
set.seed(195)
N <- 1000
p <- 1
q <- 1
total <- generacion(100*N, p, q, 2)
numeros_muestra <- sort(sample(x=1:(100*N), size=100, replace=FALSE))
A <- total[numeros_muestra,]
B <- contaminacion(N, A, numeros_muestra, p, q, 0.2, 50, 1)

#Datos que no han sido contaminados
fuzzy_orig <- data.frame(x = as.vector(t(B[A[,2]==B[,2],])), y = rep(c(0,1,1,0),80), f = rep(1:80, each=4))

#Datos que han sido contaminados
fuzzy_cont <- data.frame(x = as.vector(t(B[A[,2]!=B[,2],])), y = rep(c(0,1,1,0),20), f = rep(1:20, each=4))

#Representación gráfica
ggplot() +
  geom_path(data=fuzzy_orig, aes(x=x, y = y, group=f, color='1')) +
  geom_path(data=fuzzy_cont, aes(x=x, y = y, group=f, color='2')) +
  theme_classic() +
  scale_x_continuous(limits=c(0,1), expand = c(0, 0.0)) +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0)) +
  xlab('') +
  ylab('') +
  theme(plot.margin = unit(c(0.2,0.4,0.0,0.0),"cm"), legend.position = 'right', legend.title = element_blank()) +
  scale_color_manual('', values = c('azure4', 'red'), labels=c('No contaminado', 'Contaminado'))

#Función de densidad en el Caso 2
ggplot() +
  stat_function(fun = dbeta, args = c(10,10), n=1000, size=1) +
  xlab('') +
  ylab('') +
  theme_classic() +
  annotate(geom = 'text', parse=T, label=as.character(expression(beta(10,10))), x=0.2, y=1.25) +
  xlim(0,1) +
  ylim(0,3.75)

#Ejemplo de muestra del Caso 2
set.seed(195)
N <- 1000
p <- 10
q <- 10
total <- generacion(100*N, p, q, 1)
numeros_muestra <- sort(sample(x=1:(100*N), size=100, replace=FALSE))
A <- total[numeros_muestra,]
B <- contaminacion(N, A, numeros_muestra, p, q, 0.2, 50, 1)

#Datos que no han sido contaminados
fuzzy_orig <- data.frame(x = as.vector(t(B[A[,2]==B[,2],])), y = rep(c(0,1,1,0),80), f = rep(1:80, each=4))

#Datos que han sido contaminados
fuzzy_cont <- data.frame(x = as.vector(t(B[A[,2]!=B[,2],])), y = rep(c(0,1,1,0),20), f = rep(1:20, each=4))

#Representación gráfica
ggplot() +
  geom_path(data=fuzzy_orig, aes(x=x, y = y, group=f, color='1')) +
  geom_path(data=fuzzy_cont, aes(x=x, y = y, group=f, color='2')) +
  theme_classic() +
  scale_x_continuous(limits=c(0,1), expand = c(0, 0.0)) +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0)) +
  xlab('') +
  ylab('') +
  theme(plot.margin = unit(c(0.2,0.4,0.0,0.0),"cm"), legend.position = 'right', legend.title = element_blank()) +
  scale_color_manual('', values = c('azure4', 'red'), labels=c('No contaminado', 'Contaminado'))

#Función de densidad en el Caso 3
ggplot() +
  stat_function(fun = dbeta, args = c(10,1), n=1000, size=1) +
  xlab('') +
  ylab('') +
  theme_classic() +
  annotate(geom = 'text', parse=T, label=as.character(expression(beta(10,1))), x=0.2, y=1.25) +
  xlim(0,1) +
  ylim(0,10)

#Ejemplo de muestra del Caso 3
set.seed(195)
N <- 1000
p <- 10
q <- 1
total <- generacion(100*N, p, q, 1)
numeros_muestra <- sort(sample(x=1:(100*N), size=100, replace=FALSE))
A <- total[numeros_muestra,]
B <- contaminacion(N, A, numeros_muestra, p, q, 0.2, 50, 2)

#Datos que no han sido contaminados
fuzzy_orig <- data.frame(x = as.vector(t(B[A[,2]==B[,2],])), y = rep(c(0,1,1,0),80), f = rep(1:80, each=4))

#Datos que han sido contaminados
fuzzy_cont <- data.frame(x = as.vector(t(B[A[,2]!=B[,2],])), y = rep(c(0,1,1,0),20), f = rep(1:20, each=4))

#Representación gráfica
ggplot() +
  geom_path(data=fuzzy_orig, aes(x=x, y = y, group=f, color='1')) +
  geom_path(data=fuzzy_cont, aes(x=x, y = y, group=f, color='2')) +
  theme_classic() +
  scale_x_continuous(limits=c(0,1), expand = c(0, 0.0)) +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0)) +
  xlab('') +
  ylab('') +
  theme(plot.margin = unit(c(0.2,0.4,0.0,0.0),"cm"), legend.position = 'right', legend.title = element_blank()) +
  scale_color_manual('', values = c('azure4', 'red'), labels=c('No contaminado', 'Contaminado'))


#Generación de la muestra original y su muestra contaminada
set.seed(195)
N <- 1000
p <- 1
q <- 1
total <- generacion(100*N, p, q, 1)
numeros_muestra <- sort(sample(x=1:(100*N), size=100, replace=FALSE))
A <- total[numeros_muestra,]
B <- contaminacion(N, A, numeros_muestra, p, q, 0.2, 25, 2)

#Datos necesarios para la representación de la muestra original
fuzzy <- data.frame(x = as.vector(t(A)), y = rep(c(0,1,1,0),100), f = rep(1:100, each=4))

#Media de tipo Aumann de la muestra
media <- data.frame(x = as.vector(t(Mean(A))), y = c(0,1,1,0))

#Alfa-niveles utilizados
alfas <- seq(0,1, length.out=101)

#Mediana 1-norma de la muestra original
x_mediana <- c(Median1norm(A)[,2,1], rev(Median1norm(A)[,3,1]))
y_mediana <- c(alfas, rev(alfas))
mediana <- data.frame(x = x_mediana, y = y_mediana)

#Representación para la muestra original
ggplot() +
  geom_path(data=fuzzy, aes(x=x, y = y, group=f, color='1')) +
  geom_path(data=media, aes(x=x, y=y, color='2'), size=1.5) +
  geom_path(data=mediana, aes(x=x, y=y, color='3'), size=1.5) +
  theme_classic() +
  scale_x_continuous(limits=c(0,1), expand = c(0, 0.0)) +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0)) +
  xlab('') +
  ylab('') +
  theme(plot.margin = unit(c(0.2,0.4,0.0,0.0),"cm"), legend.position = 'right', legend.title = element_blank()) +
  scale_color_manual('', values = c('azure4', 'blue', 'green'), labels=c('Muestra', 'AUM', 'Me1'))

#Datos que no han sido contaminados
fuzzy_orig <- data.frame(x = as.vector(t(B[A[,2]==B[,2],])), y = rep(c(0,1,1,0),80), f = rep(1:80, each=4))

#Datos que han sido contaminados
fuzzy_cont <- data.frame(x = as.vector(t(B[A[,2]!=B[,2],])), y = rep(c(0,1,1,0),20), f = rep(1:20, each=4))

#Media de tipo Aumann de la muestra contaminada
media_cont <- data.frame(x = as.vector(t(Mean(B))), y = c(0,1,1,0))

#Mediana 1-norma de la muestra contaminada
x_mediana_cont <- c(Median1norm(B)[,2,1], rev(Median1norm(B)[,3,1]))
y_mediana_cont <- c(alfas, rev(alfas))
mediana_cont <- data.frame(x = x_mediana_cont, y = y_mediana_cont)

#Representación gráfica para la muestra contaminada
ggplot() +
  geom_path(data=fuzzy_orig, aes(x=x, y = y, group=f, color='1')) +
  geom_path(data=fuzzy_cont, aes(x=x, y = y, group=f, color='2')) +
  geom_path(data=media, aes(x=x, y=y, color='3'), size=1.5, linetype='dotted') +
  geom_path(data=mediana, aes(x=x, y=y, color='4'), size=1.5, linetype='dotted') +
  geom_path(data=media_cont, aes(x=x, y=y, color='5'), size=1.5) +
  geom_path(data=mediana_cont, aes(x=x, y=y, color='6'), size=1.5) +
  theme_classic() +
  scale_x_continuous(limits=c(0,1), expand = c(0, 0.0)) +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0)) +
  xlab('') +
  ylab('') +
  theme(plot.margin = unit(c(0.2,0.4,0.0,0.0),"cm"), legend.position = 'right', legend.title = element_blank()) +
  scale_color_manual('', values = c('azure4', 'red', 'blue', 'green', 'darkorchid4', 'orange'), labels=c('Original', 'Contaminado', 'AUM orig.', 'Me1 orig.', 'AUM cont.', 'Me1 cont.'))

#Cálculo de los mid del 0-nivel para ambas muestras
midA <- (A[,1]+A[,4])/2
midB <- (B[,1]+B[,4])/2

#Histograma del mid del 0-nivel para la muestra original
ggplot() +
  geom_histogram(data=data.frame(x=midA), aes(x=x), binwidth=0.1, center=0.05) +
  theme_classic() +
  xlab(expression('mid '*tilde(U)[0]*'')) +
  ylab('Frecuencia') +
  xlim(0,1) +
  theme(plot.margin = unit(c(0.2,0.4,0.0,0.0),"cm")) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle('Muestra original') +
  theme(plot.title = element_text(hjust=0.5))

#Histograma del mid del 0-nivel para la muestra contaminada
ggplot() +
  geom_histogram(data=data.frame(x=midB), aes(x=x), binwidth=0.1, center=0.05) +
  theme_classic() +
  xlab(expression('mid '*tilde(U)[0]*'')) +
  ylab('Frecuencia') +
  xlim(0,1) +
  theme(plot.margin = unit(c(0.2,0.4,0.0,0.0),"cm")) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle('Muestra contaminada') +
  theme(plot.title = element_text(hjust=0.5))
