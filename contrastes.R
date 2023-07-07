library(FuzzyStatTra)

#Números fuzzy asoaciados a la escala Likert y FLS
puntos_likert <- matrix(c(rep(0,4), rep(1/3, 4), rep(2/3, 4), rep(1, 4)), byrow=T,4,4)
puntos_fls <- data.frame(a=c(0,0,1/3,2/3), b=c(0,1/3,2/3,1), c=c(0,1/3,2/3,1), d=c(1/3,2/3,1,1))

#Test de hipótesis de igualdad de medias Aumann con muestras dependientes
#X, Y: matrices con las muestras fuzzy
#B: número de repeticiones bootstrap
testBmedia <- function(X, Y, B){
  n <- nrow(X)
  
  mediaX <- apply(X, 2, mean)
  mediaY <- apply(Y, 2, mean)
  
  num <- distTra(mediaX, mediaY, theta=1)^2
  
  dist <-rep(0,n)
  for (i in 1:n){
    dist[i] <- distTra(X[i,]+mediaY, Y[i,]+mediaX, theta=1)^2
  }
  den <- mean(dist)
  
  estadistico <- num/den
  
  #Población bootstrap
  Xboot <- matrix(nrow=n, ncol=4)
  Yboot <- matrix(nrow=n, ncol=4)
  
  for (i in 1:n){
    for (l in 1:4){
      Xboot[i,l] <- X[i,l]+mediaY[l]
      Yboot[i,l] <- Y[i,l]+mediaX[l]
    }
  }
  
  #Muestras bootstrap
  estadisticoboot <- c(length=B)
  
  for (b in 1:B){
    muestra <- sort(sample(1:n, n, replace=TRUE))
    
    muestrabootX <- Xboot[muestra,] 
    muestrabootY <- Yboot[muestra,] 
    
    mediabootX <- apply(muestrabootX, 2, mean)
    mediabootY <- apply(muestrabootY, 2, mean)
    
    
    numboot <- distTra(mediabootX, mediabootY, theta=1)^2
    
    dist <-rep(0,n)
    for (i in 1:n){
      dist[i] <- distTra(muestrabootX[i,]+mediabootY, muestrabootY[i,]+mediabootX, theta=1)^2
    }
    denboot <- mean(dist)
    estadisticoboot[b] <- numboot/denboot
    
  }
  
  pvalor <- mean(estadisticoboot>estadistico)
  return(pvalor)
}

#Aplicación del contraste de igualdad de medias a varias muestras
#n: tamaño muestral
#num_muestras: número de muestras generadas
#B: número de repeticiones bootstrap
contraste <- function(n, num_muestras, B){
  resultados <- matrix(nrow=num_muestras, ncol=3)
  total <- generacion(1000*n, 1, 1, 1)
  numeros <- 1:(1000*n)
  for(i in 1:num_muestras){
    set.seed(i+800)
    # print(i)
    numeros_muestra <- sample(x=numeros, size=n, replace=FALSE)
    numeros <- setdiff(numeros, numeros_muestra)
    FRS <- total[numeros_muestra,]
    vas <- apply(FRS, 1, mean)
    VAS <- cbind(vas, vas, vas, vas)
    lik <- apply(DthetaphiTra(FRS, puntos_likert, theta=1), 1, which.min)
    LIK <- (cbind(lik, lik, lik, lik)-1)/3
    FLS <- apply(t(apply(puntos_fls[lik,], 1, as.numeric)), 2, as.numeric)
    
    # resultados[i,] <- c(testBmedia(FRS, VAS, B), testBmedia(FRS, LIK, B), testBmedia(FRS, FLS, B))
    resultados[i,] <- c(testBmedia(FRS, VAS, B), testBmedia(FRS, LIK, B), testBmedia(FRS, FLS, B))
  }
  return(resultados)
}

#Ejemplos 
set.seed(25)
contraste(25, 20, 1000)

set.seed(50)
contraste(50, 20, 1000)