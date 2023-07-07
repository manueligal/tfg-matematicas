library(FuzzyStatTra)

#Función que genera los n*N datos siguiendo el modelo
#nN: n*N (tamaño de la muestra*número de muestras)
#p,q: parámetros de la distribución beta que determina midU1
#s: semilla para asegurar reproducibilidad
generacion <- function(nN,p,q,s){
  A <- matrix(0,nN,4)
  set.seed(s)
  
  for(i in 1:(0.05*nN)){
    A[i,] <- sort(rbeta(4,p,q))
  }
  
  for(i in (0.05*nN+1):(0.4*nN)){
    X1 <- rbeta(1,p,q)
    X2 <- runif(1,0,min(c(0.1,X1,1-X1)))
    X3 <- runif(1,0,min(c(0.2,X1-X2)))
    X4 <- runif(1,0,min(c(0.2,1-X1-X2)))
    A[i,] <- c(X1-X2-X3,X1-X2,X1+X2,X1+X2+X4)
  }
  
  for(i in (0.4*nN+1):nN){
    X1 <- rbeta(1,p,q)
    if(X1>=0.25 & X1<=0.75){
      X2 <- min(c(rexp(1,200),X1,1-X1))
    } else if(X1<0.25){
      X2 <- min(c(rexp(1,100+4*X1),X1,1-X1))
    } else{
      X2 <- min(c(rexp(1,500-4*X1),X1,1-X1))
    }
    if((X1-X2)>=0.25){
      X3 <- min(c(rgamma(1,4,100),X1-X2))
    } else{
      X3 <- min(c(rgamma(1,4,500+4*X1),X1-X2))
    }
    if((X1+X2)<=0.75){
      X4 <- min(c(rgamma(1,4,100),1-X1-X2))
    } else{
      X4 <- min(c(rgamma(1,4,500-4*X1),1-X1-X2))
    }
    A[i,] <- c(X1-X2-X3,X1-X2,X1+X2,X1+X2+X4)
  }
  return(A)
}

#Función que contamina una muestra
#N: número de muestras generadas
#muestra: muestra que se va a contaminar
#numeros_muestra: índices de las filas que se van a contaminar
#p,q: parámetros de la distribución beta
#cp: proporción de contaminación
#CD: nivel de contaminación
#: semilla para asegurar la reproducibilidad
contaminacion <- function(N, muestra, numeros_muestra, p, q, cp, CD, s){
  set.seed(s)
  cont <- muestra
  n_contaminados <- sort(sample(1:100, 100*cp, replace=FALSE))
  
  for(j in n_contaminados){
    if(numeros_muestra[j]<=0.05*N){
      cont[j,] <-sort(rbeta(4,p,q+CD))
    }else if(numeros_muestra[j]<=0.4*N){
      X1 <- rbeta(1,p,q+CD)
      X2 <- runif(1,0,min(c(0.1,X1,1-X1)))
      X3 <- runif(1,0,min(c(0.2,X1-X2)))
      X4 <- runif(1,0,min(c(0.2,1-X1-X2)))
      cont[j,] <- c(X1-X2-X3,X1-X2,X1+X2,X1+X2+X4)
    }else{
      X1 <- rbeta(1,p,q+CD)
      if(X1>=0.25 & X1<=0.75){
        X2 <- min(c(rexp(1,200),X1,1-X1))
      } else if(X1<0.25){
        X2 <- min(c(rexp(1,100+4*X1),X1,1-X1))
      } else{
        X2 <- min(c(rexp(1,500-4*X1),X1,1-X1))
      }
      if((X1-X2)>=0.25){
        X3 <- min(c(rgamma(1,4,100),X1-X2))
      } else{
        X3 <- min(c(rgamma(1,4,500+4*X1),X1-X2))
      }
      if((X1+X2)<=0.75){
        X4 <- min(c(rgamma(1,4,100),1-X1-X2))
      } else{
        X4 <- min(c(rgamma(1,4,500-4*X1),1-X1-X2))
      }
      cont[j,] <- c(X1-X2-X3,X1-X2,X1+X2,X1+X2+X4)
    }
  }
  return(cont)
}

#M-estimación con función de pérdida de Huber de la muestra fuzzy A
M_estimation <- function(A){
  est_inicial <- Median1norm(A)
  epsilon <- 10
  
  ATra <- TransfTra(A)
  
  a <- median(Dthetaphi(ATra, est_inicial, theta=1/3))
  if(a==0){
    dist <- Dthetaphi(ATra, est_inicial, theta=1/3)
    a <- min(dist[dist>0])
  }
  
  rho <- function(x){
    return(ifelse(x<=a, x^2/2, a*(x-a/2)))
  }
  
  phi <- function(x){
    return(ifelse(x<=a, 1, a/x))
  }
  
  J <- function(est){
    return(mean(rho(Dthetaphi(ATra, est, theta=1/3))))
  }
  
  estimacion <- est_inicial
  
  while(epsilon>1e-6){
    numeradores <- phi(Dthetaphi(ATra, estimacion, theta=1/3))
    w <- numeradores/sum(numeradores)
    estimacionprev <- estimacion
    estimacion[,2,] <- rowSums(sweep(TransfTra(A)[,2,],2,w,FUN="*"))
    estimacion[,3,] <- rowSums(sweep(TransfTra(A)[,3,],2,w,FUN="*"))
    epsilon <- abs(J(estimacion)-J(estimacionprev))/J(estimacionprev)
  }
  return(estimacion)
}

#M-estimación con función de pérdida de Hampel de la muestra fuzzy A
M_estimation2 <- function(A){
  est_inicial <- Median1norm(A)
  epsilon <- 10
  
  ATra <- TransfTra(A)
  
  distancias <- Dthetaphi(ATra, est_inicial, theta=1/3)
  a <- median(distancias)
  b <- as.numeric(quantile(distancias, 0.75))
  c <- as.numeric(quantile(distancias, 0.85))
  
  if(a==0){
    distancias <- distancias[distancias>0]
    a <- median(distancias)
    b <- as.numeric(quantile(distancias, 0.75))
    c <- as.numeric(quantile(distancias, 0.85))
  }
  
  rho <- function(x){
    return(ifelse(x<a, x^2/2, ifelse(x<b, a*(x-a/2), ifelse(x<c, a*(x-c)^2/(2*b-2*c) + a*(b+c-a)/2, a*(b+c-a)/2))))
  }
  
  phi <- function(x){
    return(ifelse(x<a, 1, ifelse(x<b, a/x, ifelse(x<c, a*(x-c)/(x*(b-c)), 0))))
  }
  
  J <- function(est){
    return(mean(rho(Dthetaphi(ATra, est, theta=1/3))))
  }
  
  estimacion <- est_inicial
  
  while(epsilon>1e-6){
    numeradores <- phi(Dthetaphi(ATra, estimacion, theta=1/3))
    if(sum(numeradores)==0){
      break
    }
    w <- numeradores/sum(numeradores)
    estimacionprev <- estimacion
    estimacion[,2,] <- rowSums(sweep(ATra[,2,],2,w,FUN="*"))
    estimacion[,3,] <- rowSums(sweep(ATra[,3,],2,w,FUN="*"))
    epsilon <- abs(J(estimacion)-J(estimacionprev))/J(estimacionprev)
  }
  return(estimacion)
}

#Media recortada con proporción c de la muestra fuzzy A
TrimmedMean <- function(A, c){
  return(apply(A, 2, mean, trim = c))
}

#Distancia para números trapezoidales. Versión simplificada
distTra <- function(U, V, theta){
  a1 <- U[1]
  b1 <- U[2]
  c1 <- U[3]
  d1 <- U[4]
  a2 <- V[1]
  b2 <- V[2]
  c2 <- V[3]
  d2 <- V[4]
  
  dist <- sqrt(((a1+c1)/2-(a2+c2)/2)^2+theta*((c1-a1)/2-(c2-a2)/2)^2+(1+theta)/3*(((b1-a1)/2-(b2-a2)/2)^2+
             ((d1-c1)/2-(d2-c2)/2)^2)+(1-theta)/3*(((b1-a1)/2-(b2-a2)/2)*((d1-c1)/2-(d2-c2)/2))+ ((b1-a1)/2-
             (b2-a2)/2)*(((a1+c1)/2-(a2+c2)/2)-theta*((c1-a1)/2-(c2-a2)/2))+((d1-c1)/2-(d2-c2)/2)*(((a1+c1)/2-
             (a2+c2)/2)+theta*((c1-a1)/2-(c2-a2)/2)))
  return(dist)
}
