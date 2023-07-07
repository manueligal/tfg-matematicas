#Calcula la varianza de la muestra A para números fuzzy no trapezoidales
VarNoTra <- function(A){
  media <- array(rowMeans(A, dims = 2), c(101,3,1))
  return(mean(Dthetaphi(A, media, theta=1)^2))
}

#Matrices con los números fuzzy de las escalas LIK y FLS
puntos_likert <- matrix(c(rep(0,4), rep(1/3, 4), rep(2/3, 4), rep(1, 4)), byrow=T,4,4)
puntos_fls <- data.frame(a=c(0,0,1/3,2/3), b=c(0,1/3,2/3,1), c=c(0,1/3,2/3,1), d=c(1/3,2/3,1,1))

#Robustez frente a cambios de escala usando la varianza
#nrep: número de repeticiones
#est: parámetro que determina la medida de localización que se utiliza
varianzas <- function(nrep, est){
  sumas <- 0
  
  muestra <- generacion(100*nrep, 10, 1, 2)
  numeros <- 1:(100*nrep)
  
  for(i in 1:nrep){
    set.seed(i+200)
    # print(i)
    
    numeros_muestra <- sample(x=numeros, size=100, replace=FALSE)
    numeros <- setdiff(numeros, numeros_muestra)
    
    #Transformación a las distintas escalas
    FRS <- muestra[numeros_muestra,]
    vas <- apply(FRS, 1, mean)
    VAS <- cbind(vas, vas, vas, vas)
    lik <- apply(DthetaphiTra(FRS, puntos_likert, theta=1), 1, which.min)
    LIK <- (cbind(lik, lik, lik, lik)-1)/3
    FLS <- apply(t(apply(puntos_fls[lik,], 1, as.numeric)), 2, as.numeric)
    
    if(est==1){
      estimacion <- matrix(rbind(Mean(FRS), Mean(VAS), Mean(LIK), Mean(FLS)), 4, 4)
      sumas <- sumas + Var(estimacion, theta=1)
    } else if(est==2){
      estimacion <- array(c(Median1norm(FRS)[,,1], Median1norm(VAS)[,,1], Median1norm(LIK)[,,1], Median1norm(FLS)[,,1]), c(101, 3, 4))
      sumas <- sumas + VarNoTra(estimacion)
    } else if(est==3){
      estimacion <- array(c(Medianwabl(FRS)[,,1], Medianwabl(VAS)[,,1], Medianwabl(LIK)[,,1], Medianwabl(FLS)[,,1]), c(101, 3, 4))
      sumas <- sumas + VarNoTra(estimacion)
    } else if(est==4){
      estimacion <- array(c(M_estimation(FRS)[,,1], M_estimation(VAS)[,,1], M_estimation(LIK)[,,1], M_estimation(FLS)[,,1]), c(101, 3, 4))
      sumas <- sumas + VarNoTra(estimacion)
    } else if(est==5){
      estimacion <- array(c(M_estimation2(FRS)[,,1], M_estimation2(VAS)[,,1], M_estimation2(LIK)[,,1], M_estimation2(FLS)[,,1]), c(101, 3, 4))
      sumas <- sumas + VarNoTra(estimacion)
    } else if(est==6){
      estimacion <- matrix(rbind(TrimmedMean(FRS, 0.3), TrimmedMean(VAS, 0.3), TrimmedMean(LIK, 0.3), TrimmedMean(FLS, 0.3)), 4, 4)
      sumas <- sumas + Var(estimacion, theta=1)
    } else if(est==7){
      estimacion <- matrix(rbind(TrimmedMean(FRS, 0.5), TrimmedMean(VAS, 0.5), TrimmedMean(LIK, 0.5), TrimmedMean(FLS, 0.5)), 4, 4)
      sumas <- sumas + Var(estimacion, theta=1)
    }
  }
  print('Varianza media')
  print(round(100*sumas/nrep, 4))
}

#Robustez frente a cambios de escala usando distancias medias
#nrep: número de repeticiones
#est: parámetro que determina la medida de localización que se utiliza
distescalas <- function(nrep, est){
  sumas <- 0
  
  muestra <- generacion(100*nrep, 1, 1, 3)
  numeros <- 1:(100*nrep)
  
  for(i in 1:nrep){
    set.seed(i+400)
    # print(i)
    
    numeros_muestra <- sample(x=numeros, size=100, replace=FALSE)
    numeros <- setdiff(numeros, numeros_muestra)
    
    FRS <- muestra[numeros_muestra,]
    vas <- apply(FRS, 1, mean)
    VAS <- cbind(vas, vas, vas, vas)
    lik <- apply(DthetaphiTra(FRS, puntos_likert, theta=1), 1, which.min)
    LIK <- (cbind(lik, lik, lik, lik)-1)/3
    FLS <- apply(t(apply(puntos_fls[lik,], 1, as.numeric)), 2, as.numeric)
    
    if(est==1){
      estimacion <- matrix(rbind(Mean(FRS), Mean(VAS), Mean(LIK), Mean(FLS)), 4, 4)
      sep <- DthetaphiTra(estimacion, estimacion, theta=1)
    } else if(est==2){
      estimacion <- array(c(Median1norm(FRS)[,,1], Median1norm(VAS)[,,1], Median1norm(LIK)[,,1], Median1norm(FLS)[,,1]), c(101, 3, 4))
      sep <- Dthetaphi(estimacion, estimacion, theta=1)
    } else if(est==3){
      estimacion <- array(c(Medianwabl(FRS)[,,1], Medianwabl(VAS)[,,1], Medianwabl(LIK)[,,1], Medianwabl(FLS)[,,1]), c(101, 3, 4))
      sep <- Dthetaphi(estimacion, estimacion, theta=1)
    } else if(est==4){
      estimacion <- array(c(M_estimation(FRS)[,,1], M_estimation(VAS)[,,1], M_estimation(LIK)[,,1], M_estimation(FLS)[,,1]), c(101, 3, 4))
      sep <- Dthetaphi(estimacion, estimacion, theta=1)
    } else if(est==5){
      estimacion <- array(c(M_estimation2(FRS)[,,1], M_estimation2(VAS)[,,1], M_estimation2(LIK)[,,1], M_estimation2(FLS)[,,1]), c(101, 3, 4))
      sep <- Dthetaphi(estimacion, estimacion, theta=1)
    } else if(est==6){
      estimacion <- matrix(rbind(TrimmedMean(FRS, 0.3), TrimmedMean(VAS, 0.3), TrimmedMean(LIK, 0.3), TrimmedMean(FLS, 0.3)), 4, 4)
      sep <- DthetaphiTra(estimacion, estimacion, theta=1)
    } else if(est==7){
      estimacion <- matrix(rbind(TrimmedMean(FRS, 0.5), TrimmedMean(VAS, 0.5), TrimmedMean(LIK, 0.5), TrimmedMean(FLS, 0.5)), 4, 4)
      sep <- DthetaphiTra(estimacion, estimacion, theta=1)
    }
    sumas <- sumas + c(sep[1,2], sep[1,3], sep[1,4], sep[2,3], sep[2,4], sep[3,4])
  }
  print('Distancias')
  print(round(sumas/nrep, 4))
}

#Robustez frente a cambios de escala usando parámetros de corte
#nrep: número de repeticiones
#est: parámetro que determina la medida de localización que se utiliza
corte <- function(nrep, est){
  porcentajes5 <- rep(0,6)
  porcentajes10 <- rep(0,6)
  porcentajes15 <- rep(0,6)
  
  muestra <- generacion(100*nrep, 1, 1, 4)
  numeros <- 1:(100*nrep)
  
  for(i in 1:nrep){
    set.seed(i+600)
    # print(i)
    
    numeros_muestra <- sample(x=numeros, size=100, replace=FALSE)
    numeros <- setdiff(numeros, numeros_muestra)
    
    FRS <- muestra[numeros_muestra,]
    vas <- apply(FRS, 1, mean)
    VAS <- cbind(vas, vas, vas, vas)
    lik <- apply(DthetaphiTra(FRS, puntos_likert, theta=1), 1, which.min)
    LIK <- (cbind(lik, lik, lik, lik)-1)/3
    FLS <- apply(t(apply(puntos_fls[lik,], 1, as.numeric)), 2, as.numeric)
    
    if(est==1){
      estimacion <- matrix(rbind(Mean(FRS), Mean(VAS), Mean(LIK), Mean(FLS)), 4, 4)
      sep <- DthetaphiTra(estimacion, estimacion, theta=1)
    } else if(est==2){
      estimacion <- array(c(Median1norm(FRS)[,,1], Median1norm(VAS)[,,1], Median1norm(LIK)[,,1], Median1norm(FLS)[,,1]), c(101, 3, 4))
      sep <- Dthetaphi(estimacion, estimacion, theta=1)
    } else if(est==3){
      estimacion <- array(c(Medianwabl(FRS)[,,1], Medianwabl(VAS)[,,1], Medianwabl(LIK)[,,1], Medianwabl(FLS)[,,1]), c(101, 3, 4))
      sep <- Dthetaphi(estimacion, estimacion, theta=1)
    } else if(est==4){
      estimacion <- array(c(M_estimation(FRS)[,,1], M_estimation(VAS)[,,1], M_estimation(LIK)[,,1], M_estimation(FLS)[,,1]), c(101, 3, 4))
      sep <- Dthetaphi(estimacion, estimacion, theta=1)
    } else if(est==5){
      estimacion <- array(c(M_estimation2(FRS)[,,1], M_estimation2(VAS)[,,1], M_estimation2(LIK)[,,1], M_estimation2(FLS)[,,1]), c(101, 3, 4))
      sep <- Dthetaphi(estimacion, estimacion, theta=1)
    } else if(est==6){
      estimacion <- matrix(rbind(TrimmedMean(FRS, 0.3), TrimmedMean(VAS, 0.3), TrimmedMean(LIK, 0.3), TrimmedMean(FLS, 0.3)), 4, 4)
      sep <- DthetaphiTra(estimacion, estimacion, theta=1)
    } else if(est==7){
      estimacion <- matrix(rbind(TrimmedMean(FRS, 0.5), TrimmedMean(VAS, 0.5), TrimmedMean(LIK, 0.5), TrimmedMean(FLS, 0.5)), 4, 4)
      sep <- DthetaphiTra(estimacion, estimacion, theta=1)
    }
    porcentajes5 <- porcentajes5 + (c(sep[1,2], sep[1,3], sep[1,4], sep[2,3], sep[2,4], sep[3,4])>0.05)
    porcentajes10 <- porcentajes10 + (c(sep[1,2], sep[1,3], sep[1,4], sep[2,3], sep[2,4], sep[3,4])>0.10)
    porcentajes15 <- porcentajes15 + (c(sep[1,2], sep[1,3], sep[1,4], sep[2,3], sep[2,4], sep[3,4])>0.15)
  }
  print('Porcentajes')
  print(100*porcentajes5/nrep)
  print(100*porcentajes10/nrep, 4)
  print(100*porcentajes15/nrep, 4)
}
