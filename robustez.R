#Función para determinar la robustez frente a valores atípicos con distancias
#N: número de muestras que se generan
#p,q: parámetros de la distribución beta que determina mid U1
#medida: parámetro que indica la medida de localización que se va a usar
distancias <- function(N, p, q, medida){
  sumaB <- 0
  sumaC <- 0
  sumaD <- 0
  sumaE <- 0
  total <- generacion(100*N, p, q, 1)
  numeros <- 1:(100*N)
  for(i in 1:N){
    set.seed(i)
    # print(i)
    numeros_muestra <- sort(sample(x=numeros, size=100, replace=FALSE))
    numeros <- setdiff(numeros, numeros_muestra)
    A <- total[numeros_muestra,]
    B <- contaminacion(N, A, numeros_muestra, p, q, 0.2, 25, i)
    C <- contaminacion(N, A, numeros_muestra, p, q, 0.2, 50, i)
    D <- contaminacion(N, A, numeros_muestra, p, q, 0.4, 25, i)
    E <- contaminacion(N, A, numeros_muestra, p, q, 0.4, 50, i)
    if(medida==1){
      #Cálculos con la media tipo Aumann
      sumaB <- sumaB + DthetaphiTra(Mean(A),Mean(B),theta=1)
      sumaC <- sumaC + DthetaphiTra(Mean(A),Mean(C),theta=1)
      sumaD <- sumaD + DthetaphiTra(Mean(A),Mean(D),theta=1)
      sumaE <- sumaE + DthetaphiTra(Mean(A),Mean(E),theta=1)
    } else if(medida==2){
      #Cálculos con la mediana 1-norma
      sumaB <- sumaB + Dthetaphi(Median1norm(A),Median1norm(B),theta=1)
      sumaC <- sumaC + Dthetaphi(Median1norm(A),Median1norm(C),theta=1)
      sumaD <- sumaD + Dthetaphi(Median1norm(A),Median1norm(D),theta=1)
      sumaE <- sumaE + Dthetaphi(Median1norm(A),Median1norm(E),theta=1)
    } else if(medida==3){
      #Cálculos con la mediana wabl
      sumaB <- sumaB + Dthetaphi(Medianwabl(A),Medianwabl(B),theta=1)
      sumaC <- sumaC + Dthetaphi(Medianwabl(A),Medianwabl(C),theta=1)
      sumaD <- sumaD + Dthetaphi(Medianwabl(A),Medianwabl(D),theta=1)
      sumaE <- sumaE + Dthetaphi(Medianwabl(A),Medianwabl(E),theta=1)
    }else if(medida==4){
      #Cálculos con la M-estimación de localización Huber
      sumaB <- sumaB + Dthetaphi(M_estimation(A),M_estimation(B),theta=1)
      sumaC <- sumaC + Dthetaphi(M_estimation(A),M_estimation(C),theta=1)
      sumaD <- sumaD + Dthetaphi(M_estimation(A),M_estimation(D),theta=1)
      sumaE <- sumaE + Dthetaphi(M_estimation(A),M_estimation(E),theta=1)
    }else if(medida==5){
      # Cálculos con la M-estimación de localización Hampel
      sumaB <- sumaB + Dthetaphi(M_estimation2(A),M_estimation2(B),theta=1)
      sumaC <- sumaC + Dthetaphi(M_estimation2(A),M_estimation2(C),theta=1)
      sumaD <- sumaD + Dthetaphi(M_estimation2(A),M_estimation2(D),theta=1)
      sumaE <- sumaE + Dthetaphi(M_estimation2(A),M_estimation2(E),theta=1)
    }else if(medida==6){
      #Cálculos con la media recortada 0.3
      sumaB <- sumaB + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.3))),t(as.matrix(TrimmedMean(B,0.3))),theta=1)
      sumaC <- sumaC + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.3))),t(as.matrix(TrimmedMean(C,0.3))),theta=1)
      sumaD <- sumaD + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.3))),t(as.matrix(TrimmedMean(D,0.3))),theta=1)
      sumaE <- sumaE + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.3))),t(as.matrix(TrimmedMean(E,0.3))),theta=1)
    }else if(medida==7){
      #Cálculos con la media recortada 0.5
      sumaB <- sumaB + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.5))),t(as.matrix(TrimmedMean(B,0.5))),theta=1)
      sumaC <- sumaC + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.5))),t(as.matrix(TrimmedMean(C,0.5))),theta=1)
      sumaD <- sumaD + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.5))),t(as.matrix(TrimmedMean(D,0.5))),theta=1)
      sumaE <- sumaE + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.5))),t(as.matrix(TrimmedMean(E,0.5))),theta=1)
    }
  }
  if(medida==1){
    print('Cálculo con la media tipo Aumann')
  } else if(medida==2){
    print('Cálculo con la mediana 1-norma')
  } else if(medida==3){
    print('Cálculo con la mediana wabl')
  } else if(medida==4){
    print('Cálculo con la M-estimación Huber')
  } else if(medida==5){
    print('Cálculo con la M-estimación Hampel')
  } else if(medida==6){
    print('Cálculo con la media recortada 0.3')
  } else if(medida==7){
    print('Cálculo con la media recortada 0.5')
  }
  return(c(round(sumaB/N, 4), round(sumaC/N, 4), round(sumaD/N, 4), round(sumaE/N, 4)))
}
