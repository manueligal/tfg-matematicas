#Funci�n para determinar la robustez frente a valores at�picos con distancias
#N: n�mero de muestras que se generan
#p,q: par�metros de la distribuci�n beta que determina mid U1
#medida: par�metro que indica la medida de localizaci�n que se va a usar
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
      #C�lculos con la media tipo Aumann
      sumaB <- sumaB + DthetaphiTra(Mean(A),Mean(B),theta=1)
      sumaC <- sumaC + DthetaphiTra(Mean(A),Mean(C),theta=1)
      sumaD <- sumaD + DthetaphiTra(Mean(A),Mean(D),theta=1)
      sumaE <- sumaE + DthetaphiTra(Mean(A),Mean(E),theta=1)
    } else if(medida==2){
      #C�lculos con la mediana 1-norma
      sumaB <- sumaB + Dthetaphi(Median1norm(A),Median1norm(B),theta=1)
      sumaC <- sumaC + Dthetaphi(Median1norm(A),Median1norm(C),theta=1)
      sumaD <- sumaD + Dthetaphi(Median1norm(A),Median1norm(D),theta=1)
      sumaE <- sumaE + Dthetaphi(Median1norm(A),Median1norm(E),theta=1)
    } else if(medida==3){
      #C�lculos con la mediana wabl
      sumaB <- sumaB + Dthetaphi(Medianwabl(A),Medianwabl(B),theta=1)
      sumaC <- sumaC + Dthetaphi(Medianwabl(A),Medianwabl(C),theta=1)
      sumaD <- sumaD + Dthetaphi(Medianwabl(A),Medianwabl(D),theta=1)
      sumaE <- sumaE + Dthetaphi(Medianwabl(A),Medianwabl(E),theta=1)
    }else if(medida==4){
      #C�lculos con la M-estimaci�n de localizaci�n Huber
      sumaB <- sumaB + Dthetaphi(M_estimation(A),M_estimation(B),theta=1)
      sumaC <- sumaC + Dthetaphi(M_estimation(A),M_estimation(C),theta=1)
      sumaD <- sumaD + Dthetaphi(M_estimation(A),M_estimation(D),theta=1)
      sumaE <- sumaE + Dthetaphi(M_estimation(A),M_estimation(E),theta=1)
    }else if(medida==5){
      # C�lculos con la M-estimaci�n de localizaci�n Hampel
      sumaB <- sumaB + Dthetaphi(M_estimation2(A),M_estimation2(B),theta=1)
      sumaC <- sumaC + Dthetaphi(M_estimation2(A),M_estimation2(C),theta=1)
      sumaD <- sumaD + Dthetaphi(M_estimation2(A),M_estimation2(D),theta=1)
      sumaE <- sumaE + Dthetaphi(M_estimation2(A),M_estimation2(E),theta=1)
    }else if(medida==6){
      #C�lculos con la media recortada 0.3
      sumaB <- sumaB + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.3))),t(as.matrix(TrimmedMean(B,0.3))),theta=1)
      sumaC <- sumaC + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.3))),t(as.matrix(TrimmedMean(C,0.3))),theta=1)
      sumaD <- sumaD + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.3))),t(as.matrix(TrimmedMean(D,0.3))),theta=1)
      sumaE <- sumaE + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.3))),t(as.matrix(TrimmedMean(E,0.3))),theta=1)
    }else if(medida==7){
      #C�lculos con la media recortada 0.5
      sumaB <- sumaB + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.5))),t(as.matrix(TrimmedMean(B,0.5))),theta=1)
      sumaC <- sumaC + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.5))),t(as.matrix(TrimmedMean(C,0.5))),theta=1)
      sumaD <- sumaD + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.5))),t(as.matrix(TrimmedMean(D,0.5))),theta=1)
      sumaE <- sumaE + DthetaphiTra(t(as.matrix(TrimmedMean(A,0.5))),t(as.matrix(TrimmedMean(E,0.5))),theta=1)
    }
  }
  if(medida==1){
    print('C�lculo con la media tipo Aumann')
  } else if(medida==2){
    print('C�lculo con la mediana 1-norma')
  } else if(medida==3){
    print('C�lculo con la mediana wabl')
  } else if(medida==4){
    print('C�lculo con la M-estimaci�n Huber')
  } else if(medida==5){
    print('C�lculo con la M-estimaci�n Hampel')
  } else if(medida==6){
    print('C�lculo con la media recortada 0.3')
  } else if(medida==7){
    print('C�lculo con la media recortada 0.5')
  }
  return(c(round(sumaB/N, 4), round(sumaC/N, 4), round(sumaD/N, 4), round(sumaE/N, 4)))
}
