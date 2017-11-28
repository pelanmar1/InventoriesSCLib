# Pedro Lanzagorta
library(SCperf)

resuelveWW <- function(L,c,K,h){
  encuentraSecWW <- function(T,GR){
    ordenes = c()
    n = nrow(T)
    m = ncol(T)
    limite = m
    while(limite>0){
      colAct = T[,limite]
      colMin = min(colAct,na.rm=TRUE)
      renMin = which(T[,limite] == colMin, arr.ind=TRUE)
      costoAcum = sum(GR[renMin:limite]) 
      limite = renMin -1
      ordenes <- c(ordenes,c(renMin,costoAcum))
    }
    
    #ordenes <- c(ordenes,c(NA,"?"))
    
    ordenes
    R = matrix(ordenes,ncol = 2, byrow=TRUE)
    colnames(R) <- c("M","Q*")
    for(i in 1:m){
      if(any(R[,1]==i)==FALSE){
        R = rbind(R, c(i, 0))
        
      }
    }
    R = R[order(R[,1]),]
    
    
    return( R)
    
  }
  
  calculaCostosWW <- function(R,L,c,K,h){
    n = nrow(R)
    m = ncol(R)
    costo = 0
    In=0
    G = c()
    for(i in 1:n){
      costoM = 0
      M = R[i,1]
      Q = R[i,2]
      I = (In + Q - L[M])
      Km = K[M]
      if(Q==0){
        Km=0
      }
      costoM = I*h[M] + Km + c[M]*Q
      costo = costo + costoM
      In = I
      G <- c(G,costoM)
      
    }
    R = cbind(R, G)
    R = rbind(R,T=c(NA,sum(R[,2]),sum(R[,3])))
    
    return(R)
    
    
  }
  # Algoritmo de Wagner - Whitin
  result <- WW(L,K,h,method="forward")
  T=t(result$Solution)
  R = encuentraSecWW(T,L)
  S = calculaCostosWW(R,L,c,K,h)
  return(S)
}

#L = c(1500,100,700,1200,200,1700)
#c = c(10,10.5,10,11,11,11)
#K = c(150,150,150,200,200,200)
#h = c(1,1,1,1.5,1.5,1.5)

#L=  c(60,100,10,200,120,15)
#c = c(5,5,6,6,7,7)
#K = c(150 ,150 ,150,100 ,100,100 )
#h = c(1,1,2,2,2,2)

L = c(100,100,200,100,120,80)
c = c(3.47,3.47,3.47,3.47,3.47,3.47)
K = c(50,50,50,50,50,50)
h = c(1.735,1.735,1.735,1.735,1.735,1.735)


resuelveWW(L,c,K,h)