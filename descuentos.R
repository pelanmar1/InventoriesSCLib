# Pedro Lanzagorta
optimizaQDescuentos <- function(lambda,K,I,bc,inc=FALSE){
  
  costoEOQ <- function(lambda, K, I, cj, Q){
    return(cj*lambda + K*lambda/Q + 0.5*I*cj*Q)
  }
  
  descTodasUnidades <- function(lambda, K, I, bc){
    
    Q = c()
    GB = c()
    for (j in 1:nrow(bc)){
      # Variables de cada región
      cj = bc[j,1]
      bj = bc[j,2]
      bj1 = bc[j,3]
      # Cantidad óptima por región
      Qj = sqrt(2*K*lambda/(I*cj))
      # Si solución está es factible
      if((bj<= Qj) && (Qj <= bj1)){
        Q = append(Q,Qj)
      }
    }
    
    # Cantidad óptima general
    Qo = max(Q)
    indice = which(Q==Qo)
    cm = bc[indice,1]
    
    # Costo en cantidad óptima general
    GQo = costoEOQ(lambda,K,I,cm,Qo)
    
    # Buscamos en extremos de intervalos un costo menor
    for(j in (indice+1):nrow(bc)){
      for(k in 2:3){
        ck = bc[j,1]
        Qk = bc[j,k]
        Gk = costoEOQ(lambda,K,I,ck,Qk)
        if(Gk<GQo){
          GQo = Gk
          Qo = Qk
        }
      }
    }
    x = matrix(c(Qo,GQo),byrow=TRUE,nrow=1,ncol=2)
    colnames(x) <- c("Q*","G(Q*)")
    
    return(x)
  }
  
  buscaIntervalo <- function(Q,bc){
    sigue = TRUE
    pos = -1
    i = 1
    while(sigue && i <=nrow(bc)){
      sigue = !(bc[i,2]<=Q && bc[i,3]>Q)
      if(!sigue)
        pos = i
      i = i+1
    }
    return (pos)
  }
  
  cjIncremental <- function(Q,bc){
    bi = buscaIntervalo(Q,bc)
    dQ =  Q - bc[bi,2]
    G =0
    G = G + dQ*bc[bi,1]
    Kn = 0
    if(bi>1){
      for(i in 1:(bi-1)){
        G = G+ (bc[i,3]-bc[i,2])*bc[i,1]
      }
      
      Kn = G - Q*bc[bi,1]
    }
    x = matrix(c(G,Kn),byrow=TRUE,nrow=1,ncol=2)
    colnames(x) <- c("G(Q)*","dK")
    
    return(x)
  }
  

  
  
  descIncremental <- function(lambda,K,i,bc){
    
    
    Q = c()
    for (j in 1:nrow(bc)){
      # Variables de cada región
      cj = bc[j,1]
      bj = bc[j,2]
      bj1 = bc[j,3]
      Kj = K + cjIncremental(bj,bc)[2]
      
      # Cantidad óptima por región
      Qj = sqrt(2*Kj*lambda/(i*cj))
      # Si solución está es factible
      if((bj<= Qj) && (Qj <= bj1)){
        Q = append(Q,Qj)
      }
    }
    
    # Cantidad óptima general
    
    GQ = c()
    
    for(i in 1:length(Q)){
      Qo = Q[i]
      indice = buscaIntervalo(Qo,bc)
      cjKn = cjIncremental(Qo,bc)
      cm = cjKn[1]/Qo
      Kn = cjKn[2] + K
      
      GQo = costoEOQ(lambda,Kn,I,cm,Qo)
      GQ = append(GQ,GQo)
    }
    m =which(GQ==min(GQ))
    Qo = Q[m]
    GQo = GQ[m]
    
    
    x = matrix(c(Qo,GQo),byrow=TRUE,nrow=1,ncol=2)
    colnames(x) <- c("Q*","G(Q*)")
    
    return(x)
  }
  
  # main
  
  if(inc == FALSE){
      return (descTodasUnidades(lambda,K,i,bc))
  }else{
    
      return(descIncremental(lambda,K,i,bc))
  }
  
  
}

# EJEMPLO DE USO
#lambda = 600

#K = 8
#i = 0.2

#bc = matrix(c(0.30,0,500,0.29,500,1000,0.28,1000,Inf),nrow=3,ncol=3,byrow=TRUE)

# Descuento sobre todas las unidades
#optimizaQDescuentos(lambda,K,i,bc)

# Descuento incremental
#optimizaQDescuentos(lambda,K,i,bc,inc=TRUE)

lambda = 480

K = 150
i = 0.23

bc = matrix(c(1.25,0,700,1.05,700,Inf),nrow=2,ncol=3,byrow=TRUE)

optimizaQDescuentos(lambda,K,i,bc,inc=TRUE)

