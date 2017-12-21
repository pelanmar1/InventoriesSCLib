
library(SCperf)

# MODELOS (Q,R), (s,S), (S,T) Y NIVELES DE SERVICIO

inventarios_QR_normal <- function(lambda,K,h,p,lambda_tau,sigma_tau,MAX_J=100,epsilon=0.01){
  # FUNCIONES
  
  # Funcion de perdida normal estandar
  f_perdida_norm_est = function(x){
    exp(-(x^2)/2)/(sqrt(2*pi))-x*(1-pnorm(x))
  }
  
  imprime_reporte = function(i,q,r,f){
    print("------------------------------------------")
    print(paste("Q",i," - ",q))
    print(paste("R",i," - ",r))
    print(paste("faltantes",i," - ",f))
    print("------------------------------------------")
  }
  
  
  # Optimización
  j=0
  sigue = TRUE
  # epsilon es la diferencia mínima entre los valores de t y t-1 a partir de la cual
  # se consideran como iguales
  epsilon = 0.01
  
  # El número máximo de iteraciones
  MAX_J = 100
  Qo = sqrt(2*lambda*K/h)
  z = qnorm(1-Qo*h/(p * lambda))
  Ro = lambda_tau +sigma_tau*z
  faltantes = sigma_tau*f_perdida_norm_est(z)
  
  print('Iniciando heurístico para encontrar óptimo (Q,R)')
  while(sigue || j>MAX_J){
    imprime_reporte(j,Qo,Ro,faltantes)
    Q1 = sqrt((2*lambda)*(K+p*faltantes)/h)
    z =qnorm(1-Q1*h/(p * lambda))
    R1 = lambda_tau +sigma_tau*z
    faltantes = sigma_tau*f_perdida_norm_est(z)
    
    difQ= abs(Q1 - Qo)
    difR = abs(R1 - Ro)
    
    if(difQ <= epsilon || difR <= epsilon)
      sigue = FALSE
    
    Qo = Q1
    Ro = R1
    j= j+1
    
  }
  cat(sprintf("Política optima encontrada (Q*,R*) = (%s,%s)\n", round(Qo),round(Ro)))
  cat(sprintf("Óptimo alcanzado en %f iteraciones con un criterio de epsilon = %f\n",j,epsilon))
  B = 1- faltantes/Qo
  m = matrix(c(round(Qo),round(Ro),faltantes,B),nrow = 1,ncol = 4)
  colnames(m) <-c('Q*','R*','n(R)','B')
  
  return(m)
  
}

inventarios_QR_uniforme <- function(lambda,K,h,p,a,b,MAX_J=100,epsilon=0.01){
  # FUNCIONES
  
  # Funcion de perdida uniforme
  f_perdida_uniforme <- function(x,a,b){
    return(1/(b-a)*(0.5*(b^2+x^2)-b*x))
  }
  
  f_inv_uniforme <- function(r,a,b){
    return(a+r*(b-a))
  }
  
  imprime_reporte = function(i,q,r,f){
    print("------------------------------------------")
    print(paste("Q",i," - ",q))
    print(paste("R",i," - ",r))
    print(paste("faltantes",i," - ",f))
    print("------------------------------------------")
  }
  
  
  # Optimización
  j=0
  sigue = TRUE
  # epsilon es la diferencia mínima entre los valores de t y t-1 a partir de la cual
  # se consideran como iguales
  epsilon = 0.01
  
  # El número máximo de iteraciones
  MAX_J = 100
  Qo = sqrt(2*lambda*K/h)
  Ro = f_inv_uniforme(1-Qo*h/(p * lambda),a,b)
  faltantes = f_perdida_uniforme(Ro,a,b)
  
  print('Iniciando heurístico para encontrar óptimo (Q,R)')
  while(sigue || j>MAX_J){
    imprime_reporte(j,Qo,Ro,faltantes)
    Q1 = sqrt((2*lambda)*(K+p*faltantes)/h)
    R1 = f_inv_uniforme(1-Q1*h/(p * lambda),a,b)
    faltantes = f_perdida_uniforme(R1,a,b)
    
    difQ= abs(Q1 - Qo)
    difR = abs(R1 - Ro)
    
    if(difQ <= epsilon || difR <= epsilon)
      sigue = FALSE
    
    Qo = Q1
    Ro = R1
    j= j+1
    
  }
  cat(sprintf("Política optima encontrada (Q*,R*) = (%s,%s)\n", round(Qo),round(Ro)))
  cat(sprintf("Óptimo alcanzado en %f iteraciones con un criterio de epsilon = %f\n",j,epsilon))
  
  B = 1- faltantes/Qo
  m = matrix(c(round(Qo),round(Ro),faltantes,B),nrow = 1,ncol = 4)
  colnames(m) <-c('Q*','R*','n(R)','B')
  
  return(m)
}

inventarios_QR_exponencial <- function(lambda,K,h,p,MAX_J=100,epsilon=0.01){
  # FUNCIONES
  
  # Funcion de perdida uniforme
  f_perdida_exp <- function(x,lambda){
    return(1/lambda*exp(-lambda*x))
  }
  
  f_inv_exp <- function(r,lambda){
    return(-1/lambda*log(1-r))
  }
  
  imprime_reporte = function(i,q,r,f){
    print("------------------------------------------")
    print(paste("Q",i," - ",q))
    print(paste("R",i," - ",r))
    print(paste("faltantes",i," - ",f))
    print("------------------------------------------")
  }
  
  
  # Optimización
  j=0
  sigue = TRUE
  # epsilon es la diferencia mínima entre los valores de t y t-1 a partir de la cual
  # se consideran como iguales
  epsilon = 0.01
  
  # El número máximo de iteraciones
  MAX_J = 100
  Qo = sqrt(2*lambda*K/h)
  Ro = f_inv_exp(1-Qo*h/(p * lambda),lambda)
  faltantes = f_perdida_exp(Ro,lambda)
  
  print('Iniciando heurístico para encontrar óptimo (Q,R)')
  while(sigue || j>MAX_J){
    imprime_reporte(j,Qo,Ro,faltantes)
    Q1 = sqrt((2*lambda)*(K+p*faltantes)/h)
    R1 = f_inv_exp(1-Q1*h/(p * lambda),lambda)
    faltantes = f_perdida_exp(R1,lambda)
    
    difQ= abs(Q1 - Qo)
    difR = abs(R1 - Ro)
    
    if(difQ <= epsilon || difR <= epsilon)
      sigue = FALSE
    
    Qo = Q1
    Ro = R1
    j= j+1
    
  }
  cat(sprintf("Política optima encontrada (Q*,R*) = (%s,%s)\n", round(Qo),round(Ro)))
  cat(sprintf("Óptimo alcanzado en %f iteraciones con un criterio de epsilon = %f\n",j,epsilon))
  
  B = 1- faltantes/Qo
  m = matrix(c(round(Qo),round(Ro),faltantes,B),nrow = 1,ncol = 4)
  colnames(m) <-c('Q*','R*','n(R)','B')
  
  return(m)
}

inventarios_Ss_normal <- function(lambda,K,h,p,lambda_tau,sigma_tau,MAX_J=100,epsilon=0.01){
  m = inventarios_QR_normal(lambda,K,h,p,lambda_tau,sigma_tau,MAX_J=100,epsilon=0.01)
  s = m[1,2]
  S = m[1,2]+m[1,1]
  r = matrix(c(s,S),nrow = 1,ncol = 2)
  colnames(r) <- c('s*','S*')
  return(r)
  
}

inventarios_Ss_uniforme <- function(lambda,K,h,p,a,b,MAX_J=100,epsilon=0.01){
  m = inventarios_QR_normal(lambda,K,h,p,a,b,MAX_J=100,epsilon=0.01)
  s = m[1,2]
  S = m[1,2]+m[1,1]
  r = matrix(c(s,S),nrow = 1,ncol = 2)
  colnames(r) <- c('s*','S*')
  return(r)
}

inventarios_Ss_exponencial <- function(lambda,K,h,p,MAX_J=100,epsilon=0.01){
  m = inventarios_QR_exponencial(lambda,K,h,p,MAX_J=100,epsilon=0.01)
  s = m[1,2]
  S = m[1,2]+m[1,1]
  r = matrix(c(s,S),nrow = 1,ncol = 2)
  colnames(r) <- c('s*','S*')
  return(r)
}

inventarios_ST_normal_Alfa <- function(lambda, K,h,p,tau,sigma,alfa){
  T = sqrt(2*K/(h*lambda))
  z = qnorm(alfa)
  s = z*sigma*sqrt(T+tau)
  S = lambda*(T+tau) + s
  
  m = matrix(c(s,S,T),nrow = 1,ncol = 3)
  colnames(m) <- c("s","S","T")
  return(m)
  
}

inventarios_ST_normal_Beta <- function(lambda, K,h,p,tau,sigma,beta){
  T = sqrt(2*K/(h*lambda))
  L = (1-beta)*lambda*(T+tau)/(sigma*sqrt(T+tau))
  z = qnorm(L)
  s = z*sigma*sqrt(T+tau)
  S = lambda*(T+tau) + s
  
  
  
  m = matrix(c(s,S,T),nrow = 1,ncol = 3)
  colnames(m) <- c("s","S","T")
  return(m)
  
}

inventarios_alfa_normal <- function(alfa,miu,sigma){
  return(qnorm(alfa*sigma+miu,miu,sigma))
}

inventarios_beta_normal(lambda,K,h,p,lambda_tau,sigma_tau,MAX_J=100,epsilon=0.01){
  r = inventarios_QR_normal(lambda,K,h,p,lambda_tau,sigma_tau,MAX_J=100,epsilon=0.01)
  return(r[1,4])
}

inventarios_alfa_exponencial <- function(alfa,lambda){
  return(-1/lambda*log(1-alfa))
}

inventarios_beta_exponencial(lambda,K,h,p,MAX_J=100,epsilon=0.01){
  r = inventarios_QR_exponencial(lambda,K,h,p,MAX_J=100,epsilon=0.01)
  return(r[1,4])
}

inventarios_alfa_uniforme <- function(alfa,a,b){
  return(a+alfa*(b-a))
}

inventarios_beta_uniforme(lambda,K,h,p,a,b,MAX_J=100,epsilon=0.01){
  r = inventarios_QR_uniforme(lambda,K,h,p,a,b,MAX_J=100,epsilon=0.01)
  return(r[1,4])
}

# DESCUENTOS
optimizaQDescuentos <- function(lambda,K,I,bc,inc=FALSE){
  
  costoEOQ <- function(lambda, K, I, cj, Q){
    return(cj*lambda + K*lambda/Q + 0.5*I*cj*Q)
  }
  
  descTodasUnidades <- function(lambda, K, I, bc){
    
    Q = c()
    GB = c()
    for (j in 1:nrow(bc)){
      # Variables de cada regi?n
      cj = bc[j,1]
      bj = bc[j,2]
      bj1 = bc[j,3]
      # Cantidad ?ptima por regi?n
      Qj = sqrt(2*K*lambda/(I*cj))
      # Si soluci?n est? es factible
      if((bj<= Qj) && (Qj <= bj1)){
        Q = append(Q,Qj)
      }
    }
    
    # Cantidad ?ptima general
    Qo = max(Q)
    indice = which(Q==Qo)
    cm = bc[indice,1]
    
    # Costo en cantidad ?ptima general
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
      # Variables de cada regi?n
      cj = bc[j,1]
      bj = bc[j,2]
      bj1 = bc[j,3]
      Kj = K + cjIncremental(bj,bc)[2]
      
      # Cantidad ?ptima por regi?n
      Qj = sqrt(2*Kj*lambda/(i*cj))
      # Si soluci?n est? es factible
      if((bj<= Qj) && (Qj <= bj1)){
        Q = append(Q,Qj)
      }
    }
    
    # Cantidad ?ptima general
    
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
  

  if(inc == FALSE){
    return (descTodasUnidades(lambda,K,i,bc))
  }else{
    
    return(descIncremental(lambda,K,i,bc))
  }
  
  
}

# WAGNER AUTOMATIZADO

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

# DINAMIC LOT SIZING
silverMeal <- function(D, K, h, n=1,titulo=TRUE){
  if(titulo){
    cat("An?lisis por m?todo Silver-Meal\n")
  }
  limite = length(D)
  i = 1
  suma = K
  cmTemp = suma
  if (i == limite){
    Q = D[i]
    cat(sprintf("Q%g = %f\n",n,Q))
  }else{
    repeat{
      suma = suma + h*(i)*D[i+1]
      cm = suma/(i+1)
      if(cm>=cmTemp){
        Q = sum(D[1:i])
        cat(sprintf("Q%g = %f\n",n,Q))
        nD = D[(i+1):limite]
        silverMeal(nD,K,h,n+i,titulo=FALSE)
        break
      }else if(i+1 ==limite){
        Q = sum(D[1:(i+1)])
        cat(sprintf("Q%g = %f\n",n,Q))
        break
      }
      i = i+1
      cmTemp = cm
      
      
    }
  }
}

leastUnitCost <- function(D, K, h, n=1,titulo=TRUE){
  if(titulo){
    cat("An?lisis por m?todo Least Unit Cost\n")
  }
  limite = length(D)
  i = 1
  suma = K
  cmTemp = suma
  if (i == limite){
    Q = D[i]
    cat(sprintf("Q%g = %f\n",n,Q))
  }else{
    repeat{
      suma = suma + h*(i)*D[i+1]
      cm = suma/(sum(D[1:(i+1)]))
      if(cm>=cmTemp){
        Q = sum(D[1:i])
        cat(sprintf("Q%g = %f\n",n,Q))
        nD = D[(i+1):limite]
        leastUnitCost(nD,K,h,n+i,titulo=FALSE)
        break
      }else if(i+1 ==limite){
        Q = sum(D[1:(i+1)])
        cat(sprintf("Q%g = %f\n",n,Q))
        break
      }
      i = i+1
      cmTemp = cm
      
      
    }
  }
}

partPeriodBalancing <- function(D,K,h,n=1,titulo=TRUE){
  if(titulo){
    cat("An?lisis por m?todo Part Period Balancing\n")
  }
  limite = length(D)
  i = 1
  suma = 0
  cmTemp = suma*h
  if (i == limite){
    Q = D[i]
    cat(sprintf("Q%g = %f\n",n,Q))
  }else{
    repeat{
      suma = suma + (i)*D[i+1]
      cm = h*suma
      if(cmTemp <= K && cm>=K){
        if(abs(cmTemp-K)<abs(cm-K)){
          l = 0
        }else
          l=1
        Q = sum(D[1:(i+l)])
        cat(sprintf("Q%g = %f\n",n,Q))
        nD = D[(i+1+l):limite]
        partPeriodBalancing(nD,K,h,n+i,titulo=FALSE)
        break
      }else if(i+1 ==limite){
        Q = sum(D[1:(i+1)])
        cat(sprintf("Q%g = %f\n",n,Q))
        break
      }
      i = i+1
      cmTemp = cm
      
      
    }
  }
}

petersonSilverRule <- function(D){
  # if (V < 0.25)
  #     -> use EOQ model with ?? as the demand estimate 
  # else
  #     -> use dynamic lot sizing rules
  
  V = sum(lambda^2)*length(D)/(sum(D)^2)-1
  return(V)
}

# MODELOS NEWS VENDOR
newsboy_exponencial <- function(lambda, co,cu){
  tc = cu/(cu+co)
  
  Qo = -log(1-tc)/lambda
  
  # Integrando
  #integral = function(x){
  #  lambda*exp(-lambda*x)
  #}
  
  
  #for(demanda in a:b){
  #  valor = integrate(integral, a, demanda)
  #  print(paste(demanda, valor$value))
  #}
  
  
  sobrante = function(x){
    (Qo - x)*lambda*exp(-lambda*x)
  }
  valor_sobrante = integrate(sobrante,0,Qo)
  
  faltantes = function(x){
    (x-Qo)*lambda*exp(-lambda*x)
  }
  
  valor_faltante = integrate(faltantes,Qo,Inf)
  
  costo_total = valor_sobrante$value*co + valor_faltante$value*cu 
  
  m = matrix(c(Qo, valor_faltante,valor_sobrante,costo_total),nrow = 1,ncol = 4)
  colnames(m) <- c("Q*","n(Q*)","n_(Q*)","Costo")
  return(m)
}

newsboy_normal <- function(miu,sigma, co,cu){
  tc = cu/(cu+co)
  
  Qo = qnorm(tc)
  Qo = Qo*sigma+miu
  
  sobrante = function(x){
    (Qo - x)*dnorm(x,miu,sigma)
  }
  valor_sobrante = integrate(sobrante,-Inf,Qo)
  
  faltantes = function(x){
    (x-Qo)*dnorm(x,miu,sigma)
  }
  
  valor_faltante = integrate(faltantes,Qo,Inf)
  
  costo_total = valor_sobrante$value*co + valor_faltante$value*cu 
  
  m = matrix(c(Qo, valor_faltante$value,valor_sobrante$value,costo_total),nrow = 1,ncol = 4)
  colnames(m) <- c("Q*","n(Q*)","n_(Q*)","Costo")
  return(m)
}

newsboy_uniforme <- function(a,b, co,cu){
  tc = cu/(cu+co)
  
  Qo = a+tc*(b-a)
  
  sobrante = function(x){
    (Qo - x)/(b-a)
  }
  valor_sobrante = integrate(sobrante,a,Qo)
  
  faltantes = function(x){
    (x-Qo)/(b-a)
  }
  
  valor_faltante = integrate(faltantes,Qo,b)
  
  costo_total = valor_sobrante$value*co + valor_faltante$value*cu 
  
  m = matrix(c(Qo, valor_faltante,valor_sobrante,costo_total),nrow = 1,ncol = 4)
  colnames(m) <- c("Q*","n(Q*)","n_(Q*)","Costo")
  return(m)
}




