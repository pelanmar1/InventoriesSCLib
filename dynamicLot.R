# Pedro Lanzagorta
silverMeal <- function(D, K, h, n=1,titulo=TRUE){
  if(titulo){
    cat("Análisis por método Silver-Meal\n")
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
    cat("Análisis por método Least Unit Cost\n")
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
    cat("Análisis por método Part Period Balancing\n")
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



lambda = c(55,	70,	105,	120,	115,	95,	100,	75,	120,	75,	60,	45)
K =20
h = 0.4/12

lambda = c(580,	440,	288,	202,	150,	102,	68,	50,	38,	24,	15,	12)

silverMeal(lambda, K,h)
leastUnitCost(lambda,K,h)
partPeriodBalancing(lambda,K,h)

petersonSilverRule(lambda)