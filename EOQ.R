# Pedro Lanzagorta

optimizaEOQ <- function(lambda,K,h,t=0){
  Qo = sqrt(2*K*lambda/h)
  T = Qo/lambda
  razon = t/T
  tao = floor(razon)-razon
  R = lambda*t
  
  costoInventario = h*Qo/2
  costoOrdenar = K*lambda/Qo
  return (Qo)
}