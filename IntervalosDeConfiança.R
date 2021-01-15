# Encoding: UTF-8

# Função para cálculo dos betas
ols <- function(Y,X){
  # variáveis
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  n <- nrow(Y)
  # coeficientes (X'X)^(-1)(X'Y)
  betas  <- solve(t(X) %*% X) %*% t(X) %*% Y
  # resíduos: Y-Xb
  uhat    <- as.matrix(Y-X %*% betas)
  # graus de liberdade dos resíduos (n - k - 1)
  glr    <- n - ncol(X) - 1
  # soma dos quadrados dos resíduos
  sigshat <- as.numeric((t(uhat)%*%uhat)/glr)
  Vbetahat <- sigshat * solve(t(X)%*%X)
  # erros-padrão
  se     <- as.matrix(sqrt(diag(Vbetahat)))
  return(list(betas = betas,se = se))
}

# Função principal, plota os intervalos de confiança para os dois betas.
principal <- function(n,B,perc,beta0,beta1){
  CIs = array(data = rep(0,2*2*B),
              dim = c(B,2,2))   # CIs[observacao, lower/upper, beta0/beta1]
  
  p = qnorm(1 - (1 - perc)/2)
  
  for(j in 1:B){
    x = cbind(rep(1,n),rnorm(n,10,1))
    y = beta0 + beta1*x[,2] + rnorm(n,0,10)
    
    # intercepto
    CIs[j,1,1] = ols(y,x)$betas[1]-p*ols(y,x)$se[1] # lower
    CIs[j,2,1] = ols(y,x)$betas[1]+p*ols(y,x)$se[1] # upper
    
    # beta1
    CIs[j,1,2] = ols(y,x)$betas[2]-p*ols(y,x)$se[2] # lower
    CIs[j,2,2] = ols(y,x)$betas[2]+p*ols(y,x)$se[2] # upper
  }
  
  ### plot dos CIs do intercepto
  {
  # cor
  ID = CIs[1:B,2,1] < beta0 | CIs[1:B,1,1] > beta0
  
  colors = rep(gray(0.6),B)
  colors[ID] = 'red'
  
  plot(0, 
       xlim = c(min(CIs[1:B,1,1]), max(CIs[1:B,2,1])), 
       ylim = c(1, B), 
       ylab = "Sample", 
       xlab = expression(beta), 
       main = "Confidence Intervals - Intercepto")
  
  # draw reference line
  abline(v = beta0, lty = 2)
  
  # Adiciona CIs
  for(j in 1:B) {
    lines(c(CIs[j, 1, 1], CIs[j, 2, 1]), 
          c(j, j), 
          col = colors[j], 
          lwd = 2)
  }
  }
  ### plot dos CIs do beta1
  {
  # cor
  ID = CIs[1:B,2,2] < beta1 | CIs[1:B,1,2] > beta1
  
  colors = rep(gray(0.6),B)
  colors[ID] = 'red'
  
  plot(0, 
       xlim = c(min(CIs[1:B,1,2]), max(CIs[1:B,2,2])), 
       ylim = c(1, B), 
       ylab = "Sample", 
       xlab = expression(beta), 
       main = "Confidence Intervals - Beta1")
  
  # draw reference line
  abline(v = beta1, lty = 2)
  
  # Adiciona CIs
  for(j in 1:B) {
    lines(c(CIs[j, 1, 2], CIs[j, 2, 2]), 
          c(j, j), 
          col = colors[j], 
          lwd = 2)
  }
  }
}

principal(n = 100, # Tamanho da amostra
          B = 100, # Número de amostras
          perc = 0.95, # Nível de confiança
          beta0 = 1000, # Valor real do intercepto
          beta1 = 99 # Valor real do parâmetro de inclinação
          )
