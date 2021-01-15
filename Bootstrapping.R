# Encoding: UTF-8
# Bootstrapping
rm(list=ls())

# Função para cálculo dos betas, é a mesma que está na tarefa anterior, resolvi aproveitá-la
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

ci.boot <- function(data, B, seed, alpha, beta){
  n <- nrow(data)
  estat <- rep(0, B)
  set.seed(seed)
  for(i in 1:B){
    rx <- sample(1:nrow(data), size=n, replace = TRUE)
    estat[i] <- ols(data[rx,1],data[rx,2:3])$betas[beta]
  }
  se.boot <- sd(estat)
  se.teorico <- ols(data[,1],data[,2:3])$se[beta]
  t.stat <- abs(qt(alpha/2, n-1))
  ci.se.boot <- c(mean(estat) - t.stat*se.boot, 
                  mean(estat) + t.stat*se.boot)
  ci.se.teorico <- c(ols(data[,1],data[,2:3])$betas[beta] - t.stat*se.teorico, 
                     ols(data[,1],data[,2:3])$betas[beta] + t.stat*se.teorico)
  ci.perc <- c(quantile(estat, alpha), quantile(estat, 1-alpha))
  res <- list(boot=ci.se.boot, teorico=ci.se.teorico, perc = ci.perc)
  return(res)
}

# dados para teste
n <- 100000 # número de observações na população
beta0 = 0 # beta0 da população
beta1 = 31 # beta1 da população
x = cbind(rep(1,n),rnorm(n,10,1))       # X
y = beta0 + beta1*x[,2] + rnorm(n,0,10) # Y da população
tam_amostra <- 1000 # tamanho da amostra original
data <- cbind(y,x)[sample(1:tam_amostra, replace = FALSE),] # amostra original

ci.boot(data,B = 100, seed = 777, alpha = 0.05, beta = 2) # beta = 1(para beta0) ou 2(para beta1)

