dist.amostral <- function(n = c(10,5), B = 10, media = 0, dp = 1, plot = TRUE, seed = 0){
  set.seed(seed)
  if(length(n) == 1){
    res <- rep(0,B)
    for(i in 1:B){
      amostras[[i]] <- rnorm(n, mean = media, sd = dp)
      res[i] <- mean(amostras[[i]])
    }
    if(plot){
      plot(density(res), col = 'black',
           main = 'Densidade',
           ylim = c(0,3))
    }
  } else {
    res <- matrix(0, nrow = length(n), ncol = B)
    for(i in 1:nrow(res)){
      for(j in 1:B){
        amostra <- rnorm(i, mean = media, sd = dp)
        res[i,j] <- mean(amostra)
      }
    }
    if(plot){
      plot(density(res[1,]), col = 'black',
           main = 'Densidade',
           ylim = c(0,3))
      for(i in 2:nrow(res)){
        lines(density(res[i,]), col = sample(colours(),1))
      }
    }
  }
}

dist.amostral(n = c(1:10),seed = 50, media = 10)
