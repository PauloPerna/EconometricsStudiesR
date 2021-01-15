# Encoding: UTF-8

rm(list=ls())
library(stats)
library(mfx)
library(wooldridge)

mfxBoot <- function(formula = affair ~ male + age + yrsmarr + kids + vryrel + vryhap, 
                     data = affairs, 
                     dist = 'probit',  # aceita-se probit ou logit
                     boot = 1000,
                     n_boot_sample = nrow(affairs), 
                    seed = 1312,
                    formula_nova = FALSE, # FALSE: usamos o cálculo dos efeitos marginais em aula (efeitos marginais nos x médios), 
                                          # TRUE: (efeitos marginais calculados nos y_hat médios)
                    alpha = 0.05
){
  #### 1. Rodamos o Modelo (aceita-se probit ou logit) ####
  modelo <- glm(formula, family = binomial(link = dist), data)
  
  #### 2.Calculamos os Efeitos Marginais ####
  #  Calculado de forma diferente da feita em aula, optei por essa forma para podermos rodar o código mesmo em distribuição logit
  #  o resultado é próximo mas não é o mesmo, seguem:
  # feito em aula:
  # (Intercept)         male          age      yrsmarr         kids       vryrel       vryhap 
  # -0.124616618  0.055184659 -0.006930849  0.015458645  0.055998959 -0.079915608 -0.143912126 
  # formula nova:
  # (Intercept)         male          age      yrsmarr         kids       vryrel       vryhap 
  # -0.121162428  0.053655021 -0.006738736  0.015030154  0.054446750 -0.077700464 -0.139923093 
  
  # As diferenças se devem ao ponto onde são calculadas, em aula é sobre a média dos X, abaixo é sobre a média dos y_hat,
  #  infelizmente perdi a fonte de onde retirei esta formula.
  if(formula_nova){
    pdf <- ifelse(dist == 'probit',
                  mean(dnorm(predict(modelo, type = 'link'))),
                  mean(dlogis(predict(modelo, type = 'link'))))
    efMarg <- pdf*coef(modelo)
  } else {
    x.bhat <- apply(cbind(1,modelo$model[2:ncol(modelo$model)]),2,mean)%*%modelo$coef
    if(dist == 'probit'){
      efMarg <- as.vector(dnorm(x.bhat))*modelo$coef
    } else {
      efMarg <- as.vector(dlogis(x.bhat))*modelo$coef
    }
  }
  
  #### 3. Bootstrap ####
  bootvals <- matrix(rep(NA,boot*length(coef(modelo))), nrow=boot)
  set.seed(seed)
  for(i in 1:boot){
    if(i%%100 == 0){ # debug
      print(paste0("Bootstrap Sample número ",i," de um total de ", boot))
    }
    # Bootstrapp sample
    samp1 <- data[sample(1:nrow(data),replace=T,n_boot_sample),]
    # Rodamos o modelo
    modelo1 <- glm(formula, family=binomial(link=dist),samp1)
    # Calculamos os efeitos marginais
    if(formula_nova){
      pdf1 <- ifelse(dist == 'probit',
                     mean(dnorm(predict(modelo1, type = 'link'))),
                     mean(dlogis(predict(modelo1, type = 'link'))))
      efMarg1 <- pdf1*coef(modelo1)
    } else {
      x.bhat1 <- apply(cbind(1,modelo1$model[2:ncol(modelo1$model)]),2,mean)%*%modelo1$coef
      if(dist == 'probit'){
        efMarg1 <- as.vector(dnorm(x.bhat1))*modelo1$coef
      } else {
        efMarg1 <- as.vector(dlogis(x.bhat1))*modelo1$coef
      }
    }
    bootvals[i,] <- efMarg1
  }
  
  #### 4. Calculamos os Intervalos de Confiança ####
  t.stat <- abs(qt(alpha/2, nrow(data)-1))
  sd <- apply(bootvals,2,sd)
  res <- rbind(efMarg,efMarg + sd*t.stat,efMarg - sd*t.stat)
  rownames(res) <- c('efMarg','Superior','Inferior')
  return(res)
}

### Exemplo + Plot dos intervalos de confianca calculados acima

# Calculo dos exemplos
#  exemplo 1
data("affairs")
mfxExemplo1 <- mfxBoot(formula_nova = TRUE) # valores default
mfx::probitmfx(formula = affair ~ male + age + yrsmarr + kids + vryrel + vryhap, 
               data = affairs)

#  exemplo 2
mfxExemplo2 <- mfxBoot(formula_nova = TRUE,
                       dist = 'logit') # utilizando a nova fórmula e distribuição logit
mfx::logitmfx(formula = affair ~ male + age + yrsmarr + kids + vryrel + vryhap, 
               data = affairs)

#  exemplo 3
data("alcohol")
mfxExemplo3 <- mfxBoot(formula = abuse ~ status + age + educ + beertax, # utilizando outro dataset
                       data = alcohol,
                       formula_nova = TRUE,
                       dist = 'probit',
                       boot = 500,
                       n_boot_sample = nrow(alcohol)/2,
                       alpha = 0.1)
mfx::probitmfx(formula = abuse ~ status + age + educ + beertax, 
               data = alcohol)

# Plot
library(ggplot2)
dataPlot1 <- as.data.frame(t(mfxExemplo1))
dataPlot1$V1 <- rownames(dataPlot1)

ggplot(dataPlot1, aes(efMarg, ymin = Inferior, ymax = Superior)) +
  scale_x_discrete('Variável') +
  scale_y_continuous('Efeito Marginal') +
  geom_errorbar(aes(x = V1)) +
  geom_point(aes(x = V1, y = efMarg)) +
  coord_flip()

dataPlot2 <- as.data.frame(t(mfxExemplo2))
dataPlot2$V1 <- rownames(dataPlot2)

ggplot(dataPlot2, aes(efMarg, ymin = Inferior, ymax = Superior)) +
  scale_x_discrete('Variável') +
  scale_y_continuous('Efeito Marginal') +
  geom_errorbar(aes(x = V1)) +
  geom_point(aes(x = V1, y = efMarg)) +
  coord_flip()

dataPlot3 <- as.data.frame(t(mfxExemplo3))
dataPlot3$V1 <- rownames(dataPlot3)

ggplot(dataPlot3, aes(efMarg, ymin = Inferior, ymax = Superior)) +
  scale_x_discrete('Variável') +
  scale_y_continuous('Efeito Marginal') +
  geom_errorbar(aes(x = V1)) +
  geom_point(aes(x = V1, y = efMarg)) +
  coord_flip()
