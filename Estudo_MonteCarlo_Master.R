# ==============================================================================
# Estudo_MonteCarlo_NonSpatial.R
# Extração completa de Métricas Estáticas (ESS, Bias, etc) e Resumos Dinâmicos
# ==============================================================================
rm(list = ls())
library(nimble)
library(coda)
library(parallel)
library(dplyr)
library(readr)
library(stringr)

# 1. DIRETÓRIOS E CONFIGURAÇÕES
dir_projeto <- "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Simulacao/main/main2"
setwd(dir_projeto)

n_replicas <- 50
niter      <- 50000
nburnin    <- 5000
nchains    <- 1
thin       <- 10

dir.create(file.path(dir_projeto, "resultados_MC_NonSpatial"), showWarnings = FALSE)

# ==============================================================================
# Função do Trabalhador (Worker)
# ==============================================================================
worker_fun <- function(rep_indices, T_val, w_fix, niter, nburnin, nchains, thin, dir_base) {
  setwd(dir_base) 
  library(nimble)
  library(coda)
  
  pasta_dados <- file.path(dir_base, "dados_montecarlo")
  dados_base <- readRDS(file.path(pasta_dados, sprintf("replica_T%d_rep%02d.rds", T_val, 1)))
  C <- dados_base$constants
  A_val <- C$n_regions
  p <- C$p; K <- C$K
  
  constants_nimble <- list(
    n_regions = A_val, n_times = T_val, p = p, K = K,
    h = C$h, w = w_fix, a0 = 1.0, b0 = 1.0, mu_beta = rep(0, p),
    a_unif = 0.0, b_unif = 0.1
  )
  
  # --- APENAS MODELO NÃO-ESPACIAL ---
  code_nonspatial <- nimbleCode({
    for(j in 1:p) beta[j] ~ dnorm(mu_beta[j], sd=10)
    gamma[1] ~ dunif(min=a_unif, max=b_unif)
    for(j in 2:K) gamma[j] ~ dunif(min=0, max=(1 - sum(gamma[1:(j-1)])))
    for(i in 1:n_regions) {
      epsilon[i] <- 1 - inprod(h[i,1:K],gamma[1:K])
      for(t in 1:n_times) {
        lambda[i,t] ~ dgamma(1,1)
        log(mu[i,t]) <- log(lambda[i,t])+log(E[i,t])+log(epsilon[i])+inprod(beta[1:p],x[i,t,1:p])
        Y[i,t]        ~ dpois(mu[i,t])
      }
    }
  })
  
  data_nimble <- list(Y = dados_base$Y, E = dados_base$E, x = dados_base$x)
  
  ffbs_nonspatial <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_regions<-control$n_regions; n_times<-control$n_times; p<-control$p; a0<-control$a0; b0<-control$b0; w<-control$w
      buf_size<-n_regions*(n_times+1)
      at_buf<-nimNumeric(buf_size,0); bt_buf<-nimNumeric(buf_size,0)
      calcNodes<-model$getDependencies(target,self=FALSE)
      targetNodes<-model$expandNodeNames(target)
      setupOutputs(at_buf,bt_buf)
    },
    run = function() {
      declare(i,integer()); declare(t,integer()); declare(k,integer()); declare(prod_val,double())
      declare(g_it,double()); declare(att_t,double()); declare(btt_t,double())
      declare(nu,double()); declare(idx,integer()); declare(idx_next,integer())
      for(i in 1:n_regions){
        idx<-(i-1)*(n_times+1)+1; at_buf[idx]<<-a0; bt_buf[idx]<<-b0
        for(t in 1:n_times){
          idx<-(i-1)*(n_times+1)+t; idx_next<-idx+1
          att_t<-w*at_buf[idx]; btt_t<-w*bt_buf[idx]
          prod_val<-0
          for(k in 1:p) prod_val<-prod_val+model$x[i,t,k]*model$beta[k]
          g_it<-model$E[i,t]*model$epsilon[i]*exp(prod_val) 
          at_buf[idx_next]<<-att_t+model$Y[i,t]; bt_buf[idx_next]<<-btt_t+g_it
        }
        idx<-(i-1)*(n_times+1)+n_times+1
        model$lambda[i,n_times]<<-rgamma(1,shape=at_buf[idx],rate=bt_buf[idx])
        for(t_idx in 1:(n_times-1)){
          t_back<-n_times-t_idx; idx_buf<-(i-1)*(n_times+1)+t_back+1
          nu<-rgamma(1,shape=(1-w)*at_buf[idx_buf],rate=bt_buf[idx_buf])
          model$lambda[i,t_back]<<-nu+w*model$lambda[i,t_back+1]
        }
      }
      model$calculate(calcNodes)
      copy(from=model,to=mvSaved,row=1,nodes=targetNodes,logProb=TRUE)
    },
    methods=list(reset=function(){})
  )
  
  # Compilação Única
  Rmodel <- nimbleModel(code = code_nonspatial, constants = constants_nimble, data = data_nimble)
  Cmodel <- compileNimble(Rmodel)
  conf <- configureMCMC(Rmodel, thin = thin)
  conf$removeSamplers("lambda")
  conf$addSampler(target="lambda", type=ffbs_nonspatial, control=list(
    n_regions=A_val, n_times=T_val, p=p, a0=1.0, b0=1.0, w=w_fix
  ))
  conf$removeSampler("gamma")
  conf$addSampler(target="gamma", type="AF_slice")
  conf$addMonitors(c("beta", "gamma", "lambda"))
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  # Vetores verdadeiros para os parâmetros estáticos
  beta_true_vec  <- setNames(c(-1.0, 1.0, 0.5), paste0("beta[", 1:p, "]"))
  gamma_true_vec <- setNames(c(0.05, 0.10, 0.10, 0.15), paste0("gamma[", 1:K, "]"))
  true_statics   <- c(beta_true_vec, gamma_true_vec)
  
  res_estaticos_lista <- list()
  res_dinamicos_lista <- list()
  
  for(rep in rep_indices) {
    arquivo_rep <- file.path(pasta_dados, sprintf("replica_T%d_rep%02d.rds", T_val, rep))
    dados_rep <- readRDS(arquivo_rep)
    
    Cmodel$setData(list(Y = dados_rep$Y, E = dados_rep$E, x = dados_rep$x))
    Cmcmc$run(niter = niter, nburnin = nburnin)
    amostras <- as.matrix(Cmcmc$mvSamples)
    
    # --------------------------------------------------------------------------
    # AVALIAÇÃO DE PARÂMETROS ESTÁTICOS (com ESS)
    # --------------------------------------------------------------------------
    # Calcula o ESS para toda a cadeia de uma vez (apenas colunas de beta e gamma)
    cols_estaticas <- names(true_statics)
    ess_vals <- effectiveSize(as.mcmc(amostras[, cols_estaticas]))
    
    df_estaticos_rep <- do.call(rbind, lapply(cols_estaticas, function(nm) {
      samps <- amostras[, nm]
      tv <- true_statics[nm]
      me <- mean(samps)
      
      # Quantis para HPD rápido
      hpd <- HPDinterval(as.mcmc(samps), prob=0.95) 
      
      data.frame(
        Replica = rep,
        Parametro = nm,
        True = tv,
        Mean = me,
        Bias = me - tv,
        MSE = (me - tv)^2 + var(samps),
        Cov = as.integer(tv >= hpd[1] & tv <= hpd[2]),
        ESS = ess_vals[nm],
        stringsAsFactors = FALSE
      )
    }))
    res_estaticos_lista[[as.character(rep)]] <- df_estaticos_rep
    
    # --------------------------------------------------------------------------
    # AVALIAÇÃO GLOBAL DE VARIÁVEIS DINÂMICAS (Lambda)
    # --------------------------------------------------------------------------
    lambda_cols <- grep("^lambda\\[", colnames(amostras), value = TRUE)
    
    # Para economizar memória, ao invés de HPD de 3625 variáveis, usamos 
    # quantis empíricos que dão resultados 99% idênticos em cadeias bem comportadas
    lambda_means <- colMeans(amostras[, lambda_cols])
    lambda_lower <- apply(amostras[, lambda_cols], 2, quantile, probs = 0.025)
    lambda_upper <- apply(amostras[, lambda_cols], 2, quantile, probs = 0.975)
    
    # Extrair índices reais do array lambda_true salvo nos dados
    lt_true_vec <- as.vector(t(dados_rep$lambda_true)) # Alinhado por região e tempo
    
    # Construir tabela temporária para calcular agregados de forma fácil
    df_dyn <- data.frame(
      Mean = lambda_means, Lower = lambda_lower, Upper = lambda_upper, True = lt_true_vec
    )
    df_dyn$Bias <- df_dyn$Mean - df_dyn$True
    df_dyn$MSE <- df_dyn$Bias^2
    df_dyn$Cov <- ifelse(df_dyn$True >= df_dyn$Lower & df_dyn$True <= df_dyn$Upper, 1, 0)
    
    # Salvar apenas o resumo global e temporal para esta réplica
    res_dinamicos_lista[[as.character(rep)]] <- data.frame(
      Replica = rep,
      Lambda_Global_Bias = mean(df_dyn$Bias),
      Lambda_Global_MSE = mean(df_dyn$MSE),
      Lambda_Global_Cov = mean(df_dyn$Cov)
    )
    
    Cmcmc$mvSamples$resize(0) 
  }
  
  return(list(
    Estaticos = bind_rows(res_estaticos_lista),
    Dinamicos = bind_rows(res_dinamicos_lista)
  ))
}

# ==============================================================================
# Execução da Grade
# ==============================================================================
# Apenas modelo Non_spatial, variando W
scenarios <- expand.grid(
  T_val = c(25), # Ajuste para c(25, 100) se quiser rodar ambos
  w_fix = c(0.7, 0.9),
  stringsAsFactors = FALSE
)

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)

for(i in 1:nrow(scenarios)) {
  S <- scenarios[i, ]
  cat(sprintf("\n=== RODANDO CENÁRIO: T=%d | w=%.1f ===\n", S$T_val, S$w_fix))
  
  rep_chunks <- split(1:n_replicas, cut(1:n_replicas, n_cores, labels = FALSE))
  
  res_lista <- parLapply(cl, rep_chunks, worker_fun, 
                         T_val = S$T_val, w_fix = S$w_fix, 
                         niter = niter, nburnin = nburnin, nchains = nchains, thin = thin, 
                         dir_base = dir_projeto) 
  
  # Juntando as tabelas retornadas por todos os cores
  df_estaticos <- bind_rows(lapply(res_lista, function(x) x$Estaticos))
  df_dinamicos <- bind_rows(lapply(res_lista, function(x) x$Dinamicos))
  
  df_estaticos$T_val <- S$T_val
  df_estaticos$w_fix <- S$w_fix
  
  df_dinamicos$T_val <- S$T_val
  df_dinamicos$w_fix <- S$w_fix
  
  write_csv(df_estaticos, file.path(dir_projeto, "resultados_MC_NonSpatial", sprintf("Estaticos_T%d_w%.1f.csv", S$T_val, S$w_fix)))
  write_csv(df_dinamicos, file.path(dir_projeto, "resultados_MC_NonSpatial", sprintf("DinamicosGlobais_T%d_w%.1f.csv", S$T_val, S$w_fix)))
}

stopCluster(cl)
cat("\nSimulação concluída com sucesso!\n")