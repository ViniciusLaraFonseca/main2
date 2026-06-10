# ==============================================================================
# Estudo_MC_CustoDinamica_OTIMIZADO_COM_LAMBDA.R
# ==============================================================================

library(nimble)
library(coda)
library(parallel)
library(dplyr)
library(readr)
library(stringr)

# ==============================================================================
# 1. CONFIGURAÇÕES
# ==============================================================================
dir_projeto <- "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Simulacao/main/main2"
setwd(dir_projeto)

n_replicas <- 50
niter      <- 50000
nburnin    <- 5000
thin       <- 10
w_fit      <- 0.8

# ==============================================================================
# 2. FUNÇÕES AUXILIARES
# ==============================================================================

compute_waic <- function(loglik_matrix) {
  lppd   <- sum(log(colMeans(exp(loglik_matrix))))
  p_waic <- sum(apply(loglik_matrix, 2, var))
  waic   <- -2 * (lppd - p_waic)
  list(lppd = lppd, p_waic = p_waic, waic = waic)
}

# ==============================================================================
# 3. WORKER
# ==============================================================================
worker_fun <- function(rep_indices, T_val, dgp, w_fit,
                       niter, nburnin, thin, dir_base) {
  
  library(nimble); library(coda); library(dplyr); library(stringr)
  
  pasta_dados <- file.path(dir_base, sprintf("dados_MC_Custo_DGP%d", dgp))
  
  dados_tmpl <- readRDS(
    file.path(pasta_dados, sprintf("replica_T%d_rep%02d.rds", T_val, rep_indices[1]))
  )
  
  C     <- dados_tmpl$constants
  A_val <- C$n_regions
  p     <- C$p
  K     <- C$K
  
  # ================= MODELOS =================
  
  code_ffbs <- nimbleCode({
    for(j in 1:p) beta[j] ~ dnorm(0, sd = 10)
    
    gamma[1] ~ dunif(0, 0.1)
    for(j in 2:K) gamma[j] ~ dunif(0, 1 - sum(gamma[1:(j-1)]))
    
    for(i in 1:n_regions) {
      epsilon[i] <- 1 - inprod(h[i,1:K], gamma[1:K])
      
      for(t in 1:n_times) {
        lambda[i,t] ~ dgamma(1,1)
        
        log(mu[i,t]) <- log(lambda[i,t]) +
          log(E[i,t]) +
          log(epsilon[i]) +
          inprod(beta[1:p], x[i,t,1:p])
        
        Y[i,t] ~ dpois(mu[i,t])
      }
    }
  })
  
  code_fixed <- nimbleCode({
    for(j in 1:p) beta[j] ~ dnorm(0, sd = 10)
    
    gamma[1] ~ dunif(0, 0.1)
    for(j in 2:K) gamma[j] ~ dunif(0, 1 - sum(gamma[1:(j-1)]))
    
    for(i in 1:n_regions) {
      alpha[i] ~ dnorm(0,1)
      
      epsilon[i] <- 1 - inprod(h[i,1:K], gamma[1:K])
      
      for(t in 1:n_times) {
        log(mu[i,t]) <- alpha[i] +
          log(E[i,t]) +
          log(epsilon[i]) +
          inprod(beta[1:p], x[i,t,1:p])
        
        Y[i,t] ~ dpois(mu[i,t])
      }
    }
  })
  
  const_ffbs  <- list(n_regions=A_val, n_times=T_val, p=p, K=K, h=C$h, w=w_fit, a0=1, b0=1)
  const_fixed <- list(n_regions=A_val, n_times=T_val, p=p, K=K, h=C$h)
  
  data_tmpl <- list(Y=dados_tmpl$Y, E=dados_tmpl$E, x=dados_tmpl$x)
  
  Cf  <- compileNimble(nimbleModel(code_ffbs, constants=const_ffbs, data=data_tmpl))
  Cfx <- compileNimble(nimbleModel(code_fixed, constants=const_fixed, data=data_tmpl))
  
  for(rep in rep_indices) {
    
    dados <- readRDS(
      file.path(pasta_dados, sprintf("replica_T%d_rep%02d.rds", T_val, rep))
    )
    
    # ================= FFBS =================
    
    Cf$setData(list(Y=dados$Y, E=dados$E, x=dados$x))
    
    mcmc_f <- nimbleMCMC(
      Cf,
      niter=niter, nburnin=nburnin, thin=thin,
      monitors=c("beta","gamma","lambda"),
      progressBar=FALSE
    )
    
    # ================= LOG-LIK FFBS =================
    
    n_samp <- nrow(mcmc_f)
    n_obs  <- A_val * T_val
    
    loglik_f <- matrix(NA, n_samp, n_obs)
    
    for(s in 1:n_samp) {
      
      beta_s  <- mcmc_f[s, grep("^beta", colnames(mcmc_f))]
      gamma_s <- mcmc_f[s, grep("^gamma", colnames(mcmc_f))]
      
      idx <- 1
      
      for(i in 1:A_val) {
        
        epsilon_i <- 1 - sum(C$h[i,] * gamma_s)
        
        for(t in 1:T_val) {
          
          lambda_st <- mcmc_f[s, paste0("lambda[",i,",",t,"]")]
          xb <- sum(dados$x[i,t,] * beta_s)
          
          mu <- lambda_st * dados$E[i,t] * epsilon_i * exp(xb)
          
          loglik_f[s, idx] <- dpois(dados$Y[i,t], mu, log=TRUE)
          
          idx <- idx + 1
        }
      }
    }
    
    waic_f <- compute_waic(loglik_f)
    
    # ================= LAMBDA FFBS =================
    
    lambda_cols <- grep("^lambda", colnames(mcmc_f))
    
    lambda_df_f <- lapply(lambda_cols, function(col) {
      
      samples <- mcmc_f[, col]
      hpd <- HPDinterval(as.mcmc(samples), prob=0.95)
      
      idx <- str_extract_all(colnames(mcmc_f)[col], "\\d+")[[1]]
      
      data.frame(
        Replica = rep,
        i = as.integer(idx[1]),
        t = as.integer(idx[2]),
        Mean  = mean(samples),
        Lower = hpd[1],
        Upper = hpd[2],
        Model = "FFBS"
      )
    }) %>% bind_rows()
    
    # ================= FIXED =================
    
    Cfx$setData(list(Y=dados$Y, E=dados$E, x=dados$x))
    
    mcmc_x <- nimbleMCMC(
      Cfx,
      niter=niter, nburnin=nburnin, thin=thin,
      monitors=c("beta","gamma","alpha"),
      progressBar=FALSE
    )
    
    # ================= LOG-LIK FIXED =================
    
    loglik_x <- matrix(NA, n_samp, n_obs)
    
    for(s in 1:n_samp) {
      
      beta_s  <- mcmc_x[s, grep("^beta", colnames(mcmc_x))]
      gamma_s <- mcmc_x[s, grep("^gamma", colnames(mcmc_x))]
      
      idx <- 1
      
      for(i in 1:A_val) {
        
        epsilon_i <- 1 - sum(C$h[i,] * gamma_s)
        alpha_i   <- mcmc_x[s, paste0("alpha[",i,"]")]
        
        for(t in 1:T_val) {
          
          xb <- sum(dados$x[i,t,] * beta_s)
          
          mu <- exp(alpha_i + xb) * dados$E[i,t] * epsilon_i
          
          loglik_x[s, idx] <- dpois(dados$Y[i,t], mu, log=TRUE)
          
          idx <- idx + 1
        }
      }
    }
    
    waic_x <- compute_waic(loglik_x)
    
    # ================= LAMBDA FIXED =================
    
    lambda_fixed <- lapply(1:A_val, function(i) {
      
      alpha_samples <- mcmc_x[, paste0("alpha[",i,"]")]
      
      lapply(1:T_val, function(t) {
        
        samples <- exp(alpha_samples)
        hpd <- HPDinterval(as.mcmc(samples), prob=0.95)
        
        data.frame(
          Replica = rep,
          i = i,
          t = t,
          Mean  = mean(samples),
          Lower = hpd[1],
          Upper = hpd[2],
          Model = "Static"
        )
      }) %>% bind_rows()
    }) %>% bind_rows()
    
    # ================= SALVAR =================
    
    lambda_all <- bind_rows(lambda_df_f, lambda_fixed)
    
    lambda_all$T_val <- T_val
    lambda_all$DGP   <- dgp
    
    saveRDS(
      list(
        lambda = lambda_all,
        loglik_ffbs = loglik_f,
        loglik_fixed = loglik_x,
        waic_ffbs = waic_f,
        waic_fixed = waic_x
      ),
      file = file.path(dir_base,
                       sprintf("lambda_rep_T%d_DGP%d_rep%02d.rds",
                               T_val, dgp, rep))
    )
  }
  
  return(NULL)
}

# ==============================================================================
# 4. EXECUÇÃO PARALELA
# ==============================================================================

scenarios <- expand.grid(T_val=c(25,100), dgp=1:3)

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)

clusterExport(cl, ls())
clusterEvalQ(cl, {
  library(nimble); library(coda); library(dplyr); library(stringr)
})

for(i in 1:nrow(scenarios)) {
  
  S <- scenarios[i,]
  
  cat("\nRodando DGP", S$dgp, "T=", S$T_val, "\n")
  
  chunks <- split(1:n_replicas,
                  cut(1:n_replicas, n_cores, labels=FALSE))
  
  parLapply(cl, chunks, worker_fun,
            T_val=S$T_val,
            dgp=S$dgp,
            w_fit=w_fit,
            niter=niter,
            nburnin=nburnin,
            thin=thin,
            dir_base=dir_projeto)
}

stopCluster(cl)

cat("\n=== FINALIZADO COM SUCESSO ===\n")