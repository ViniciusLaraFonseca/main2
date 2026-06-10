# ==============================================================================
# Geracao_Dados_MC_CustoDinamica.R
#
# Gera 50 réplicas para 3 DGPs e T ∈ {25, 100} com A = 75
#
# DGP 1 — Intercepto regional estático
#   λ_i ~ Normal(0, 1), sorteado uma vez por réplica/região, constante no tempo.
#   Representa o cenário mais simples: heterogeneidade espacial pura, sem dinâmica.
#
# DGP 2 — Intercepto dinâmico por cluster (FFBS cluster-level)
#   λ_{k,t} evolui via FFBS agregado por cluster k.
#   Regiões no mesmo cluster compartilham o mesmo λ_{k,t}.
#   Dinâmica temporal existe, mas é compartilhada dentro de cada grupo.
#
# DGP 3 — Intercepto dinâmico espaço-temporal individual (FFBS completo)
#   λ_{i,t} evolui de forma independente por região e tempo.
#   Este é o processo que o modelo FFBS proposto foi desenhado para capturar.
#

# ==============================================================================
rm(list = ls())
library(spdep)

set.seed(2024)

dir_projeto <- "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Simulacao/main/main2"
setwd(dir_projeto)

# ==============================================================================
# 1. ESTRUTURA ESPACIAL — Mapa Minas Gerais (A = 75 regiões)
# ==============================================================================



source("_dataCaseStudy.r")
A_val       <- data$N
adj_vec     <- as.integer(data$adj)
num_vec     <- as.integer(data$num)
n_adj_val   <- as.integer(data$sumNumNeigh)
weights_vec <- rep(1.0, n_adj_val)

# ==============================================================================
# 2. PARÂMETROS VERDADEIROS FIXOS
# ==============================================================================
K          <- 4
p          <- 3
beta_true  <- c(-1.0,  1.0,  0.5)
gamma_true <- c(0.05, 0.10, 0.10, 0.15)
a0_true    <- 1.0
b0_true    <- 1.0
w_gen      <- 0.8   # Persistência usada nos DGPs 2 e 3

# Atribuição de clusters: K blocos consecutivos de regiões
cluster_ids <- as.integer(cut(seq_len(A_val), breaks = K, labels = FALSE))

# Matriz cumulativa h (A × K): h[i, k] = 1 se região i pertence a cluster <= k
h_cumulativo <- matrix(0, nrow = A_val, ncol = K)
for(i in seq_len(A_val)) h_cumulativo[i, seq_len(cluster_ids[i])] <- 1

# ε_i = 1 - h_i' γ  (frações de exposição por cluster)
epsilon_true <- as.numeric(1 - h_cumulativo %*% gamma_true)

# Lista de constantes salva em cada réplica
make_constants <- function(T_val) {
  list(
    n_regions = A_val, n_times = T_val, p = p, K = K,
    h         = h_cumulativo,
    adj       = adj_vec, num = num_vec,
    weights   = weights_vec, n_adj = n_adj_val,
    cluster_ids = cluster_ids
  )
}

# Criação das pastas de saída
for(dgp in 1:3) dir.create(sprintf("dados_MC_Custo_DGP%d", dgp), showWarnings = FALSE)

# ==============================================================================
# 3. FUNÇÕES DE GERAÇÃO POR DGP
# ==============================================================================

# ------------------------------------------------------------------------------
# DGP 1: λ_i estático (um único valor Gamma por região, constante no tempo)
# ------------------------------------------------------------------------------
gerar_DGP1 <- function(T_val) {
  x <- array(rnorm(A_val * T_val * p), dim = c(A_val, T_val, p))
  
  E <- {
    E_raw <- matrix(runif(A_val * T_val, 150, 250), A_val)
    E_raw / mean(E_raw)
  }
  
  lambda_i <- rnorm(A_val, 0, 1)  # intercepto no log
  
  Y <- matrix(NA_integer_, A_val, T_val)
  mu <- matrix(NA_real_, A_val, T_val)
  
  for(i in seq_len(A_val))
    for(t in seq_len(T_val)) {
      
      eta_it <- lambda_i[i] + sum(x[i, t, ] * beta_true)
      
      mu[i, t] <- exp(eta_it) * epsilon_true[i] * E[i, t]
      
      Y[i, t] <- rpois(1, mu[i, t])
    }
  
  list(Y = Y, E = E, x = x, mu = mu,
       lambda_true = lambda_i,
       constants = make_constants(T_val))
}


gerar_DGP2 <- function(T_val) {
  x <- array(rnorm(A_val * T_val * p), dim = c(A_val, T_val, p))
  
  E <- {
    E_raw <- matrix(runif(A_val * T_val, 150, 250), A_val)
    E_raw / mean(E_raw)
  }
  
  # -----------------------------
  # Intercepto dinâmico por cluster (i.i.d. no tempo)
  # λ_{k,t} ~ N(0,1)
  # -----------------------------
  lambda_kt <- matrix(rnorm(K * T_val, 0, 1), nrow = K, ncol = T_val)
  
  # Expandir para regiões
  lambda_true <- matrix(NA_real_, A_val, T_val)
  for(i in seq_len(A_val)) {
    k <- cluster_ids[i]
    lambda_true[i, ] <- lambda_kt[k, ]
  }
  
  # -----------------------------
  # Geração dos dados
  # -----------------------------
  Y  <- matrix(NA_integer_, A_val, T_val)
  mu <- matrix(NA_real_, A_val, T_val)
  
  for(i in seq_len(A_val))
    for(t in seq_len(T_val)) {
      
      eta_it <- lambda_true[i, t] + sum(x[i, t, ] * beta_true)
      
      mu[i, t] <- exp(eta_it) * epsilon_true[i] * E[i, t]
      
      Y[i, t] <- rpois(1, mu[i, t])
    }
  
  list(
    Y = Y,
    E = E,
    x = x,
    mu = mu,
    lambda_true = lambda_true,
    constants = make_constants(T_val)
  )
}
# ------------------------------------------------------------------------------
# DGP 3: λ_{i,t} dinâmico individual via FFBS completo
#
# Forward filter por região i:  a_{i,0} = a0, b_{i,0} = b0
# Preditiva em t:  λ_{i,t} ~ Gamma(w·a_{i,t-1}, w·b_{i,t-1})
# Y_{i,t} ~ Poisson(λ_{i,t} · g_{i,t})
# Atualização:     a_{i,t} = w·a_{i,t-1} + Y_{i,t}
#                  b_{i,t} = w·b_{i,t-1} + g_{i,t}
# ------------------------------------------------------------------------------
gerar_DGP3 <- function(T_val) {
  x   <- array(rnorm(A_val * T_val * p), dim = c(A_val, T_val, p))
  E   <- { E_raw <- matrix(runif(A_val * T_val, 150, 250), A_val)
            E_raw / mean(E_raw) }

  lambda_true <- matrix(NA_real_, A_val, T_val)
  Y           <- matrix(NA_integer_, A_val, T_val)
  at          <- matrix(a0_true, A_val, T_val + 1)
  bt          <- matrix(b0_true, A_val, T_val + 1)

  for(i in seq_len(A_val))
    for(t in seq_len(T_val)) {
      att <- w_gen * at[i, t]
      btt <- w_gen * bt[i, t]
      lambda_true[i, t] <- rgamma(1, shape = att, rate = btt)
      g_it   <- E[i, t] * epsilon_true[i] * exp(sum(x[i, t, ] * beta_true))
      Y[i, t] <- rpois(1, lambda_true[i, t] * g_it)
      at[i, t + 1] <- att + Y[i, t]
      bt[i, t + 1] <- btt + g_it
    }

  list(Y = Y, E = E, x = x, lambda_true = lambda_true,
       w_gen = w_gen, constants = make_constants(T_val))
}

# ==============================================================================
# 4. LOOP DE GERAÇÃO
# ==============================================================================
n_replicas <- 50

for(T_val in c(25, 100)) {
  cat(sprintf("\n--- T = %d ---\n", T_val))
  for(rep in seq_len(n_replicas)) {

    saveRDS(gerar_DGP1(T_val),
            sprintf("dados_MC_Custo_DGP1/replica_T%d_rep%02d.rds", T_val, rep))

    saveRDS(gerar_DGP2(T_val),
            sprintf("dados_MC_Custo_DGP2/replica_T%d_rep%02d.rds", T_val, rep))

    saveRDS(gerar_DGP3(T_val),
            sprintf("dados_MC_Custo_DGP3/replica_T%d_rep%02d.rds", T_val, rep))

    if(rep %% 10 == 0) cat(sprintf("  Réplica %d/%d\n", rep, n_replicas))
  }
}

cat("\n=== Geração de dados concluída! ===\n")
cat(sprintf("Estrutura: A=%d regiões | K=%d clusters | p=%d covariáveis\n", A_val, K, p))
cat(sprintf("DGPs 2 e 3 gerados com w_gen = %.2f\n", w_gen))
cat(sprintf("Parâmetros verdadeiros:\n"))
cat(sprintf("  beta  = (%s)\n", paste(beta_true,  collapse = ", ")))
cat(sprintf("  gamma = (%s)\n", paste(gamma_true, collapse = ", ")))
