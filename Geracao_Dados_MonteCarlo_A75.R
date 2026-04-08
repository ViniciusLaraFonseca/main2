# ==============================================================================
# Geracao_Dados_MonteCarlo_A75.R
# Gera 50 réplicas para T=25 e T=100, com w_i ~ Unif(0.7, 0.99) e A=75
# Usa _dataCaseStudy.R para a estrutura de matriz e vizinhança.
# ==============================================================================
rm(list=ls())
library(nimble)

# --- 1. Carregar Dados do Estudo de Caso (A=75) ---
# Certifique-se de que _dataCaseStudy.R (ou _dataCaseStudy.r) está no seu diretório
source("_dataCaseStudy.r") 

A_val <- data$N # 75
K <- 4
p <- 3

# Extrair matrizes espaciais estáticas
h_cumulativo <- data$hAI
adj_vec      <- as.integer(data$adj)
num_vec      <- as.integer(data$num)
n_adj_val    <- as.integer(data$sumNumNeigh)
weights_vec  <- rep(1.0, n_adj_val)

# Parâmetros verdadeiros
beta_true  <- c(-1.0, 1.0, 0.5)
gamma_true <- c(0.05, 0.10, 0.10, 0.15)
a0_true <- 1.0
b0_true <- 1.0

# Epsilon fixo por região
epsilon_true <- numeric(A_val)
for(i in 1:A_val) epsilon_true[i] <- 1 - sum(h_cumulativo[i, ] * gamma_true)

# Criar pasta dedicada
dir.create("dados_montecarlo_A75", showWarnings = FALSE)

# --- 2. Função de Geração de Réplicas ---
gerar_replicas <- function(T_val, num_replicas = 50) {
  cat(sprintf("\nGerando %d réplicas para A=%d e T=%d...\n", num_replicas, A_val, T_val))
  
  for(rep in 1:num_replicas) {
    # w sorteado novo para cada réplica
    w_true <- runif(A_val, 0.7, 0.99) 
    
    x_true <- array(rnorm(A_val * T_val * p), dim = c(A_val, T_val, p))
    E_raw  <- matrix(runif(A_val * T_val, 150, 250), nrow = A_val)
    E_true <- E_raw / mean(E_raw)
    
    g_it_true <- array(NA, dim = c(A_val, T_val))
    for(i in 1:A_val) {
      for(t in 1:T_val) {
        prod_val <- sum(x_true[i, t, ] * beta_true)
        # Gerado sem efeito espacial (s=0)
        g_it_true[i, t] <- E_true[i, t] * epsilon_true[i] * exp(prod_val) 
      }
    }
    
    lambda_true <- matrix(NA, nrow = A_val, ncol = T_val)
    Y_ini       <- matrix(NA, nrow = A_val, ncol = T_val)
    at_true     <- matrix(NA, nrow = A_val, ncol = T_val + 1)
    bt_true     <- matrix(NA, nrow = A_val, ncol = T_val + 1)
    
    for(i in 1:A_val) {
      at_true[i, 1] <- a0_true
      bt_true[i, 1] <- b0_true
      for(t in 1:T_val) {
        att_true_val      <- w_true[i] * at_true[i, t]
        btt_true_val      <- w_true[i] * bt_true[i, t]
        lambda_true[i, t] <- rgamma(1, shape = att_true_val, rate = btt_true_val)
        mu_it             <- lambda_true[i, t] * g_it_true[i, t]
        Y_ini[i, t]       <- rpois(1, mu_it)
        at_true[i, t+1]   <- att_true_val + Y_ini[i, t]
        bt_true[i, t+1]   <- btt_true_val + g_it_true[i, t]
      }
    }
    
    # Salvar
    replica_data <- list(
      Y = Y_ini, E = E_true, x = x_true, 
      lambda_true = lambda_true, w_true = w_true,
      constants = list(
        n_regions = A_val, n_times = T_val, p = p, K = K, h = h_cumulativo,
        adj = adj_vec, num = num_vec, weights = weights_vec, n_adj = n_adj_val
      )
    )
    
    nome_arquivo <- sprintf("dados_montecarlo_A75/replica_T%d_rep%02d.rds", T_val, rep)
    saveRDS(replica_data, nome_arquivo)
  }
}

set.seed(999) # Seed fixa para reprodutibilidade
gerar_replicas(T_val = 25, num_replicas = 50)
gerar_replicas(T_val = 100, num_replicas = 50)
cat("\nGeração concluída! Dados salvos na pasta 'dados_montecarlo_A75/'.\n")