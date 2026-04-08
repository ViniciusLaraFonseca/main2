# ==============================================================================
# Geracao_Dados_MonteCarlo_A145.R
# Gera 50 réplicas para T=25 e T=100, com w_i ~ Unif(0.7, 0.99) e A=145
# Usa geobr e spdep para criar a estrutura espacial real.
# ==============================================================================
library(geobr)
library(spdep)
library(sf)
library(dplyr)

set.seed(12345)
dir.create("dados_montecarlo", showWarnings = FALSE)

# 1. Obter o Mapa das Microrregiões de Saúde do Sudeste
cat("Baixando mapa das microrregiões do Sudeste...\n")
mapa_sudeste <- read_health_region(year = 2013, showProgress = FALSE) %>%
  filter(abbrev_state %in% c("MG", "ES", "RJ", "SP"))

A_val <- nrow(mapa_sudeste)
cat("Número de regiões encontradas:", A_val, "\n") # Geralmente muito próximo ou exato 145

# 2. Construir Matriz de Vizinhança
nb <- poly2nb(mapa_sudeste, queen = TRUE)

# Tratar possíveis ilhas (regiões sem vizinhos) conectando ao mais próximo
coords <- st_coordinates(st_centroid(mapa_sudeste))
if(any(card(nb) == 0)) {
  ilhas <- which(card(nb) == 0)
  for(i in ilhas) {
    distancias <- spDistsN1(coords, coords[i,])
    distancias[i] <- Inf # Ignorar a si mesmo
    vizinho_mais_proximo <- which.min(distancias)
    nb[[i]] <- as.integer(vizinho_mais_proximo)
    nb[[vizinho_mais_proximo]] <- sort(unique(c(nb[[vizinho_mais_proximo]], i)))
  }
}
wb <- nb2WB(nb)

# 3. Criar Matriz Cumulativa de Clusters (h) - 4 clusters aleatórios
K <- 4
cluster_ids <- sample(1:K, A_val, replace = TRUE)
h_cumulativo <- matrix(0, nrow = A_val, ncol = K)
for(i in 1:A_val) {
  h_cumulativo[i, 1:cluster_ids[i]] <- 1
}

# 4. Parâmetros Verdadeiros Fixos
p <- 3
beta_true <- c(-1.0, 1.0, 0.5)
gamma_true <- c(0.05, 0.10, 0.10, 0.15)
a0_true <- 1.0
b0_true <- 1.0

epsilon_true <- numeric(A_val)
for(i in 1:A_val) epsilon_true[i] <- 1 - sum(h_cumulativo[i, ] * gamma_true)

# 5. Função para gerar réplicas para um dado T
gerar_replicas <- function(T_val, num_replicas = 50) {
  cat(sprintf("\nGerando %d réplicas para T=%d...\n", num_replicas, T_val))
  
  for(rep in 1:num_replicas) {
    # w_i diferente para cada região, sorteado na réplica
    w_true <- runif(A_val, 0.7, 0.99) 
    
    x_true <- array(rnorm(A_val * T_val * p), dim = c(A_val, T_val, p))
    E_raw  <- matrix(runif(A_val * T_val, 150, 250), nrow = A_val)
    E_true <- E_raw / mean(E_raw)
    
    g_it_true <- array(NA, dim = c(A_val, T_val))
    for(i in 1:A_val) {
      for(t in 1:T_val) {
        prod_val <- sum(x_true[i, t, ] * beta_true)
        g_it_true[i, t] <- E_true[i, t] * epsilon_true[i] * exp(prod_val) # Gerado SEM efeito espacial (s=0)
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
    
    # Salvar tudo necessário num arquivo .RData para a réplica
    replica_data <- list(
      Y = Y_ini, E = E_true, x = x_true, 
      lambda_true = lambda_true, w_true = w_true,
      constants = list(
        n_regions = A_val, n_times = T_val, p = p, K = K, h = h_cumulativo,
        adj = wb$adj, num = wb$num, weights = wb$weights, n_adj = length(wb$adj)
      )
    )
    
    nome_arquivo <- sprintf("dados_montecarlo/replica_T%d_rep%02d.rds", T_val, rep)
    saveRDS(replica_data, nome_arquivo)
  }
}

gerar_replicas(T_val = 25, num_replicas = 50)
gerar_replicas(T_val = 100, num_replicas = 50)
cat("\nGeração concluída com sucesso. Dados em 'dados_montecarlo/'.\n")