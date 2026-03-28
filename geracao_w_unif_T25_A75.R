# ==============================================================================
# Geracao_w_unif_T_25_A_75.R  — ROBUSTO (sem efeito espacial, w~Unif)
# ==============================================================================

library(nimble)
cat("--- Iniciando Simulação SEM Efeito Espacial ICAR ---\n")

# ==============================================================================
# 1. CARREGAR DADOS DE FORMA SEGURA (SEM attach / SEM data GLOBAL)
# ==============================================================================
env_data <- new.env(parent = baseenv())
source("data.r", local = env_data)

if(!exists("data_list", envir = env_data))
  stop("ERRO: data.r deve criar objeto 'data_list'")

data_list <- env_data$data_list
rm(env_data)

# ==============================================================================
# 2. DIMENSÕES E PARÂMETROS
# ==============================================================================
N_regions <- data_list$N
n_times   <- 25
p         <- 3
K         <- 4

beta_true <- c(-1.0, 1.0, 0.5)
w_true    <- runif(N_regions, 0.7, 0.99)
a0_true   <- 1.0
b0_true   <- 1.0

set.seed(1)

# ==============================================================================
# 3. ESTRUTURA DE CLUSTERS
# ==============================================================================
h_cumulativo <- data_list$hAI
gamma_true   <- c(0.05, 0.10, 0.10, 0.15)

epsilon_true <- 1 - as.vector(h_cumulativo %*% gamma_true)

# ==============================================================================
# 4. COVARIÁVEIS E EXPOSIÇÃO
# ==============================================================================
x_true <- array(rnorm(N_regions * n_times * p),
                dim = c(N_regions, n_times, p))

E_raw  <- matrix(runif(N_regions * n_times, 150, 250),
                 nrow = N_regions)

E_true <- E_raw / mean(E_raw)

# ==============================================================================
# 5. COMPONENTE g_it (SEM efeito espacial)
# ==============================================================================
g_it_true <- array(0, dim = c(N_regions, n_times))

for(i in 1:N_regions) {
  for(t in 1:n_times) {
    prod_val <- sum(x_true[i, t, ] * beta_true)
    g_it_true[i, t] <- E_true[i, t] * epsilon_true[i] * exp(prod_val)
  }
}

# ==============================================================================
# 6. SIMULAÇÃO FFBS-LIKE
# ==============================================================================
lambda_true <- matrix(NA, N_regions, n_times)
Y_ini       <- matrix(NA, N_regions, n_times)

at_true <- matrix(NA, N_regions, n_times + 1)
bt_true <- matrix(NA, N_regions, n_times + 1)

for(i in 1:N_regions) {
  
  at_true[i, 1] <- a0_true
  bt_true[i, 1] <- b0_true
  
  for(t in 1:n_times) {
    
    att <- w_true[i] * at_true[i, t]
    btt <- w_true[i] * bt_true[i, t]
    
    lambda_true[i, t] <- rgamma(1, shape = att, rate = btt)
    
    mu_it <- lambda_true[i, t] * g_it_true[i, t]
    Y_ini[i, t] <- rpois(1, mu_it)
    
    at_true[i, t + 1] <- att + Y_ini[i, t]
    bt_true[i, t + 1] <- btt + g_it_true[i, t]
  }
}

# ==============================================================================
# 7. ESTRUTURA ESPACIAL (usada no modelo spatial)
# ==============================================================================
adj_vec     <- as.integer(data_list$adj)
num_vec     <- as.integer(data_list$num)
n_adj_val   <- as.integer(data_list$sumNumNeigh)
weights_vec <- rep(1.0, n_adj_val)

stopifnot(sum(num_vec) == n_adj_val)

cat("Verificação espacial OK: sum(num) = length(adj) =", n_adj_val, "\n")

# ==============================================================================
# 8. OBJETOS PARA NIMBLE
# ==============================================================================
constants_nimble <- list(
  n_regions = N_regions,
  n_times   = n_times,
  p         = p,
  K         = K,
  h         = h_cumulativo,
  mu_beta   = rep(0, p),
  a_unif    = 0.0,
  b_unif    = 0.1,
  a0        = a0_true,
  b0        = b0_true,
  w         = 0.9,
  adj       = adj_vec,
  num       = num_vec,
  weights   = weights_vec,
  n_adj     = n_adj_val
)

data_nimble <- list(
  Y = Y_ini,
  E = E_true,
  x = x_true
)

# ==============================================================================
# 9. VALIDAÇÃO (ESSENCIAL)
# ==============================================================================
stopifnot(
  all(dim(Y_ini) == c(N_regions, n_times)),
  all(dim(E_true) == c(N_regions, n_times)),
  all(dim(x_true) == c(N_regions, n_times, p))
)

# ==============================================================================
# 10. OUTPUT FINAL
# ==============================================================================
cat("--- Geração de Dados Concluída (ROBUSTA) ---\n")