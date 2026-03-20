# ==============================================================================
# Geracao_w_unif_T_25_A_75.R  —  SEM efeito espacial, w~Unif(0.7,0.99), T=25, A=75
# ==============================================================================
library(nimble)
cat("--- Iniciando Simulação SEM Efeito Espacial ICAR ---\n")

source("_dataCaseStudy.r")   # carrega objeto 'data' e faz attach()

N_regions <- data$N          # 75
n_times   <- 25
p         <- 3
K         <- 4

beta_true <- c(-1.0, 1.0, 0.5)
w_true    <- runif(N_regions, 0.7, 0.99)   # vetor por região
a0_true   <- 1.0
b0_true   <- 1.0

set.seed(1)

# Estrutura real de clusters
h_cumulativo <- data$hAI
gamma_true   <- c(0.05, 0.10, 0.10, 0.15)

epsilon_true <- numeric(N_regions)
for(i in 1:N_regions)
  epsilon_true[i] <- 1 - sum(h_cumulativo[i, ] * gamma_true)

x_true <- array(rnorm(N_regions * n_times * p), dim = c(N_regions, n_times, p))
E_raw  <- matrix(runif(N_regions * n_times, 150, 250), nrow = N_regions)
E_true <- E_raw / mean(E_raw)

# g_it sem efeito espacial
g_it_true <- array(NA, dim = c(N_regions, n_times))
for(i in 1:N_regions)
  for(t in 1:n_times) {
    prod_val <- sum(x_true[i, t, ] * beta_true)
    g_it_true[i, t] <- E_true[i, t] * epsilon_true[i] * exp(prod_val)
  }

lambda_true <- matrix(NA, nrow = N_regions, ncol = n_times)
Y_ini       <- matrix(NA, nrow = N_regions, ncol = n_times)
at_true     <- matrix(NA, nrow = N_regions, ncol = n_times + 1)
bt_true     <- matrix(NA, nrow = N_regions, ncol = n_times + 1)

for(i in 1:N_regions) {
  at_true[i, 1] <- a0_true
  bt_true[i, 1] <- b0_true
  for(t in 1:n_times) {
    att_true_val      <- w_true[i] * at_true[i, t]   # usa w_true[i]
    btt_true_val      <- w_true[i] * bt_true[i, t]
    lambda_true[i, t] <- rgamma(1, shape = att_true_val, rate = btt_true_val)
    mu_it             <- lambda_true[i, t] * g_it_true[i, t]
    Y_ini[i, t]       <- rpois(1, mu_it)
    at_true[i, t+1]   <- att_true_val + Y_ini[i, t]
    bt_true[i, t+1]   <- btt_true_val + g_it_true[i, t]
  }
}

adj_vec     <- as.integer(data$adj)
num_vec     <- as.integer(data$num)
n_adj_val   <- as.integer(data$sumNumNeigh)
weights_vec <- rep(1.0, n_adj_val)
stopifnot("adj/num inconsistentes" = sum(num_vec) == n_adj_val)
cat("Verificação espacial OK: sum(num) = length(adj) =", n_adj_val, "\n")

constants_nimble <- list(
  n_regions = N_regions, n_times = n_times, p = p, K = K,
  h = h_cumulativo, mu_beta = rep(0, p),
  a_unif = 0.0, b_unif = 0.1,
  a0 = a0_true, b0 = b0_true,
  w = 0.9,
  adj = adj_vec, num = num_vec, weights = weights_vec, n_adj = n_adj_val
)

data_nimble <- list(Y = Y_ini, E = E_true, x = x_true)

inits_nimble_1 <- list(beta = beta_true, gamma = gamma_true, lambda = lambda_true)
inits_nimble_2 <- list(
  beta   = rnorm(p, 0, 0.5),
  gamma  = gamma_true * 0.9,
  lambda = matrix(rgamma(N_regions * n_times, 1, 1), nrow = N_regions)
)
inits_list_nimble <- list(inits_nimble_1, inits_nimble_2)

cat("--- Geração de Dados Concluída ---\n")
