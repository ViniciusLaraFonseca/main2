# ==============================================================================
# Geração_w_0_7_T_100_A_75.R
# ==============================================================================
library(nimble)
cat("--- Iniciando Simulação SEM Efeito Espacial ICAR ---\n")

# --- PASSO 1: CARREGAR DADOS REAIS (vizinhança + clusters + covariáveis) ---
source("_dataCaseStudy.r")   # carrega objeto 'data' e faz attach()

N_regions <- data$N          # 75
n_times   <- 100
p         <- 3
K         <- 4

# Parâmetros de Regressão e Dinâmicos
beta_true <- c(-1.0, 1.0, 0.5)
w_true    <- 0.7
a0_true   <- 1.0
b0_true   <- 1.0

set.seed(1)

# --- PASSO 2: CLUSTERS E GAMMAS (estrutura real do estudo de caso) ---
# Usa hAI e clAI do arquivo em vez de amostragem aleatória
h_cumulativo <- data$hAI    # matriz 75x4 já pronta

gamma_true <- c(0.05, 0.10, 0.10, 0.15)

# --- PASSO 3: EPSILON E PREDITOR LINEAR ---
epsilon_true <- numeric(N_regions)
for(i in 1:N_regions) {
  epsilon_true[i] <- 1 - sum(h_cumulativo[i, ] * gamma_true)
}

# Simulando Covariáveis (x) e Offset (E)
x_true <- array(rnorm(N_regions * n_times * p), dim = c(N_regions, n_times, p))
E_raw  <- matrix(runif(N_regions * n_times, 150, 250), nrow = N_regions)
E_true <- E_raw / mean(E_raw)

# --- PASSO 4: CÁLCULO DE g_it SEM EFEITO ESPACIAL ---
g_it_true <- array(NA, dim = c(N_regions, n_times))
for(i in 1:N_regions) {
  for(t in 1:n_times) {
    prod_val <- sum(x_true[i, t, ] * beta_true)
    g_it_true[i, t] <- E_true[i, t] * epsilon_true[i] * exp(prod_val)
  }
}

# --- PASSO 5: SIMULAR LAMBDA E Y ---
lambda_true <- matrix(NA, nrow = N_regions, ncol = n_times)
Y_ini       <- matrix(NA, nrow = N_regions, ncol = n_times)
at_true     <- matrix(NA, nrow = N_regions, ncol = n_times + 1)
bt_true     <- matrix(NA, nrow = N_regions, ncol = n_times + 1)

for(i in 1:N_regions) {
  at_true[i, 1] <- a0_true
  bt_true[i, 1] <- b0_true
  
  for(t in 1:n_times) {
    att_true_val <- w_true * at_true[i, t]
    btt_true_val <- w_true * bt_true[i, t]
    
    lambda_true[i, t] <- rgamma(1, shape = att_true_val, rate = btt_true_val)
    
    mu_it <- lambda_true[i, t] * g_it_true[i, t]
    Y_ini[i, t] <- rpois(1, mu_it)
    
    at_true[i, t+1] <- att_true_val + Y_ini[i, t]
    bt_true[i, t+1] <- btt_true_val + g_it_true[i, t]
  }
}

# --- PASSO 6: PREPARAR OBJETOS PARA O NIMBLE ---

# Estrutura de vizinhança vinda do arquivo — garantidamente consistente
adj_vec     <- as.integer(data$adj)
num_vec     <- as.integer(data$num)
n_adj_val   <- as.integer(data$sumNumNeigh)   # 386
weights_vec <- rep(1.0, n_adj_val)

# Verificação de segurança
stopifnot("adj/num inconsistentes" = sum(num_vec) == n_adj_val)
cat("Verificação espacial OK: sum(num) = length(adj) =", n_adj_val, "\n")

constants_nimble <- list(
  n_regions = N_regions,
  n_times   = n_times,
  p         = p,
  K         = K,
  h         = h_cumulativo,
  mu_beta   = rep(0, p),
  a_unif    = 0.0, b_unif = 0.1,
  a0        = a0_true, b0 = b0_true,
  w         = 0.9,
  # Estrutura espacial — disponível para o modelo espacial no script principal
  adj     = adj_vec,
  num     = num_vec,
  weights = weights_vec,
  n_adj   = n_adj_val
)

data_nimble <- list(
  Y = Y_ini,
  E = E_true,
  x = x_true
)

inits_nimble_1 <- list(
  beta   = beta_true,
  gamma  = gamma_true,
  lambda = lambda_true
)
inits_nimble_2 <- list(
  beta   = rnorm(p, 0, 0.5),
  gamma  = gamma_true * 0.9,
  lambda = matrix(rgamma(N_regions * n_times, 1, 1), nrow = N_regions)
)
inits_list_nimble <- list(inits_nimble_1, inits_nimble_2)

cat("--- Geração de Dados Concluída ---\n")