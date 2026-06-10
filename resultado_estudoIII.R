library(dplyr)
library(purrr)
library(readr)

# ==============================================================================
# 1. Função geral
# ==============================================================================

processa_cenario <- function(T_val, A_val, w_vals, base_path) {
  
  arquivos <- paste0(
    base_path,
    "/Estaticos_T", T_val, "_w", w_vals, ".csv"
  )
  
  dados <- map2_dfr(arquivos, w_vals, ~ {
    read_csv(.x) %>%
      mutate(w_fix = .y)
  })
  
  dados %>%
    mutate(
      Parametro = as.factor(Parametro),
      T_val = T_val,
      A_val = A_val
    ) %>%
    group_by(Parametro, w_fix, T_val, A_val) %>%
    summarise(
      Mean_bias = mean(Bias, na.rm = TRUE),
      Mean_MSE  = mean(MSE, na.rm = TRUE),
      Mean_cov  = mean(Cov, na.rm = TRUE),
      Mean_ESS  = mean(ESS, na.rm = TRUE),
      .groups = "drop"
    )
}

# ==============================================================================
# 2. Grid com caminhos
# ==============================================================================

cenarios <- tibble(
  T_val = c(25, 100, 25, 100),
  A_val = c(75, 75, 145, 145),
  base_path = c(
    "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Simulacao/main/main2/resultados_MC_NonSpatial_A75",
    "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Simulacao/main/main2/resultados_MC_NonSpatial_A75",
    "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Simulacao/main/main2/resultados_MC_NonSpatial",
    "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Simulacao/main/main2/resultados_MC_NonSpatial"
  )
)

w_vals <- c(0.7, 0.9)

# ==============================================================================
# 3. Executa tudo
# ==============================================================================

resultado_final <- pmap_dfr(
  cenarios,
  ~ processa_cenario(..1, ..2, w_vals, ..3)
)

print(resultado_final, n = Inf)