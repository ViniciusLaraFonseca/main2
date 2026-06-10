# ==============================================================================
# pos_processamento_MC_CustoDinamica.R
# Pós-processamento dos resultados do Estudo_MC_CustoDinamica.R
# ==============================================================================
# Objetivo: Consolidar resultados das 50 réplicas, calcular métricas de performance
#           (Viés, MSE, Coverage, ESS) e gerar visualizações comparativas 
#           entre modelos FFBS e Static para os 3 DGPs
# ==============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(stringr)
library(coda)

# ==============================================================================
# CONFIGURAÇÃO
# ==============================================================================

rm(list = ls())

dir_projeto <- "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Simulacao/main/main2"
setwd(dir_projeto)

# Cenários rodados: T ∈ {25, 100}, DGP ∈ {1, 2, 3}
scenarios <- expand.grid(T_val = c(25, 100), dgp = 1:3)
n_replicas <- 50

# Cores para os modelos
cores_modelos <- c("FFBS" = "#2166AC", "Static" = "#D6604D")
cores_dgp     <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")

# Criar pasta de saída para pós-processamento
dir.create("pos_processamento", showWarnings = FALSE)

# ==============================================================================
# 1. FUNÇÃO DE CONSOLIDAÇÃO
# ==============================================================================

consolidate_results <- function(T_val, dgp, n_replicas) {
  
  cat(sprintf("\n--- Consolidando T=%d, DGP=%d ---\n", T_val, dgp))
  
  lambda_list <- list()
  loglik_list <- list()
  waic_list   <- list()
  
  for(rep in 1:n_replicas) {
    
    fname <- sprintf("lambda_rep_T%d_DGP%d_rep%02d.rds", T_val, dgp, rep)
    
    if(!file.exists(fname)) {
      cat(sprintf("  AVISO: Arquivo não encontrado: %s\n", fname))
      next
    }
    
    res <- readRDS(fname)
    
    # ─ Lambda resumos (já inclui região i, tempo t, True, Mean, Lower, Upper, Model)
    lambda_list[[rep]] <- res$lambda %>%
      mutate(Replica_ID = rep)
    
    # ─ Log-verossimilhança
    loglik_list[[rep]] <- list(
      ffbs = res$loglik_ffbs,
      static = res$loglik_fixed
    )
    
    # ─ WAIC
    waic_list[[rep]] <- data.frame(
      Replica = rep,
      Model = c("FFBS", "Static"),
      WAIC = c(res$waic_ffbs$waic, res$waic_fixed$waic),
      lppd = c(res$waic_ffbs$lppd, res$waic_fixed$lppd),
      p_waic = c(res$waic_ffbs$p_waic, res$waic_fixed$p_waic),
      stringsAsFactors = FALSE
    )
  }
  
  lambda_all <- bind_rows(lambda_list)
  waic_all   <- bind_rows(waic_list)
  
  return(list(
    lambda = lambda_all,
    waic = waic_all,
    loglik = loglik_list
  ))
}

# Consolidar todos os cenários
cat("\n════════════════════════════════════════════════════════════════\n")
cat("  CONSOLIDANDO RESULTADOS...\n")
cat("════════════════════════════════════════════════════════════════\n")

all_results <- list()

for(i in 1:nrow(scenarios)) {
  S <- scenarios[i, ]
  key <- sprintf("T%d_DGP%d", S$T_val, S$dgp)
  all_results[[key]] <- consolidate_results(S$T_val, S$dgp, n_replicas)
}

# ==============================================================================
# 2. ANÁLISE DE LAMBDA: VIÉS, MSE E COVERAGE
# ==============================================================================

cat("\n--- Calculando métricas de Lambda ---\n")

lambda_metrics_all <- bind_rows(lapply(
  names(all_results),
  function(key) {
    results_list <- all_results[[key]]
    lambda_df <- results_list$lambda
    
    # Resumir por região, tempo e modelo
    metrics <- lambda_df %>%
      group_by(i, t, Model, DGP, T_val) %>%
      summarise(
        n_rep = n_distinct(Replica_ID),
        Mean_est = mean(Mean, na.rm = TRUE),
        SD_est = sd(Mean, na.rm = TRUE),
        MSE = mean((Mean - True)^2, na.rm = TRUE),
        Bias = mean(Mean - True, na.rm = TRUE),
        Coverage_95 = mean((Lower <= True & True <= Upper), na.rm = TRUE),
        .groups = "drop"
      )
    
    return(metrics)
  }
))

# Resumo por cenário
lambda_scenario_summary <- lambda_metrics_all %>%
  group_by(T_val, DGP, Model) %>%
  summarise(
    N_lambda_params = n(),
    MSE_medio = mean(MSE, na.rm = TRUE),
    Bias_medio = mean(Bias, na.rm = TRUE),
    Coverage_95_pct = 100 * mean(Coverage_95, na.rm = TRUE),
    SD_est_medio = mean(SD_est, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(T_val, DGP, Model)

write_csv(lambda_metrics_all, "pos_processamento/lambda_metrics_detailed.csv")
write_csv(lambda_scenario_summary, "pos_processamento/lambda_scenario_summary.csv")

cat("\n=== RESUMO: MSE por Cenário ===\n")
print(lambda_scenario_summary)

# ==============================================================================
# 3. ANÁLISE WAIC: COMPARAÇÃO ENTRE MODELOS
# ==============================================================================

cat("\n--- Analisando WAIC ---\n")

waic_summary <- bind_rows(lapply(
  names(all_results),
  function(key) {
    parts <- str_match(key, "T(\\d+)_DGP(\\d+)")
    T_val <- as.integer(parts[2])
    dgp <- as.integer(parts[3])
    
    all_results[[key]]$waic %>%
      mutate(T_val = T_val, DGP = dgp, .before = 1)
  }
))

waic_by_scenario <- waic_summary %>%
  group_by(T_val, DGP, Model) %>%
  summarise(
    WAIC_mean = mean(WAIC, na.rm = TRUE),
    WAIC_sd = sd(WAIC, na.rm = TRUE),
    WAIC_min = min(WAIC, na.rm = TRUE),
    WAIC_max = max(WAIC, na.rm = TRUE),
    lppd_mean = mean(lppd, na.rm = TRUE),
    p_waic_mean = mean(p_waic, na.rm = TRUE),
    .groups = "drop"
  )

# Calcular WAIC relativo (diferença do melhor modelo)
waic_by_scenario <- waic_by_scenario %>%
  group_by(T_val, DGP) %>%
  mutate(
    WAIC_best = min(WAIC_mean),
    Delta_WAIC = WAIC_mean - WAIC_best,
    WAIC_rel_ranking = rank(WAIC_mean)
  ) %>%
  ungroup() %>%
  select(-WAIC_best)

write_csv(waic_summary, "pos_processamento/waic_all_replicas.csv")
write_csv(waic_by_scenario, "pos_processamento/waic_scenario_summary.csv")

cat("\n=== RESUMO WAIC ===\n")
print(waic_by_scenario %>% select(T_val, DGP, Model, WAIC_mean, Delta_WAIC, WAIC_rel_ranking))

# ==============================================================================
# 4. GRÁFICOS: VIÉS E MSE
# ==============================================================================

cat("\n--- Gerando gráficos de desempenho ---\n")

# Viés e MSE por cenário
bias_mse_plot <- lambda_scenario_summary %>%
  pivot_longer(cols = c(Bias_medio, MSE_medio),
               names_to = "Metrica", values_to = "Valor") %>%
  mutate(Metrica = factor(Metrica, levels = c("Bias_medio", "MSE_medio"),
                          labels = c("Viés Médio", "MSE Médio")))

ggsave("pos_processamento/fig_01_bias_mse_comparison.png",
  ggplot(bias_mse_plot, aes(x = factor(DGP), y = Valor, fill = Model)) +
    geom_col(position = "dodge", width = 0.7) +
    facet_grid(Metrica ~ T_val, 
               labeller = labeller(T_val = c("25" = "T=25", "100" = "T=100"),
                                   Metrica = c("Bias_medio" = "Viés Médio",
                                              "MSE_medio" = "MSE Médio")),
               scales = "free_y") +
    scale_fill_manual(values = cores_modelos) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom") +
    labs(title = "Viés e MSE de λ̂ por Modelo e Cenário",
         x = "DGP", y = "Valor",
         fill = "Modelo"),
  width = 11, height = 7)

# ==============================================================================
# 5. COVERAGE: HPD 95%
# ==============================================================================

coverage_summary <- lambda_scenario_summary %>%
  select(T_val, DGP, Model, Coverage_95_pct)

ggsave("pos_processamento/fig_02_coverage_hpd95.png",
  ggplot(coverage_summary, aes(x = factor(DGP), y = Coverage_95_pct, fill = Model)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red", linewidth = 1) +
    facet_wrap(~T_val, labeller = labeller(T_val = c("25" = "T=25", "100" = "T=100"))) +
    scale_y_continuous(limits = c(80, 100), breaks = seq(80, 100, 5)) +
    scale_fill_manual(values = cores_modelos) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom") +
    labs(title = "Cobertura do HPD 95% para λ̂",
         x = "DGP", y = "Cobertura (%)",
         fill = "Modelo"),
  width = 10, height = 5)

# ==============================================================================
# 6. WAIC: COMPARAÇÃO VISUAL
# ==============================================================================

ggsave("pos_processamento/fig_03_waic_by_scenario.png",
  ggplot(waic_by_scenario, aes(x = factor(DGP), y = WAIC_mean, fill = Model)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = WAIC_mean - WAIC_sd, ymax = WAIC_mean + WAIC_sd),
                  position = position_dodge(0.7), width = 0.3) +
    facet_wrap(~T_val, labeller = labeller(T_val = c("25" = "T=25", "100" = "T=100"))) +
    scale_fill_manual(values = cores_modelos) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom") +
    labs(title = "WAIC por Modelo e Cenário",
         x = "DGP", y = "WAIC",
         fill = "Modelo"),
  width = 10, height = 5)

# Diferença relativa de WAIC (ΔW)
ggsave("pos_processamento/fig_04_waic_delta.png",
  ggplot(waic_by_scenario, aes(x = factor(DGP), y = Delta_WAIC, fill = Model)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    facet_wrap(~T_val, labeller = labeller(T_val = c("25" = "T=25", "100" = "T=100"))) +
    scale_fill_manual(values = cores_modelos) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom") +
    labs(title = "ΔW = WAIC(modelo) - WAIC(melhor) [menor é melhor]",
         x = "DGP", y = "Diferença de WAIC",
         fill = "Modelo"),
  width = 10, height = 5)

# ==============================================================================
# 7. HEATMAP: MSE por DGP e Modelo
# ==============================================================================

mse_matrix <- lambda_scenario_summary %>%
  select(T_val, DGP, Model, MSE_medio) %>%
  mutate(Cenario = paste0("T=", T_val, " DGP=", DGP)) %>%
  select(Cenario, Model, MSE_medio) %>%
  pivot_wider(names_from = Model, values_from = MSE_medio)

ggsave("pos_processamento/fig_05_mse_heatmap.png",
  ggplot(lambda_scenario_summary, 
         aes(x = factor(DGP), y = factor(T_val), fill = MSE_medio)) +
    geom_tile(width = 0.8, height = 0.8) +
    facet_wrap(~Model) +
    scale_fill_gradient(low = "white", high = "darkred", name = "MSE") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right") +
    labs(title = "MSE de λ̂: Comparação de Padrões",
         x = "DGP", y = "Tamanho Amostral (T)"),
  width = 10, height = 5)

# ==============================================================================
# 8. LOG-VEROSSIMILHANÇA: DISTRIBUIÇÃO POR RÉPLICA
# ==============================================================================

cat("\n--- Processando Log-Verossimilhança ---\n")

loglik_df <- bind_rows(lapply(
  names(all_results),
  function(key) {
    parts <- str_match(key, "T(\\d+)_DGP(\\d+)")
    T_val <- as.integer(parts[2])
    dgp <- as.integer(parts[3])
    
    loglik_list <- all_results[[key]]$loglik
    
    bind_rows(lapply(
      seq_along(loglik_list),
      function(rep) {
        if(length(loglik_list[[rep]]) == 0) return(NULL)
        
        data.frame(
          T_val = T_val,
          DGP = dgp,
          Replica = rep,
          Model = c("FFBS", "Static"),
          LogLik_mean = c(
            mean(colMeans(loglik_list[[rep]]$ffbs)),
            mean(colMeans(loglik_list[[rep]]$static))
          ),
          stringsAsFactors = FALSE
        )
      }
    ))
  }
))

loglik_summary <- loglik_df %>%
  group_by(T_val, DGP, Model) %>%
  summarise(
    LogLik_mean = mean(LogLik_mean, na.rm = TRUE),
    LogLik_sd = sd(LogLik_mean, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(loglik_df, "pos_processamento/loglik_by_replica.csv")

ggsave("pos_processamento/fig_06_loglik_distribution.png",
  ggplot(loglik_df, aes(x = Model, y = LogLik_mean, fill = Model)) +
    geom_boxplot(alpha = 0.7) +
    facet_grid(T_val ~ DGP, 
               labeller = labeller(
                 T_val = c("25" = "T=25", "100" = "T=100"),
                 DGP = c("1" = "DGP=1", "2" = "DGP=2", "3" = "DGP=3"))) +
    scale_fill_manual(values = cores_modelos) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom") +
    labs(title = "Log-Verossimilhança Média por Réplica",
         x = "Modelo", y = "Log-Verossimilhança Média",
         fill = "Modelo"),
  width = 12, height = 7)

# ==============================================================================
# 9. COMPARAÇÃO MSE vs COVERAGE
# ==============================================================================

ggsave("pos_processamento/fig_07_mse_vs_coverage.png",
  ggplot(lambda_scenario_summary, 
         aes(x = MSE_medio, y = Coverage_95_pct, color = Model, shape = factor(DGP))) +
    geom_point(size = 4, alpha = 0.7) +
    facet_wrap(~T_val, labeller = labeller(T_val = c("25" = "T=25", "100" = "T=100"))) +
    scale_color_manual(values = cores_modelos) +
    scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 18), name = "DGP") +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red", alpha = 0.5) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom") +
    labs(title = "Relação entre MSE e Cobertura do HPD",
         x = "MSE Médio", y = "Coverage 95% (%)",
         color = "Modelo", shape = "DGP"),
  width = 11, height = 6)

# ==============================================================================
# 10. TABELA RESUMO GERAL
# ==============================================================================

tabela_resumo <- lambda_scenario_summary %>%
  left_join(
    waic_by_scenario %>% select(T_val, DGP, Model, WAIC_mean, Delta_WAIC),
    by = c("T_val", "DGP", "Model")
  ) %>%
  arrange(T_val, DGP, Model) %>%
  select(T_val, DGP, Model, N_lambda_params,
         MSE_medio, Bias_medio, Coverage_95_pct,
         WAIC_mean, Delta_WAIC)

colnames(tabela_resumo) <- c("T", "DGP", "Modelo", "N(λ)",
                              "MSE", "Viés", "Cov.95%", "WAIC", "ΔW")

write_csv(tabela_resumo, "pos_processamento/tabela_resumo_geral.csv")

cat("\n")
print(tabela_resumo)

# ==============================================================================
# 11. RANKING DOS MODELOS
# ==============================================================================

cat("\n=== RANKING POR CRITÉRIO ===\n")

ranking <- lambda_scenario_summary %>%
  select(T_val, DGP, Model, MSE_medio) %>%
  group_by(T_val, DGP) %>%
  mutate(Ranking_MSE = rank(MSE_medio)) %>%
  ungroup() %>%
  left_join(
    waic_by_scenario %>% select(T_val, DGP, Model, WAIC_rel_ranking),
    by = c("T_val", "DGP", "Model")
  ) %>%
  rename(Ranking_WAIC = WAIC_rel_ranking) %>%
  mutate(Score = Ranking_MSE + Ranking_WAIC) %>%
  arrange(T_val, DGP, Score)

write_csv(ranking, "pos_processamento/ranking_modelos.csv")
print(ranking)

# ==============================================================================
# 12. ANÁLISE DE SENSIBILIDADE: DGP vs PERFORMANCE
# ==============================================================================

cat("\n--- Gerando visualizações de sensibilidade ---\n")

# Por DGP: como cada modelo se comporta nos 3 DGPs
dgp_analysis <- lambda_scenario_summary %>%
  pivot_longer(cols = c(MSE_medio, Coverage_95_pct),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, 
                        levels = c("MSE_medio", "Coverage_95_pct"),
                        labels = c("MSE", "Coverage (%)")))

ggsave("pos_processamento/fig_08_dgp_sensitivity.png",
  ggplot(dgp_analysis, aes(x = factor(DGP), y = Value, color = Model, group = Model)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    facet_grid(Metric ~ T_val, scales = "free_y",
               labeller = labeller(T_val = c("25" = "T=25", "100" = "T=100"))) +
    scale_color_manual(values = cores_modelos) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom") +
    labs(title = "Sensibilidade dos Modelos aos DGPs",
         x = "DGP", y = "Valor",
         color = "Modelo"),
  width = 11, height = 7)

# ==============================================================================
# RESUMO FINAL
# ==============================================================================

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("  PÓS-PROCESSAMENTO CONCLUÍDO COM SUCESSO\n")
cat("════════════════════════════════════════════════════════════════\n")
cat(sprintf("Total de cenários: %d (T×DGP)\n", nrow(scenarios)))
cat(sprintf("Total de réplicas por cenário: %d\n", n_replicas))
cat(sprintf("Total de parâmetros λ analisados: %d\n", 
            nrow(lambda_metrics_all)))
cat("\n📊 ARQUIVOS GERADOS:\n")
cat("\n📁 Dados em CSV:\n")
cat("  • lambda_metrics_detailed.csv — Métricas de cada λ̂\n")
cat("  • lambda_scenario_summary.csv — Resumo por cenário\n")
cat("  • waic_all_replicas.csv — WAIC de todas as réplicas\n")
cat("  • waic_scenario_summary.csv — WAIC por cenário\n")
cat("  • loglik_by_replica.csv — Log-verossimilhança por réplica\n")
cat("  • tabela_resumo_geral.csv — Tabela consolidada\n")
cat("  • ranking_modelos.csv — Ranking dos modelos\n")
cat("\n📈 Gráficos em PNG:\n")
cat("  • fig_01_bias_mse_comparison.png\n")
cat("  • fig_02_coverage_hpd95.png\n")
cat("  • fig_03_waic_by_scenario.png\n")
cat("  • fig_04_waic_delta.png\n")
cat("  • fig_05_mse_heatmap.png\n")
cat("  • fig_06_loglik_distribution.png\n")
cat("  • fig_07_mse_vs_coverage.png\n")
cat("  • fig_08_dgp_sensitivity.png\n")
cat("\n🎯 Local: pos_processamento/\n")
cat("════════════════════════════════════════════════════════════════\n")
