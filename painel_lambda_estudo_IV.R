# ==============================================================================
# PAINEL FINAL - NIVEL REVISTA
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)

# ==============================================================================
# CONFIGURAÇÃO DE LOCALE (acentos)
# ==============================================================================
Sys.setlocale("LC_ALL", "Portuguese_Brazil.utf8")

# ==============================================================================
# CORES (mais elegantes e equilibradas)
# ==============================================================================
cores <- c(
  "FFBS"  = "#1b9e77",
  "Fixed" = "#d95f02"
)

# ==============================================================================
# PLOT FINAL
# ==============================================================================

p <- ggplot(df_plot, aes(x = t)) +
  
  # λ verdadeiro (linha referência)
  geom_line(aes(y = lambda_true),
            color = "black",
            linewidth = 0.9,
            linetype = "dashed") +
  
  # HPD (sem legenda duplicada)
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = model),
              alpha = 0.15,
              color = NA,
              show.legend = FALSE) +
  
  # Linhas principais
  geom_line(aes(y = Mean, color = model),
            linewidth = 1) +
  
  facet_grid(Cluster ~ DGP,
             labeller = labeller(
               DGP = function(x) paste("DGP", x),
               Cluster = function(x) paste("Cluster", x)
             )) +
  
  scale_color_manual(values = cores) +
  scale_fill_manual(values = cores) +
  
  labs(
    title = "Trajetória de λ_it: modelo dinâmico vs estático",
    subtitle = "Linha tracejada: valor verdadeiro | Linhas: médias a posteriori | Faixas: intervalos HPD",
    x = "Tempo",
    y = "λ",
    color = "Modelo"
  ) +
  
  theme_minimal(base_size = 14, base_family = "serif") +
  
  theme(
    # legenda
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    
    # facet
    strip.text = element_text(face = "bold", size = 12),
    
    # grade
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    # eixos
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    
    # títulos
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    
    # margens
    plot.margin = margin(12, 18, 12, 18)
  )

# ==============================================================================
# SALVAR
# ==============================================================================

ggsave("painel_lambda_nivel_revista.png", p,
       width = 12,
       height = 10,
       dpi = 400)

print(p)