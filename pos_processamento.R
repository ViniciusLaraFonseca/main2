# ==============================================================================
# R/pos_processamento.R
# Função worker COMPARTILHADA para todos os cenários
# Melhorias implementadas:
#   1. niter automático: 50000 para T>=100 ou A>=510, 20000 caso contrário
#   2. Painel de mu_{it} = lambda_{it} * E_{it} * epsilon_i * exp(x*beta)
#   3. Epsilon médio posterior exibido nos painéis de lambda / theta / mu
#   4. Traceplots com cadeias separadas (sobrepostas) + média ergódica
#   5. Diagnóstico de ACF: correlação da cadeia e lag efetivo
#
# Uso em cada run_*.R:
#   source("R/pos_processamento.R")     # define run_model no ambiente global
#   clusterExport(cl, c(..., "run_model"))  # exporta ao cluster
# ==============================================================================

run_model <- function(model_type, output_dir) {
  
  library(nimble); library(coda)
  library(dplyr);  library(ggplot2)
  library(readr);  library(stringr)
  
  # ── MUDANÇA 1: niter automático ────────────────────────────────────────────
  is_heavy <- (n_times >= 100 || n_regions >= 510)
  niter   <- if (is_heavy) 50000 else 20000
  nburnin <- if (is_heavy) 10000 else  5000
  thin    <- if (n_regions >= 510 && n_times >= 100) 10 else 1
  cat(sprintf("\n[%s] niter=%d | nburnin=%d | thin=%d\n",
              model_type, niter, nburnin, thin))
  
  # ── FFBS ESPACIAL ──────────────────────────────────────────────────────────
  ffbs_spatial <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_regions <- control$n_regions;  n_times <- control$n_times
      p  <- control$p;  a0 <- control$a0;  b0 <- control$b0;  w <- control$w
      buf_size    <- n_regions * (n_times + 1)
      at_buf      <- nimNumeric(buf_size, 0)
      bt_buf      <- nimNumeric(buf_size, 0)
      calcNodes   <- model$getDependencies(target, self = FALSE)
      targetNodes <- model$expandNodeNames(target)
      setupOutputs(at_buf, bt_buf)
    },
    run = function() {
      declare(i, integer()); declare(t, integer()); declare(k, integer())
      declare(prod_val, double()); declare(g_it, double())
      declare(att_t,  double()); declare(btt_t,  double())
      declare(shape_tmp, double()); declare(rate_tmp, double())
      declare(lambda_futuro, double()); declare(nu, double())
      declare(idx, integer()); declare(idx_next, integer())
      for (i in 1:n_regions) {
        idx <- (i - 1) * (n_times + 1) + 1
        at_buf[idx] <<- a0;  bt_buf[idx] <<- b0
        for (t in 1:n_times) {
          idx      <- (i - 1) * (n_times + 1) + t
          idx_next <- idx + 1
          att_t <- w * at_buf[idx];  btt_t <- w * bt_buf[idx]
          prod_val <- 0
          for (k in 1:p) prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
          g_it <- model$E[i, t] * model$epsilon[i] * exp(prod_val + model$s[i])
          at_buf[idx_next] <<- att_t + model$Y[i, t]
          bt_buf[idx_next] <<- btt_t + g_it
        }
        idx       <- (i - 1) * (n_times + 1) + n_times + 1
        model$lambda[i, n_times] <<- rgamma(1, shape = at_buf[idx], rate = bt_buf[idx])
        for (t_idx in 1:(n_times - 1)) {
          t_back  <- n_times - t_idx
          idx_buf <- (i - 1) * (n_times + 1) + t_back + 1
          nu <- rgamma(1, shape = (1 - w) * at_buf[idx_buf], rate = bt_buf[idx_buf])
          model$lambda[i, t_back] <<- nu + w * model$lambda[i, t_back + 1]
        }
      }
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # ── FFBS NÃO-ESPACIAL ──────────────────────────────────────────────────────
  ffbs_nonspatial <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_regions <- control$n_regions;  n_times <- control$n_times
      p  <- control$p;  a0 <- control$a0;  b0 <- control$b0;  w <- control$w
      buf_size    <- n_regions * (n_times + 1)
      at_buf      <- nimNumeric(buf_size, 0)
      bt_buf      <- nimNumeric(buf_size, 0)
      calcNodes   <- model$getDependencies(target, self = FALSE)
      targetNodes <- model$expandNodeNames(target)
      setupOutputs(at_buf, bt_buf)
    },
    run = function() {
      declare(i, integer()); declare(t, integer()); declare(k, integer())
      declare(prod_val, double()); declare(g_it, double())
      declare(att_t,  double()); declare(btt_t,  double())
      declare(shape_tmp, double()); declare(rate_tmp, double())
      declare(lambda_futuro, double()); declare(nu, double())
      declare(idx, integer()); declare(idx_next, integer())
      for (i in 1:n_regions) {
        idx <- (i - 1) * (n_times + 1) + 1
        at_buf[idx] <<- a0;  bt_buf[idx] <<- b0
        for (t in 1:n_times) {
          idx      <- (i - 1) * (n_times + 1) + t
          idx_next <- idx + 1
          att_t <- w * at_buf[idx];  btt_t <- w * bt_buf[idx]
          prod_val <- 0
          for (k in 1:p) prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
          g_it <- model$E[i, t] * model$epsilon[i] * exp(prod_val)   # sem s[i]
          at_buf[idx_next] <<- att_t + model$Y[i, t]
          bt_buf[idx_next] <<- btt_t + g_it
        }
        idx       <- (i - 1) * (n_times + 1) + n_times + 1
        model$lambda[i, n_times] <<- rgamma(1, shape = at_buf[idx], rate = bt_buf[idx])
        for (t_idx in 1:(n_times - 1)) {
          t_back  <- n_times - t_idx
          idx_buf <- (i - 1) * (n_times + 1) + t_back + 1
          nu <- rgamma(1, shape = (1 - w) * at_buf[idx_buf], rate = bt_buf[idx_buf])
          model$lambda[i, t_back] <<- nu + w * model$lambda[i, t_back + 1]
        }
      }
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # ── Seleção do modelo ───────────────────────────────────────────────────────
  is_spatial <- (model_type == "spatial")
  model_code <- if (is_spatial) code_spatial      else code_nonspatial
  constants  <- if (is_spatial) constants_spatial  else constants_nonspatial
  inits_list <- if (is_spatial) inits_list_spatial else inits_list_nonspatial
  ffbs_fn    <- if (is_spatial) ffbs_spatial       else ffbs_nonspatial
  
  cat("\n--- Iniciando modelo:", model_type, "---\n")
  scenario_dir <- file.path(output_dir, model_type)
  dir.create(file.path(scenario_dir, "lambdas"),    recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(scenario_dir, "traceplots"), recursive = TRUE, showWarnings = FALSE)
  
  # ── Build e run MCMC ───────────────────────────────────────────────────────
  model  <- nimbleModel(code = model_code, constants = constants,
                        data = data_nimble, inits = inits_list[[1]], check = FALSE)
  Cmodel <- compileNimble(model)
  conf   <- configureMCMC(model)
  conf$removeSamplers("lambda")
  conf$addSampler(target = "lambda", type = ffbs_fn,
                  control = list(n_regions = n_regions, n_times = n_times, p = p,
                                 a0 = constants$a0, b0 = constants$b0, w = constants$w))
  conf$removeSampler("gamma")
  conf$addSampler(target = "gamma", type = "AF_slice")
  
  # Monitoring: seletivo para A>=510, completo para A<510
  monitors_base <- c("beta", "gamma", "logLik_Y")
  if (is_spatial) monitors_base <- c(monitors_base, "s", "sigma_s", "tau_s")
  conf$addMonitors(monitors_base)
  
  if (exists("LAMBDA_MONITORS")) {
    conf$addMonitors(LAMBDA_MONITORS)   # A=510: monitoramento seletivo
  } else {
    conf$addMonitors("lambda")          # A=75:  monitora tudo
  }
  conf$printSamplers()
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  nchains <- 2
  samples <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin, nchains = nchains,
                     thin = thin, inits = inits_list,
                     samplesAsCodaMCMC = TRUE, summary = FALSE, WAIC = FALSE)
  saveRDS(samples, file.path(scenario_dir, "samples.rds"))
  
  samples_mat    <- as.matrix(samples)
  mcmc_list_full <- mcmc.list(lapply(seq_len(nchains), function(ch) as.mcmc(samples[[ch]])))
  n_saved        <- nrow(as.matrix(samples[[1]]))    # iterações por cadeia (após thin)
  
  # libera memória no caso A=510
  if (n_regions >= 510) { rm(samples); gc() }
  
  # ── Helpers ─────────────────────────────────────────────────────────────────
  compute_metrics <- function(sv, tv) {
    if (var(sv) < 1e-12)
      return(data.frame(Mean = mean(sv), SD = sd(sv),
                        HPD_Lower = NA, HPD_Upper = NA,
                        Bias = mean(sv) - tv,
                        MSE  = (mean(sv) - tv)^2, Coverage = NA))
    hpd <- HPDinterval(as.mcmc(sv), prob = 0.95)
    me  <- mean(sv);  sd_e <- sd(sv)
    data.frame(Mean = me, SD = sd_e,
               HPD_Lower = hpd[1], HPD_Upper = hpd[2],
               Bias = me - tv,
               MSE  = (me - tv)^2 + sd_e^2,
               Coverage = as.integer(tv >= hpd[1] & tv <= hpd[2]))
  }
  
  safe_gelman <- function(obj) {
    tryCatch(gelman.diag(obj, autoburnin = FALSE)$psrf[, 1],
             error = function(e) rep(NA, nvar(obj)))
  }
  
  # ── Nomes ───────────────────────────────────────────────────────────────────
  beta_names  <- paste0("beta[",  seq_len(p), "]")
  gamma_names <- paste0("gamma[", seq_len(K), "]")
  h_mat       <- constants$h    # N x K
  
  # ── MUDANÇA 3: epsilon posterior por região ─────────────────────────────────
  # epsilon_i = 1 - h[i,] %*% gamma,  calculado para cada draw
  gamma_draws     <- samples_mat[, gamma_names, drop = FALSE]       # n_draw × K
  epsilon_draws   <- 1 - tcrossprod(gamma_draws, t(h_mat))          # n_draw × N  →  transposto abaixo
  # tcrossprod(A, B) = A %*% t(B), queremos [n_draw × N]:
  # epsilon_draws[d, i] = 1 - sum_k gamma_draws[d,k] * h_mat[i,k]
  epsilon_draws   <- 1 - gamma_draws %*% t(h_mat)                   # n_draw × N
  epsilon_mean    <- colMeans(epsilon_draws)                         # vetor de tamanho N
  epsilon_hpd     <- apply(epsilon_draws, 2,
                           function(x) HPDinterval(as.mcmc(x), prob = 0.95))
  epsilon_summary <- data.frame(
    Region    = seq_len(n_regions),
    Eps_Mean  = epsilon_mean,
    Eps_Lower = epsilon_hpd[1, ],
    Eps_Upper = epsilon_hpd[2, ]
  )
  write_csv(epsilon_summary, file.path(scenario_dir, "epsilon_summary.csv"))
  
  # Função para rótulo de facet com epsilon
  make_region_label <- function(regions_vec) {
    setNames(
      sprintf("Região %d\nε̂=%.3f", regions_vec, epsilon_mean[regions_vec]),
      as.character(regions_vec)
    )
  }
  
  # ── beta e gamma ────────────────────────────────────────────────────────────
  beta_metrics <- cbind(
    Parameter = beta_names,
    do.call(rbind, lapply(seq_len(p), function(j)
      compute_metrics(samples_mat[, beta_names[j]], beta_true[j])))
  )
  beta_metrics$ESS  <- effectiveSize(mcmc_list_full[, beta_names])
  beta_metrics$Rhat <- safe_gelman(mcmc_list_full[, beta_names])
  write_csv(beta_metrics, file.path(scenario_dir, "beta_metrics.csv"))
  
  gamma_metrics <- cbind(
    Parameter = gamma_names,
    do.call(rbind, lapply(seq_len(K), function(k)
      compute_metrics(samples_mat[, gamma_names[k]], gamma_true[k])))
  )
  gamma_metrics$ESS  <- effectiveSize(mcmc_list_full[, gamma_names])
  gamma_metrics$Rhat <- safe_gelman(mcmc_list_full[, gamma_names])
  write_csv(gamma_metrics, file.path(scenario_dir, "gamma_metrics.csv"))
  
  # ── Efeito espacial s (apenas spatial) ─────────────────────────────────────
  ESS_s_mean <- NA; ESS_tau <- NA; corr_s <- NA
  
  if (is_spatial) {
    s_names   <- paste0("s[", seq_len(n_regions), "]")
    idx_show  <- seq_len(min(100, n_regions))
    s_metrics <- do.call(rbind, lapply(seq_len(n_regions), function(i) {
      samp <- samples_mat[, s_names[i]]
      hpd  <- HPDinterval(as.mcmc(samp), prob = 0.95)
      data.frame(Mean = mean(samp), SD = sd(samp),
                 HPD_Lower = hpd[1], HPD_Upper = hpd[2],
                 Bias = NA, MSE = NA, Coverage = NA)
    }))
    s_metrics <- cbind(Region = seq_len(n_regions), s_metrics)
    s_metrics$ESS <- effectiveSize(mcmc_list_full[, s_names])
    write_csv(s_metrics, file.path(scenario_dir, "s_metrics.csv"))
    ESS_s_mean <- mean(s_metrics$ESS, na.rm = TRUE)
    
    tau_samp    <- samples_mat[, "tau_s"]
    hpd_tau     <- HPDinterval(as.mcmc(tau_samp), prob = 0.95)
    tau_metrics <- data.frame(Parameter = "tau_s",
                              Mean = mean(tau_samp), SD = sd(tau_samp),
                              HPD_Lower = hpd_tau[1], HPD_Upper = hpd_tau[2],
                              Bias = NA, MSE = NA, Coverage = NA)
    tau_metrics$ESS  <- effectiveSize(mcmc_list_full[, "tau_s"])
    tau_metrics$Rhat <- safe_gelman(mcmc_list_full[, "tau_s"])
    write_csv(tau_metrics, file.path(scenario_dir, "tau_metrics.csv"))
    ESS_tau <- tau_metrics$ESS
    
    ggsave(file.path(scenario_dir, "s_posterior.png"),
           ggplot(data.frame(Region  = idx_show,
                             Mean    = s_metrics$Mean[idx_show],
                             Lower   = s_metrics$HPD_Lower[idx_show],
                             Upper   = s_metrics$HPD_Upper[idx_show]),
                  aes(x = Region, y = Mean)) +
             geom_point(size = 0.8) +
             geom_errorbar(aes(ymin = Lower, ymax = Upper),
                           width = 0.3, linewidth = 0.3) +
             geom_hline(yintercept = 0, linetype = "dashed") +
             theme_bw() +
             labs(title = sprintf("Efeito espacial s (%s primeiras regiões): média posterior e HPD 95%%",
                                  length(idx_show)),
                  y = "s[i]", x = "Região"),
           width = 10, height = 5)
  }
  
  # ── Lambdas: regiões de interesse ──────────────────────────────────────────
  lambda_names_all <- grep("^lambda\\[", colnames(samples_mat), value = TRUE)
  
  # Para A<510, filtramos por regions_interest; para A>=510 já vêm filtrados
  if (n_regions < 200) {
    regions_interest <- c(1, 8, 15, 19, 22, 31, 34, 40, 46, 55, 65, n_regions)
    lambda_names_used <- lambda_names_all[
      vapply(lambda_names_all, function(nm) {
        i <- as.integer(str_match(nm, "lambda\\[(\\d+),")[, 2])
        i %in% regions_interest
      }, logical(1))
    ]
  } else {
    lambda_names_used <- lambda_names_all   # monitoramento já foi seletivo
  }
  
  ESS_lambda_mean <- mean(effectiveSize(mcmc_list_full[, lambda_names_used]),
                          na.rm = TRUE)
  
  beta_cols <- grep("^beta\\[", colnames(samples_mat), value = TRUE)
  
  # ── Lambda summary ──────────────────────────────────────────────────────────
  lambda_summary <- do.call(rbind, lapply(lambda_names_used, function(nm) {
    idx <- str_match(nm, "lambda\\[(\\d+),\\s*(\\d+)\\]")
    i   <- as.integer(idx[2]);  t <- as.integer(idx[3])
    hpd <- HPDinterval(as.mcmc(samples_mat[, nm]))
    data.frame(Region = i, Time = t,
               True   = lambda_true[i, t],
               Mean   = mean(samples_mat[, nm]),
               Lower  = hpd[1], Upper = hpd[2],
               model  = model_type)
  }))
  write_csv(lambda_summary, file.path(scenario_dir, "lambda_selected.csv"))
  
  # ── Theta summary ───────────────────────────────────────────────────────────
  theta_summary <- do.call(rbind, lapply(lambda_names_used, function(nm) {
    idx <- str_match(nm, "lambda\\[(\\d+),\\s*(\\d+)\\]")
    i   <- as.integer(idx[2]);  t <- as.integer(idx[3])
    if (t > n_times) return(NULL)
    ldraws   <- samples_mat[, nm]
    bdraws   <- samples_mat[, beta_cols, drop = FALSE]
    x_it     <- data_nimble$x[i, t, ]
    lin_pred <- as.vector(bdraws %*% x_it)
    theta    <- ldraws * exp(lin_pred)
    hpd      <- HPDinterval(as.mcmc(theta))
    tv       <- lambda_true[i, t] * exp(sum(x_it * beta_true))
    data.frame(Region = i, Time = t, True = tv,
               Mean = mean(theta), Lower = hpd[1], Upper = hpd[2],
               model = model_type)
  }))
  write_csv(theta_summary, file.path(scenario_dir, "theta_selected.csv"))
  
  # ── MUDANÇA 2: Mu summary ───────────────────────────────────────────────────
  # mu_{it} = lambda_{it} * E_{it} * epsilon_i * exp(x_{it}^T beta)
  #         = theta_{it} * E_{it} * epsilon_i
  mu_summary <- do.call(rbind, lapply(lambda_names_used, function(nm) {
    idx <- str_match(nm, "lambda\\[(\\d+),\\s*(\\d+)\\]")
    i   <- as.integer(idx[2]);  t <- as.integer(idx[3])
    if (t > n_times) return(NULL)
    ldraws         <- samples_mat[, nm]
    bdraws         <- samples_mat[, beta_cols, drop = FALSE]
    x_it           <- data_nimble$x[i, t, ]
    lin_pred       <- as.vector(bdraws %*% x_it)
    epsilon_i_draw <- epsilon_draws[, i]         # já calculado acima
    mu_draws       <- ldraws * exp(lin_pred) * data_nimble$E[i, t] * epsilon_i_draw
    hpd            <- HPDinterval(as.mcmc(mu_draws))
    mu_true_val    <- lambda_true[i, t] * exp(sum(x_it * beta_true)) *
      data_nimble$E[i, t] * (1 - sum(h_mat[i, ] * gamma_true))
    data.frame(Region = i, Time = t, True = mu_true_val,
               Mean = mean(mu_draws), Lower = hpd[1], Upper = hpd[2],
               model = model_type)
  }))
  write_csv(mu_summary, file.path(scenario_dir, "mu_selected.csv"))
  
  # ── WAIC / LPML ─────────────────────────────────────────────────────────────
  loglik_names <- grep("logLik_Y", colnames(samples_mat), value = TRUE)
  waic <- NA;  LPML <- NA
  if (length(loglik_names) > 0) {
    lm    <- samples_mat[, loglik_names, drop = FALSE]
    lppd  <- sum(apply(lm, 2, function(x) { mx <- max(x); mx + log(mean(exp(x - mx))) }))
    p_waic <- sum(apply(lm, 2, var))
    waic  <- -2 * (lppd - p_waic)
    LPML  <- sum(log(1 / apply(lm, 2, function(x) mean(exp(-x)))))
    write_csv(data.frame(WAIC = waic, LPML = LPML, lppd = lppd, pWAIC = p_waic),
              file.path(scenario_dir, "criteria.csv"))
  }
  
  # ── ESS / Rhat estruturais ───────────────────────────────────────────────────
  params_struct <- c(beta_names, gamma_names)
  if (is_spatial) params_struct <- c(params_struct, "tau_s")
  ESS_struct  <- effectiveSize(mcmc_list_full[, params_struct])
  Rhat_struct <- safe_gelman(mcmc_list_full[, params_struct])
  
  # ── MUDANÇA 5: Diagnóstico de ACF ───────────────────────────────────────────
  acf_results <- do.call(rbind, lapply(params_struct, function(nm) {
    ac  <- acf(samples_mat[, nm], lag.max = 200, plot = FALSE)
    lags <- as.vector(ac$lag[-1])
    acfs <- as.vector(ac$acf[-1])
    # Lag onde |acf| cai abaixo dos limiares
    lag_01  <- lags[which(abs(acfs) < 0.10)[1]]
    lag_005 <- lags[which(abs(acfs) < 0.05)[1]]
    data.frame(
      Parameter = nm,
      ESS       = ESS_struct[nm],
      Rhat      = Rhat_struct[nm],
      lag_0.10  = if (is.na(lag_01))  Inf else lag_01,
      lag_0.05  = if (is.na(lag_005)) Inf else lag_005,
      acf_lag1  = acfs[1],
      acf_lag10 = acfs[min(10, length(acfs))]
    )
  }))
  write_csv(acf_results, file.path(scenario_dir, "acf_diagnostics.csv"))
  cat(sprintf("[%s] Maior lag(|acf|<0.10): %g | lag(|acf|<0.05): %g\n",
              model_type,
              max(acf_results$lag_0.10, na.rm = TRUE),
              max(acf_results$lag_0.05, na.rm = TRUE)))
  
  # ACF plot dos parâmetros estruturais
  acf_df <- do.call(rbind, lapply(params_struct, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    data.frame(Parameter = nm, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(file.path(scenario_dir, "acf_params.png"),
         ggplot(acf_df, aes(x = Lag, y = ACF)) +
           geom_col(width = 0.6, fill = "grey50") +
           geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed",
                      color = "blue", linewidth = 0.5) +
           geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted",
                      color = "red",  linewidth = 0.5) +
           facet_wrap(~Parameter, scales = "free_y") +
           theme_bw(base_size = 11) +
           labs(title  = paste("ACF dos parâmetros estruturais (", model_type, ")"),
                subtitle = "Azul tracejado: |0.10| | Vermelho pontilhado: |0.05|"),
         width = 10, height = 6)
  
  # ── Gráficos: densidades beta / gamma ───────────────────────────────────────
  ggsave(file.path(scenario_dir, "dens_beta.png"),
         ggplot(data.frame(Value     = as.vector(samples_mat[, beta_names]),
                           Parameter = rep(beta_names, each = nrow(samples_mat))),
                aes(x = Value)) +
           geom_density(fill = "grey70") +
           facet_wrap(~Parameter, scales = "free") +
           theme_bw() +
           labs(title = paste("Posterior densities — beta (", model_type, ")")),
         width = 8, height = 6)
  
  ggsave(file.path(scenario_dir, "dens_gamma.png"),
         ggplot(data.frame(Value     = as.vector(samples_mat[, gamma_names]),
                           Parameter = rep(gamma_names, each = nrow(samples_mat))),
                aes(x = Value)) +
           geom_density(fill = "grey70") +
           facet_wrap(~Parameter, scales = "free") +
           theme_bw() +
           labs(title = paste("Posterior densities — gamma (", model_type, ")")),
         width = 8, height = 6)
  
  # ── MUDANÇA 3+4: Lambda painel com epsilon nos rótulos ─────────────────────
  lambda_regs <- sort(unique(lambda_summary$Region))
  ggsave(file.path(scenario_dir, "painel_lambda.png"),
         ggplot(lambda_summary, aes(x = Time)) +
           geom_ribbon(aes(ymin = Lower, ymax = Upper),
                       fill = "grey70", alpha = 0.5) +
           geom_line(aes(y = Mean),  color = "black") +
           geom_line(aes(y = True),  color = "red", linetype = "dashed") +
           facet_wrap(~Region, scales = "free_y", ncol = 3,
                      labeller = labeller(Region = make_region_label(lambda_regs))) +
           theme_bw(base_size = 10) +
           labs(title    = paste("λ estimado (", model_type, ")"),
                subtitle = "Vermelho tracejado: valor verdadeiro | ε̂ no título da célula",
                x = "Tempo", y = expression(lambda[i * t])),
         width = 14, height = 10)
  
  # ── Theta painel com epsilon nos rótulos ───────────────────────────────────
  theta_regs <- sort(unique(theta_summary$Region))
  ggsave(file.path(scenario_dir, "painel_theta.png"),
         ggplot(theta_summary, aes(x = Time)) +
           geom_ribbon(aes(ymin = Lower, ymax = Upper),
                       fill = "steelblue", alpha = 0.3) +
           geom_line(aes(y = Mean),  color = "steelblue") +
           geom_line(aes(y = True),  color = "red", linetype = "dashed") +
           facet_wrap(~Region, scales = "free_y", ncol = 3,
                      labeller = labeller(Region = make_region_label(theta_regs))) +
           theme_bw(base_size = 10) +
           labs(title    = paste("θ estimado (", model_type, ")"),
                subtitle = "Vermelho tracejado: valor verdadeiro",
                x = "Tempo", y = expression(theta[i * t])),
         width = 14, height = 10)
  
  # ── MUDANÇA 2: Mu painel ────────────────────────────────────────────────────
  mu_regs <- sort(unique(mu_summary$Region))
  ggsave(file.path(scenario_dir, "painel_mu.png"),
         ggplot(mu_summary, aes(x = Time)) +
           geom_ribbon(aes(ymin = Lower, ymax = Upper),
                       fill = "darkorange", alpha = 0.3) +
           geom_line(aes(y = Mean),  color = "darkorange") +
           geom_line(aes(y = True),  color = "red", linetype = "dashed") +
           facet_wrap(~Region, scales = "free_y", ncol = 3,
                      labeller = labeller(Region = make_region_label(mu_regs))) +
           theme_bw(base_size = 10) +
           labs(title    = paste("μ estimado (", model_type, ")"),
                subtitle = expression(paste(mu[it] == lambda[it] %.% E[it] %.% epsilon[i] %.% e^{x[it]^T * beta})),
                x = "Tempo", y = expression(mu[i * t])),
         width = 14, height = 10)
  
  # ── MUDANÇA 4: Traceplots com cadeias separadas + média ergódica ─────────
  trace_params <- c(beta_names, gamma_names)
  if (is_spatial) trace_params <- c(trace_params, "tau_s")
  
  df_trace <- do.call(rbind, lapply(seq_len(nchains), function(ch) {
    cm <- as.matrix(mcmc_list_full[[ch]])
    do.call(rbind, lapply(trace_params, function(nm) {
      vals     <- cm[, nm]
      erg_mean <- cumsum(vals) / seq_along(vals)
      data.frame(
        Iter      = seq_along(vals),
        Value     = vals,
        ErgMedia  = erg_mean,
        Parameter = nm,
        Cadeia    = paste0("Cadeia ", ch),
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  cores_cadeia <- c("Cadeia 1" = "#2166AC", "Cadeia 2" = "#D6604D")
  
  ggsave(file.path(scenario_dir, "traceplots.png"),
         ggplot(df_trace, aes(x = Iter, color = Cadeia, fill = Cadeia)) +
           # traço fino para a cadeia completa
           geom_line(aes(y = Value),    alpha = 0.30, linewidth = 0.25) +
           # linha grossa para a média ergódica
           geom_line(aes(y = ErgMedia), alpha = 0.90, linewidth = 0.80,
                     linetype = "solid") +
           scale_color_manual(values = cores_cadeia) +
           scale_fill_manual(values  = cores_cadeia) +
           facet_wrap(~Parameter, scales = "free_y") +
           theme_bw(base_size = 11) +
           theme(legend.position = "bottom") +
           labs(title    = paste("Traceplots + Média Ergódica (", model_type, ")"),
                subtitle = "Linha grossa = média ergódica | Linha fina = cadeia",
                x = "Iteração (pós-burnin)", y = "Valor",
                color = "Cadeia", fill = "Cadeia"),
         width = 12, height = max(6, 3 * ceiling(length(trace_params) / 3)))
  
  # ── Epsilon painel ───────────────────────────────────────────────────────────
  eps_show <- min(100, n_regions)
  ggsave(file.path(scenario_dir, "epsilon_posterior.png"),
         ggplot(epsilon_summary[seq_len(eps_show), ],
                aes(x = Region, y = Eps_Mean)) +
           geom_point(size = 0.8) +
           geom_errorbar(aes(ymin = Eps_Lower, ymax = Eps_Upper),
                         width = 0.3, linewidth = 0.3) +
           geom_hline(yintercept = 1 - sum(gamma_true * colMeans(h_mat)),
                      color = "red", linetype = "dashed") +
           theme_bw() +
           labs(title    = sprintf("ε posterior (%d primeiras regiões)", eps_show),
                subtitle = "Vermelho tracejado: valor médio verdadeiro",
                y = expression(epsilon[i]), x = "Região"),
         width = 10, height = 4)
  
  # ── Retorno ──────────────────────────────────────────────────────────────────
  data.frame(
    model           = model_type,
    niter           = niter,
    nburnin         = nburnin,
    thin            = thin,
    WAIC            = waic,
    LPML            = LPML,
    corr_s          = corr_s,
    MSE_beta        = mean(beta_metrics$MSE,  na.rm = TRUE),
    MSE_gamma       = mean(gamma_metrics$MSE, na.rm = TRUE),
    Coverage_beta   = mean(beta_metrics$Coverage,  na.rm = TRUE),
    Coverage_gamma  = mean(gamma_metrics$Coverage, na.rm = TRUE),
    ESS_beta_min    = min(beta_metrics$ESS,  na.rm = TRUE),
    ESS_gamma_min   = min(gamma_metrics$ESS, na.rm = TRUE),
    ESS_tau         = ESS_tau,
    ESS_s_mean      = ESS_s_mean,
    ESS_lambda_mean = ESS_lambda_mean,
    ESS_global_min  = min(ESS_struct, na.rm = TRUE),
    Rhat_max        = max(Rhat_struct, na.rm = TRUE),
    lag_max_0.10    = max(acf_results$lag_0.10, na.rm = TRUE),
    lag_max_0.05    = max(acf_results$lag_0.05, na.rm = TRUE)
  )
}
