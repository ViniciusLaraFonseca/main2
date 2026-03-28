# ==============================================================================
# run_spatial_vs_nonspatial_w_unif_T_100_A_510_parallel.R
# A=510: seleção de 12 regiões (3 por cluster) e monitoramento seletivo
# ==============================================================================
rm(list = ls())
inicio_global <- Sys.time()
setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Simulacao/main/main2")

pkgs <- c("nimble","coda","parallel","dplyr","ggplot2","tidyr","readr","stringr")
for(pkg in pkgs) { if(!require(pkg,character.only=TRUE)) install.packages(pkg); library(pkg,character.only=TRUE) }
Sys.setenv(OMP_NUM_THREADS="1"); Sys.setenv(MKL_NUM_THREADS="1")
if(requireNamespace("RhpcBLASctl",quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(1)

# ------ 1. Dados ------
source("Geracao_w_unif_T_100_A_510.R")

n_regions <- constants_nimble$n_regions; n_times <- constants_nimble$n_times
p <- constants_nimble$p;                  K       <- constants_nimble$K
beta_true   <- get("beta_true",   envir=.GlobalEnv)
gamma_true  <- get("gamma_true",  envir=.GlobalEnv)
lambda_true <- get("lambda_true", envir=.GlobalEnv)
cat("Dados carregados. n_regions =",n_regions,", n_times =",n_times,"\n")

# --- Seleção de 12 regiões: 3 de cada um dos K clusters (K = 4) ---
h_mat <- constants_nimble$h
cluster_ids <- apply(h_mat, 1, sum)   # soma dos uns = número do cluster
set.seed(123)
regions_per_cluster <- lapply(1:K, function(cl) {
  regs <- which(cluster_ids == cl)
  if(length(regs) >= 3) sample(regs, 3) else regs
})
REGIONS_INTEREST <- unlist(regions_per_cluster[1:K])
cat("Regiões selecionadas:\n"); print(sort(REGIONS_INTEREST))

LAMBDA_MONITORS  <- unlist(lapply(REGIONS_INTEREST,
                                  function(r) paste0("lambda[",r,", ",seq_len(n_times),"]")))

# ------ 2a. Modelo ESPACIAL ------
code_spatial <- nimbleCode({
  for(j in 1:p) beta[j] ~ dnorm(mu_beta[j], sd=10)
  gamma[1] ~ dunif(min=a_unif, max=b_unif)
  for(j in 2:K) gamma[j] ~ dunif(min=0, max=(1 - sum(gamma[1:(j-1)])))
  sigma_s ~ T(dt(0,1,1),0,)
  tau_s   <- 1/(sigma_s^2)
  s[1:n_regions] ~ dcar_normal(adj[1:n_adj],weights[1:n_adj],num[1:n_regions],tau_s,zero_mean=1)
  for(i in 1:n_regions) {
    epsilon[i] <- 1 - inprod(h[i,1:K],gamma[1:K])
    for(t in 1:n_times) {
      lambda[i,t] ~ dgamma(1,1)
      log(mu[i,t]) <- log(lambda[i,t])+log(E[i,t])+log(epsilon[i])+inprod(beta[1:p],x[i,t,1:p])+s[i]
      Y[i,t]        ~ dpois(mu[i,t])
      logLik_Y[i,t] <- dpois(Y[i,t],mu[i,t],log=TRUE)
    }
  }
})

# ------ 2b. Modelo NÃO-ESPACIAL ------
code_nonspatial <- nimbleCode({
  for(j in 1:p) beta[j] ~ dnorm(mu_beta[j], sd=10)
  gamma[1] ~ dunif(min=a_unif, max=b_unif)
  for(j in 2:K) gamma[j] ~ dunif(min=0, max=(1 - sum(gamma[1:(j-1)])))
  for(i in 1:n_regions) {
    epsilon[i] <- 1 - inprod(h[i,1:K],gamma[1:K])
    for(t in 1:n_times) {
      lambda[i,t] ~ dgamma(1,1)
      log(mu[i,t]) <- log(lambda[i,t])+log(E[i,t])+log(epsilon[i])+inprod(beta[1:p],x[i,t,1:p])
      Y[i,t]        ~ dpois(mu[i,t])
      logLik_Y[i,t] <- dpois(Y[i,t],mu[i,t],log=TRUE)
    }
  }
})

# ------ 3. Constantes e inits ------
constants_spatial    <- constants_nimble
constants_nonspatial <- constants_nimble[setdiff(names(constants_nimble),c("adj","num","weights","n_adj"))]

inits_list_spatial <- list(
  list(beta=beta_true,gamma=gamma_true,lambda=lambda_true,sigma_s=0.5,s=rep(0,n_regions)),
  list(beta=rnorm(p,0,0.5),gamma=gamma_true*0.9,
       lambda=matrix(rgamma(n_regions*n_times,1,1),nrow=n_regions),sigma_s=1,s=rep(0,n_regions))
)
inits_list_nonspatial <- list(
  list(beta=beta_true,gamma=gamma_true,lambda=lambda_true),
  list(beta=rnorm(p,0,0.5),gamma=gamma_true*0.9,
       lambda=matrix(rgamma(n_regions*n_times,1,1),nrow=n_regions))
)

# ------ 4. Função worker ------
run_model <- function(model_type, output_dir) {
  library(nimble); library(coda); library(dplyr)
  library(ggplot2); library(readr); library(stringr)
  
  ffbs_spatial <- nimbleFunction(
    contains=sampler_BASE,
    setup=function(model,mvSaved,target,control){
      n_regions<-control$n_regions; n_times<-control$n_times; p<-control$p
      a0<-control$a0; b0<-control$b0; w<-control$w
      buf_size<-n_regions*(n_times+1)
      at_buf<-nimNumeric(buf_size,0); bt_buf<-nimNumeric(buf_size,0)
      calcNodes<-model$getDependencies(target,self=FALSE)
      targetNodes<-model$expandNodeNames(target)
      setupOutputs(at_buf,bt_buf)
    },
    run=function(){
      declare(i,integer()); declare(t,integer()); declare(k,integer())
      declare(prod_val,double()); declare(g_it,double())
      declare(att_t,double()); declare(btt_t,double())
      declare(shape_tmp,double()); declare(rate_tmp,double())
      declare(lambda_futuro,double()); declare(nu,double())
      declare(idx,integer()); declare(idx_next,integer())
      for(i in 1:n_regions){
        idx<-(i-1)*(n_times+1)+1; at_buf[idx]<<-a0; bt_buf[idx]<<-b0
        for(t in 1:n_times){
          idx<-(i-1)*(n_times+1)+t; idx_next<-idx+1
          att_t<-w*at_buf[idx]; btt_t<-w*bt_buf[idx]
          prod_val<-0
          for(k in 1:p) prod_val<-prod_val+model$x[i,t,k]*model$beta[k]
          g_it<-model$E[i,t]*model$epsilon[i]*exp(prod_val+model$s[i])
          at_buf[idx_next]<<-att_t+model$Y[i,t]; bt_buf[idx_next]<<-btt_t+g_it
        }
        idx<-(i-1)*(n_times+1)+n_times+1
        model$lambda[i,n_times]<<-rgamma(1,shape=at_buf[idx],rate=bt_buf[idx])
        for(t_idx in 1:(n_times-1)){
          t_back<-n_times-t_idx; idx_buf<-(i-1)*(n_times+1)+t_back+1
          nu<-rgamma(1,shape=(1-w)*at_buf[idx_buf],rate=bt_buf[idx_buf])
          model$lambda[i,t_back]<<-nu+w*model$lambda[i,t_back+1]
        }
      }
      model$calculate(calcNodes)
      copy(from=model,to=mvSaved,row=1,nodes=targetNodes,logProb=TRUE)
    },
    methods=list(reset=function(){})
  )
  
  ffbs_nonspatial <- nimbleFunction(
    contains=sampler_BASE,
    setup=function(model,mvSaved,target,control){
      n_regions<-control$n_regions; n_times<-control$n_times; p<-control$p
      a0<-control$a0; b0<-control$b0; w<-control$w
      buf_size<-n_regions*(n_times+1)
      at_buf<-nimNumeric(buf_size,0); bt_buf<-nimNumeric(buf_size,0)
      calcNodes<-model$getDependencies(target,self=FALSE)
      targetNodes<-model$expandNodeNames(target)
      setupOutputs(at_buf,bt_buf)
    },
    run=function(){
      declare(i,integer()); declare(t,integer()); declare(k,integer())
      declare(prod_val,double()); declare(g_it,double())
      declare(att_t,double()); declare(btt_t,double())
      declare(shape_tmp,double()); declare(rate_tmp,double())
      declare(lambda_futuro,double()); declare(nu,double())
      declare(idx,integer()); declare(idx_next,integer())
      for(i in 1:n_regions){
        idx<-(i-1)*(n_times+1)+1; at_buf[idx]<<-a0; bt_buf[idx]<<-b0
        for(t in 1:n_times){
          idx<-(i-1)*(n_times+1)+t; idx_next<-idx+1
          att_t<-w*at_buf[idx]; btt_t<-w*bt_buf[idx]
          prod_val<-0
          for(k in 1:p) prod_val<-prod_val+model$x[i,t,k]*model$beta[k]
          g_it<-model$E[i,t]*model$epsilon[i]*exp(prod_val)   # sem s[i]
          at_buf[idx_next]<<-att_t+model$Y[i,t]; bt_buf[idx_next]<<-btt_t+g_it
        }
        idx<-(i-1)*(n_times+1)+n_times+1
        model$lambda[i,n_times]<<-rgamma(1,shape=at_buf[idx],rate=bt_buf[idx])
        for(t_idx in 1:(n_times-1)){
          t_back<-n_times-t_idx; idx_buf<-(i-1)*(n_times+1)+t_back+1
          nu<-rgamma(1,shape=(1-w)*at_buf[idx_buf],rate=bt_buf[idx_buf])
          model$lambda[i,t_back]<<-nu+w*model$lambda[i,t_back+1]
        }
      }
      model$calculate(calcNodes)
      copy(from=model,to=mvSaved,row=1,nodes=targetNodes,logProb=TRUE)
    },
    methods=list(reset=function(){})
  )
  
  is_spatial <- (model_type=="spatial")
  model_code <- if(is_spatial) code_spatial     else code_nonspatial
  constants  <- if(is_spatial) constants_spatial else constants_nonspatial
  inits_list <- if(is_spatial) inits_list_spatial else inits_list_nonspatial
  ffbs_fn    <- if(is_spatial) ffbs_spatial      else ffbs_nonspatial
  
  cat("\n--- Iniciando modelo:",model_type,"---\n")
  scenario_dir <- file.path(output_dir,model_type)
  dir.create(file.path(scenario_dir,"lambdas"),   recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path(scenario_dir,"traceplots"),recursive=TRUE,showWarnings=FALSE)
  
  model  <- nimbleModel(code=model_code,constants=constants,data=data_nimble,
                        inits=inits_list[[1]],check=FALSE)
  Cmodel <- compileNimble(model)
  conf   <- configureMCMC(model)
  conf$removeSamplers("lambda")
  conf$addSampler(target="lambda",type=ffbs_fn,
                  control=list(n_regions=n_regions,n_times=n_times,p=p,
                               a0=constants$a0,b0=constants$b0,w=constants$w))
  conf$removeSampler("gamma")
  conf$addSampler(target="gamma",type="AF_slice")
  
  # Monitoring seletivo — apenas 12 regiões × n_times lambdas
  monitors_base <- c("beta","gamma","logLik_Y")
  if(is_spatial) monitors_base <- c(monitors_base,"s","sigma_s","tau_s")
  conf$addMonitors(monitors_base)
  conf$addMonitors(LAMBDA_MONITORS)
  conf$printSamplers()
  
  Rmcmc <- buildMCMC(conf); Cmcmc <- compileNimble(Rmcmc,project=model)
  # ── Parâmetros para T=100 ───────────────────────────────────────────────
  niter   <- 50000
  nburnin <- 10000
  nchains <- 2
  thin    <- 10
  cat(sprintf("[%s] niter=%d | nburnin=%d | thin=%d\n", model_type, niter, nburnin, thin))
  
  samples <- runMCMC(Cmcmc,niter=niter,nburnin=nburnin,nchains=nchains,thin=thin,
                     inits=inits_list,samplesAsCodaMCMC=TRUE,summary=FALSE,WAIC=FALSE)
  saveRDS(samples,file.path(scenario_dir,"samples.rds"))
  
  samples_mat    <- as.matrix(samples)
  mcmc_list_full <- mcmc.list(lapply(1:nchains,function(ch) as.mcmc(samples[[ch]])))
  rm(samples); gc()
  
  compute_metrics <- function(sv,tv){
    if(var(sv)<1e-12) return(data.frame(Mean=mean(sv),SD=sd(sv),HPD_Lower=NA,HPD_Upper=NA,
                                        Bias=mean(sv)-tv,MSE=(mean(sv)-tv)^2,Coverage=NA))
    hpd<-HPDinterval(as.mcmc(sv),prob=0.95); me<-mean(sv); sd_e<-sd(sv)
    data.frame(Mean=me,SD=sd_e,HPD_Lower=hpd[1],HPD_Upper=hpd[2],
               Bias=me-tv,MSE=(me-tv)^2+sd_e^2,Coverage=as.integer(tv>=hpd[1]&tv<=hpd[2]))
  }
  safe_gelman <- function(obj) tryCatch(gelman.diag(obj,autoburnin=FALSE)$psrf[,1],
                                        error=function(e) rep(NA,nvar(obj)))
  
  beta_names  <- paste0("beta[",1:p,"]")
  gamma_names <- paste0("gamma[",1:K,"]")
  
  # epsilon posterior por cluster (garantindo igualdade dentro do cluster)
  h_mat         <- constants$h                          # N x K
  cluster_ids   <- apply(h_mat, 1, sum)                 # vetor de clusters (1..K)
  gamma_draws   <- samples_mat[,gamma_names,drop=FALSE] # n_draw x K
  epsilon_draws <- 1 - gamma_draws %*% t(h_mat)         # n_draw x N
  
  # Média posterior por cluster
  epsilon_by_cluster <- tapply(1:n_regions, cluster_ids, function(regs) {
    colMeans(epsilon_draws[, regs, drop = FALSE])
  })
  epsilon_cluster_mean <- sapply(epsilon_by_cluster, function(x) x[1])
  epsilon_mean <- epsilon_cluster_mean[cluster_ids]
  
  # HPD por cluster (usando todas as regiões do cluster)
  epsilon_hpd <- matrix(NA, nrow = 2, ncol = n_regions)
  for(cl in 1:K) {
    regs <- which(cluster_ids == cl)
    if(length(regs) > 0) {
      combined_draws <- as.vector(epsilon_draws[, regs])
      hpd_int <- HPDinterval(as.mcmc(combined_draws), prob = 0.95)
      epsilon_hpd[, regs] <- hpd_int
    }
  }
  
  epsilon_summary <- data.frame(Region = 1:n_regions,
                                Cluster = cluster_ids,
                                Eps_Mean = epsilon_mean,
                                Eps_Lower = epsilon_hpd[1, ],
                                Eps_Upper = epsilon_hpd[2, ])
  write_csv(epsilon_summary, file.path(scenario_dir, "epsilon_summary.csv"))
  
  # Rótulo de facet com epsilon médio e cluster
  make_region_label <- function(regions_vec) {
    setNames(sprintf("Região %d (C%d)\nε̂=%.3f", 
                     regions_vec, 
                     cluster_ids[regions_vec], 
                     epsilon_mean[regions_vec]),
             as.character(regions_vec))
  }
  
  beta_metrics  <- cbind(Parameter=beta_names,
                         do.call(rbind,lapply(1:p,function(j) compute_metrics(samples_mat[,beta_names[j]],beta_true[j]))))
  beta_metrics$ESS  <- effectiveSize(mcmc_list_full[,beta_names])
  beta_metrics$Rhat <- safe_gelman(mcmc_list_full[,beta_names])
  write_csv(beta_metrics,file.path(scenario_dir,"beta_metrics.csv"))
  
  gamma_metrics <- cbind(Parameter=gamma_names,
                         do.call(rbind,lapply(1:K,function(k) compute_metrics(samples_mat[,gamma_names[k]],gamma_true[k]))))
  gamma_metrics$ESS  <- effectiveSize(mcmc_list_full[,gamma_names])
  gamma_metrics$Rhat <- safe_gelman(mcmc_list_full[,gamma_names])
  write_csv(gamma_metrics,file.path(scenario_dir,"gamma_metrics.csv"))
  
  ESS_s_mean<-NA; ESS_tau<-NA; corr_s<-NA
  if(is_spatial){
    s_names   <- paste0("s[",1:n_regions,"]")
    s_metrics <- do.call(rbind,lapply(1:n_regions,function(i){
      samp<-samples_mat[,s_names[i]]; hpd<-HPDinterval(as.mcmc(samp),prob=0.95)
      data.frame(Mean=mean(samp),SD=sd(samp),HPD_Lower=hpd[1],HPD_Upper=hpd[2],Bias=NA,MSE=NA,Coverage=NA)
    }))
    s_metrics     <- cbind(Region=1:n_regions,s_metrics)
    s_metrics$ESS <- effectiveSize(mcmc_list_full[,s_names])
    write_csv(s_metrics,file.path(scenario_dir,"s_metrics.csv"))
    ESS_s_mean <- mean(s_metrics$ESS,na.rm=TRUE)
    tau_samp   <- samples_mat[,"tau_s"]; hpd_t<-HPDinterval(as.mcmc(tau_samp),prob=0.95)
    tau_m      <- data.frame(Parameter="tau_s",Mean=mean(tau_samp),SD=sd(tau_samp),
                             HPD_Lower=hpd_t[1],HPD_Upper=hpd_t[2],Bias=NA,MSE=NA,Coverage=NA)
    tau_m$ESS  <- effectiveSize(mcmc_list_full[,"tau_s"])
    tau_m$Rhat <- safe_gelman(mcmc_list_full[,"tau_s"])
    write_csv(tau_m,file.path(scenario_dir,"tau_metrics.csv"))
    ESS_tau <- tau_m$ESS
    idx100 <- seq_len(min(100,n_regions))
    ggsave(file.path(scenario_dir,"s_posterior.png"),
           ggplot(data.frame(Region=idx100,Mean=s_metrics$Mean[idx100],
                             Lower=s_metrics$HPD_Lower[idx100],Upper=s_metrics$HPD_Upper[idx100]),
                  aes(x=Region,y=Mean))+
             geom_point(size=0.8)+geom_errorbar(aes(ymin=Lower,ymax=Upper),width=0.3,linewidth=0.3)+
             geom_hline(yintercept=0,linetype="dashed")+theme_bw()+
             labs(title="Efeito espacial s (primeiras 100 regiões): média posterior e HPD 95%",
                  y="s[i]",x="Região"),width=10,height=5)
  }
  
  # Lambdas — usar LAMBDA_MONITORS
  lambda_names    <- LAMBDA_MONITORS
  ESS_lambda_mean <- mean(effectiveSize(mcmc_list_full[,lambda_names]),na.rm=TRUE)
  
  beta_cols <- grep("^beta\\[",colnames(samples_mat),value=TRUE)
  
  lambda_summary <- data.frame()
  for(nm in lambda_names){
    idx<-str_match(nm,"lambda\\[(\\d+),\\s*(\\d+)\\]"); i<-as.numeric(idx[2]); t<-as.numeric(idx[3])
    hpd<-HPDinterval(as.mcmc(samples_mat[,nm]))
    lambda_summary<-rbind(lambda_summary,data.frame(Region=i,Time=t,True=lambda_true[i,t],
                                                    Mean=mean(samples_mat[,nm]),Lower=hpd[1],Upper=hpd[2],model=model_type))
  }
  write_csv(lambda_summary,file.path(scenario_dir,"lambda_selected.csv"))
  
  theta_summary <- data.frame()
  for(nm in lambda_names){
    idx<-str_match(nm,"lambda\\[(\\d+),\\s*(\\d+)\\]"); i<-as.numeric(idx[2]); t<-as.numeric(idx[3])
    ldraws<-samples_mat[,nm]; bdraws<-samples_mat[,beta_cols,drop=FALSE]
    x_it<-data_nimble$x[i,t,]; lin<-as.vector(bdraws %*% x_it); theta<-ldraws*exp(lin)
    hpd<-HPDinterval(as.mcmc(theta)); tv<-lambda_true[i,t]*exp(sum(x_it*beta_true))
    theta_summary<-rbind(theta_summary,data.frame(Region=i,Time=t,True=tv,
                                                  Mean=mean(theta),Lower=hpd[1],Upper=hpd[2],model=model_type))
  }
  write_csv(theta_summary,file.path(scenario_dir,"theta_selected.csv"))
  
  # Mu posterior
  mu_summary <- data.frame()
  for(nm in lambda_names){
    idx<-str_match(nm,"lambda\\[(\\d+),\\s*(\\d+)\\]"); i<-as.numeric(idx[2]); t<-as.numeric(idx[3])
    ldraws         <- samples_mat[,nm]
    bdraws         <- samples_mat[,beta_cols,drop=FALSE]
    x_it           <- data_nimble$x[i,t,]
    lin            <- as.vector(bdraws %*% x_it)
    epsilon_i_draw <- epsilon_draws[,i]
    mu_draws       <- ldraws * exp(lin) * data_nimble$E[i,t] * epsilon_i_draw
    hpd            <- HPDinterval(as.mcmc(mu_draws))
    mu_true_val    <- lambda_true[i,t] * exp(sum(x_it*beta_true)) *
      data_nimble$E[i,t] * (1 - sum(h_mat[i,]*gamma_true))
    mu_summary <- rbind(mu_summary, data.frame(Region=i,Time=t,True=mu_true_val,
                                               Mean=mean(mu_draws),Lower=hpd[1],Upper=hpd[2],
                                               model=model_type))
  }
  write_csv(mu_summary,file.path(scenario_dir,"mu_selected.csv"))
  
  loglik_names<-grep("logLik_Y",colnames(samples_mat),value=TRUE); waic<-NA; LPML<-NA
  if(length(loglik_names)>0){
    lm<-samples_mat[,loglik_names,drop=FALSE]
    lppd  <-sum(apply(lm,2,function(x){mx<-max(x); mx+log(mean(exp(x-mx)))}))
    p_waic<-sum(apply(lm,2,var)); waic<- -2*(lppd-p_waic)
    LPML  <-sum(log(1/apply(lm,2,function(x) mean(exp(-x)))))
    write_csv(data.frame(WAIC=waic,LPML=LPML,lppd=lppd,pWAIC=p_waic),file.path(scenario_dir,"criteria.csv"))
  }
  
  params_struct<-c(beta_names,gamma_names); if(is_spatial) params_struct<-c(params_struct,"tau_s")
  ESS_struct <-effectiveSize(mcmc_list_full[,params_struct])
  Rhat_struct<-safe_gelman(mcmc_list_full[,params_struct])
  
  # ACF diagnostics
  acf_results <- do.call(rbind,lapply(params_struct,function(nm){
    ac   <- acf(samples_mat[,nm],lag.max=200,plot=FALSE)
    lags <- as.vector(ac$lag[-1]); acfs <- as.vector(ac$acf[-1])
    lag_01  <- lags[which(abs(acfs)<0.10)[1]]
    lag_005 <- lags[which(abs(acfs)<0.05)[1]]
    data.frame(Parameter=nm,
               ESS      =ESS_struct[nm],
               Rhat     =Rhat_struct[nm],
               lag_0.10 =if(is.na(lag_01))  Inf else lag_01,
               lag_0.05 =if(is.na(lag_005)) Inf else lag_005,
               acf_lag1 =acfs[1],
               acf_lag10=acfs[min(10,length(acfs))])
  }))
  write_csv(acf_results,file.path(scenario_dir,"acf_diagnostics.csv"))
  cat(sprintf("[%s] Maior lag(|acf|<0.10): %g | lag(|acf|<0.05): %g\n",
              model_type,
              max(acf_results$lag_0.10,na.rm=TRUE),
              max(acf_results$lag_0.05,na.rm=TRUE)))
  
  acf_df <- do.call(rbind,lapply(params_struct,function(nm){
    ac<-acf(samples_mat[,nm],lag.max=100,plot=FALSE)
    data.frame(Parameter=nm,Lag=as.vector(ac$lag[-1]),ACF=as.vector(ac$acf[-1]))
  }))
  ggsave(file.path(scenario_dir,"acf_params.png"),
         ggplot(acf_df,aes(x=Lag,y=ACF))+
           geom_col(width=0.6,fill="grey50")+
           geom_hline(yintercept=c(-0.10,0.10),linetype="dashed",color="blue",linewidth=0.5)+
           geom_hline(yintercept=c(-0.05,0.05),linetype="dotted",color="red",linewidth=0.5)+
           facet_wrap(~Parameter,scales="free_y")+theme_bw(base_size=11)+
           labs(title=paste("ACF dos parâmetros estruturais (",model_type,")"),
                subtitle="Azul tracejado: |0.10| | Vermelho pontilhado: |0.05|"),
         width=10,height=6)
  
  # Densidades
  ggsave(file.path(scenario_dir,"dens_beta.png"),
         ggplot(data.frame(Value=as.vector(samples_mat[,beta_names]),
                           Parameter=rep(beta_names,each=nrow(samples_mat))),aes(x=Value))+
           geom_density(fill="grey70")+facet_wrap(~Parameter,scales="free")+theme_bw()+
           labs(title=paste("Posterior densities - beta (",model_type,")")),width=8,height=6)
  ggsave(file.path(scenario_dir,"dens_gamma.png"),
         ggplot(data.frame(Value=as.vector(samples_mat[,gamma_names]),
                           Parameter=rep(gamma_names,each=nrow(samples_mat))),aes(x=Value))+
           geom_density(fill="grey70")+facet_wrap(~Parameter,scales="free")+theme_bw()+
           labs(title=paste("Posterior densities - gamma (",model_type,")")),width=8,height=6)
  
  # Painéis com rótulos de epsilon médio e cluster
  lambda_regs <- sort(unique(lambda_summary$Region))
  ggsave(file.path(scenario_dir,"painel_lambda.png"),
         ggplot(lambda_summary,aes(x=Time))+
           geom_ribbon(aes(ymin=Lower,ymax=Upper),fill="grey70",alpha=0.5)+
           geom_line(aes(y=Mean),color="black")+
           geom_line(aes(y=True),color="red",linetype="dashed")+
           facet_wrap(~Region,scales="free_y",ncol=3,
                      labeller=labeller(Region=make_region_label(lambda_regs)))+
           theme_bw(base_size=10)+
           labs(title=paste("Lambda estimado (",model_type,")"),
                subtitle="Vermelho tracejado: valor verdadeiro | ε̂ no título da célula",
                x="Tempo",y=expression(lambda[i*t])),
         width=12,height=10)
  
  theta_regs <- sort(unique(theta_summary$Region))
  ggsave(file.path(scenario_dir,"painel_theta.png"),
         ggplot(theta_summary,aes(x=Time))+
           geom_ribbon(aes(ymin=Lower,ymax=Upper),fill="steelblue",alpha=0.3)+
           geom_line(aes(y=Mean),color="steelblue")+
           geom_line(aes(y=True),color="red",linetype="dashed")+
           facet_wrap(~Region,scales="free_y",ncol=3,
                      labeller=labeller(Region=make_region_label(theta_regs)))+
           theme_bw(base_size=10)+
           labs(title=paste("Theta estimado (",model_type,")"),
                subtitle="Vermelho tracejado: valor verdadeiro",
                x="Tempo",y=expression(theta[i*t])),
         width=12,height=10)
  
  mu_regs <- sort(unique(mu_summary$Region))
  ggsave(file.path(scenario_dir,"painel_mu.png"),
         ggplot(mu_summary,aes(x=Time))+
           geom_ribbon(aes(ymin=Lower,ymax=Upper),fill="darkorange",alpha=0.3)+
           geom_line(aes(y=Mean),color="darkorange")+
           geom_line(aes(y=True),color="red",linetype="dashed")+
           facet_wrap(~Region,scales="free_y",ncol=3,
                      labeller=labeller(Region=make_region_label(mu_regs)))+
           theme_bw(base_size=10)+
           labs(title=paste("Mu estimado (",model_type,")"),
                subtitle=expression(paste(mu[it]==lambda[it]%.%E[it]%.%epsilon[i]%.%e^{x[it]^T*beta})),
                x="Tempo",y=expression(mu[i*t])),
         width=12,height=10)
  
  # Epsilon posterior com pontos dos valores verdadeiros e cores por cluster
  true_epsilon <- 1 - apply(constants$h * gamma_true, 1, sum)
  eps_show <- min(100,n_regions)
  eps_df <- epsilon_summary[seq_len(eps_show), ]
  eps_df$True <- true_epsilon[seq_len(eps_show)]
  
  ggsave(file.path(scenario_dir,"epsilon_posterior.png"),
         ggplot(eps_df, aes(x=Region, y=Eps_Mean, color=factor(Cluster))) +
           geom_point(size=0.8) +
           geom_errorbar(aes(ymin=Eps_Lower, ymax=Eps_Upper), width=0.3, linewidth=0.3) +
           geom_point(aes(y=True), color="red", size=1.2, shape=4, show.legend=FALSE) +
           scale_color_discrete(name="Cluster") +
           theme_bw() +
           labs(title=sprintf("Epsilon posterior (%d primeiras regiões)", eps_show),
                subtitle="Cruz vermelha: valor verdadeiro",
                y=expression(epsilon[i]), x="Região"),
         width=10, height=4)
  
  # Traceplots com média ergódica
  trace_p <- c(beta_names,gamma_names)
  if(is_spatial) trace_p <- c(trace_p,"tau_s")
  
  cores_cadeia <- c("Cadeia 1"="#2166AC","Cadeia 2"="#D6604D")
  
  df_trace <- do.call(rbind,lapply(seq_len(nchains),function(ch){
    cm <- as.matrix(mcmc_list_full[[ch]])
    do.call(rbind,lapply(trace_p,function(nm){
      vals     <- cm[,nm]
      erg_mean <- cumsum(vals)/seq_along(vals)
      data.frame(Iter=seq_along(vals),Value=vals,ErgMedia=erg_mean,
                 Parameter=nm,Cadeia=paste0("Cadeia ",ch),
                 stringsAsFactors=FALSE)
    }))
  }))
  
  ggsave(file.path(scenario_dir,"traceplots.png"),
         ggplot(df_trace,aes(x=Iter,color=Cadeia,fill=Cadeia))+
           geom_line(aes(y=Value),alpha=0.25,linewidth=0.20)+
           geom_line(aes(y=ErgMedia),alpha=0.90,linewidth=0.75,linetype="solid")+
           scale_color_manual(values=cores_cadeia)+
           scale_fill_manual(values=cores_cadeia)+
           facet_wrap(~Parameter,scales="free_y")+
           theme_bw(base_size=11)+theme(legend.position="bottom")+
           labs(title=paste("Traceplots + Média Ergódica (",model_type,")"),
                subtitle="Linha grossa = média ergódica | Linha fina = cadeia",
                x="Iteração (pós-burnin)",y="Valor",color="Cadeia",fill="Cadeia"),
         width=10,height=max(6,3*ceiling(length(trace_p)/3)))
  
  data.frame(model=model_type,niter=niter,nburnin=nburnin,thin=thin,
             WAIC=waic,LPML=LPML,corr_s=corr_s,
             MSE_beta=mean(beta_metrics$MSE),MSE_gamma=mean(gamma_metrics$MSE),
             Coverage_beta=mean(beta_metrics$Coverage,na.rm=TRUE),
             Coverage_gamma=mean(gamma_metrics$Coverage,na.rm=TRUE),
             ESS_beta_min=min(beta_metrics$ESS,na.rm=TRUE),
             ESS_gamma_min=min(gamma_metrics$ESS,na.rm=TRUE),
             ESS_tau=ESS_tau,ESS_s_mean=ESS_s_mean,ESS_lambda_mean=ESS_lambda_mean,
             ESS_global_min=min(ESS_struct,na.rm=TRUE),Rhat_max=max(Rhat_struct,na.rm=TRUE),
             lag_max_0.10=max(acf_results$lag_0.10,na.rm=TRUE),
             lag_max_0.05=max(acf_results$lag_0.05,na.rm=TRUE))
}

# ------ 5. Paralelo ------
model_types <- c("spatial","non_spatial")
n_cores     <- min(length(model_types),parallel::detectCores()-1); if(n_cores<1) n_cores<-1

output_dir <- "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Simulacao/main/resultados_spatial_vs_nonspatial_w_unif_T_100_A_510"
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

cl <- makeCluster(n_cores)
clusterExport(cl,c("constants_spatial","constants_nonspatial","data_nimble",
                   "inits_list_spatial","inits_list_nonspatial",
                   "n_regions","n_times","p","K","beta_true","gamma_true","lambda_true",
                   "code_spatial","code_nonspatial","run_model","output_dir",
                   "REGIONS_INTEREST","LAMBDA_MONITORS"))
clusterEvalQ(cl,{
  library(nimble); library(coda); library(dplyr); library(ggplot2); library(readr)
  setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Simulacao/main/main2")
  Sys.setenv(OMP_NUM_THREADS="1"); Sys.setenv(MKL_NUM_THREADS="1")
  if(requireNamespace("RhpcBLASctl",quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(1)
})
resultados <- parLapply(cl,model_types,function(m) run_model(m,output_dir))
stopCluster(cl)

# ------ 6. Consolidação ------
library(dplyr); library(readr); library(ggplot2); library(stringr)
resumo <- bind_rows(resultados)
write_csv(resumo,file.path(output_dir,"resumo_comparativo.csv"))
print(resumo)

# Função para gráficos comparativos com rótulos incluindo cluster e ε dos dois modelos
build_compare_plot <- function(tipo, titulo_expr) {
  all_df <- lapply(model_types, function(m) {
    path <- file.path(output_dir, m, paste0(tipo, "_selected.csv"))
    if (!file.exists(path)) {
      warning(paste("Arquivo não encontrado:", path))
      return(NULL)
    }
    read_csv(path, show_col_types = FALSE)
  }) %>% bind_rows()
  
  if (is.null(all_df) || nrow(all_df) == 0) {
    stop("Nenhum dado carregado para o tipo ", tipo)
  }
  
  eps_list <- lapply(model_types, function(m) {
    path <- file.path(output_dir, m, "epsilon_summary.csv")
    if (!file.exists(path)) {
      warning(paste("Arquivo epsilon_summary não encontrado:", path))
      return(NULL)
    }
    eps <- read_csv(path, show_col_types = FALSE)
    eps$model <- m
    eps
  }) %>% bind_rows()
  
  if (is.null(eps_list) || nrow(eps_list) == 0) {
    warning("Nenhum epsilon_summary carregado. Os rótulos não conterão ε̂.")
    region_labels <- setNames(as.character(unique(all_df$Region)), unique(all_df$Region))
  } else {
    regions <- unique(all_df$Region)
    region_info <- eps_list %>%
      filter(Region %in% regions) %>%
      select(Region, Cluster, Eps_Mean, model) %>%
      distinct() %>%
      pivot_wider(names_from = model, values_from = Eps_Mean, names_prefix = "Eps_")
    
    region_labels <- setNames(lapply(regions, function(r) {
      info <- region_info %>% filter(Region == r)
      if(nrow(info) == 0) {
        return(sprintf("Região %d\nε_sp = NA\nε_ns = NA", r))
      }
      cl <- info$Cluster[1]
      eps_sp <- if("Eps_spatial" %in% names(info)) info$Eps_spatial else NA
      eps_ns <- if("Eps_non_spatial" %in% names(info)) info$Eps_non_spatial else NA
      sprintf("Região %d (C%d)\nε_sp = %.3f\nε_ns = %.3f", r, cl, eps_sp, eps_ns)
    }), as.character(regions))
    cat("\n--- Rótulos das facetas para", tipo, "---\n")
    print(head(region_labels))
  }
  
  ggplot(all_df, aes(x = Time, y = Mean, color = model, fill = model, group = model)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_line(data = all_df %>% distinct(Region, Time, True),
              aes(x = Time, y = True, group = Region), inherit.aes = FALSE,
              color = "black", linetype = "dashed", linewidth = 0.7) +
    facet_wrap(~ Region, scales = "free_y", ncol = 3,
               labeller = labeller(Region = function(x) region_labels[x])) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom") +
    labs(title = titulo_expr,
         subtitle = "Linha preta tracejada: valor verdadeiro",
         x = "Tempo", color = "Modelo", fill = "Modelo")
}

ggsave(file.path(output_dir,"lambda_comparativo.png"),
       build_compare_plot("lambda", expression(paste("Comparação de ",lambda[i,t],": spatial vs non_spatial"))),
       width = 14, height = 10, dpi = 300)
ggsave(file.path(output_dir,"theta_comparativo.png"),
       build_compare_plot("theta", expression(paste("Comparação de ",theta[i,t],": spatial vs non_spatial"))),
       width = 14, height = 10, dpi = 300)
ggsave(file.path(output_dir,"mu_comparativo.png"),
       build_compare_plot("mu", expression(paste("Comparação de ",mu[i,t],": spatial vs non_spatial"))),
       width = 14, height = 10, dpi = 300)

acf_comparativo <- lapply(model_types, function(m)
  read_csv(file.path(output_dir, m, "acf_diagnostics.csv"), show_col_types = FALSE) %>%
    mutate(model = m)) %>% bind_rows()
write_csv(acf_comparativo, file.path(output_dir, "acf_comparativo.csv"))
cat("\n--- ACF Diagnóstico por modelo ---\n"); print(acf_comparativo)

cat("\nTempo total:\n"); print(Sys.time() - inicio_global)