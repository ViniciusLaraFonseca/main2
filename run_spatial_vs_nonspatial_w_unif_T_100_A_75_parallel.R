# ==============================================================================
# run_spatial_vs_nonspatial_w_unif_T_100_A_75_parallel.R
# ==============================================================================
rm(list = ls())
inicio_global <- Sys.time()
setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Simulacao/main/main2")

pkgs <- c("nimble","coda","parallel","dplyr","ggplot2","tidyr","readr","stringr")
for(pkg in pkgs) { if(!require(pkg,character.only=TRUE)) install.packages(pkg); library(pkg,character.only=TRUE) }

Sys.setenv(OMP_NUM_THREADS="1"); Sys.setenv(MKL_NUM_THREADS="1")
if(requireNamespace("RhpcBLASctl",quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(1)

# ------ 1. Dados ------
source("Geracao_w_unif_T_100_A_75.R")

n_regions <- constants_nimble$n_regions; n_times <- constants_nimble$n_times
p <- constants_nimble$p;                  K       <- constants_nimble$K
beta_true   <- get("beta_true",   envir=.GlobalEnv)
gamma_true  <- get("gamma_true",  envir=.GlobalEnv)
lambda_true <- get("lambda_true", envir=.GlobalEnv)
cat("Dados carregados. n_regions =",n_regions,", n_times =",n_times,"\n")

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
  list(beta=beta_true,  gamma=gamma_true,  lambda=lambda_true,
       sigma_s=0.5, s=rep(0,n_regions)),
  list(beta=rnorm(p,0,0.5), gamma=gamma_true*0.9,
       lambda=matrix(rgamma(n_regions*n_times,1,1),nrow=n_regions),
       sigma_s=1, s=rep(0,n_regions))
)
inits_list_nonspatial <- list(
  list(beta=beta_true,  gamma=gamma_true,  lambda=lambda_true),
  list(beta=rnorm(p,0,0.5), gamma=gamma_true*0.9,
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

  monitors <- c("beta","gamma","lambda","logLik_Y")
  if(is_spatial) monitors <- c(monitors,"s","sigma_s","tau_s")
  conf$addMonitors(monitors)
  conf$printSamplers()

  Rmcmc <- buildMCMC(conf); Cmcmc <- compileNimble(Rmcmc,project=model)
  niter<-20000; nburnin<-5000; nchains<-2

  samples <- runMCMC(Cmcmc,niter=niter,nburnin=nburnin,nchains=nchains,
                     inits=inits_list,samplesAsCodaMCMC=TRUE,summary=FALSE,WAIC=FALSE)
  saveRDS(samples,file.path(scenario_dir,"samples.rds"))

  samples_mat    <- as.matrix(samples)
  mcmc_list_full <- mcmc.list(lapply(1:nchains,function(ch) as.mcmc(samples[[ch]])))

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
    ggsave(file.path(scenario_dir,"s_posterior.png"),
           ggplot(data.frame(Region=1:n_regions,Mean=s_metrics$Mean,
                             Lower=s_metrics$HPD_Lower,Upper=s_metrics$HPD_Upper),
                  aes(x=Region,y=Mean))+
             geom_point()+geom_errorbar(aes(ymin=Lower,ymax=Upper),width=0.3)+
             geom_hline(yintercept=0,linetype="dashed")+theme_bw()+
             labs(title="Efeito espacial s: média posterior e HPD 95%",y="s[i]",x="Região"),
           width=8,height=5)
  }

  lambda_names    <- grep("^lambda\\[",colnames(samples_mat),value=TRUE)
  ESS_lambda_mean <- mean(effectiveSize(mcmc_list_full[,lambda_names]),na.rm=TRUE)
  regions_interest <- c(1,8,15,19,22,31,34,40,46,55,65,75)

  lambda_summary <- data.frame()
  for(nm in lambda_names){
    idx<-str_match(nm,"lambda\\[(\\d+),\\s*(\\d+)\\]"); i<-as.numeric(idx[2]); t<-as.numeric(idx[3])
    if(i %in% regions_interest){
      hpd<-HPDinterval(as.mcmc(samples_mat[,nm]))
      lambda_summary<-rbind(lambda_summary,data.frame(Region=i,Time=t,True=lambda_true[i,t],
                                                       Mean=mean(samples_mat[,nm]),Lower=hpd[1],Upper=hpd[2],model=model_type))
    }
  }
  write_csv(lambda_summary,file.path(scenario_dir,"lambda_selected.csv"))

  beta_cols   <- grep("^beta\\[",  colnames(samples_mat),value=TRUE)
  lambda_cols <- grep("^lambda\\[",colnames(samples_mat),value=TRUE)
  theta_summary <- data.frame()
  for(nm in lambda_cols){
    idx<-str_match(nm,"lambda\\[(\\d+),\\s*(\\d+)\\]"); i<-as.numeric(idx[2]); t<-as.numeric(idx[3])
    if(i %in% regions_interest && t<=n_times){
      ldraws<-samples_mat[,nm]; bdraws<-samples_mat[,beta_cols,drop=FALSE]
      x_it<-data_nimble$x[i,t,]; lin<-as.vector(bdraws %*% x_it); theta<-ldraws*exp(lin)
      hpd<-HPDinterval(as.mcmc(theta)); tv<-lambda_true[i,t]*exp(sum(x_it*beta_true))
      theta_summary<-rbind(theta_summary,data.frame(Region=i,Time=t,True=tv,
                                                      Mean=mean(theta),Lower=hpd[1],Upper=hpd[2],model=model_type))
    }
  }
  write_csv(theta_summary,file.path(scenario_dir,"theta_selected.csv"))

  loglik_names<-grep("logLik_Y",colnames(samples_mat),value=TRUE); waic<-NA; LPML<-NA
  if(length(loglik_names)>0){
    lm<-samples_mat[,loglik_names,drop=FALSE]
    lppd  <-sum(apply(lm,2,function(x){mx<-max(x); mx+log(mean(exp(x-mx)))}))
    p_waic<-sum(apply(lm,2,var)); waic<- -2*(lppd-p_waic)
    LPML  <-sum(log(1/apply(lm,2,function(x) mean(exp(-x)))))
    write_csv(data.frame(WAIC=waic,LPML=LPML,lppd=lppd,pWAIC=p_waic),file.path(scenario_dir,"criteria.csv"))
  }

  params_struct<-c(beta_names,gamma_names); if(is_spatial) params_struct<-c(params_struct,"tau_s")
  ESS_struct<-effectiveSize(mcmc_list_full[,params_struct])
  Rhat_struct<-safe_gelman(mcmc_list_full[,params_struct])

  # Gráficos
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
  ggsave(file.path(scenario_dir,"painel_lambda.png"),
         ggplot(lambda_summary,aes(x=Time))+
           geom_ribbon(aes(ymin=Lower,ymax=Upper),fill="grey70",alpha=0.5)+
           geom_line(aes(y=Mean),color="black")+geom_line(aes(y=True),color="red",linetype="dashed")+
           facet_wrap(~Region,scales="free_y",ncol=3)+theme_bw()+
           labs(title=paste("Lambda estimado (",model_type,")")),width=12,height=10)
  trace_p<-c(beta_names,gamma_names); if(is_spatial) trace_p<-c(trace_p,"tau_s")
  ggsave(file.path(scenario_dir,"traceplots.png"),
         ggplot(data.frame(Iter=rep(1:nrow(samples_mat),times=length(trace_p)),
                           Value=as.vector(samples_mat[,trace_p]),
                           Parameter=rep(trace_p,each=nrow(samples_mat))),aes(x=Iter,y=Value))+
           geom_line(alpha=0.3)+facet_wrap(~Parameter,scales="free_y")+theme_bw()+
           labs(title=paste("Traceplots (",model_type,")")),width=10,height=6)

  data.frame(model=model_type,WAIC=waic,LPML=LPML,corr_s=corr_s,
             MSE_beta=mean(beta_metrics$MSE),MSE_gamma=mean(gamma_metrics$MSE),
             Coverage_beta=mean(beta_metrics$Coverage,na.rm=TRUE),
             Coverage_gamma=mean(gamma_metrics$Coverage,na.rm=TRUE),
             ESS_beta_min=min(beta_metrics$ESS,na.rm=TRUE),
             ESS_gamma_min=min(gamma_metrics$ESS,na.rm=TRUE),
             ESS_tau=ESS_tau,ESS_s_mean=ESS_s_mean,ESS_lambda_mean=ESS_lambda_mean,
             ESS_global_min=min(ESS_struct,na.rm=TRUE),Rhat_max=max(Rhat_struct,na.rm=TRUE))
}

# ------ 5. Paralelo ------
model_types <- c("spatial","non_spatial")
n_cores     <- min(length(model_types),parallel::detectCores()-1); if(n_cores<1) n_cores<-1

output_dir <- "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Simulacao/main/resultados_spatial_vs_nonspatial_w_unif_T_100_A_75"
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

cl <- makeCluster(n_cores)
clusterExport(cl,c("constants_spatial","constants_nonspatial","data_nimble",
                   "inits_list_spatial","inits_list_nonspatial",
                   "n_regions","n_times","p","K","beta_true","gamma_true","lambda_true",
                   "code_spatial","code_nonspatial","run_model","output_dir"))
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

build_compare_plot <- function(tipo, titulo_expr) {
  all_df <- lapply(model_types,function(m)
    read_csv(file.path(output_dir,m,paste0(tipo,"_selected.csv")),show_col_types=FALSE)) %>% bind_rows()
  ggplot(all_df,aes(x=Time,y=Mean,color=model,fill=model,group=model))+
    geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.15,color=NA)+geom_line(linewidth=0.8)+
    geom_line(data=all_df %>% distinct(Region,Time,True),
              aes(x=Time,y=True,group=Region),inherit.aes=FALSE,
              color="black",linetype="dashed",linewidth=0.7)+
    facet_wrap(~Region,scales="free_y",ncol=3)+theme_bw(base_size=12)+
    theme(legend.position="bottom")+
    labs(title=titulo_expr,subtitle="Linha preta tracejada: valor verdadeiro",
         x="Tempo",color="Modelo",fill="Modelo")
}
ggsave(file.path(output_dir,"lambda_comparativo.png"),
       build_compare_plot("lambda",expression(paste("Comparação de ",lambda[i,t],": spatial vs non_spatial"))),
       width=14,height=10,dpi=300)
ggsave(file.path(output_dir,"theta_comparativo.png"),
       build_compare_plot("theta",expression(paste("Comparação de ",theta[i,t],": spatial vs non_spatial"))),
       width=14,height=10,dpi=300)

cat("\nTempo total:\n"); print(Sys.time()-inicio_global)
