library(atmle)
library(purrr)
library(xgboost)
library(tmle)
library(foreach)
library(doParallel)
library(parallel)

source("~/atmle-fusionRCT/sim_data.R")
source("~/atmle-fusionRCT/summary_stats.R")
B <- 100
n <- 5000
p_S <- 0.2
psi0 <- get_psi0(p_S)
psi0S <- -10
psi0S_interaction <- -15

RNGkind(sample.kind = "Rejection")
set.seed(1517)
seeds <- sample(1:2^15, B)
save(seeds, file = "~/atmle-fusionRCT/output/seeds.Rdata")

registerDoParallel(cores = detectCores()-2)
getDoParWorkers()
result_list <- list()
for(i in 1:B){
  print(paste0("Starting simulation ", i))
  
  set.seed(seeds[i])
  d <- sim_data(p_S, n, AS_interaction = TRUE)
  
  ############### subgroup-specific data fusion effect estimates ###############
  start <- proc.time()
  set.seed(seeds[i])
  atmle_res <- atmle(data = d,
                     S = 1,
                     W = c(2, 3),
                     A = 4,
                     Y = 5,
                     controls_only = FALSE,
                     family = "gaussian",
                     theta_method = "sl3",
                     Pi_method = "sl3",
                     g_method = "sl3",
                     Q_method = "sl3",
                     theta_tilde_method = "sl3",
                     bias_working_model = "HAL",
                     pooled_working_model = "HAL",
                     enumerate_basis_args = list(num_knots = 10),
                     verbose = FALSE)
  res_atmle <- data.frame(
    psi = atmle_res$est, 
    se = (atmle_res$upper-atmle_res$est)/1.96
  )
  print(paste0("atmle_res for simulation ", i, " time: ", as.numeric((proc.time()-start)["elapsed"])/60))
  
  start <- proc.time()
  set.seed(seeds[i])
  atmle_res_glmnet <- atmle(data = d,
                            S = 1,
                            W = c(2, 3),
                            A = 4,
                            Y = 5,
                            controls_only = FALSE,
                            family = "gaussian",
                            verbose = FALSE)
  res_atmle_glmnet <- data.frame(
    psi = atmle_res_glmnet$est, 
    se = (atmle_res_glmnet$upper-atmle_res_glmnet$est)/1.96
  )
  print(paste0("atmle_res_glmnet for simulation ", i, " time: ", as.numeric((proc.time()-start)["elapsed"])/60))
  
  start <- proc.time()
  set.seed(seeds[i])
  atmle_res_hal <- atmle(data = d,
                         S = 1,
                         W = c(2, 3),
                         A = 4,
                         Y = 5,
                         controls_only = FALSE,
                         family = "gaussian",
                         theta_method = "glm",
                         Pi_method = "HAL",
                         g_method = "glm",
                         Q_method = "HAL",
                         theta_tilde_method = "glm",
                         bias_working_model = "HAL",
                         pooled_working_model = "HAL",
                         enumerate_basis_args = list(num_knots = 10),
                         verbose = FALSE)
  res_atmle_hal <- data.frame(
    psi = atmle_res_hal$est, 
    se = (atmle_res_hal$upper-atmle_res_hal$est)/1.96
  )
  print(paste0("atmle_res_hal for simulation ", i, " time: ", as.numeric((proc.time()-start)["elapsed"])/60))
  
  #################### subgroup-specific effect estimates ######################
  d_S <- d[d$S == 1,]
  
  # TMLE ATE estimator w GLM estimation for Q
  start <- proc.time()
  set.seed(seeds[i])
  tmle_res_glmQ <- tmle(
    Y = d_S$Y, A = d_S$A, W = d_S[, c("W1", "W2")], 
    g1W = rep(0.5, nrow(d_S)), Q.SL.library = c("SL.glm"),
    family = "gaussian"
  )
  res_tmle_glmQ <- data.frame(psi=tmle_res_glmQ$estimates$ATE$psi, 
                              se=sqrt(tmle_res_glmQ$estimates$ATE$var.psi))
  print(paste0("tmle_res_glmQ for sim ", i, " time: ", as.numeric((proc.time()-start)["elapsed"])/60))
  
  # TMLE ATE estimator w SL estimation for Q
  # start <- proc.time()
  # set.seed(seeds[i])
  # tmle_res_slQ <- tmle(
  #   Y = d_S$Y, A = d_S$A, W = d_S[, c("W1", "W2")], 
  #   g1W = tmle_res_glmQ$g$g1W,
  #   family = "gaussian"
  # )
  # res_tmle_slQ <- data.frame(psi=tmle_res_slQ$estimates$ATE$psi, 
  #                            se=sqrt(tmle_res_slQ$estimates$ATE$var.psi))
  # print(paste0("tmle_res_slQ for sim ", i, " time: ", as.numeric((proc.time()-start)["elapsed"])/60))
  # 
  # TMLE ATE estimator w GLM estimation for Q + GLM estimation for g
  start <- proc.time()
  set.seed(seeds[i])
  tmle_res_glmQ_glmG <- tmle(
    Y = d_S$Y, A = d_S$A, W = d_S[, c("W1", "W2")], 
    g.SL.library = c("SL.glm"), Q = tmle_res_glmQ$Qinit$Q,
    family = "gaussian"
  )
  res_tmle_glmQ_glmG <- data.frame(psi=tmle_res_glmQ_glmG$estimates$ATE$psi, 
                                   se=sqrt(tmle_res_glmQ_glmG$estimates$ATE$var.psi))
  print(paste0("tmle_res_glmQ_glmG for sim ", i, " time: ", as.numeric((proc.time()-start)["elapsed"])/60))
  
  # TMLE ATE estimator w SL estimation for Q + GLM estimation for g
  # start <- proc.time()
  # set.seed(seeds[i])
  # tmle_res_slQ_glmG <- tmle(
  #   Y = d_S$Y, A = d_S$A, W = d_S[, c("W1", "W2")], 
  #   g1W = tmle_res_glmQ_glmG$g$g1W, Q = tmle_res_slQ$Qinit$Q,
  #   family = "gaussian"
  # )
  # res_tmle_slQ_glmG <- data.frame(psi=tmle_res_slQ_glmG$estimates$ATE$psi, 
  #                                 se=sqrt(tmle_res_slQ_glmG$estimates$ATE$var.psi))
  # print(paste0("tmle_res_glmQ_glmG for sim ", i, " time: ", as.numeric((proc.time()-start)["elapsed"])/60))
  # 
  # unadjusted ATE estimator
  mod_unadj <- lm(d_S$Y ~ d_S$A)
  res_unadj <- data.frame(psi=summary(mod_unadj)$coefficients[2,1], 
                          se=summary(mod_unadj)$coefficients[2,2])
  
  res_list <- list(
    res_atmle = res_atmle, res_atmle_glmnet = res_atmle_glmnet, 
    res_atmle_hal = res_atmle_hal,
    res_tmle_glmQ = res_tmle_glmQ, 
    # res_tmle_slQ = res_tmle_slQ, 
    res_tmle_glmQ_glmG = res_tmle_glmQ_glmG, 
    # res_tmle_slQ_glmG = res_tmle_slQ_glmG, 
    res_unadj = res_unadj
  )
  start <- proc.time()
  save(res_list, file = paste0("~/atmle-fusionRCT/output/ASinter_res",i,".Rdata"))
  print(paste0("saving sim ", i, " time: ", as.numeric((proc.time()-start)["elapsed"])/60))
  
  result_list[[i]] <- res_list
}

result_list <- list()
for (i in 1:B) {
  print(i)
  load(paste0("~/atmle-fusionRCT/output/ASinter_res",i,".Rdata"))
  result_list[[i]] <- res_list
}

atmle_res <- do.call(rbind, lapply(result_list, "[[", "res_atmle"))
atmle_res_glmnet <- do.call(rbind, lapply(result_list, "[[", "res_atmle_glmnet"))
atmle_res_hal <- do.call(rbind, lapply(result_list, "[[", "res_atmle_hal"))
tmle_glmQ_res <- do.call(rbind, lapply(result_list, "[[", "res_tmle_glmQ"))
#tmle_slQ_res <- do.call(rbind, lapply(result_list, "[[", "res_tmle_slQ"))
tmle_glmQ_glmG_res <- do.call(rbind, lapply(result_list, "[[", "res_tmle_glmQ_glmG"))
#tmle_slQ_glmG_res <- do.call(rbind, lapply(result_list, "[[", "res_tmle_slQ_glmG"))
unadj_res <- do.call(rbind, lapply(result_list, "[[", "res_unadj"))

atmle_stats <- summary_stats(estimator_name = "ATMLE", psi = atmle_res$psi, 
                             se = atmle_res$se, psi_0 = psi0S_interaction)
atmle_glmnet_stats <- summary_stats(estimator_name = "ATMLE_glmnet", psi = atmle_res_glmnet$psi, 
                                    se = atmle_res_glmnet$se, psi_0 = psi0S_interaction)
atmle_hal_stats <- summary_stats(estimator_name = "ATMLE_HAL", psi = atmle_res_hal$psi, 
                                 se = atmle_res_hal$se, psi_0 = psi0S_interaction)
tmle_glmQ_stats <- summary_stats(estimator_name = "TMLE_glmQ", psi = tmle_glmQ_res$psi, 
                                 se = tmle_glmQ_res$se, psi_0 = psi0S_interaction)
# tmle_slQ_stats <- summary_stats(estimator_name = "TMLE_slQ", psi = tmle_slQ_res$psi, 
#                                 se = tmle_slQ_res$se, psi_0 = psi0S_interaction)
tmle_glmQ_glmG_stats <- summary_stats(estimator_name = "TMLE_glmQ_glmG", psi = tmle_glmQ_glmG_res$psi, 
                                      se = tmle_glmQ_glmG_res$se, psi_0 = psi0S_interaction)
# tmle_slQ_glmG_stats <- summary_stats(estimator_name = "TMLE_slQ_glmG", psi = tmle_slQ_glmG_res$psi, 
#                                      se = tmle_slQ_glmG_res$se, psi_0 = psi0S_interaction)
unadj_stats <- summary_stats(estimator_name = "Unadjusted", psi = unadj_res$psi, 
                             se = unadj_res$se, psi_0 = psi0S_interaction)

res <- data.frame(
  rbind(atmle_stats, atmle_glmnet_stats, atmle_hal_stats, 
        tmle_glmQ_stats, 
        # tmle_slQ_stats, 
        tmle_glmQ_glmG_stats, 
        #tmle_slQ_glmG_stats, 
        unadj_stats)
)
write.csv(res, file = "~/atmle-fusionRCT/results_ASinteraction.csv", row.names = F)
