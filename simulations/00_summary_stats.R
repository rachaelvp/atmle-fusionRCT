# calculates summary statistics using truth and vector of estimates
summary_stats <- function(estimator_name, psi, se = NULL, psi_0 = 1.685959,
                          var = NULL){
  Bias <- mean(psi - psi_0)
  SE_bootstrap <- sqrt(var(psi))
  Ratio_BiasSE_bootstrap <- Bias/SE_bootstrap
  bootstrap_CIlower <- psi - 1.96*SE_bootstrap
  bootstrap_CIupper <- psi + 1.96*SE_bootstrap
  Coverage_bootstrapSE <- mean(bootstrap_CIlower <= psi_0 & 
                                 bootstrap_CIupper >= psi_0)
  MSE <- mean((psi - psi_0)^2)
  
  if(is.null(var) & is.null(se)) {
    Coverage_analyticSE <- NA
    Ratio_BiasSE_analytic <- NA
    avg_analyticSE <- NA
  } else if(!is.null(se)) {
    Ratio_BiasSE_analytic <- mean((psi - psi_0)/se)
    ci_lower <- psi - 1.96*se
    ci_upper <- psi + 1.96*se
    Coverage_analyticSE <- mean(ci_lower <= psi_0 & ci_upper >= psi_0)
    avg_analyticSE <- mean(se)
  } else if(!is.null(var)) {
    Ratio_BiasSE_analytic <- mean((psi[i] - psi_0)/sqrt(var))
    ci_lower <- psi - 1.96*sqrt(var)
    ci_upper <- psi + 1.96*sqrt(var)
    Coverage_analyticSE <- mean(ci_lower <= psi_0 & ci_upper >= psi_0)
    avg_analyticSE <- mean(sqrt(var))
  }
  
  return_vec <- data.frame(Estimator = estimator_name, Bias = Bias, MSE = MSE, 
                           SE_bootstrap = SE_bootstrap, 
                           avg_analyticSE = avg_analyticSE, 
                           Ratio_BiasSE_bootstrap = Ratio_BiasSE_bootstrap,
                           Ratio_BiasSE_analytic = Ratio_BiasSE_analytic,
                           Coverage_bootstrapSE = Coverage_bootstrapSE, 
                           Coverage_analyticSE = Coverage_analyticSE)
  return(return_vec)
}