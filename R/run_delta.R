##' Run the model with the specified configurations
##' @param data$tmb_data data input for the TMB model run
##' @param data$tmb_params parameter input for the TMB model run
##' @param conf  configurations file for the model and data
##' @details the model is run in two steps to make model parameters estimable
##' @return different outputs needed for summarizing and plotting the results at later stage
##' @useDynLib mackerel_ST
##' @export fit_mackerelST
##' 


fit_mackerelST_delta <- function(data, conf, conf_extra=NULL) {

  version <- paste0(getwd(), "/src/spatial_SDM_RE_delta"); dllversion = "spatial_SDM_RE_delta"
  
  #dyn.unload(dynlib(version))
  # if (.Platform$OS.type == "unix") compile(paste0(version, ".cpp"),"-O0 -g")
  # if (.Platform$OS.type == "windows") compile(paste0(version, ".cpp"), "-O1 -g", DLLFLAGS="")
  dyn.load(dynlib(version))


# for the delta model, starting to estimate only the positive values

Map_phase1 <- list()
Map_phase1$beta_absc <- factor(rep(NA, length(data$tmb_params$beta_absc)))

Map_phase1$s50_absc <- factor(NA)
Map_phase1$logslope_absc <- factor(NA)
Map_phase1$s50 <- factor(NA)
Map_phase1$logslope <- factor(NA)

Map_phase1$transf_rho_age_absc <- factor(rep(NA, length(data$tmb_params$transf_rho_age_absc)))
Map_phase1$transf_rho_age <- factor(rep(NA, length(data$tmb_params$transf_rho_age)))
Map_phase1$transf_rho_absc <- factor(rep(NA, length(data$tmb_params$transf_rho_absc)))
Map_phase1$logKappa <- factor(rep(NA, length(data$tmb_params$logKappa)))
Map_phase1$logKappa_absc <- factor(rep(NA, length(data$tmb_params$logKappa_absc)))
Map_phase1$omega_absc <- factor(rep(NA, length(data$tmb_params$omega_absc)))
Map_phase1$omega <- factor(rep(NA, length(data$tmb_params$omega)))
Map_phase1$logTauO_absc <- factor(rep(NA, length(data$tmb_params$logTauO_absc)))
Map_phase1$logTauO <- factor(rep(NA, length(data$tmb_params$logTauO)))
Map_phase1$epsilon_st_absc <- factor(rep(NA, length(data$tmb_params$epsilon_st_absc)))
Map_phase1$logTauE_absc <- factor(rep(NA, length(data$tmb_params$logTauE_absc)))
Map_phase1$logsds <- factor(rep(NA, length(data$tmb_params$logsds)))
Map_phase1$logTauE <- factor(rep(NA, length(data$tmb_params$logTauE)))

Map_phase1$ln_tau_G_absc <- factor(rep(NA, length(data$tmb_params$ln_tau_G_absc)))# not estimating the random vessel effect
Map_phase1$RE_absc <- factor(rep(NA,length(data$tmb_params$RE_absc)))# not estimating the random vessel effect

if (conf$family != 0) Map_phase1$thetaf = factor(rep(NA, length(data$tmb_params$thetaf)))

# better use a stepwise optimization 
Map_phase1$epsilon_st <- factor(rep(NA, length(data$tmb_params$epsilon_st)))
Map_phase1$logKappa <- factor(rep(NA, length(data$tmb_params$logKappa)))


tmb_obj_phase1 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase1,
                                 random = c("RE"), DLL = dllversion, silent = TRUE)

startTime = Sys.time()
opt_phase1 <- fit_tmb(tmb_obj_phase1, lower=-15, upper=15, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))          
endTime = Sys.time()
endTime-startTime

# 
# rm(tmb_obj_phase2)
# # for the delta model, add the zero inflation part
# 
Map_phase2 <- list()
Map_phase2$beta_absc <- factor(rep(NA, length(data$tmb_params$beta_absc)))

Map_phase2$s50_absc <- factor(NA)
Map_phase2$logslope_absc <- factor(NA)
Map_phase2$s50 <- factor(NA)
Map_phase2$logslope <- factor(NA)

Map_phase2$transf_rho_age_absc <- factor(rep(NA, length(data$tmb_params$transf_rho_age_absc)))
Map_phase2$transf_rho_age <- factor(rep(NA, length(data$tmb_params$transf_rho_age)))

Map_phase2$transf_rho_absc <- factor(rep(NA, length(data$tmb_params$transf_rho_absc)))
if(conf$ARorIID[2] != 0) Map_phase2$transf_rho <- factor(rep(NA, length(data$tmb_params$transf_rho)))
if(conf$ARorIID[2] == 1) Map_phase2$transf_rho <- factor(rep(300000, length(data$tmb_params$transf_rho)))

Map_phase2$logKappa_absc <- factor(rep(NA, length(data$tmb_params$logKappa_absc)))
Map_phase2$logKappa <- factor(rep(123456, length(data$tmb_params$logKappa)))
Map_phase2$omega_absc <- factor(rep(NA, length(data$tmb_params$omega_absc)))
Map_phase2$omega <- factor(rep(NA, length(data$tmb_params$omega)))
Map_phase2$logTauO_absc <- factor(rep(NA, length(data$tmb_params$logTauO_absc)))
Map_phase2$logTauO <- factor(rep(NA, length(data$tmb_params$logTauO)))
Map_phase2$epsilon_st_absc <- factor(rep(NA, length(data$tmb_params$epsilon_st_absc)))
Map_phase2$logTauE_absc <- factor(rep(NA, length(data$tmb_params$logTauE_absc)))
if(conf$ARorIID[2] != 2) Map_phase2$logsds <- factor(rep(NA, length(data$tmb_params$logsds)))

Map_phase2$ln_tau_G_absc <- factor(rep(NA, length(data$tmb_params$ln_tau_G_absc)))# not estimating the random vessel effect
Map_phase2$RE_absc <- factor(rep(NA,length(data$tmb_params$RE_absc)))# not estimating the random vessel effect

if (conf$family != 0) Map_phase2$thetaf = factor(rep(NA, length(data$tmb_params$thetaf)))

old_par <- set_par_value(opt_phase1$par, "beta")
Map_phase2$beta <- matrix(1:(ncol(data$tmb_data$X)*data$tmb_data$Nage), nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)
Map_phase2$beta <- factor(Map_phase2$beta)
map_beta <- Map_phase2$beta
map_beta <- factor(map_beta, labels = seq(1, sum(!is.na(map_beta))))
data$tmb_params$beta <- matrix(old_par[map_beta], nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)
data$tmb_params$ln_phi <- set_par_value(opt_phase1$par, "ln_phi")
data$tmb_params$ln_tau_G <- matrix(set_par_value(opt_phase1$par, "ln_tau_G"), ncol= data$tmb_data$Nage, byrow=F)
data$tmb_params$RE <- matrix(set_par_value(tmb_obj_phase1$env$last.par.best, "RE"), nrow=data$tmb_data$nobs_RE, ncol= data$tmb_data$Nage, byrow=F)

tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
                                 random = c("RE", "epsilon_st"), DLL = dllversion, silent = TRUE)
startTime = Sys.time()
opt_phase2 <- fit_tmb(tmb_obj_phase2, lower=-15, upper=15, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))
endTime = Sys.time()
endTime-startTime


#### Activating all parameter estimation including especially the ZI part
Map_phase3 <- list()
Map_phase3$beta <- matrix(1:(ncol(data$tmb_data$X)*data$tmb_data$Nage), nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)  
Map_phase3$beta <- factor(Map_phase3$beta)         
if(conf$ARorIID[2] != 2) Map_phase3$logsds <- factor(rep(NA, length(data$tmb_params$logsds)))

if (conf$density_dependence == FALSE){
  Map_phase3$s50_absc <- factor(NA)
  Map_phase3$s50 <- factor(NA)
  Map_phase3$logslope_absc <- factor(NA)
  Map_phase3$logslope <- factor(NA)
}

if (conf$mixture_model == 3) {
  if(conf$ARorIID[1] == 0) Map_phase3$transf_rho_absc <- factor(rep(NA, length(data$tmb_params$transf_rho_absc)))
  if(conf$ARorIID[1] == 1) Map_phase3$transf_rho_absc <- factor(rep(20000, length(data$tmb_params$transf_rho_absc)))
  if(conf$ARorIID[2] != 0) Map_phase3$transf_rho <- factor(rep(NA, length(data$tmb_params$transf_rho)))
  if(conf$ARorIID[2] == 1) Map_phase3$transf_rho <- factor(rep(300000, length(data$tmb_params$transf_rho)))
  
	Map_phase3$transf_rho_age_absc <- factor(rep(NA, length(data$tmb_params$transf_rho_age_absc)))
	Map_phase3$transf_rho_age <- factor(rep(NA, length(data$tmb_params$transf_rho_age)))

  # Map_phase3$transf_rho <- factor(rep(NA, length(data$tmb_params$transf_rho)))
  # Map_phase3$logTauE <- factor(rep(NA, length(data$tmb_params$logTauE)))
	Map_phase3$logKappa <- factor(rep(123456, length(data$tmb_params$logKappa)))
	Map_phase3$logKappa_absc <- factor(rep(1234567, length(data$tmb_params$logKappa_absc)))
  
  Map_phase3$ln_tau_G_absc <- factor(rep(NA, length(data$tmb_params$ln_tau_G_absc)))
  Map_phase3$RE_absc <- factor(rep(NA,length(data$tmb_params$RE_absc)))
  
  Map_phase3$logTauE_absc <- factor(rep(30000, length(data$tmb_params$logTauE_absc)))
  Map_phase3$logTauE_absc <- factor(rep(30000, length(data$tmb_params$logTauE_absc)))
  # Map_phase3$logKappa <- factor(rep(40000, length(data$tmb_params$logKappa)))
  # Map_phase3$logTauE <- factor(rep(40000, length(data$tmb_params$logKappa)))
}

if (conf$family != 0) Map_phase3$thetaf = factor(rep(NA, length(data$tmb_params$thetaf)))

if (conf$keep_omega == FALSE) {
  Map_phase3$omega_absc <- factor(rep(NA, length(data$tmb_params$omega_absc)))
  Map_phase3$omega <- factor(rep(NA, length(data$tmb_params$omega)))
  Map_phase3$logTauO_absc <- factor(rep(NA, length(data$tmb_params$logTauO_absc)))
  Map_phase3$logTauO <- factor(rep(NA, length(data$tmb_params$logTauO)))
}  
  # use the phase 2 parameter estimates with the map setting
  # positive model init from phase 2
  old_par <- set_par_value(opt_phase1$par, "beta")
  map_beta <- Map_phase3$beta
  map_beta <- factor(map_beta, labels = seq(1, sum(!is.na(map_beta))))		
  data$tmb_params$beta <- matrix(old_par[map_beta], nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)   
  data$tmb_params$ln_tau_G <- matrix(set_par_value(opt_phase2$par, "ln_tau_G"), nrow=length(data$tmb_data$nobs_RE), ncol= data$tmb_data$Nage, byrow=F)  
  data$tmb_params$RE <- matrix(set_par_value(tmb_obj_phase2$env$last.par.best, "RE"), nrow=data$tmb_data$nobs_RE, ncol= data$tmb_data$Nage, byrow=F) 
  data$tmb_params$logTauE <- set_par_value(opt_phase2$par, "logTauE")
  if(conf$ARorIID[2] == 1) data$tmb_params$transf_rho <- rep(set_par_value(opt_phase2$par, "transf_rho"), length(data$tmb_params$transf_rho))
  data$tmb_params$logKappa <- rep(set_par_value(opt_phase2$par, "logKappa"), length(data$tmb_params$logKappa))
  data$tmb_params$epsilon_st <- array(set_par_value(tmb_obj_phase2$env$last.par.best, "epsilon_st"), dim=c(data$tmb_data$Nmesh, data$tmb_data$Nyear, data$tmb_data$Nage))
  # 
  
  
  tmb_obj_phase3 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase3,
                                   random = c("RE", "epsilon_st", "epsilon_st_absc"), DLL = dllversion, silent = TRUE)
  startTime = Sys.time()
  opt_phase3 <- fit_tmb(tmb_obj_phase3, lower=-15, upper=15, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))          
  endTime = Sys.time()
  endTime-startTime
 
  
  sdrep <- NA #sdreport(tmb_obj_phase2)
  # 
  dens_pred <- NA
  epsilon_pred <- NA
  
  if (conf$Do_predict == 1) {
    if (conf$bias_correct == TRUE) sdrep <- sdreport(tmb_obj_phase3, bias.correct =TRUE, bias.correct.control = list(sd= TRUE))
    if (conf$bias_correct == FALSE) sdrep <- sdreport(tmb_obj_phase3, bias.correct =FALSE)
    # extract the density estimates on the prediction grid
    pred_grid = data$pred_grid
    mu_proj2_new <- tmb_obj_phase3$report()$mu_proj
    dens_pred <- data.table::as.data.table(mu_proj2_new)
    colnames(dens_pred) <- c("LOC", "YEAR", "AGE", "Pred")
    dens_pred$YEAR <- as.factor(dens_pred$YEAR)
    dens_pred$YEAR <- as.factor(as.numeric(as.character(factor(dens_pred$YEAR, labels=sort(unique(data$data$YEAR))))))
    pred_grid$YEAR <- as.factor(as.numeric(as.character(factor(pred_grid$YEAR, labels=sort(unique(data$data$YEAR))))))
    pred_grid$LOC <- rep(1:(nrow(pred_grid)/conf$Nyear), conf$Nyear)
    dens_pred <- dens_pred %>% left_join(pred_grid)
    dens_pred$logPred <- log(dens_pred$Pred)
    
    # extract estimated spatio-temporal field on the prediction grid
    epsilon2_new <- tmb_obj_phase3$report()$epsilon_st_A_proj
    epsilon_pred <- data.table::as.data.table(epsilon2_new)
    colnames(epsilon_pred) <- c("LOC", "YEAR", "AGE", "Pred")
    epsilon_pred$YEAR <- as.factor(epsilon_pred$YEAR)
    epsilon_pred$YEAR <- as.factor(as.numeric(as.character(factor(epsilon_pred$YEAR, labels=sort(unique(data$data$YEAR))))))
    pred_grid$LOC <- rep(1:(nrow(pred_grid)/conf$Nyear), conf$Nyear)
    epsilon_pred <- epsilon_pred %>% left_join(pred_grid)
  }
  
  
  output_run <- list(obj = tmb_obj_phase3, opt = opt_phase3, data = data, rep = rep, sdrep = sdrep, map= Map_phase3, dens_pred=dens_pred, epsilon_pred=epsilon_pred)
  
  return(output_run)
} 
  