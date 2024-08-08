##' Run the model with the specified configurations
##' @param data$tmb_data data input for the TMB model run
##' @param data$tmb_params parameter input for the TMB model run
##' @param conf  configurations file for the model and data
##' @details the model is run in two steps to make model parameters estimable
##' @return different outputs needed for summarizing and plotting the results at later stage
##' @useDynLib mackerel_ST
##' @export fit_mackerelST
##' 


fit_mackerelST <- function(data, conf, conf_extra=NULL) {
  
  ### Step1: (not estimating the spatio-temporal component)
 
    if (conf$mixture_model == 0 & conf$density_dependence == FALSE) { version <- paste0(getwd(), "/src/spatial_SDM_RE"); dllversion = "spatial_SDM_RE"} 
    if (conf$mixture_model == 1) { version <- paste0(getwd(), "/src/spatial_SDM_mixture"); dllversion = "spatial_SDM_mixture"}
    if (conf$mixture_model == 2) { version <- paste0(getwd(), "/src/spatial_SDM_variablebarrier"); dllversion = "spatial_SDM_variablebarrier"}
    if (conf$mixture_model == 3) { version <- paste0(getwd(), "/src/spatial_SDM_RE_delta"); dllversion = "spatial_SDM_RE_delta"}
    if (conf$mixture_model == 0 & conf$density_dependence == TRUE) { version <- paste0(getwd(), "/src/spatial_SDM_RE_dd"); dllversion = "spatial_SDM_RE_dd"}
    if (!is.null(conf$cohort)) { version <- paste0(getwd(), "/src/spatial_SDM_cohort"); dllversion = "spatial_SDM_cohort"}
    if (conf$mixture_model == 4) { version <- paste0(getwd(), "/src/spatial_SDM_condlogit"); dllversion = "spatial_SDM_condlogit"}

    # dyn.unload(dynlib(version))
    # if (.Platform$OS.type == "unix") compile(paste0(version, ".cpp"),"-O0 -g")
    # if (.Platform$OS.type == "windows") compile(paste0(version, ".cpp"))#, "-O1 -g", DLLFLAGS="")
    dyn.load(dynlib(version))
    
    Map_phase1 <- list()
    Map_phase1$beta <- matrix(1:(ncol(data$tmb_data$X)*data$tmb_data$Nage), nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)  
    Map_phase1$beta <- factor(Map_phase1$beta)        
    Map_phase1$epsilon_st <- factor(rep(NA, length(data$tmb_params$epsilon_st)))
    Map_phase1$s50 <- factor(NA)
    Map_phase1$logslope <- factor(NA)
    Map_phase1$transf_rho <- factor(rep(NA, length(data$tmb_params$transf_rho)))
    Map_phase1$logTauE <- factor(rep(NA, length(data$tmb_params$logTauE)))
    Map_phase1$logKappa <- factor(rep(NA, length(data$tmb_params$logKappa)))
    # Map_phase1$logKappa <- factor(rep(4321, length(data$tmb_params$logKappa)))
    Map_phase1$ln_tau_G <- factor(rep(1234,length(data$tmb_params$ln_tau_G)))
    
    # the SVC effect to include the density dependent (population size) effect on trend in local density changes
    if (conf$density_dependence == TRUE)
    {
      Map_phase1$teta_svc <- factor(rep(NA,length(data$tmb_params$teta_svc)))
      Map_phase1$logTauSVC <- factor(NA)
    }
    
		if(conf$keep_epsilon == FALSE) {
		  Map_phase1$logsds <- factor(rep(NA, length(data$tmb_params$logsds)))
		  Map_phase1$logKappa <- factor(rep(NA, length(data$tmb_params$logKappa)))
		}
    if(conf$corr_str != 1 | conf$ARorIID != 0) {
      Map_phase1$logsds <- factor(rep(NA, length(data$tmb_params$logsds)))
    }    
    if (conf$corr_str != 3)  
    {
      Map_phase1$transf_rho_age <- factor(NA)
    }
    if (conf$add_nugget != 1)  
    {
      Map_phase1$transf_rho_nugget <- factor(NA)
      Map_phase1$logTau_nugget <- factor(NA)
      Map_phase1$nugget_effect <- factor(rep(NA, length(data$tmb_params$nugget_effect)))
    }
    
    if (conf$mixture_model == 3){
      Map_phase1$epsilon_st_absc <- factor(rep(NA, length(data$tmb_params$epsilon_st_absc)))
      Map_phase1$s50_absc <- factor(NA)
      Map_phase1$logslope_absc <- factor(NA)
      Map_phase1$transf_rho_absc <- factor(rep(NA, length(data$tmb_params$transf_rho_absc)))
      Map_phase1$logTauE_absc <- factor(rep(NA, length(data$tmb_params$logTauE_absc)))
      Map_phase1$logKappa_absc <- factor(rep(NA, length(data$tmb_params$logKappa_absc)))
      Map_phase1$omega_absc <- factor(rep(NA, length(data$tmb_params$omega_absc)))
      Map_phase1$logTauO_absc <- factor(rep(NA, length(data$tmb_params$logTauO_absc)))
      
      Map_phase1$ln_tau_G_absc <- factor(rep(NA, length(data$tmb_params$ln_tau_G_absc)))# not estimating the random vessel effect
      Map_phase1$RE_absc <- factor(rep(NA,length(data$tmb_params$RE_absc)))# not estimating the random vessel effect
      
      Map_phase1$logsds <- factor(rep(NA,length(data$tmb_params$logsds)))# not estimating the random vessel effect
      if (conf$corr_str != 3) Map_phase1$transf_rho_age_absc  <- factor(NA)
    }
    
    if (conf$family != 0) Map_phase1$thetaf = factor(rep(NA, length(data$tmb_params$thetaf)))
      
    if (conf$mixture_model == 1){
      data$tmb_params$beta_mix = data$tmb_params$beta_mix - 14
    }
    
		if (conf$cohort == 0){
			Map_phase1$init_age <- factor(rep(NA, length(data$tmb_params$init_age)))
			Map_phase1$init_year <- factor(rep(NA, length(data$tmb_params$init_year)))
			Map_phase1$RE_diff <- factor(rep(NA, length(data$tmb_params$RE_diff)))
			Map_phase1$logcohort_SD <- factor(NA)
	}
		
		
    if (!is.null(conf$conf_model)) Map_phase1 <- append(Map_phase1, conf$conf_model)
    
    if (FALSE %in% is.na(attributes(conf)$RE_effects)) random_var <- c("RE") 
    if (TRUE %in% is.na(attributes(conf)$RE_effects)) random_var <- NULL
    if (conf$cohort == 1) random_var = c(random_var, "RE_diff")
    if (conf$add_nugget == 1) random_var = c(random_var, "nugget_effect")

     if (conf$keep_omega == FALSE) {
      Map_phase1$omega <- factor(rep(NA, length(data$tmb_params$omega)))
      Map_phase1$logTauO <- factor(rep(NA, length(data$tmb_params$logTauO)))
      if (FALSE %in% is.na(attributes(conf)$RE_effects)) {
        if (conf$mixture_model %in% c(0,2)) tmb_obj_phase1 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase1,
                                                       random = random_var, DLL = dllversion, silent = TRUE)
        if (conf$cohort == 1) tmb_obj_phase1 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase1,
                                                       random = random_var, DLL = dllversion, silent = TRUE)
        if (conf$mixture_model == 3) tmb_obj_phase1 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase1,
                                                       random = random_var, DLL = dllversion, silent = TRUE)
        if (conf$mixture_model == 1) tmb_obj_phase1 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase1,
                                                       random = random_var, DLL = dllversion, silent = TRUE)
      }
      if (TRUE %in% is.na(attributes(conf)$RE_effects)) {
        Map_phase1$RE <- factor(rep(NA,length(data$tmb_params$RE)))
        Map_phase1$ln_tau_G <- factor(rep(NA,length(data$tmb_params$ln_tau_G)))
        if (conf$mixture_model != 1) tmb_obj_phase1 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase1,
                                                       random = NULL, DLL = dllversion, silent = TRUE)
        if (conf$mixture_model == 1) tmb_obj_phase1 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase1,
                                                       random = NULL, DLL = dllversion, silent = TRUE)
      }
      
    }
    # if (conf$keep_omega == TRUE) {
    #   tmb_obj_phase1 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase1,
    #                                                    random = "omega", DLL = "spatial_SDM_mixture", silent = TRUE)
    # 
    # }
    
    opt_phase1 <- fit_tmb(tmb_obj_phase1, lower=-15, upper=15, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))          
		# run1 = list(opt=opt_phase1, obj=tmb_obj_phase1, data = data, conf=conf)
    # save(run1, file="aaa.Rdata")
		# fixed1 <- tmb_obj_phase1$report(tmb_obj_phase1$env$last.par.best)$fixed
		# beta1 <- opt_phase1$par
		# fixed2 <- tmb_obj_phase1$report(tmb_obj_phase1$env$last.par.best)$fixed
		# beta2 <- opt_phase1$par
		# fixed3 <- tmb_obj_phase1$report(tmb_obj_phase1$env$last.par.best)$fixed_noRE
    # opt_phase1 <- nlminb(start=tmb_obj_phase1$par, objective = tmb_obj_phase1$fn, gradient = tmb_obj_phase1$gr, lower=-15, upper=15)

  # detecting identifiability issue (whether label switching has been taken care of)
  # check_estimability(tmb_obj_phase1)
  
  # Check if anything hits the bounds
  opt_phase1$diagnostics[which(abs(opt_phase1$diagnostics$MLE) > 14),]
  
  
  ### Step2: Putting the spatio-temporal component into the estimation (starting at the value from step1)
  if (conf$keep_epsilon == TRUE) {

    if (conf$mixture_model != 3) {
    Map_phase2 <- list()
    Map_phase2$beta <- matrix(1:(ncol(data$tmb_data$X)*data$tmb_data$Nage), nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)  
    Map_phase2$beta <- factor(Map_phase2$beta)         

    Map_phase2$ln_tau_G <- factor(rep(1234,length(data$tmb_params$ln_tau_G)))
    Map_phase2$logKappa <- factor(rep(4321, length(data$tmb_params$logKappa)))
    
    Map_phase2$s50 <- factor(NA)
    Map_phase2$logslope <- factor(NA)

    Map_phase2$transf_rho <- factor(rep(NA, length(data$tmb_params$transf_rho)))
    if (conf$mixture_model == 3) {
      Map_phase2$beta_absc <- matrix(1:(ncol(data$tmb_data$X)*data$tmb_data$Nage), nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)  
      Map_phase2$beta_absc <- factor(Map_phase2$beta_absc) 
      Map_phase2$logTauE <- factor(rep(NA, length(data$tmb_params$logTauE)))

      Map_phase2$s50_absc <- factor(NA)
      Map_phase2$logslope_absc <- factor(NA)

      Map_phase2$transf_rho_absc <- factor(rep(NA, length(data$tmb_params$transf_rho_absc)))
      if (conf$corr_str != 3)  
      {
        Map_phase2$transf_rho_age_absc  <- factor(NA)
      }
      if (conf$keep_omega == FALSE) {
        Map_phase2$omega_absc <- factor(rep(NA, length(data$tmb_params$omega_absc)))
        Map_phase2$logTauO_absc <- factor(rep(NA, length(data$tmb_params$logTauO_absc)))
      }
      old_par_absc <- set_par_value(opt_phase1$par, "beta_absc")
      map_beta_absc <- Map_phase2$beta_absc
      map_beta_absc <- factor(map_beta_absc, labels = seq(1, sum(!is.na(map_beta_absc))))		
      data$tmb_params$beta_absc <- matrix(old_par_absc[map_beta_absc], nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)   
      
    }
    
    if (!is.null(conf$conf_model)) Map_phase2 <- append(Map_phase2, conf$conf_model)
    
		if (!is.null(conf_extra$map)) Map_phase2$beta <- factor(conf_extra$beta)
    
    if(conf$corr_str != 1 | conf$ARorIID != 0) {
      Map_phase2$logsds <- factor(rep(NA, length(data$tmb_params$logsds)))
    } 
    
    if (2 %in% conf$ARorIID)  
    {
      Map_phase2$logTauE <- factor(rep(NA, length(data$tmb_params$logTauE)))
      # Map_phase2$logKappa <- factor(rep(1234, length(data$tmb_params$logKappa)))
      # Map_phase2$logsds <- factor(rep(NA, length(data$tmb_params$logsds)))
    }
    
    if (conf$corr_str == 0)  Map_phase2$transf_rho <- factor(rep(NA, length(Map_phase2$transf_rho)))
    if (conf$corr_str != 3)  
    {
      Map_phase2$transf_rho_age <- factor(NA)
    }
 
 		if (conf$cohort == 0){
			Map_phase2$init_age <- factor(rep(NA, length(data$tmb_params$init_age)))
			Map_phase2$init_year <- factor(rep(NA, length(data$tmb_params$init_year)))
			Map_phase2$RE_diff <- factor(rep(NA, length(data$tmb_params$RE_diff)))
			Map_phase2$logcohort_SD <- factor(NA)
		}
		
    if (conf$add_nugget != 1)  
    {
      Map_phase2$transf_rho_nugget <- factor(NA)
      Map_phase2$logTau_nugget <- factor(NA)
      Map_phase2$nugget_effect <- factor(rep(NA, length(data$tmb_params$nugget_effect)))
    }
       
    if (conf$family != 0) Map_phase2$thetaf = factor(rep(NA, length(data$tmb_params$thetaf)))
    
    if (conf$corr_str == 3)  {
      Map_phase2$transf_rho <- factor(rep(20000, length(data$tmb_params$transf_rho)))
      Map_phase2$logTauE <- factor(rep(30000, length(data$tmb_params$logTauE)))
      Map_phase2$logKappa <- factor(rep(40000, length(data$tmb_params$logKappa)))
    }
    
    if (2 %in% conf$ARorIID & conf$mixture_model == 3)  
    {
      Map_phase2$transf_rho <- factor(rep(NA, length(data$tmb_params$transf_rho)))
    }
    
    if (conf$keep_omega == FALSE) {
      Map_phase2$omega <- factor(rep(NA, length(data$tmb_params$omega)))
      Map_phase2$logTauO <- factor(rep(NA, length(data$tmb_params$logTauO)))
      
			# use the phase 1 parameter estimates with the map setting
			old_par <- set_par_value(opt_phase1$par, "beta")
			map_beta <- Map_phase2$beta
			map_beta <- factor(map_beta, labels = seq(1, sum(!is.na(map_beta))))		
			data$tmb_params$beta <- matrix(old_par[map_beta], nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)   

			data$tmb_params$ln_phi <- set_par_value(opt_phase1$par, "ln_phi")   
			if (conf$cohort == 1) data$tmb_params$init_year <- set_par_value(opt_phase1$par, "init_year")   
			if (conf$cohort == 1) data$tmb_params$init_age <- set_par_value(opt_phase1$par, "init_age")   
			if (conf$cohort == 1) data$tmb_params$logcohort_SD <- set_par_value(opt_phase1$par, "logcohort_SD")   
			# data$tmb_params$ln_tau_G <- matrix(set_par_value(opt_phase1$par, "ln_tau_G"), ncol= data$tmb_data$Nage, byrow=F)   
			data$tmb_params$ln_tau_G <- matrix(rep(set_par_value(opt_phase1$par, "ln_tau_G"),data$tmb_data$Nage), ncol= data$tmb_data$Nage, byrow=F)   
			data$tmb_params$RE <- matrix(set_par_value(tmb_obj_phase1$env$last.par.best, "RE"), nrow=data$tmb_data$nobs_RE, ncol= data$tmb_data$Nage, byrow=F)    

			if (conf$mixture_model == 3) {
			  # these are all the parameters for the positive model
  			Map_phase2$beta <- factor(rep(NA, length(data$tmb_params$beta)))       
  			Map_phase2$epsilon_st <- factor(rep(NA, length(data$tmb_params$epsilon_st)))
  			Map_phase2$transf_rho <- factor(rep(NA, length(data$tmb_params$transf_rho)))
  			Map_phase2$transf_rho_age <- factor(NA)
  			Map_phase2$logKappa <- factor(rep(NA, length(data$tmb_params$logKappa)))
  			Map_phase2$logTauE <- factor(rep(NA, length(data$tmb_params$logTauE)))
  			Map_phase2$ln_tau_G <- factor(rep(NA, length(data$tmb_params$ln_tau_G)))
  			Map_phase2$RE <- factor(rep(NA, length(data$tmb_params$RE)))
  			Map_phase2$ln_phi <- factor(rep(NA, length(data$tmb_params$ln_phi)))
			}
			
			if (FALSE %in% is.na(attributes(conf)$RE_effects)) random_var <- c("RE", "epsilon_st") 
			if (TRUE %in% is.na(attributes(conf)$RE_effects)) random_var <- c("epsilon_st")
			if (conf$cohort == 1) random_var = c(random_var, "RE_diff")
			if (conf$add_nugget == 1) random_var = c(random_var, "nugget_effect")
			
      if (FALSE %in% is.na(attributes(conf)$RE_effects)) {
        if (conf$mixture_model %in% c(0,2) & conf$density_dependence == FALSE) tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
                                                       random = c("RE", "epsilon_st"), DLL = dllversion, silent = TRUE)
        if (conf$cohort == 1) tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
                                                       random = random_var, DLL = dllversion, silent = TRUE)
        if (conf$mixture_model %in% c(0,2) & conf$density_dependence == TRUE) tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
                                                       random = c("RE", "epsilon_st","teta_svc"), DLL = dllversion, silent = TRUE)
        if (conf$mixture_model == 3) tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
                                                       random = c("RE_absc", "RE", "epsilon_st_absc"), DLL = dllversion, silent = TRUE)
        if (conf$mixture_model == 1) tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
                                                       random = c("RE", "epsilon_st"), DLL = dllversion, silent = TRUE)
      }
      if (TRUE %in% is.na(attributes(conf)$RE_effects)) {
        Map_phase2$RE <- factor(rep(NA,length(data$tmb_params$RE)))
        Map_phase2$ln_tau_G <- factor(rep(NA,length(data$tmb_params$ln_tau_G)))
        
        if (conf$mixture_model != 1) tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
                                                       random = c("epsilon_st"), DLL = dllversion, silent = TRUE)
        if (conf$mixture_model == 1) tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
                                                       random = c("epsilon_st"), DLL = dllversion, silent = TRUE)
      }
    }
    # if (conf$keep_omega == TRUE) {
    #   data$tmb_params$beta <- matrix(set_par_value(opt_phase1, "beta"), nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)   # making sure to keep at 0 the age1 groups fixed effects
    #   data$tmb_params$logTauO <- set_par_value(opt_phase1, "logTauO")
    #   data$tmb_params$logKappa <- set_par_value(opt_phase1, "logKappa")
    #   
    #   tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
    #                                                    random = c("omega", "epsilon_st"), DLL = "spatial_SDM_mixture", silent = FALSE)
    # }
    
    opt_phase2 <- fit_tmb(tmb_obj_phase2, lower=-15, upper=15, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))          
    }
    
    if (conf$mixture == 3) {
      
      Map_phase2 <- list()
      Map_phase2$beta <- matrix(1:(ncol(data$tmb_data$X)*data$tmb_data$Nage), nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)  
      Map_phase2$beta <- factor(Map_phase2$beta)        
      # Map_phase2$epsilon_st <- factor(rep(NA, length(data$tmb_params$epsilon_st)))
      Map_phase2$s50 <- factor(NA)
      Map_phase2$logslope <- factor(NA)
      Map_phase2$transf_rho <- factor(rep(NA, length(data$tmb_params$transf_rho)))
      # Map_phase2$logTauE <- factor(rep(NA, length(data$tmb_params$logTauE)))
      # Map_phase2$logKappa <- factor(rep(NA, length(data$tmb_params$logKappa)))
      
      if(conf$keep_epsilon == FALSE) Map_phase2$logsds <- factor(rep(NA, length(data$tmb_params$logsds)))
      
      if (conf$corr_str != 3)  
      {
        Map_phase2$transf_rho_age <- factor(NA)
      }
      
      if (conf$mixture_model == 3){
        Map_phase2$epsilon_st_absc <- factor(rep(NA, length(data$tmb_params$epsilon_st_absc)))
        Map_phase2$s50_absc <- factor(NA)
        Map_phase2$logslope_absc <- factor(NA)
        Map_phase2$transf_rho_absc <- factor(rep(NA, length(data$tmb_params$transf_rho_absc)))
        Map_phase2$logTauE_absc <- factor(rep(NA, length(data$tmb_params$logTauE_absc)))
        Map_phase2$logKappa_absc <- factor(rep(NA, length(data$tmb_params$logKappa_absc)))
        Map_phase2$omega_absc <- factor(rep(NA, length(data$tmb_params$omega_absc)))
        Map_phase2$logTauO_absc <- factor(rep(NA, length(data$tmb_params$logTauO_absc)))
        
        Map_phase2$ln_tau_G_absc <- factor(rep(NA, length(data$tmb_params$ln_tau_G_absc)))# not estimating the random vessel effect
        Map_phase2$RE_absc <- factor(rep(NA,length(data$tmb_params$RE_absc)))# not estimating the random vessel effect
        
        Map_phase2$logsds <- factor(rep(NA,length(data$tmb_params$logsds)))# not estimating the random vessel effect
        if (conf$corr_str != 3) Map_phase2$transf_rho_age_absc  <- factor(NA)
      }
      
      if (conf$family != 0) Map_phase2$thetaf = factor(rep(NA, length(data$tmb_params$thetaf)))
      
      if (conf$mixture_model == 1){
        data$tmb_params$beta_mix = data$tmb_params$beta_mix - 14
      }
      
      if (!is.null(conf$conf_model)) Map_phase2 <- append(Map_phase2, conf$conf_model)
      
      if (!is.null(conf_extra$map)) Map_phase2$beta <- factor(conf_extra$beta)
      
      if (conf$keep_omega == FALSE) {
        Map_phase2$omega <- factor(rep(NA, length(data$tmb_params$omega)))
        Map_phase2$logTauO <- factor(rep(NA, length(data$tmb_params$logTauO)))
        
        # use the phase 1 parameter estimates with the map setting
        old_par <- set_par_value(opt_phase1$par, "beta")
        map_beta <- Map_phase2$beta
        map_beta <- factor(map_beta, labels = seq(1, sum(!is.na(map_beta))))		
        data$tmb_params$beta <- matrix(old_par[map_beta], nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)   
        
        data$tmb_params$ln_phi <- set_par_value(opt_phase1$par, "ln_phi")   
        data$tmb_params$ln_tau_G <- matrix(set_par_value(opt_phase1$par, "ln_tau_G"), ncol= data$tmb_data$Nage, byrow=F)   
        data$tmb_params$RE <- matrix(set_par_value(tmb_obj_phase1$env$last.par.best, "RE"), nrow=data$tmb_data$nobs_RE, ncol= data$tmb_data$Nage, byrow=F)    
        
        
        if (FALSE %in% is.na(attributes(conf)$RE_effects)) {
          if (conf$mixture_model %in% c(0,2)) tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
                                                                               random = c("RE"), DLL = dllversion, silent = TRUE)
          if (conf$mixture_model == 3) tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
                                                                        random = c("RE", "epsilon_st"), DLL = dllversion, silent = TRUE)
          if (conf$mixture_model == 1) tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
                                                                        random = c("RE"), DLL = dllversion, silent = TRUE)
        }
        if (TRUE %in% is.na(attributes(conf)$RE_effects)) {
          Map_phase2$RE <- factor(rep(NA,length(data$tmb_params$RE)))
          Map_phase2$ln_tau_G <- factor(rep(NA,length(data$tmb_params$ln_tau_G)))
          if (conf$mixture_model != 1) tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
                                                                        random = NULL, DLL = dllversion, silent = TRUE)
          if (conf$mixture_model == 1) tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
                                                                        random = NULL, DLL = dllversion, silent = TRUE)
        }
        
      }
      # if (conf$keep_omega == TRUE) {
      #   tmb_obj_phase2 <- TMB::MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase2,
      #                                                    random = "omega", DLL = "spatial_SDM_mixture", silent = TRUE)
      # 
      # }
      
      opt_phase2 <- fit_tmb(tmb_obj_phase2, lower=-15, upper=15, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))          
    }
    # Check if anything hits the bounds
    opt_phase2$diagnostics[which(abs(opt_phase2$diagnostics$MLE) > 14),]

  
    ### Step3: Activating all parameters for delta model
    if (conf$mixture_model == 3) {
      Map_phase3 <- list()
      Map_phase3$beta <- matrix(1:(ncol(data$tmb_data$X)*data$tmb_data$Nage), nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)  
      Map_phase3$beta <- factor(Map_phase3$beta)         
      # Map_phase3$beta_absc <- matrix(1:(ncol(data$tmb_data$X_pres)*data$tmb_data$Nage), nrow=ncol(data$tmb_data$X_pres), ncol= data$tmb_data$Nage, byrow=F)  
      # Map_phase3$beta_absc <- factor(Map_phase3$beta_absc)         
      Map_phase3$beta_absc <- factor(rep(NA, length(data$tmb_params$beta_absc)))
      
        Map_phase3$s50_absc <- factor(NA)
        Map_phase3$s50 <- factor(NA)
        Map_phase3$logslope_absc <- factor(NA)
        Map_phase3$logslope <- factor(NA)

      if (conf$mixture_model != 3) {
        Map_phase3$beta_absc <- factor(rep(NA, length(data$tmb_params$beta_absc)))
        Map_phase3$epsilon_st_absc <- factor(rep(NA, length(data$tmb_params$epsilon_st_absc)))
        Map_phase3$logKappa_absc <- factor(rep(NA, length(data$tmb_params$logKappa_absc)))
        Map_phase3$logTauE_absc <- factor(rep(NA, length(data$tmb_params$logTauE_absc)))
        Map_phase3$ln_tau_G_absc <- factor(rep(NA, length(data$tmb_params$ln_tau_G_absc)))
        Map_phase3$RE_absc <- factor(rep(NA, length(data$tmb_params$RE_absc)))
        Map_phase3$transf_rho_absc <- factor(rep(NA, length(data$tmb_params$transf_rho_absc)))
      }
        
      if (conf$mixture_model == 3) {
        if(conf$ARorIID[1] == 0) Map_phase3$transf_rho_absc <- factor(rep(NA, length(data$tmb_params$transf_rho_absc)))
        if(conf$ARorIID[1] == 1) Map_phase3$transf_rho_absc <- factor(rep(20000, length(data$tmb_params$transf_rho_absc)))
        if(conf$ARorIID[2] == 0)  Map_phase3$transf_rho <- factor(rep(NA, length(data$tmb_params$transf_rho)))
        if(conf$ARorIID[2] == 2) Map_phase3$logTauE <- factor(rep(50000, length(data$tmb_params$logTauE)))
        Map_phase3$ln_tau_G_absc <- factor(rep(NA, length(data$tmb_params$ln_tau_G_absc)))
        Map_phase3$RE_absc <- factor(rep(NA,length(data$tmb_params$RE_absc)))
        
        Map_phase3$logTauE_absc <- factor(rep(30000, length(data$tmb_params$logTauE_absc)))
        Map_phase3$logKappa_absc <- factor(rep(40000, length(data$tmb_params$logKappa_absc)))
        
      }
      
      if (!is.null(conf$conf_model)) Map_phase3 <- append(Map_phase3, conf$conf_model)
      
      if (!is.null(conf_extra$map)) Map_phase3$beta <- factor(conf_extra$beta)
      
      if (conf$corr_str == 0)  Map_phase3$transf_rho <- factor(rep(NA, length(Map_phase3$transf_rho)))
      if (conf$corr_str != 3)  {
        Map_phase3$transf_rho_age_absc <- factor(NA)
        Map_phase3$transf_rho_age <- factor(NA)
      }
      
      if (conf$family != 0) Map_phase3$thetaf = factor(rep(NA, length(data$tmb_params$thetaf)))
      
      if (conf$corr_str == 3)  {
        Map_phase3$transf_rho <- factor(rep(20000, length(data$tmb_params$transf_rho)))
        Map_phase3$logTauE <- factor(rep(30000, length(data$tmb_params$logTauE)))
        Map_phase3$logKappa <- factor(rep(40000, length(data$tmb_params$logKappa)))
      }
      
      if (conf$add_nugget != 1)  
      {
        Map_phase3$transf_rho_nugget <- factor(NA)
        Map_phase3$logTau_nugget <- factor(NA)
        Map_phase3$nugget_effect <- factor(rep(NA, length(data$tmb_params$nugget_effect)))
      }
        
              
      if (conf$keep_omega == FALSE) {
        Map_phase3$omega_absc <- factor(rep(NA, length(data$tmb_params$omega_absc)))
        Map_phase3$omega <- factor(rep(NA, length(data$tmb_params$omega)))
        Map_phase3$logTauO_absc <- factor(rep(NA, length(data$tmb_params$logTauO_absc)))
        Map_phase3$logTauO <- factor(rep(NA, length(data$tmb_params$logTauO)))
        
        # use the phase 1 parameter estimates with the map setting
        # positive model init from phase 1
        old_par <- set_par_value(opt_phase1$par, "beta")
        map_beta <- Map_phase3$beta
        map_beta <- factor(map_beta, labels = seq(1, sum(!is.na(map_beta))))		
        data$tmb_params$beta <- matrix(old_par[map_beta], nrow=ncol(data$tmb_data$X), ncol= data$tmb_data$Nage, byrow=F)   
        data$tmb_params$ln_tau_G <- matrix(set_par_value(opt_phase1$par, "ln_tau_G"), nrow=length(data$tmb_data$nobs_RE), ncol= data$tmb_data$Nage, byrow=F)  
        data$tmb_params$RE <- matrix(set_par_value(tmb_obj_phase1$env$last.par.best, "RE"), nrow=data$tmb_data$nobs_RE, ncol= data$tmb_data$Nage, byrow=F) 
        
        # pos from phase 2
          # old_par_absc <- set_par_value(opt_phase2$par, "beta_absc")
          # map_beta_absc <- Map_phase3$beta_absc
          # map_beta_absc <- factor(map_beta_absc, labels = seq(1, sum(!is.na(map_beta_absc))))		
          # data$tmb_params$beta_absc <- matrix(old_par_absc[map_beta_absc], nrow=ncol(data$tmb_data$X_pres), ncol= data$tmb_data$Nage, byrow=F)   
          # data$tmb_params$logTauE_absc <- set_par_value(opt_phase2$par, "logTauE_absc")   
          # data$tmb_params$logKappa_absc <- set_par_value(opt_phase2$par, "logKappa_absc")   
          data$tmb_params$epsilon_st <- array(set_par_value(tmb_obj_phase2$env$last.par.best, "epsilon_st"), dim=c(data$tmb_data$Nmesh, data$tmb_data$Nyear, data$tmb_data$Nage))    
  
        # Map_phase3$beta_absc <- factor(rep(NA, length(data$tmb_params$beta_absc)))
        # Map_phase3$epsilon_st_absc <- factor(rep(NA, length(data$tmb_params$epsilon_st_absc)))
        # Map_phase3$logKappa_absc <- factor(rep(NA, length(data$tmb_params$logKappa_absc)))
        # Map_phase3$logTauE_absc <- factor(rep(NA, length(data$tmb_params$logTauE_absc)))
        # Map_phase3$ln_tau_G_absc <- factor(rep(NA, length(data$tmb_params$ln_tau_G_absc)))
        # Map_phase3$RE_absc <- factor(rep(NA, length(data$tmb_params$RE_absc)))
  
        if (!is.na(attributes(conf)$RE_effects)) {
          if (conf$mixture_model %in% c(0,2)) tmb_obj_phase3 <- MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase3,
                                                                               random = c("RE", "epsilon_st"), DLL = dllversion, silent = TRUE)
          if (conf$mixture_model == 3) tmb_obj_phase3 <- MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase3,
                                                                        random = c("RE", "epsilon_st_absc", "epsilon_st"), DLL = dllversion, silent = TRUE)
          if (conf$mixture_model == 1) tmb_obj_phase3 <- MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase3,
                                                                        random = c("RE", "epsilon_st"), DLL = dllversion, silent = TRUE)
        }
        if (is.na(attributes(conf)$RE_effects)) {
          Map_phase3$RE <- factor(rep(NA,length(data$tmb_params$RE)))
          Map_phase3$ln_tau_G <- factor(rep(NA,length(data$tmb_params$ln_tau_G)))
          
          if (conf$mixture_model != 1) tmb_obj_phase3 <- MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase3,
                                                                        random = c("epsilon_st"), DLL = dllversion, silent = TRUE)
          if (conf$mixture_model == 1) tmb_obj_phase3 <- MakeADFun(data = data$tmb_data, parameters = data$tmb_params, map = Map_phase3,
                                                                        random = c("epsilon_st"), DLL = dllversion, silent = TRUE)
        }
      }
  
      opt_phase3 <- fit_tmb(tmb_obj_phase3, lower=-15, upper=15, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))          
    }
    
    if (conf$mixture_model != 3) {
      opt_phase3 = opt_phase2
      tmb_obj_phase3 = tmb_obj_phase2
      Map_phase3 = Map_phase2
    }
  }
  
  if (conf$keep_epsilon == FALSE) {
    opt_phase3 = opt_phase1
    tmb_obj_phase3 = tmb_obj_phase1
    Map_phase3 = Map_phase1
  }
  
    opt_phase3$diagnostics
  

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
    # if (conf$keep_epsilon == TRUE) {
    #   epsilon2_new <- tmb_obj_phase3$report()$epsilon_st_A_proj
    #   epsilon_pred <- data.table::as.data.table(epsilon2_new)
    #   colnames(epsilon_pred) <- c("LOC", "YEAR", "AGE", "Pred")
    #   epsilon_pred$YEAR <- as.factor(epsilon_pred$YEAR)
    #   epsilon_pred$YEAR <- as.factor(as.numeric(as.character(factor(epsilon_pred$YEAR, labels=sort(unique(data$data$YEAR))))))
    #   pred_grid$LOC <- rep(1:(nrow(pred_grid)/conf$Nyear), conf$Nyear)
    #   epsilon_pred <- epsilon_pred %>% left_join(pred_grid)
    # }
    
  }
    
  
  output_run <- list(obj = tmb_obj_phase3, opt = opt_phase3, data = data, rep = rep, sdrep = sdrep, map= Map_phase3, dens_pred=dens_pred, epsilon_pred=epsilon_pred)

  return(output_run)
}



##' Run the model with the specified configurations (only used for the combined models with 1. aggregated CPUE, 2. P(age)
##' @param data$tmb_data data input for the TMB model run
##' @param data$tmb_params parameter input for the TMB model run
##' @param conf  configurations file for the model and data
##' @details the model is run in two steps to make model parameters estimable
##' @return different outputs needed for summarizing and plotting the results at later stage
##' @useDynLib mackerel_ST
##' @export fit_mackerelST
##' 

# a simplified version that might work better for the cond logit model
fit_mackerelST_condlogit <- function(data_object, conf, conf_extra=NULL) {
  
  ### Step1: (not estimating the spatio-temporal component)
    version <- paste0(getwd(), "/src/spatial_SDM_condlogit_v2"); dllversion = "spatial_SDM_condlogit_v2"
    
    # dyn.unload(dynlib(version))
    # if (.Platform$OS.type == "unix") compile(paste0(version, ".cpp"),"-O0 -g")
    # if (.Platform$OS.type == "windows") compile(paste0(version, ".cpp"))#, "-O1 -g", DLLFLAGS="")
    dyn.load(dynlib(version))
    

      Map_phase2 <- list()
       if (!is.null(attributes(conf)$fixed_effect)){
         Map_phase2$beta <- matrix(1:(ncol(data_object$tmb_data$X)*(data_object$tmb_data$Nage-1)), nrow=ncol(data_object$tmb_data$X), ncol= data_object$tmb_data$Nage-1, byrow=F)  
         Map_phase2$beta <- factor(Map_phase2$beta)         
       }
      
      # If I want to make sure to aggregate all parameters
      Map_phase2$logTauE <- factor(rep(30000, length(data_object$tmb_params$logTauE)))  # to simplify something
      Map_phase2$logKappa <- factor(rep(4321, length(data_object$tmb_params$logKappa)))
      # Map_phase2$logKappa_cpue <- factor(rep(4321, length(data_object$tmb_params$logKappa_cpue))) # putting the same as Kappa
      
      Map_phase2$s50 <- factor(NA)
      Map_phase2$logslope <- factor(NA)
      
      if (conf$ARorIID_cpue != 1)  Map_phase2$transf_rho_cpue <- factor(rep(NA, length(data_object$tmb_params$transf_rho_cpue)))
      if (conf$ARorIID != 1) Map_phase2$transf_rho <- factor(rep(NA, length(data_object$tmb_params$transf_rho)))
        
      if(conf$corr_str != 1 | conf$ARorIID != 0) {
        Map_phase2$logsds <- factor(rep(NA, length(data_object$tmb_params$logsds)))
      } 
      
      if (2 %in% conf$ARorIID)  
      {
        Map_phase2$logTauE <- factor(rep(NA, length(data_object$tmb_params$logTauE)))
        # Map_phase2$logKappa <- factor(rep(1234, length(data_object$tmb_params$logKappa)))
        # Map_phase2$logsds <- factor(rep(NA, length(data_object$tmb_params$logsds)))
      }
      
      if (conf$corr_str == 0)  Map_phase2$transf_rho <- factor(rep(NA, length(Map_phase2$transf_rho)))
      if (conf$corr_str != 3)  
      {
        Map_phase2$transf_rho_age <- factor(NA)
      }
      
      if (conf$corr_str == 3)  {
        Map_phase2$transf_rho <- factor(rep(20000, length(data_object$tmb_params$transf_rho)))
        Map_phase2$logTauE <- factor(rep(30000, length(data_object$tmb_params$logTauE)))
        Map_phase2$logKappa <- factor(rep(40000, length(data_object$tmb_params$logKappa)))
      }

      if(conf$keep_omega == TRUE) cli_warn(c("keep_omega = FALSE is the only option that is available for now"))
      
      if (conf$cohort == 0){
        Map_phase2$init_age <- factor(rep(NA, length(data_object$tmb_params$init_age)))
        Map_phase2$init_year <- factor(rep(NA, length(data_object$tmb_params$init_year)))
        Map_phase2$RE_diff <- factor(rep(NA, length(data_object$tmb_params$RE_diff)))
        Map_phase2$logcohort_SD <- factor(NA)
      }
      
      random_var = c("epsilon_st", "epsilon_st_cpue")
      if (conf$cohort == 1) random_var = c(random_var, "RE_diff")
      if (conf$cohort == 0) random_var = random_var
      # random_var = c("RE_diff","init_year","init_age", "epsilon_st", "epsilon_st_cpue")
      
      if (conf$keep_omega == FALSE) {

        Map_phase2$omega <- factor(rep(NA, length(data_object$tmb_params$omega)))
        Map_phase2$omega_cpue <- factor(rep(NA, length(data_object$tmb_params$omega_cpue)))
        
        if (FALSE %in% is.na(attributes(conf)$RE_effects))  random_var = c(random_var, "RE")
        if (FALSE %in% is.na(attributes(conf)$RE_effects_cpue))  random_var = c(random_var, "RE_cpue")
        if (TRUE %in% is.na(attributes(conf)$RE_effects)) {
          Map_phase2$RE <- factor(rep(NA,length(data_object$tmb_params$RE)))
          Map_phase2$ln_tau_G <- factor(rep(NA,length(data_object$tmb_params$ln_tau_G)))
        }  
        if (TRUE %in% is.na(attributes(conf)$RE_effects_cpue)) {
          Map_phase2$RE_cpue <- factor(rep(NA,length(data_object$tmb_params$RE_cpue)))
          Map_phase2$ln_tau_G_cpue <- factor(rep(NA,length(data_object$tmb_params$ln_tau_G_cpue)))
        }  
      
       
        tmb_obj_phase2 <- TMB::MakeADFun(data = data_object$tmb_data, parameters = data_object$tmb_params, map = Map_phase2,
                                                                 random = random_var, DLL = dllversion, silent = TRUE)

      }
 
      opt_phase2 <- fit_tmb(tmb_obj_phase2, lower=-15, upper=15, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))          
    
   # Check if anything hits the bounds
    opt_phase2$diagnostics[which(abs(opt_phase2$diagnostics$MLE) > 14),]
    
    
    opt_phase3 = opt_phase2
    tmb_obj_phase3 = tmb_obj_phase2
    Map_phase3 = Map_phase2
  

  if (conf$keep_epsilon == FALSE) {
    opt_phase3 = opt_phase1
    tmb_obj_phase3 = tmb_obj_phase1
    Map_phase3 = Map_phase1
  }
  
  # opt_phase3$diagnostics
  
  
  sdrep <- NA #sdreport(tmb_obj_phase2)
  # 
  dens_pred <- NA
  epsilon_pred <- NA
  
  if (conf$Do_predict == 1) {
    if (conf$bias_correct == TRUE) sdrep <- sdreport(tmb_obj_phase3, bias.correct =TRUE, bias.correct.control = list(sd= TRUE))
    if (conf$bias_correct == FALSE) sdrep <- sdreport(tmb_obj_phase3, bias.correct =FALSE)
    # extract the density estimates on the prediction grid
    pred_grid = data_object$pred_grid
    mu_proj2_new <- tmb_obj_phase3$report()$mu_proj
    dens_pred <- data.table::as.data.table(mu_proj2_new)
    colnames(dens_pred) <- c("LOC", "YEAR", "AGE", "Pred")
    dens_pred$YEAR <- as.factor(dens_pred$YEAR)
    dens_pred$YEAR <- as.factor(as.numeric(as.character(factor(dens_pred$YEAR, labels=sort(unique(data_object$data$YEAR))))))
    pred_grid$YEAR <- as.factor(as.numeric(as.character(factor(pred_grid$YEAR, labels=sort(unique(data_object$data$YEAR))))))
    pred_grid$LOC <- rep(1:(nrow(pred_grid)/conf$Nyear), conf$Nyear)
    dens_pred <- dens_pred %>% left_join(pred_grid)
    dens_pred$logPred <- log(dens_pred$Pred)
    
    # extract estimated spatio-temporal field on the prediction grid
    # if (conf$keep_epsilon == TRUE) {
    #   epsilon2_new <- tmb_obj_phase3$report()$epsilon_st_A_proj
    #   epsilon_pred <- data.table::as.data.table(epsilon2_new)
    #   colnames(epsilon_pred) <- c("LOC", "YEAR", "AGE", "Pred")
    #   epsilon_pred$YEAR <- as.factor(epsilon_pred$YEAR)
    #   epsilon_pred$YEAR <- as.factor(as.numeric(as.character(factor(epsilon_pred$YEAR, labels=sort(unique(data_object$data_object$YEAR))))))
    #   pred_grid$LOC <- rep(1:(nrow(pred_grid)/conf$Nyear), conf$Nyear)
    #   epsilon_pred <- epsilon_pred %>% left_join(pred_grid)
    # }
    
  }
  
  
  output_run <- list(obj = tmb_obj_phase3, opt = opt_phase3, data = data, rep = rep, sdrep = sdrep, map= Map_phase3, dens_pred=dens_pred, epsilon_pred=epsilon_pred)
  
  return(output_run)
}





