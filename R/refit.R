##' Refit the model in case prediction was not activated in the model
##' @param run model run outputs obtained via the 'fit_mackerelST" function 
##' @param pred_data data object for prediction. used only with predict_type == 2 
##' @param predict_type prediction type? 1 = predict at real scale, 2 = predict at link scale. 
##' @param conf configuration file 
##' @details rerun the model based on previous model run to activate the prediction
##' @return returns the obj, opt, report, sdreport, density prediction on grid, st_RE on grid
##' @export


refit <- function(run=run, obj=NULL, map=NULL, oldpar=NULL, pred_data=NULL, predict_type = 1, conf, bias_correct=FALSE, ...) {

  # oldpar <- run$obj$env$last.par.best
  # if (!is.null(oldpar)) oldpar = oldpar
  if (!is.null(obj)) run$obj = obj
  if (!is.null(map)) run$map = map
 
  tmb_data_new <- run$data$tmb_data
   if (!is.null(pred_data)) {  
     X_proj = mgcv::predict.gam(run$data$mgcv_mod, type = "lpmatrix", newdata = pred_data)
     Npred = nrow(X_proj)
     tmb_data_new$X_proj = X_proj
     tmb_data_new$Npred = Npred
     tmb_data_new$do_predict = predict_type
   } else { 
     pred_grid <- run$data$pred_grid
   }
  
   if (predict_type == 2) {
    X_proj = mgcv::predict.gam(run$data$mgcv_mod, type = "lpmatrix", newdata = pred_data)
    Npred = nrow(X_proj)
    tmb_data_new$X_proj = X_proj
    tmb_data_new$Npred = Npred
    tmb_data_new$do_predict = predict_type
   }
  
  if (predict_type == 1) {
    tmb_data_new$do_predict = 1
  }
  
  if (conf$mixture_model == 0) { version <- paste0(getwd(), "/src/spatial_SDM_RE"); dllversion = "spatial_SDM_RE"} 
  if (conf$mixture_model == 1) { version <- paste0(getwd(), "/src/spatial_SDM_mixture"); dllversion = "spatial_SDM_mixture"}
  if (conf$mixture_model == 2) { version <- paste0(getwd(), "/src/spatial_SDM_variablebarrier"); dllversion = "spatial_SDM_variablebarrier"}
  if (conf$cohort == 1) { version <- paste0(getwd(), "/src/spatial_SDM_condlogit_v2"); dllversion = "spatial_SDM_condlogit_v2" }
  if (conf$mixture_model == 3) { version <- paste0(getwd(), "/src/spatial_SDM_RE_delta"); dllversion = "spatial_SDM_RE_delta"}
  dyn.load(dynlib(version))
  
  if (conf$mixture_model != 3) { 
    new_tmb_obj <- TMB::MakeADFun(
      data = tmb_data_new,
      parameters = run$data$tmb_params,
      map = run$map,
      random = c("RE", "epsilon_st"),
      DLL = dllversion,
      silent = TRUE
    )
  }
  if (conf$mixture_model == 3) {
    new_tmb_obj <- TMB::MakeADFun(
      data = tmb_data_new,
      parameters = run$data$tmb_params,
      map = run$map,
      random = c("RE", "epsilon_st_absc", "epsilon_st"),
      DLL = dllversion,
      silent = TRUE
    )
  }
      
    # need to initialize the new TMB object once:
    # new_tmb_obj$fn(oldpar)  
    new_tmb_obj$fn(run$opt$par)  

    # preparing the output
    dens_pred <- NULL
    epsilon_pred <- NULL
    rep <- NULL
    
    if (predict_type == 2) {
      dens_pred <- new_tmb_obj$report()$fixed_proj1
    } 
    
    if (predict_type == 1) {
			if (bias_correct == TRUE) rep <- sdreport(new_tmb_obj, bias.correct =TRUE, bias.correct.control = list(sd= TRUE))
			if (bias_correct == FALSE) rep <- sdreport(new_tmb_obj, bias.correct =FALSE)
 
			# extract the density estimates on the prediction grid
      mu_proj2_new <- new_tmb_obj$report()$mu_proj
      dens_pred <- data.table::as.data.table(mu_proj2_new)
      colnames(dens_pred) <- c("LOC", "YEAR", "AGE", "Pred")
      dens_pred$YEAR <- as.factor(dens_pred$YEAR)
      dens_pred$YEAR <- as.factor(as.numeric(as.character(factor(dens_pred$YEAR, labels=sort(unique(run$data$data$YEAR))))))
      pred_grid$YEAR <- as.factor(as.numeric(as.character(factor(pred_grid$YEAR, labels=sort(unique(run$data$data$YEAR))))))
      dens_pred$AGE <- factor(dens_pred$AGE, levels=conf$ages)
      pred_grid$LOC <- rep(1:(nrow(pred_grid)/conf$Nyear), conf$Nyear)
      dens_pred <- dens_pred %>% left_join(pred_grid)
      dens_pred$logPred <- log(dens_pred$Pred)
      
      # # extract estimated spatio-temporal field on the prediction grid
      # epsilon2_new <- new_tmb_obj$report()$epsilon_st_A_proj
      # epsilon_pred <- data.table::as.data.table(epsilon2_new)
      # colnames(epsilon_pred) <- c("LOC", "YEAR", "AGE", "Pred")
      # epsilon_pred$YEAR <- as.factor(epsilon_pred$YEAR)
      # epsilon_pred$YEAR <- as.factor(as.numeric(as.character(factor(epsilon_pred$YEAR, labels=sort(unique(data$data$YEAR))))))
      # epsilon_pred$AGE <- factor(epsilon_pred$AGE, levels=conf$ages)
      # pred_grid$LOC <- rep(1:(nrow(pred_grid)/conf$Nyear), conf$Nyear)
      # epsilon_pred <- epsilon_pred %>% left_join(pred_grid)
    }
    
    
    return(list(obj = new_tmb_obj, opt=run$opt, rep = new_tmb_obj$report(), sdrep = rep, dens_pred=dens_pred, epsilon_pred=epsilon_pred))
  }



##' Perform a cross validation using the model
##' @param run model run outputs obtained via the 'fit_mackerelST" function 
##' @param data the data file that was generated using the prepare_data function. 
##' @param k_fold=5 the number of fold for the CV. 
##' @param seed the seed to replicate the cross-validation work. 
##' @param conf configuration file 
##' @details rerun the model based on previous model run to activate the prediction
##' @return returns the obj, opt, report, sdreport, density prediction on grid, st_RE on grid
##' @export


cv <- function(run=run, data, k_fold = 5, conf, seed=123, ...) {

  set.seed(seed)
	bla <- data$data
	bla$id <- 1:nrow(bla)
	
	# function to create the k-fold split
	kfold_split <- function(x, k_fold){
		temp = split(x$id, sort(x$id %% k_fold))           # perfoming the split
		folds = as.numeric(rep(names(temp), sapply(temp, length)))+1  # adding 1 so that groups starts at 1
		folds = sample(folds, length(folds), replace=F)   # resample
		return(folds)
	}
			
	# performing the k-fold split by grouping variables		
	qqq <- bla %>% group_by(YEAR, Area) %>% group_map( ~ kfold_split(x=.x, k_fold=k_fold)) %>% map_dfr(~ as.data.frame(.))
	names(qqq) <- "folds"
	bla$folds = qqq
	
	# now need to rerun the model (no change to it) by leaving out the specific fold 
	SSE =  c()
	MRE = c()
	LL =  c()
	ELPD =  c()
	  
	  
	for (cv_fold in 1:k_fold){
	keep = matrix(0, nrow=nrow(bla), ncol=conf$Nage)
	keep[which(bla$folds != cv_fold),] = 1
	
  tmb_data_new <- data$tmb_data
	tmb_data_new$to_keep = keep  

	if (conf$mixture_model == 0& conf$cohort != 1) { dllversion = "spatial_SDM_RE"} 
	if (conf$mixture_model == 1& conf$cohort != 1) { dllversion = "spatial_SDM_mixture"}
	if (conf$mixture_model == 2& conf$cohort != 1) { dllversion = "spatial_SDM_variablebarrier"}
	if (conf$mixture_model == 3& conf$cohort != 1) { dllversion = "spatial_SDM_RE_delta"}
	if (conf$cohort == 1) dllversion = "spatial_SDM_cohort"
	
	if (conf$mixture_model != 3 & conf$cohort != 1) { 
	  new_tmb_obj <- TMB::MakeADFun(
	    data = tmb_data_new,
	    parameters = get_pars(run$obj),
	    map = run$map,
	    random = c("RE", "epsilon_st"),
	    DLL = dllversion,
	    silent = TRUE
	  )
	}
	if (conf$cohort == 1) { 
	  new_tmb_obj <- TMB::MakeADFun(
	    data = tmb_data_new,
	    parameters = get_pars(run$obj),
	    map = run$map,
	    random = c("RE", "RE_diff", "epsilon_st"),
	    DLL = dllversion,
	    silent = TRUE
	  )
	}
	if (conf$mixture_model == 3) {
	  new_tmb_obj <- TMB::MakeADFun(
	    data = tmb_data_new,
	    parameters = get_pars(run$obj),
	    map = run$map,
	    random = c("RE", "epsilon_st_absc", "epsilon_st"),
	    DLL = dllversion,
	    silent = TRUE
	  )
	}      
      
    old_par <- run$opt$par
    # need to initialize the new TMB object once:
    new_tmb_obj$fn(old_par)  

		# optimizing it
    opt_phase3 <- fit_tmb(new_tmb_obj, lower=-15, upper=15, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))          

    # preparing the output
    response <- data$tmb_data$yobs[which(bla$folds == cv_fold),]
    mu <- run$obj$report()$mu[which(bla$folds == cv_fold),]
    
    sse <- apply((response - mu)^2, 2, sum, na.rm=T)
    
    NLL = sapply(1:ncol(response), function(age) ll_tweedie(opt_phase3, withheld_y=response[!is.na(response[,1]),age], withheld_mu=mu[!is.na(response[,1]), age]))
    ll <- apply(NLL, 2, sum)
    elpd <- apply(NLL, 2, function(.x) log_sum_exp(.x) - log(length(.x)))
    
    SSE = rbind(SSE, sse)
    LL = rbind(LL, ll)
    ELPD = rbind(ELPD, elpd)
	}    
    
    return(list(SSE=SSE, LL=LL, ELPD=ELPD, yobs = response, mu=mu))
  }

