#### Main code to examine works done in 
# Effect of environmental drivers on the spatiotemporal distribution of mackerel at age in the Nordic Seas during 2010−20
# by: Ono K. et al. 
# in: ICES 2024 
#
# The current code lacks the underlying survey data but can be requested to the main author after obtaining the approval from PGNAPES.


rm(list=ls())
load(".RData")

### libraries
	library(Matrix)
	library(tidyverse)
	library(sf)
	library(sfheaders)
	library(raster)
	library(sp)
	library(assertthat)
	library(gridExtra)
	library(sdmTMB)
	library(sdmTMBextra)									 
	library(ggplot2)
	library(ggnewscale)
	library(ggpmisc)
	library(ggrepel)
	library(TMB)
	library(TMBhelper)
	library(GGally)
	library(MASS)
	library(DHARMa)
  library(ggpubr)
  library(ggConvexHull) #devtools::install_github("cmartin/ggConvexHull")


### Setting projection method and geographical extent of the study area 
	projection <- 3035 # 32630 or 3035
	if (projection == 32630) projection_km <- "+proj=utm +zone=30 +ellps=WGS84 +units=km +no_defs"
	if (projection == 3035) projection_km <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=km +no_defs"

### Loading data & functions
  source("R/utils.R")                 # bunch of utility functions
  source("R/data_load.R")             # for loading and "checking" the data
  source("R/conf.R")                  # for setting the data/model configuration
  source("R/data.R")                  # for preparing the data
  source("R/run.R")                   # for running the model with the above settings
  source("R/refit.R")                 # for rerunning the model in case prediction was disactivated 
  source("R/residual_diagnostics.R")	# for residual diagnostics
  source("R/Plotting_funcs.R")        # for various plotting functions

### A bit more tidying of the data for the analysis to follow: 
### Below, "dat" is a data frame in wide format with all information on covariates and CPUE 
	dat$SSTfront_scl <- dat$distance_SST_fronts_ostia_scl
	dat$VESSEL <- as.factor(dat$VESSEL)
	
### Adding some polygons info to the data
	polygons_df <- st_as_sf(x=All_strata %>% filter(ID != 13), coords = c("lon", "lat")) %>% st_set_crs(4326) %>% group_by(ID) %>% 
	  summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") %>% st_transform(projection_km)
	## calculating distance of each observation to polygon because not all points are within polygons
  	bla <- st_distance(polygons_df, st_as_sf(x=dat, coords = c("LON", "LAT")) %>% st_set_crs(4326) %>% st_transform(projection_km) %>% st_cast("POINT"))
  	dat$Area = apply(bla, 2, which.min)
  	dat <- dat %>% st_drop_geometry()
  	dat$YearArea = apply(dat[,c('YEAR', 'Area')], 1, function(x) paste(x, collapse="_"))
  	
	
### The data frame names "bio" contains all information on the biological samples along with all other useful covariates
	head(bio)
	bio_sf <- st_as_sf(bio %>% filter(!is.na(LON), !is.na(LAT)), coords = c("LON", "LAT"), crs=4326) %>% st_transform(projection_km)
	bio$X <- st_coordinates(bio_sf)[,1]/1000  # in km
	bio$Y <- st_coordinates(bio_sf)[,2]/1000  # in km
	
	
#======================================================================================================#
#===========  Fitting a spatio-temporal model to the age-based CPUE data ==============================#
#======================================================================================================#
	
### Define configurations (both data and model) for the TMB model run
### inside "formula", one write codes as in other glm(m) packages and it can include smoother as implemented in mgcv (borrows its functionality) 
	conf = conf_TMB(Model_type = "Annual_barrier", # deprecated (just keep it as is)
									keep_omega = FALSE, # keep the fixed spatial field or not
									keep_epsilon = TRUE,# if false here and keep_omega=FALSE, then it is a non-spatial model 
									include_age_correlation = "allinone",     # deprecated. "none", "extra", "with_epsilon", "allinone" (keep it as is)
									plotting = FALSE,
									Do_predict = 1,  # do prediction? FALSE saves some time
									mixture_model = 0,  # 0 = spatial_SDM_RE (default), 1 = spatial_SDM_RE_mixture, 2=spatial_SDM_RE_barriereffect, 3=spatial_SDM_RE_delta
									density_dependence = FALSE,   # include density dependent effect 
									years = 2010:2020, # if including mixed_depth stop at 2019
									plusgroup = 16,  
									ages = 3:10,
									family = 0, # 0 is tweedie, 2 = gamma, 3 = student t
									link = 1, # 1 is log link, 4 is boxcox
									corr_str = 1,  # 1 = AR1 correlation in spatial distribution between age (default), for each year, 2 = AR1 correlation in spatial distribution between year, for each age
									ARorIID = 0, # c(1,2), # 0 = IID, 1 = AR1, 2 = user specified corr matrix
									formula_pres = formula(CPUE ~ 1),  # only used with delta model 
									# formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_0m_scl,k=3) + s(REPNRTCHL_0m_scl, k=3) + s(OMLT_0m_scl, k=3) + (1|VESSEL)),
									formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_0m_scl,k=3) + s(OMLT_0m_scl, k=3) + (1|VESSEL)),
									formula_mix = formula(CPUE ~ 1),   # only used when using a mixture model (deprecated)
									mesh_type = "cutoff", # or "kmeans" "cutoff", or "own" (which modifies to the mackerel case)
									knots_mesh = 250, # default 250
									bias_correct = FALSE								
									)


	# corr_list <- list()
	# for (yr in seq_along(conf$years)) corr_list[[yr]] <- cor(subset(datdat$data, YEAR == conf$years[yr])[,98:106], use = "complete.obs")
	# sapply(corr_list, function(x) eigen(x)$values)
	# par(mfrow=c(4,3), mar=c(0,0,0,0))
	# for (i in 1:11) image(corr_list[[i]])
	# corrmat <- nearPD(corrmat, corr=TRUE)$mat


### Set up data & Define parameters
	datdat = Prepare_data(bio=bio, catch=dat, conf=conf, rerun_stack=FALSE)
	
### Calculating the values of the boxcox transformation
	boxcox_pow = rep(0, length(conf$ages))
	do_boxcox = FALSE 
	if (do_boxcox == TRUE) {
		boxcox_pow = NULL
		for (i in conf$ages){
			selcol <-  which(colnames(datdat$data) == i)
			tempdat <- data.frame( x = datdat$data[which(datdat$data[,selcol]>0),selcol] )
			lm1 <-  lm(x ~ 1, data=tempdat)
			asd <-  boxcox(lm1, lambda = seq(-1,1,by=0.0001))
			boxcox_pow <- c(boxcox_pow,asd$x[which.max(asd$y)])
		} 
	}

	datdat$tmb_data$boxcox_pow <- boxcox_pow


### Fitting the model
	startTime = Sys.time()
	run = fit_mackerelST(data=datdat, conf, conf_extra = NULL)
	endTime = Sys.time()
	endTime-startTime

	View(run$opt$diagnostics)
	(run$opt$max_gradient)
	pl <- as.list(run$sdrep, "Estimate")

		
### Doing plotting of the results 

    sfolder <- paste0(getwd(), "/plots/", conf$save_folder, "/", paste(attributes(conf)$fixed_effect[-1], collapse="_"), "_", attributes(conf)$RE_effects, conf$knots_mesh, "_corr",conf$corr_str, "_", conf$ARorIID)
    if(file.exists(sfolder) == FALSE) dir.create(sfolder)
    
    ## save model output
      ALL <- list()
      ALL$conf <- conf
      ALL$data <- datdat
      ALL$run <- run
      saveRDS(ALL, file=paste0(sfolder, "/run.rds"))

    ## Load the null model 
      sfolder_nullmodel <- paste0(getwd(), "/plots/New_run/nullmodel250")
      ALL_nullmodel <- readRDS(paste0(sfolder_nullmodel, "/run.rds"))
      
          
    ## Model diagnostics
      # AIC
        run$opt$AIC 
      
      # Check if all parameters are estimable  
        fixed_est <- summary(run$sdrep, "fixed")
        random_est <- summary(run$sdrep, "random")
      
      # qqlot (only works for non delta model)
        if (conf$mixture_model != 3) out <- Plot_residuals(run=run, data=datdat, conf=conf, name="QQplot", folder=sfolder)

      # cross validation 
        cross_val <- cv(run=run, data=datdat, k_fold = 10, conf, seed=123)
        saveRDS(cross_val, file=paste0(sfolder, "/cross_val.rds"))
        
      # Self-testing of the model 
        Nsim = 30 # originally run at 10 but should do 30 for final iteration
        ysim <- lapply(seq_len(Nsim), function(x) run$obj$simulate(par= run$obj$env$last.par.best, complete=TRUE))
        dat_fake <- datdat
        vals <- matrix(NA, nrow=Nsim, ncol=length(run$opt$par))
        for (i in 1:Nsim) {
          startTime = Sys.time()
          dat_fake$tmb_data <- ysim[[i]][1:33]
          conf1 <- conf
          conf1$Do_predict = 0
          run1 = fit_mackerelST(data=dat_fake, conf1, conf_extra = NULL)
          vals[i,] = run1$opt$par
          endTime = Sys.time()
          endTime-startTime 
          rm(run1)
          gc()
        }
        saveRDS(vals, file=paste0(sfolder, "/simulation_test.rds"))
        
        jpeg(paste0(sfolder, "/simulation_test.jpeg"), width= 15, height=8, units="in", res=300)
        if(conf$keep_epsilon == FALSE) {par(mfrow=c(2,3), mar=c(2,2,2,2), oma=c(3,2,1,1)); nrows = 6}
        if(conf$keep_epsilon == TRUE) {par(mfrow=c(2,2), mar=c(2,2,2,2), oma=c(3,2,1,1)); nrows = 3}
        for (i in 1:nrows){
          if (i < nrows) {
            boxplot(names=names(run$opt$par)[(i-1)*50+1:50], vals[,(i-1)*50+1:50], las=2, cex.axis=0.9)
            points(1:50,run$opt$par[(i-1)*50+1:50], col="red", pch=16)
          }
          if (i == nrows) {
            boxplot(names=names(run$opt$par)[(i-1)*50+1:(length(run$opt$par)%%50)], vals[,(i-1)*50+1:(length(run$opt$par)%%50)], 
                    las=2, cex.axis=0.9)
            points(1:(length(run$opt$par)%%50),run$opt$par[(i-1)*50+1:(length(run$opt$par)%%50)], col="red", pch=16)
          }
        }
        dev.off()
        
                
      # DHARma-like residual test (can be used for any kind of models) - code from sdmTMB for the one_sample_posterior to 
			# get a MVN sample from the random effects while leaving the fixed effects at their MLE
        Nsim = 500
				one_sample_posterior <- function(object, seed=123) {
					set.seed(seed) 
				  tmp <- object$obj$env$MC(n = 1L, keep = TRUE, antithetic = FALSE)
					re_samp <- as.vector(attr(tmp, "samples"))
					lp <- object$obj$env$last.par.best
					p <- numeric(length(lp))
					fe <- object$obj$env$lfixed()
					re <- object$obj$env$lrandom()
					p[re] <- re_samp
					p[fe] <- lp[fe]
					p
				}

				newpars <- one_sample_posterior(run, seed=123)  # this is to simulate a sample of the random effects from a MVN distribution
        s_nb2 <- sapply(seq_len(Nsim), function(x) run$obj$simulate(par= newpars, complete=FALSE)$yobs, simplify = 'array')
        pred_orig = run$obj$report()$mu
      
     # compare the number of zeros in the data with the number of zeros in the prediction
        grep_val = paste0("^", conf$ages[1], "$")
        col_sel = grep(grep_val, colnames(datdat$data))+c(1:(conf$Nage))-1
        nb_zero <- apply(datdat$data[,col_sel], 2, function(x) sum(x==0, na.rm=T)/length(x))
        pred_zero <- sapply(seq_len(Nsim), function(x) apply(s_nb2[,,x], 2, function(xx) sum(xx==0, na.rm=T)/length(xx)))
        par(mfrow=c(4,3), mar=c(2,2,2,2), oma=c(3,2,1,1))
        for (age in 1:(conf$Nage)){
          hist(pred_zero[age,])
          abline(v=nb_zero[age], col="red")
        }
        nb_zero; apply(pred_zero,1,mean); # quite similar!
        write.table(data.frame(obs=nb_zero, pred=pred_zero), file=paste0(sfolder, "/P_zero_04_2024a.txt"))
        
        for (age in 1:(conf$Nage)){
          to_keep = which(is.na(datdat$data[,col_sel[age]])==FALSE)
          r_nb2 <- DHARMa::createDHARMa(
            simulatedResponse = s_nb2[to_keep,age,],
            observedResponse = datdat$data[to_keep,col_sel[age]],
            fittedPredictedResponse = pred_orig[to_keep,age]
          )
          png(paste0(sfolder, "/resid_age", conf$ages[age], "_04_2024a.png"), width= 10, height=6, units="in", res=300)
          plot(r_nb2)
          dev.off()
          png(paste0(sfolder, "/resid_age", conf$ages[age], "v2_04_2024a.png"), width= 6, height=6, units="in", res=300)
          plotQQunif(simulationOutput = r_nb2)
          dev.off()
          png(paste0(sfolder, "/resid_age", conf$ages[age], "_04_30_2024.png"), width= 6, height=6, units="in", res=300)
          plotQQunif(simulationOutput = r_nb2, testUniformity = T, testOutliers = F,
                     testDispersion = F)
          dev.off()
        }
				
      # We can also use the mcmc approach mentioned by Rufener to perform the QQplot analysis
      # Idea is to take 1 sample of the random effect (using stan or tmb directly) because laplace approx might have failed 
      # so QQplot might look perform even if model is OK (and vice)
        # if stan:
        # library(tmbstan)
        # obj <- TMB::MakeADFun(run$obj$env$data, pl, map = map, DLL = "spatial_SDM_RE")
        # m_stan <- tmbstan(obj, iter = 10000, warmup = 9000, chains = 1)
        # post <- extract_mcmc(m_stan)
        
        # if TMB (we can take 10 random sample and see whether it changes the plot or not)
        # random <- unique(names(run$obj$env$par[run$obj$env$random]))
        # pl <- as.list(run$sdrep, "Estimate")
        # fixed <- !(names(pl) %in% random)
        # map_new <- lapply(pl[fixed], function(x) factor(rep(NA, length(x))))
        # C0 <- solve(run$obj$env$spHess(random=TRUE))    ## Covariance matrix of random effects        
        # Xr0 <- MASS::mvrnorm(10,random_est[,1],C0)       ## Generate one sample of random effects
        # oldpar = run$obj$env$last.par.best
        
        # for (it in 3:10) {
          # oldpar[c(grep("epsilon_st",names(oldpar)),grep("RE",names(oldpar)))] <- Xr0[it,]
          # # pred <- run$obj$report()$mu
          # # pred1 <- run$obj$report(oldpar)$mu
          # out <- Plot_residuals(run=run, data=datdat, parval= oldpar, conf=conf, name=paste0("QQplot_it",it), folder=sfolder, do_all=FALSE)
        # }
        
					
    ## Maps of spatial residuals
       var.select = "SST_0m"
       truncate = FALSE
       spatial_res <- datdat$tmb_data$yobs - run$obj$report(run$obj$env$last.par.best)$mu
       colnames(spatial_res) <- paste0("res",colnames(spatial_res))
       spatial_res <- cbind(datdat$data, spatial_res)
       to_remove <- as.numeric(which(is.na(datdat$tmb_data$yobs[,1]) == TRUE))
       for (i in seq_along(conf$ages)){
         ageval = conf$ages[i]
         model_output <- list()
         model_output$response <- datdat$tmb_data$yobs[-to_remove,i]
         model_output$mu <- run$obj$report(run$obj$env$last.par.best)$mu[-to_remove,i]
         model_output$family <- "tweedie"
         if (!is.null(conf$conf_model$ln_phi)) model_output$ln_phi <- run$opt$par[grep("ln_phi", names(run$opt$par))][as.numeric(conf$conf_model$ln_phi)[i]]
         if (is.null(conf$conf_model$ln_phi)) model_output$ln_phi <- run$opt$par[grep("ln_phi", names(run$opt$par))][i]
         if (!is.null(conf$conf_model$thetaf)) model_output$thetaf <- run$opt$par[grep("thetaf", names(run$opt$par))][as.numeric(conf$conf_model$thetaf)[i]]
         if (is.null(conf$conf_model$thetaf)) model_output$thetaf <- run$opt$par[grep("thetaf", names(run$opt$par))][i]
         model_output$obj <- run$obj
         model_output$tmb_data <- datdat$tmb_data
         
         d <- datdat$data[-to_remove,]
       
         d$residuals <- residuals.sdmTMB(model_output)
         spatial_res$resp = spatial_res[,grep(paste0("res",ageval), colnames(spatial_res))[1]]

         pred_dat = run$dens_pred %>% filter(AGE == 3)
         pred_dat$resp <- pred_dat[,eval(parse(text=var.select))]
         if (var.select == "SST_0m") pred_dat$resp <- pred_dat$resp-273.15
         pred_dat <- pred_dat %>% mutate(resp_trunc = ifelse(resp > quantile(resp, 0.99), quantile(resp, 0.99), ifelse(resp < quantile(resp, 0.01), quantile(resp, 0.01), resp)))
         if (truncate==TRUE) pred_dat$resp_val = pred_dat$resp_trunc
         if (truncate==FALSE) pred_dat$resp_val = pred_dat$resp
          
         bad_temp <- c(7.5, 7.2, 6.2, 5.6, 5.2, rep(4.4, 3))
         for (yr in 2010:2020){
           pred_dat1 = pred_dat %>% filter(YEAR == yr) %>% mutate(resp_val_bin = ifelse(resp_val < bad_temp[i], paste0("<", bad_temp[i]), paste0("≥", bad_temp[i])))
           d1 = d %>% filter(YEAR == yr)
           p1 <- ggplot(Atlantic_proj) + geom_sf(color ="grey27", size = .2) +
             geom_raster(data = pred_dat1, aes(x=X1000, y=Y1000, fill = factor(resp_val_bin))) +
             xlab("") + ylab("") + 
             coord_sf(xlim = c(1000,5000), ylim = c(4000, 6200),
                      datum = projection_km, expand = FALSE) + theme_bw() + 
             # facet_wrap(~YEAR) + 
             scale_fill_manual(values = c("lightyellow","white"), name = "Temperature") + 
             new_scale_color() + 
             geom_point(data=d1, aes(x=X1000, y=Y1000, col = residuals), size=1) +
             geom_point(data=spatial_res %>% filter(is.na(resp)), aes(x=X1000, y=Y1000), col = "black", size=0.7, shape=4) +
             scale_color_gradient2(low="blue", mid="snow", high="red2", midpoint=0) +
             ggtitle(paste0("Residual distribution")) + 
             theme(axis.title = element_text(size=15),
                   axis.text = element_text(size=12), 
                   axis.text.x = element_text(angle = 90), 
                   plot.title = element_text(hjust=0.5, size=18),
                   strip.text = element_text(size=14))
           ggsave(p1, filename = paste0(sfolder, "/Residuals_age", ageval, "year", yr, ".jpeg"), dpi ="retina", width = 12, height = 8, device = "jpeg")
          }
       }
       
       
    ## Maps of spatial distribution of residuals along with SST (the most important variable)
       var.select = "SST_0m"
       truncate = FALSE
       spatial_res <- datdat$tmb_data$yobs - run$obj$report(run$obj$env$last.par.best)$mu
       colnames(spatial_res) <- paste0("res",colnames(spatial_res))
       spatial_res <- cbind(datdat$data, spatial_res)
       to_remove <- as.numeric(which(is.na(datdat$tmb_data$yobs[,1]) == TRUE))
     
       for (i in seq_along(conf$ages)){
         ageval = conf$ages[i]
         model_output <- list()
         model_output$response <- datdat$tmb_data$yobs[-to_remove,i]
         model_output$mu <- run$obj$report(run$obj$env$last.par.best)$mu[-to_remove,i]
         model_output$family <- "tweedie"
         if (!is.null(conf$conf_model$ln_phi)) model_output$ln_phi <- run$opt$par[grep("ln_phi", names(run$opt$par))][as.numeric(conf$conf_model$ln_phi)[i]]
         if (is.null(conf$conf_model$ln_phi)) model_output$ln_phi <- run$opt$par[grep("ln_phi", names(run$opt$par))][i]
         if (!is.null(conf$conf_model$thetaf)) model_output$thetaf <- run$opt$par[grep("thetaf", names(run$opt$par))][as.numeric(conf$conf_model$thetaf)[i]]
         if (is.null(conf$conf_model$thetaf)) model_output$thetaf <- run$opt$par[grep("thetaf", names(run$opt$par))][i]
         model_output$obj <- run$obj
         model_output$tmb_data <- datdat$tmb_data
         
         d <- datdat$data[-to_remove,]
         
         d$residuals <- residuals.sdmTMB(model_output)
         spatial_res$resp = spatial_res[,grep(paste0("res",ageval), colnames(spatial_res))[1]]

         pred_dat = run$dens_pred %>% mutate(AGE = AGE + min(conf$ages) - min(AGE))
         pred_dat = pred_dat %>% filter(AGE == ageval)
         pred_dat$resp <- pred_dat[,eval(parse(text=var.select))]
         if (var.select == "SST_0m") pred_dat$resp <- pred_dat$resp-273.15
         pred_dat <- pred_dat %>% mutate(resp_trunc = ifelse(resp > quantile(resp, 0.99), quantile(resp, 0.99), ifelse(resp < quantile(resp, 0.01), quantile(resp, 0.01), resp)))
         if (truncate==TRUE) pred_dat$resp_val = pred_dat$resp_trunc
         if (truncate==FALSE) pred_dat$resp_val = pred_dat$resp
          
         bad_temp <- c(7.5, 7.2, 6.2, 5.6, 5.2, rep(4.4, 3))
         for (yr in 2010:2020){
           pred_dat1 = pred_dat %>% filter(YEAR == yr) %>% mutate(resp_val_bin = ifelse(resp_val < bad_temp[i], paste0("<", bad_temp[i]), paste0("≥", bad_temp[i])))
           d1 = d %>% filter(YEAR == yr)
           p1 <- ggplot(Atlantic_proj) + geom_sf(color ="grey27", size = .2) +
             geom_raster(data = pred_dat1, aes(x=X1000, y=Y1000, fill = factor(resp_val_bin))) +
             xlab("") + ylab("") + 
             coord_sf(xlim = c(1000,5000), ylim = c(4000, 6200),
                      datum = projection_km, expand = FALSE) + theme_bw() + 
             # facet_wrap(~YEAR) + 
             scale_fill_manual(values = c("lightyellow","white"), name = "Temperature") + 
             new_scale_color() + 
             geom_point(data=d1, aes(x=X1000, y=Y1000, col = residuals), size=1) +
             geom_point(data=spatial_res %>% filter(is.na(resp)), aes(x=X1000, y=Y1000), col = "black", size=0.7, shape=4) +
             scale_color_gradient2(low="blue", mid="snow", high="red2", midpoint=0) +
             ggtitle(paste0("Residual distribution")) + 
             theme(axis.title = element_text(size=15),
                   axis.text = element_text(size=12), 
                   axis.text.x = element_text(angle = 90), 
                   plot.title = element_text(hjust=0.5, size=18),
                   strip.text = element_text(size=14))
           ggsave(p1, filename = paste0(sfolder, "/Residuals_age", ageval, "year", yr, ".jpeg"), dpi ="retina", width = 12, height = 8, device = "jpeg")
      
           p2 <- ggplot(Atlantic_proj) + geom_sf(color ="grey27", size = .2) +
             geom_raster(data = pred_dat1, aes(x=X1000, y=Y1000, fill = factor(resp_val_bin))) +
             scale_fill_manual(values = c("black","NA"), name = "Temperature") + 
             new_scale_fill() + 
             geom_raster(data = pred_dat1 %>% mutate(logPred_trunc = ifelse(logPred<1, 1, logPred)), aes(x=X1000, y=Y1000, fill = logPred_trunc), alpha=0.8) +
             scale_fill_viridis_c(option = "plasma") +
             xlab("") + ylab("") + 
             coord_sf(xlim = c(1000,5000), ylim = c(4000, 6200),
                      datum = projection_km, expand = FALSE) + theme_bw() + 
            ggtitle(paste0("Residual distribution")) + 
             theme(axis.title = element_text(size=15),
                   axis.text = element_text(size=12), 
                   axis.text.x = element_text(angle = 90), 
                   plot.title = element_text(hjust=0.5, size=18),
                   strip.text = element_text(size=14))
           ggsave(p2, filename = paste0(sfolder, "/predict_lowtemp", ageval, "year", yr, ".jpeg"), dpi ="retina", width = 12, height = 8, device = "jpeg")
          }
       }
       
       
    ## Maps of spatial distribution of observations along with SST (the most important variable)
       var.select = "SST_0m"
       truncate = FALSE
       to_remove <- as.numeric(which(is.na(datdat$tmb_data$yobs[,1]) == TRUE))
     
       for (i in seq_along(conf$ages)){
         ageval = conf$ages[i]
         d <- datdat$data[-to_remove,]
         d$obs <- d[,"3"]
         
         
         pred_dat = run$dens_pred %>% mutate(AGE = AGE + min(conf$ages) - min(AGE))
         pred_dat = pred_dat %>% filter(AGE == ageval)
         pred_dat$resp <- pred_dat[,eval(parse(text=var.select))]
         if (var.select == "SST_0m") pred_dat$resp <- pred_dat$resp-273.15
         pred_dat <- pred_dat %>% mutate(resp_trunc = ifelse(resp > quantile(resp, 0.99), quantile(resp, 0.99), ifelse(resp < quantile(resp, 0.01), quantile(resp, 0.01), resp)))
         if (truncate==TRUE) pred_dat$resp_val = pred_dat$resp_trunc
         if (truncate==FALSE) pred_dat$resp_val = pred_dat$resp
          
         bad_temp <- c(7.5, 7.2, 6.2, 5.6, 5.2, rep(4.4, 3))
         for (yr in 2010:2020){
           pred_dat1 = pred_dat %>% filter(YEAR == yr) %>% mutate(resp_val_bin = ifelse(resp_val < bad_temp[i], paste0("<", bad_temp[i]), paste0("≥", bad_temp[i])))
           d1 = d %>% filter(YEAR == yr)
           p1 <- ggplot(Atlantic_proj) + geom_sf(color ="grey27", size = .2) +
             geom_raster(data = pred_dat1, aes(x=X1000, y=Y1000, fill = factor(resp_val_bin))) +
             xlab("") + ylab("") + 
             coord_sf(xlim = c(1000,5000), ylim = c(4000, 6200),
                      datum = projection_km, expand = FALSE) + theme_bw() + 
             # facet_wrap(~YEAR) + 
             scale_fill_manual(values = c("lightyellow","white"), name = "Temperature") + 
             geom_point(data=d1, aes(x=X1000, y=Y1000, col = obs<1), size=1) +
             # scale_color_gradient2(low="blue", mid="snow", high="red2", midpoint=1) +
             ggtitle(paste0("Residual distribution")) + 
             theme(axis.title = element_text(size=15),
                   axis.text = element_text(size=12), 
                   axis.text.x = element_text(angle = 90), 
                   plot.title = element_text(hjust=0.5, size=18),
                   strip.text = element_text(size=14))
           ggsave(p1, filename = paste0(sfolder, "/Residuals_age", ageval, "year", yr, ".jpeg"), dpi ="retina", width = 12, height = 8, device = "jpeg")
      
           p2 <- ggplot(Atlantic_proj) + geom_sf(color ="grey27", size = .2) +
             geom_raster(data = pred_dat1, aes(x=X1000, y=Y1000, fill = factor(resp_val_bin))) +
             scale_fill_manual(values = c("black","NA"), name = "Temperature") + 
             new_scale_fill() + 
             geom_raster(data = pred_dat1 %>% mutate(logPred_trunc = ifelse(logPred<1, 1, logPred)), aes(x=X1000, y=Y1000, fill = logPred_trunc), alpha=0.8) +
             scale_fill_viridis_c(option = "plasma") +
             xlab("") + ylab("") + 
             coord_sf(xlim = c(1000,5000), ylim = c(4000, 6200),
                      datum = projection_km, expand = FALSE) + theme_bw() + 
            ggtitle(paste0("Residual distribution")) + 
             theme(axis.title = element_text(size=15),
                   axis.text = element_text(size=12), 
                   axis.text.x = element_text(angle = 90), 
                   plot.title = element_text(hjust=0.5, size=18),
                   strip.text = element_text(size=14))
           ggsave(p2, filename = paste0(sfolder, "/predict_lowtemp", ageval, "year", yr, ".jpeg"), dpi ="retina", width = 12, height = 8, device = "jpeg")
          }
       }
       
       
    ## Calculating the spatial range parameter
      sqrt(8)/exp(run$opt$par[grep("Kappa", names(run$opt$par))])
   
    ## Looking at the covariance matrix
      covar <- run$obj$report()$covar
      
    ## Calculating the conditional R2 (in link space) for each age group
        fixed <- run$obj$report(run$obj$env$last.par.best)$fixed_noRE
        spatiotemporal <- run$obj$report(run$obj$env$last.par.best)$epsilon_st_A_mat
        nullmodel_method = "obs"  # "obs" or "Nakagawa"
        
        R2_cond <- c()
        fixed_contribution <- c()
        ST_contribution <- c()
        Vessel_effect <- c()
        for (i in seq_along(conf$ages)){
          age = conf$ages[i]
          var_fixed = var(fixed[,i])
          var_stRE = 0
          if (!is.null(spatiotemporal)) var_stRE = var(spatiotemporal[,i])
          if (length(grep("tau_G", names(run$opt$par))) ==  length(conf$ages)) var_random = (exp(run$opt$par[grep("tau_G", names(run$opt$par))])[i])^2
          if (length(grep("tau_G", names(run$opt$par))) == 1) var_random = (exp(run$opt$par[grep("tau_G", names(run$opt$par))])[1])^2
          if (nullmodel_method == "obs") mu_null = mean(datdat$data[,grep(paste0("^",conf$ages[i],"$"), colnames(datdat$data))], na.rm=T)  # this is the fixed effect of the null model (just exp(intercept))
          if (nullmodel_method == "Nakagawa") mu_null = exp(ALL_nullmodel$run$opt$par[grep("beta", names(ALL_nullmodel$run$opt$par))][i] + 0.5*var_random)
          s1 = boot::inv.logit(run$opt$par[grep("thetaf", names(run$opt$par))[i]]) + 1
          phi = exp(run$opt$par[grep("ln_phi", names(run$opt$par))[i]])
          var_obs = mu_null^(s1-2)* phi #- 1/4*phi^2*mu_null^(2*s1-4))
          R2_cond[i] <- (var_fixed + var_stRE + var_random)/(var_fixed + var_stRE + var_random + var_obs)
          fixed_contribution[i] <- (var_fixed)/(var_fixed + var_stRE + var_random + var_obs)
          ST_contribution[i] <- (var_stRE)/(var_fixed + var_stRE + var_random + var_obs)
          Vessel_effect[i] <- (var_random)/(var_fixed + var_stRE + var_random + var_obs)
        }
        R2_cond
        fixed_contribution/R2_cond
        Vessel_effect/R2_cond
        
        write.table(data.frame(R2=R2_cond, fixed_contribution=fixed_contribution/R2_cond), file=paste0(sfolder, "/R2_contribution.txt"))
        
    ## Calculating the conditional R2 (in link space) for each age group
    ## But this time, we divide by spatial strata
        fixed <- run$obj$report(run$obj$env$last.par.best)$fixed_noRE
        beta_est <- run$obj$report(run$obj$env$last.par.best)$beta
        X1 <- datdat$tmb_data$X; X1[, which((1:ncol(X1) %in% grep("YEAR", colnames(X1))) == FALSE)] <- 0; Year_est <- X1 %*% beta_est
        X1 <- datdat$tmb_data$X; X1[, which((1:ncol(X1) %in% grep("SST", colnames(X1))) == FALSE)] <- 0; SST_est <- X1 %*% beta_est
        X1 <- datdat$tmb_data$X; X1[, which((1:ncol(X1) %in% grep("OMLT", colnames(X1))) == FALSE)] <- 0; OMLT_est <- X1 %*% beta_est
        spatiotemporal <- run$obj$report(run$obj$env$last.par.best)$epsilon_st_A_mat
        nullmodel_method = "obs"  # "obs" or "Nakagawa"
        
        bla <- st_distance(polygons_df, st_as_sf(x=datdat$data, coords = c("LON", "LAT")) %>% st_set_crs(4326) %>% st_transform(projection_km) %>% st_cast("POINT"))
        datdat$data$Strata = apply(bla, 2, which.min)
        strata_available = polygons_df %>% st_drop_geometry() %>% unlist() %>% as.numeric()
        datdat$data$Strata = strata_available[datdat$data$Strata]
          
        R2_cond <- matrix(NA, 11, conf$Nage)   # R2 per age and strata
        fixed_contribution <- matrix(NA, 11, conf$Nage)
        ST_contribution <- matrix(NA, 11, conf$Nage)
        fixed_SST_contribution <- matrix(NA, 11, conf$Nage)
        fixed_Year_contribution <- matrix(NA, 11, conf$Nage)
        fixed_OMLT_contribution <- matrix(NA, 11, conf$Nage)
        Vessel_effect <- matrix(NA, 11, conf$Nage)
        for (i in seq_along(conf$ages)){
          for (stratum in 1:11){
            selected_row = which(datdat$data$Strata == strata_available[stratum])
            age = conf$ages[i]
            var_fixed = var(fixed[selected_row,i])
            var_fixed_year = var(Year_est[selected_row,i])
            var_fixed_SST = var(SST_est[selected_row,i])
            var_fixed_OMLT = var(OMLT_est[selected_row,i])
            var_stRE = 0
            if (!is.null(spatiotemporal)) var_stRE = var(spatiotemporal[selected_row,i])
            if (length(grep("tau_G", names(run$opt$par))) ==  length(conf$ages)) var_random = (exp(run$opt$par[grep("tau_G", names(run$opt$par))])[i])^2
            if (length(grep("tau_G", names(run$opt$par))) == 1) var_random = (exp(run$opt$par[grep("tau_G", names(run$opt$par))])[1])^2
            if (nullmodel_method == "obs") mu_null = mean(datdat$data[selected_row,grep(paste0("^",conf$ages[i],"$"), colnames(datdat$data))], na.rm=T)  # this is the fixed effect of the null model (just exp(intercept))
            if (nullmodel_method == "Nakagawa") mu_null = exp(ALL_nullmodel$run$opt$par[grep("beta", names(ALL_nullmodel$run$opt$par))][i] + 0.5*var_random)
            s1 = boot::inv.logit(run$opt$par[grep("thetaf", names(run$opt$par))[i]]) + 1
            phi = exp(run$opt$par[grep("ln_phi", names(run$opt$par))[i]])
            var_obs = mu_null^(s1-2)* phi #- 1/4*phi^2*mu_null^(2*s1-4))
            R2_cond[stratum,i] <- (var_fixed + var_stRE + var_random)/(var_fixed + var_stRE + var_random + var_obs)
            fixed_contribution[stratum,i] <- (var_fixed)/(var_fixed + var_stRE + var_random + var_obs)
            fixed_SST_contribution[stratum,i] <- (var_fixed_SST)/(var_fixed_SST + var_fixed_year + var_fixed_OMLT + var_stRE + var_random + var_obs)
            fixed_Year_contribution[stratum,i] <- (var_fixed_year)/(var_fixed_SST + var_fixed_year + var_fixed_OMLT + var_stRE + var_random + var_obs)
            fixed_OMLT_contribution[stratum,i] <- (var_fixed_OMLT)/(var_fixed_SST + var_fixed_year + var_fixed_OMLT + var_stRE + var_random + var_obs)
            ST_contribution[stratum,i] <- (var_stRE)/(var_fixed + var_stRE + var_random + var_obs)
            Vessel_effect[stratum,i] <- (var_random)/(var_fixed + var_stRE + var_random + var_obs)
          }
        }
        rownames(R2_cond) = strata_available; colnames(R2_cond) = conf$ages; R2_cond
        fixed_perc <- fixed_contribution/R2_cond
        rownames(fixed_perc) = strata_available; colnames(fixed_perc) = conf$ages; fixed_perc
        Vessel_effect/R2_cond
        rownames(fixed_contribution) = strata_available; colnames(fixed_contribution) = conf$ages; fixed_contribution
        rownames(fixed_SST_contribution) = strata_available; colnames(fixed_SST_contribution) = conf$ages; fixed_SST_contribution
        rownames(fixed_Year_contribution) = strata_available; colnames(fixed_Year_contribution) = conf$ages; fixed_Year_contribution
        rownames(fixed_OMLT_contribution) = strata_available; colnames(fixed_OMLT_contribution) = conf$ages; fixed_OMLT_contribution
        Vessel_effect/R2_cond
        
        datdat$data %>% group_by(YEAR, Strata) %>% mutate(obs = `3`) %>% summarize(zeros = length(obs == 0)) %>% View()
        
        write.table(data.frame(R2=R2_cond, 
                               fixed_contribution=fixed_contribution/R2_cond,
                               fixed_SST_contribution=fixed_SST_contribution/R2_cond,
                               fixed_Year_contribution=fixed_Year_contribution/R2_cond,
                               fixed_OMLT_contribution=fixed_OMLT_contribution/R2_cond
                               ), file=paste0(sfolder, "/R2_contribution_bystrata.txt"))


  ## Case when we want to reactivate prediction
      if (load_result == TRUE){
        sfolder <- paste0(getwd(), "/plots/New_run/SST_0m_scl_OMLT_0m_scl_VESSELnewmap/")
        # sfolder <- paste0(getwd(), "/plots/New_run/SST_0m_scl_OMLT_0m_scl_VESSEL250/")
        # sfolder <- paste0(getwd(), "/plots/New_run/OMLT_0m_scl_VESSELnewmap/")
        ALL <- readRDS(paste0(sfolder, "/run.rds"))
        conf <- ALL$conf
        datdat <- ALL$data
        run <- ALL$run   
        version <- paste0(getwd(), "/src/spatial_SDM_RE")
        dyn.load(dynlib(version))
        run$obj$retape()
        run$obj$fn(run$opt$par)
      }
        
    ## Producing the marginal effect plots (done within the model for more flexibility & complex covaraite structure)
      ## Marginal SST effect 
      SST_dat <- data.frame(expand.grid(SST_0m_scl = seq(min(datdat$data$SST_0m_scl),max(datdat$data$SST_0m_scl),length.out=50),
			                                  REPNRTCHL_0m_scl = 0, OMLT_0m_scl = 0, Area = 2,  
																				YEAR = 2019, X1000 = mean(datdat$data$X1000), Y1000=mean(datdat$data$Y1000), CPUE=1))
			Marginal_SST <- refit(run=run, data=datdat, pred_data=SST_dat, predict_type = 2, conf=conf) 
			SST_marg <- as.data.frame(Marginal_SST$dens_pred)
			colnames(SST_marg) <- paste0("Age",conf$ages)
			SST_dat <- cbind(SST_dat, SST_marg)
			SST_dat$SST_0m <- SST_dat$SST_0m_scl * sd(dat$SST_0m) + mean(dat$SST_0m) - 273.15 
			SST_dat_df <- pivot_longer(SST_dat, cols =  starts_with("Age"), names_to="Age")
			fake = datdat$data
			fake$value = 1
			fake$Age = "Age3"
			fake$Age_cont = "3"
			fake$Age_bin = factor("0")
			fake$SST_0m = fake$SST_0m - 273.15 
			SST_dat_df <- SST_dat_df %>% mutate(Age_bin = letters[(as.numeric(substring(Age, 4)) %% 3) + 1], 
			                                    Age_cont = as.numeric(substring(Age, 4)), 
			                                    Age = as.factor(Age), Age = factor(Age, levels=paste0("Age",conf$ages))
			                                    )
			p1a <- ggplot(SST_dat_df, aes(x=SST_0m, y = value, col=Age, linetype=Age)) + 
			  geom_line(linewidth=1.3) + 
				scale_color_viridis_d() + theme_bw() + labs(y="Marginal effect (link scale)", x="SST (°C)") + 
				geom_rug(data=fake, col="black")
      ggsave(p1a, filename = paste0(sfolder, "/Marginal_SST_filter.jpeg"), dpi ="retina", width = 6, height = 4, device = "jpeg")
 
		## Marginal CHL effect 
			CHL_dat <- data.frame(expand.grid(SST_0m_scl = 0,
																				REPNRTCHL_0m_scl = seq(min(datdat$data$REPNRTCHL_0m_scl),max(datdat$data$REPNRTCHL_0m_scl),length.out=50),
																				YEAR = 2019, OMLT_0m_scl = 0, Area = 2,
																				X1000 = mean(datdat$data$X1000), Y1000=mean(datdat$data$Y1000), CPUE=1))
			Marginal_CHL <- refit(run=run, data=datdat, pred_data=CHL_dat, predict_type = 2, conf=conf) 
			CHL_marg <- as.data.frame(Marginal_CHL$dens_pred)
			colnames(CHL_marg) <- paste0("Age",conf$ages)
			CHL_dat <- cbind(CHL_dat, CHL_marg)
			CHL_dat$REPNRTCHL_0m <- CHL_dat$REPNRTCHL_0m_scl * sd(dat$REPNRTCHL_0m) + mean(dat$REPNRTCHL_0m)
			CHL_dat_df <- pivot_longer(CHL_dat, cols =  starts_with("Age"), names_to="Age")
			fake = datdat$data
			fake$value = 1
			fake$Age = "Age3"
			CHL_dat_df <- CHL_dat_df %>% mutate(Age_bin = letters[(as.numeric(substring(Age, 4)) %% 3) + 1], 
			                                       Age_cont = as.numeric(substring(Age, 4)), 
			                                       Age = as.factor(Age), Age = factor(Age, levels=paste0("Age",conf$ages))
			)
			p1b <- ggplot(CHL_dat_df, aes(x=REPNRTCHL_0m, y = value, col=Age, linetype=Age)) + 
			  geom_line(linewidth=1.3) +
				scale_color_viridis_d() + theme_bw() + labs(y="Marginal effect (link scale)", x="Chlorophyl concentration (mg.m-3)") + 
				geom_rug(data=fake, col="black")
      ggsave(p1b, filename = paste0(sfolder, "/Marginal_CHL_filter.jpeg"), dpi ="retina", width = 6, height = 4, device = "jpeg")
      
		## Marginal OMLT effect 
			OMLT_dat <- data.frame(expand.grid(SST_0m_scl = 0,
																				REPNRTCHL_0m_scl = 0,
																				OMLT_0m_scl = seq(min(datdat$data$OMLT_0m_scl),max(datdat$data$OMLT_0m_scl),length.out=50),
																				YEAR = 2019, Area = 2,
																				X1000 = mean(datdat$data$X1000), Y1000=mean(datdat$data$Y1000), CPUE=1))
			
			Marginal_OMLT <- refit(run=run, data=datdat, pred_data=OMLT_dat, predict_type = 2, conf=conf) 
			OMLT_marg <- as.data.frame(Marginal_OMLT$dens_pred)
			colnames(OMLT_marg) <- paste0("Age",conf$ages)
			OMLT_dat <- cbind(OMLT_dat, OMLT_marg)
			OMLT_dat$OMLT_0m <- OMLT_dat$OMLT_0m_scl * sd(dat$OMLT_0m, na.rm=T) + mean(dat$OMLT_0m, na.rm=T)
			OMLT_dat_df <- pivot_longer(OMLT_dat, cols =  starts_with("Age"), names_to="Age")
			fake = datdat$data
			fake$value = 1
			fake$Age = "Age3"
			OMLT_dat_df <- OMLT_dat_df %>% mutate(Age_bin = letters[(as.numeric(substring(Age, 4)) %% 3) + 1], 
			                                      Age_cont = as.numeric(substring(Age, 4)), 
			                                      Age = as.factor(Age), Age = factor(Age, levels=paste0("Age",conf$ages))
			)
			p1c <- ggplot(OMLT_dat_df , aes(x=OMLT_0m, y = value, col=Age, linetype=Age)) + 
			  geom_line(linewidth=1.3) +
				scale_color_viridis_d() + theme_bw() + labs(y="Marginal effect (link scale)", x="Mixed layer depth (m)")  
			p1c1 <- p1c +	geom_rug(data=fake, col="black")
      ggsave(p1c, filename = paste0(sfolder, "/Marginal_OMLT_filter.jpeg"), dpi ="retina", width = 6, height = 4, device = "jpeg")
      
      ## Combine the plots
      p_legend <- get_legend(p1c + theme(legend.title = element_text(size = 15),legend.text = element_text(size = 14)), position = NULL) %>% as_ggplot()
      
      p1a = p1a + ggtitle("(a)") + theme(legend.position="none")
      p1c1 = p1c1 + ggtitle("(b)") + theme(legend.position="none")
      p1 <- grid.arrange(p1a,p1c1,p_legend, nrow=1, layout_matrix=matrix(c(1,1,1,1,2,2,2,2,3),nrow=1,byrow=T))
      ggsave(p1, filename = paste0(sfolder, "/Marginal_all_filter.jpeg"), dpi ="retina", width = 10, height = 4, device = "jpeg")
        
		## SST effect with OMLT
			SST_dat <- data.frame(expand.grid(SST_0m_scl = seq(min(datdat$data$SST_0m_scl),max(datdat$data$SST_0m_scl),length.out=50),
																				OMLT_0m_scl = seq(min(datdat$data$OMLT_0m_scl),max(datdat$data$OMLT_0m_scl),length.out=50),
                                          REPNRTCHL_0m_scl = 0,
																				YEAR = 2020, X1000 = mean(datdat$data$X1000), Y1000=mean(datdat$data$Y1000), CPUE=1))
			Marginal_SST_OMLT <- refit(run=run, data=datdat, pred_data=SST_dat, predict_type = 2, conf=conf)
			SST_OMLT_marg <- as.data.frame(Marginal_SST_OMLT$dens_pred)
			colnames(SST_OMLT_marg) <- paste0("Age",conf$ages)
			SST_OMLT_dat <- cbind(SST_dat, SST_OMLT_marg)
			SST_OMLT_dat$OMLT_0m <- SST_OMLT_dat$OMLT_0m_scl * sd(dat$OMLT_0m, na.rm=T) + mean(dat$OMLT_0m, na.rm=T)
			SST_OMLT_dat$SST_0m <- SST_OMLT_dat$SST_0m_scl * sd(dat$SST_0m) + mean(dat$SST_0m) - 273.15 
			SST_OMLT_dat_df <- pivot_longer(SST_OMLT_dat, cols =  starts_with("Age"), names_to="Age")
			fake = datdat$data
			fake <- fake %>% dplyr::select(SST_0m, OMLT_0m)
			fake$value = 1
			fake$SST_0m = fake$SST_0m - 273.15 
			fake <- do.call("rbind", replicate(8, fake, simplify=FALSE))
        fake$Age <- rep(paste0("Age",3:10), each = nrow(datdat$data))
        fake$Age <- as.factor(fake$Age)
        fake$Age <- factor(fake$Age, levels=paste0("Age",3:10))
			SST_OMLT_dat_df <- SST_OMLT_dat_df %>% mutate(Age = as.factor(Age), Age = factor(Age, levels=paste0("Age",1:10)))
			p1d <- ggplot(SST_OMLT_dat_df %>% filter(as.numeric(Age) > 2), aes(x=SST_0m, y = OMLT_0m, col=value, fill=value)) + 
			  geom_tile() + facet_wrap(. ~ Age) + scale_color_viridis_c(option="viridis", name="Effect") +
			  theme_bw() + scale_fill_viridis_c(option="viridis", name="Effect") +
			  geom_point(data= fake, aes(x=SST_0m, y = OMLT_0m), col=grey(0.3, alpha=0.8), fill=NA, cex=0.01, alpha = 0.2) +
			  labs(y="Mixed layer depth (m)", x="SST (°C)")
			ggsave(p1d, filename = paste0(sfolder, "/Marginal_SST_OMLT_filter.jpeg"), dpi ="retina", width = 8, height = 6, device = "jpeg")
			p1 <- grid.arrange(p1d,p1b, nrow=2, layout_matrix=matrix(c(1,1,2), ncol=1))
			ggsave(p1, filename = paste0(sfolder, "/Marginal_all_filter_bivar.jpeg"), dpi ="retina", width = 6.5, height = 8, device = "jpeg")
			p1d <- ggplot(SST_OMLT_dat_df %>% filter(as.numeric(Age) > 2), aes(x=SST_0m, y = OMLT_0m, col=value, fill=value)) + 
			  geom_tile() + facet_wrap(. ~ Age) + scale_color_viridis_c(option="viridis", name="Effect") +
			  theme_bw() + scale_fill_viridis_c(option="viridis", name="Effect") + 
			  labs(y="Mixed layer depth (m)", x="SST (°C)")
			ggsave(p1d, filename = paste0(sfolder, "/Marginal_SST_OMLT_filter_nodots.jpeg"), dpi ="retina", width = 8, height = 6, device = "jpeg")
			p1 <- grid.arrange(p1d,p1b, nrow=2, layout_matrix=matrix(c(1,1,2), ncol=1))
			ggsave(p1, filename = paste0(sfolder, "/Marginal_all_filter_bivar_nodots.jpeg"), dpi ="retina", width = 6.5, height = 8, device = "jpeg")
			
		
       
    ## Calculating the adreport variables: activating the prediction (if not already activated) 
      Mod_pred= run
			
    ## Index of abundance
      IA_IESSNS <- readRDS(paste0(getwd(), "/data/new_data/WGWIDE_IESSN_Biomass_billionTon.rds"))
      IA_IESSNS[1,9]<- 0.01
      ttt <- as.matrix(IA_IESSNS[conf$ages,conf$years-2007])
      qqq <- reshape2::melt(ttt)
      colnames(qqq) <- c("Age", "Year", "IA_IESSNS")
      qqq <- qqq %>% mutate(Year = as.numeric(Year), Age = as.numeric(as.character(factor(Age, labels=conf$ages))), IA_IESSNS = as.numeric(IA_IESSNS))
      
      IA2 <- as.data.frame(summary(Mod_pred$sdrep, "report"))
      if (ncol(IA2) ==2 ) colnames(IA2) <- c("IA", "IA_sd")
      if (ncol(IA2) ==4 ) colnames(IA2) <- c("IA_old", "IA_sd_old", "IA", "IA_sd")
      IA2$Year <- rep(1:conf$Nyear, conf$Nage)
      IA2$Age <- rep(conf$ages, each=conf$Nyear)
      IA2$IA_type <- rep(c("IA", "IA_log", "IA_corrected", "IA_log_corrected"), each=nrow(IA2)/4)
      write.table(IA2, file=paste0(sfolder, "/IA.txt"))
      
      IAs <- IA2 %>% mutate(Year = 2009+Year) %>% filter(IA_type == "IA_log_corrected")
      write.table(IAs, file=paste0(sfolder, "/Index_log.txt"))
      
      p12 <- ggplot(IAs %>% filter(IA_type == "IA_corrected"), aes(x = Year, y = IA)) + theme_bw() + facet_wrap(.~Age, scale = "free_y") +
        geom_ribbon(aes(ymin = IA - 1.96*IA_sd, ymax=IA + 1.96*IA_sd), fill=gray(0.7, 0.5)) + geom_line(size=1.2) 
      ggsave(p12, filename = paste0(sfolder, "IA_model.png"), dpi ="retina", width = 8, height = 7, device = "png")
      
      
      plot_IA(mod = Mod_pred, conf = conf, IA_assessment = qqq, name = "IAnew", folder = sfolder)
      CG_plots(run = Mod_pred, name="IA", data=datdat, no_fish_zone=FALSE, folder = sfolder, conf=conf)
      CG_plots(run = Mod_pred, name="IA", data=datdat, no_fish_zone=TRUE, folder = sfolder, conf=conf)
      
      Plot_maps(run=Mod_pred, data=datdat, name="Best", conf=conf, no_fish_zone=FALSE, folder=sfolder)
 
      ## Plotting the env covariates
      Plot_env(run = Mod_pred, var = "SST_0m", truncate=TRUE, folder=sfolder)
      Plot_env(run = Mod_pred, var = "REPNRTCHL_0m", truncate=TRUE, folder=sfolder)
      Plot_env(run = Mod_pred, var = "OMLT_0m", truncate=TRUE, folder=sfolder)
     
      
      
 
      