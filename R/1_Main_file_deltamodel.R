#### 
# Include cohort-varying env effect: should be expressed in terms of random effect (to save some degree of freedom)

rm(list=ls())

load(".Rdata")
### libraries

library(tidyverse)
library(Matrix)
library(sf)
library(sfheaders)
library(raster)
library(sp)
library(assertthat)
library(gridExtra)
library(sdmTMB)
library(leaflet)
library(leafsync)
library(ggplot2)
library(ggnewscale)
library(mapview)
library(ggpmisc)
library(ggrepel)
library(TMB)
library(TMBhelper)
library(GGally)
library(DHARMa)



### Setting projection method and geographical extent of the study area 
projection <- 3035 # 32630 or 3035
# projection_km <- "+proj=utm +zone=30 +ellps=WGS84 +units=km +no_defs"
if (projection == 32630) projection_km <- "+proj=utm +zone=30 +ellps=WGS84 +units=km +no_defs"
if (projection == 3035) projection_km <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=km +no_defs"

### Loading data & functions
  source("R/utils.R")
  source("R/data_load.R")
  source("R/conf.R")                  # for setting the data/model configuration
  source("R/data.R")                  # for preparing the data
  source("R/run.R")                   # for running the model with the above settings
  source("R/run_delta.R")             # for running the model with the above settings
  source("R/refit.R")                 # for rerunning the model in case prediction was disactivated 
  source("R/residual_diagnostics.R")	# for residual diagnostics
  source("R/Plotting_funcs.R")        # for various plotting functions

#======================================================================================================#
#===========  Fitting a spatio-temporal model to the age-disaggregated CPUE data ================#
#======================================================================================================#

dat$SSTfront_scl <- dat$distance_SST_fronts_ostia_scl
dat$VESSEL <- as.factor(dat$VESSEL)
# dat$YEAR_fct <- as.factor(dat$YEAR)

### Adding some polygons info to the data
polygons_df <- st_as_sf(x=All_strata %>% filter(ID != 13), coords = c("lon", "lat")) %>% st_set_crs(4326) %>% group_by(ID) %>% 
  summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") %>% st_transform(projection_km)
# calculating distance of each observation to polygon because not all points are within polygons
bla <- st_distance(polygons_df, st_as_sf(x=dat, coords = c("LON", "LAT")) %>% st_set_crs(4326) %>% st_transform(projection_km) %>% st_cast("POINT"))
dat$Area = apply(bla, 2, which.min)
dat <- dat %>% st_drop_geometry()
dat$YearArea = apply(dat[,c('YEAR', 'Area')], 1, function(x) paste(x, collapse="_"))


## Define configurations (both data and model) for the TMB model run
conf = conf_TMB(Model_type = "Annual_barrier",
                keep_omega = FALSE,
                keep_epsilon = TRUE,# if false here and keep_omega=FALSE, then it is a non-spatial model 
                include_age_correlation = "allinone",     # "none", "extra", "with_epsilon", "allinone"
                plotting = FALSE,
                Do_predict = 1,
                mixture_model = 3,  # 0 = spatial_SDM_RE (default), 1 = spatial_SDM_RE_mixture, 2=spatial_SDM_RE_barriereffect, 3=spatial_SDM_RE_delta
                density_dependence = FALSE,
                years = 2010:2020, # if including mixed_depth stop at 2019
                plusgroup = 16,  
                ages = 3:10,
                family = 3, # 0 is tweedie, 1=gaussian, 2=gamma, 3=student-t, 4=lognormal
                link = 1, # log link
                corr_str = 1,  # 1 = AR1 correlation in spatial distribution between age (default), for each year, 2 = AR1 correlation in spatial distribution between year, for each age
                ARorIID = c(1,2), # 0 = IID, 1 = AR1, 2 = user specified corr matrix
                # conf_model = list(transf_rho = factor(rep(NA, 10))),
                # conf_model = NULL,
                # conf_model = list(ln_phi     = factor(c(1001:1010)),# this is the number of age groups
                #                   # ln_phi     = factor(c(1001:1010)),
                #                   # transf_rho = factor(rep(3000, 11)), # this is the number of year if corr_str = 0
                #                   transf_rho = factor(rep(3001, 10)), # this is the number of ages if corr_str = 1
                #                   # transf_rho1 = factor(rep(5000, 11)), # this is the number of year
                #                   # logKappa = factor(rep(4000, 11)), # this is the number of year if corr_str = 0
                #                   logKappa = factor(rep(4001, 10)), # this is the number of ages if corr_str = 1
                #                   # logKappa1 = factor(rep(6000, 4)), # this is the number of year
                #                   # thetaf     = factor(c(2000,2001,rep(2002,2),rep(2003,6)))),# this is the number of age groups
                #                   thetaf     = factor(c(2001:2010))),# this is the number of age groups
                #                   # thetaf     = factor(c(2001:2010))),# this is the number of age groups
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_0m_scl,k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_0m_scl,k=3) + s(REPNRTCHL_0m_scl,k=3)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_0m_scl,k=3) + s(REPNRTCHL_0m_scl,k=3) + s(SSTfront_scl, k=3)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_clim_0m_scl,k=3) + s(CHL_clim_0m_scl,k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_clim_0m_scl,k=3, by=YEAR_fct) + s(CHL_clim_0m_scl,k=3, by=YEAR_fct) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_0m_scl,k=3) + s(REPNRTCHL_0m_scl,k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) ),
                formula_pres = formula(CPUE ~ 1),
                # formula_pres = formula(CPUE ~ -1 + as.factor(YEAR)),
                # formula_pres = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_0m_scl,k=3) + s(REPNRTCHL_0m_scl, k=3)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR)+ (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SSTfront_scl,k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(CHL_clim_0m_scl,k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_clim_0m_scl,k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_clim_0m_scl,k=3) + s(CHL_clim_0m_scl, k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + s(SST_0m_scl,k=3) + s(REPNRTCHL_0m_scl, k=3) + (1|VESSEL) + (1|YearArea)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_0m_scl,k=3) + s(REPNRTCHL_0m_scl, k=3) + s(OMLT_0m_scl, k=3) + (1|VESSEL)),
                formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_0m_scl,k=3) + s(OMLT_0m_scl, k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) +  te(SST_clim_0m_scl, CHL_clim_0m_scl, k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + ti(SST_clim_0m_scl,k=3) + ti(CHL_clim_0m_scl, k=3) + ti(SST_clim_0m_scl, CHL_clim_0m_scl, k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_clim_0m_scl,k=3) + s(CHL_clim_0m_scl, k=3) + s(CHL_clim_1m_scl, k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_clim_0m_scl,k=3) + s(CHL_clim_0m_scl, k=3) + s(SST_clim_1m_scl,k=3) + s(CHL_clim_1m_scl, k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(CHL_clim_0m_scl,k=3) + s(SSTfront_scl, k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_clim_0m_scl,k=3) + s(CHL_clim_0m_scl, k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_0m_scl,k=3) + s(REPNRTCHL_0m_scl,k=3) + s(SSTfront_scl, k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_clim_0m_scl,k=3) + s(REPNRTCHL_0m_scl,k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_clim_0m_scl,k=3) + s(CHL_clim_0m_scl,k=3) + s(SST_clim_1m_scl,k=3) + s(CHL_clim_1m_scl,k=3) + s(SSTfront_scl, k=3) + (1|VESSEL)),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) ),
                formula_mix = formula(CPUE ~ 1),
                mesh_type = "cutoff", # or "kmeans" "cutoff", or "own" (which modifies to the mackerel case)
                knots_mesh = 250,
								bias_correct = FALSE								
                )

# library(GGally)
# dat[,c("SST_clim_0m_scl", "SST_clim_1m_scl", "CHL_clim_0m_scl", "CHL_clim_1m_scl", "SSTfront_scl")] %>% st_drop_geometry() %>% 
# ggcorr( )

# corr_list <- list()
# for (yr in seq_along(conf$years)) corr_list[[yr]] <- cor(subset(datdat$data, YEAR == conf$years[yr])[,82:91], use = "complete.obs")
# sapply(corr_list, function(x) eigen(x)$values)
# par(mfrow=c(4,3), mar=c(0,0,0,0))
# for (i in 1:11) image(corr_list[[i]])
# corrmat <- nearPD(corrmat, corr=TRUE)$mat


## Set up data & Define parameters
datdat = Prepare_data(bio=bio, catch=dat, conf=conf, rerun_stack=FALSE)

saveRDS(datdat, file="delta_data.rds")
saveRDS(conf, file="delta_conf.rds")

# datdat$tmb_params$beta = matrix(1, nrow=ncol(datdat$tmb_data$X), ncol=conf$Nage)
## If you want to fix some fixed effect parameters after the fact (to simplify the model)
# aaa <- matrix(1:(ncol(datdat$tmb_data$X)*datdat$tmb_data$Nage), nrow=ncol(datdat$tmb_data$X), ncol= datdat$tmb_data$Nage, byrow=F)  
# colnames(datdat$tmb_data$X)
# aaa[16:19,1] = NA
# conf_extra <- list(beta = aaa)

## Fit model
startTime = Sys.time()
run = fit_mackerelST_delta(data=datdat, conf, conf_extra = NULL)
endTime = Sys.time()
endTime-startTime

betas <- matrix(run$opt$par[grep("beta", names(run$opt$par))], nrow=ncol(datdat$tmb_data$X), ncol= datdat$tmb_data$Nage, byrow=F)  

## Extract output
# output <- extract_output(run, datdat)


gg_control <- theme(axis.title = element_text(size=15),
                    axis.text = element_text(size=14), 
                    plot.title = element_text(hjust=0.5, size=18),
                    strip.text = element_text(size=14))  

### Adding some polygons info to the data
# 	polygons_df <- st_as_sf(x=All_strata %>% filter(ID != 13), coords = c("lon", "lat")) %>% st_set_crs(4326) %>% group_by(ID) %>% 
#   summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") %>% st_transform(projection_km)
# 	# calculating distance of each observation to polygon because not all points are within polygons
# 	bla <- st_distance(polygons_df, st_as_sf(x=datdat$data, coords = c("sdm_x", "sdm_y")) %>% st_set_crs(projection_km) %>% st_cast("POINT"))
# 	datdat$data$Area = apply(bla, 2, which.min)
	
### Doing plotting of the results 

    sfolder <- paste0(getwd(), "/plots/", conf$save_folder, "/Delta_", paste(attributes(conf)$fixed_effect[-1], collapse="_"), "_", attributes(conf)$RE_effects, "_student")
   
    ## save model output
      ALL <- list()
      ALL$conf <- conf
      ALL$data <- datdat
      ALL$run <- run
      saveRDS(ALL, file=paste0(sfolder, "/run.rds"))
    
    ## Model diagnostics
      # AIC
        run$opt$AIC 
      
      ## Check if all parameters are estimable  
        fixed_est <- summary(run$sdrep, "fixed")
        random_est <- summary(run$sdrep, "random")
      
      # qqlot (only works for non delta model)
        if (conf$mixture_model != 3) out <- Plot_residuals(run=run, data=datdat, conf=conf, name="QQplot", folder=sfolder)

      # cross validation 
        cross_val <- cv(run=run, data=datdat, k_fold = 10, conf, seed=123)
        saveRDS(cross_val, file=paste0(sfolder, "/cross_val.rds"))
        
      # DHARma-like residual test (can be used for any kind of models) 
        Nsim = 500
        s_nb2 <- sapply(seq_len(Nsim), function(x) run$obj$simulate(par= run$obj$env$last.par.best, complete=FALSE)$yobs, simplify = 'array')
        
      # compare the number of zeros in the data with the number of zeros in the prediction
        grep_val = paste0("^", conf$ages[1], "$")
        col_sel = grep(grep_val, colnames(datdat$data))+c(1:(conf$Nage))-1
        nb_zero <- apply(datdat$data[,col_sel], 2, function(x) sum(x==0, na.rm=T)/length(x))
        pred_zero <- apply(sapply(seq_len(Nsim), function(x) apply(s_nb2[,,x], 2, function(xx) sum(xx==0, na.rm=T)/length(xx))),1,mean)
        nb_zero; pred_zero; # quite similar!
        write.table(data.frame(obs=nb_zero, pred=pred_zero), file=paste0(sfolder, "/P_zero.txt"))
        
        for (age in 1:(conf$Nage)){
          to_keep = which(is.na(datdat$data[,col_sel[age]])==FALSE)
          r_nb2 <- DHARMa::createDHARMa(
            simulatedResponse = s_nb2[to_keep,age,],
            observedResponse = datdat$data[to_keep,col_sel[age]],
            fittedPredictedResponse = fixed[to_keep,age]
          )
          png(paste0(sfolder, "/resid_age", age, ".png"), width= 10, height=6, units="in", res=300)
          plot(r_nb2)
          dev.off()
        }
        
      # one step ahead residual that avoids the issue of pit residuals for state-space model
  			pred <- oneStepPredict(run$obj, observation.name ="yobs", data.term.indicator = "keep", discrete=FALSE)
  			qqnorm(pred$residual[sel[which(!is.na(pred$observation[sel]) == TRUE)]]); abline(0,1)
  			qqnorm(pred$residual[1:datdat$tmb_data$Nobs]); abline(0,1)
				qqnorm(pred$residual[datdat$tmb_data$Nobs+1:2*datdat$tmb_data$Nobs]); abline(0,1)
			
      # check consistency
       # chk <- checkConsistency(run$obj, n=250)
       # summary(chk)
       
       # p-value=0, but the simulation sample size is small so this is logical for a very complex model. 
       # but bias is small so it is not an issue. 4.584637e-07 to 8.881380e-02
       # but marginal bias is high...  0.002215884 - 21.966392837
       
       
               
    ## Calculating the spatial range parameter
      sqrt(8)/exp(run$opt$par[grep("Kappa", names(run$opt$par))])
   
    ## Looking at the covariance matrix
      covar <- run$obj$report()$covar
      
    ## Calculating the conditional R2 (in link space) for each age group
        fixed <- run$obj$report(run$obj$env$last.par.best)$fixed
        spatiotemporal <- run$obj$report(run$obj$env$last.par.best)$epsilon_st_A_mat

        R2_cond <- c()
        fixed_contribution <- c()
        ST_contribution <- c()
        for (i in seq_along(conf$ages)){
          age = conf$ages[i]
					var_fixed = var(fixed[,i])
          var_stRE = var(spatiotemporal[,i])
          var_random = exp(run$opt$par[grep("tau_G", names(run$opt$par))])[i]
          mu_null = mean(datdat$data[,grep(paste0("^",conf$ages[i],"$"), colnames(datdat$data))], na.rm=T)  # this is the fixed effect of the null model (just exp(intercept))
          s1 = boot::inv.logit(run$opt$par[grep("thetaf", names(run$opt$par))[i]]) + 1
          phi = exp(run$opt$par[grep("ln_phi", names(run$opt$par))[i]])
          var_obs = mu_null^(s1-2)* phi #- 1/4*phi^2*mu_null^(2*s1-4))
          R2_cond[i] <- (var_fixed + var_stRE + var_random)/(var_fixed + var_stRE + var_random + var_obs)
          fixed_contribution[i] <- (var_fixed)/(var_fixed + var_stRE + var_random + var_obs)
          ST_contribution[i] <- (var_stRE)/(var_fixed + var_stRE + var_random + var_obs)
        }
				R2_cond
        fixed_contribution/R2_cond
        write.table(data.frame(R2=R2_cond, fixed_contribution=fixed_contribution/R2_cond), file=paste0(sfolder, "/R2_contribution.txt"))
        
        
#### If we want to refit something      

  ## Case when we want to reactivate prediction
      if (load_result == TRUE){
        ALL <- readRDS(paste0(getwd(), "\\plots\\Annual_barrier_omegaFALSE_agecorallinone\\SST_clim_0m_scl_CHL_clim_0m_scl_VESSEL\\run.rds"))
        ALL <- readRDS(paste0(getwd(), "\\plots\\New_run\\Env_strata\\run.rds"))
        conf <- ALL$conf
        datdat <- ALL$data
        run <- ALL$run   
        sfolder <- paste0(getwd(), "/plots/", conf$save_folder, "/", paste(attributes(conf)$fixed_effect[-1], collapse="_"), "_", attributes(conf)$RE_effects)
        sfolder <- paste0(getwd(), "\\plots\\New_run\\Env_strata\\")
        version <- paste0(getwd(), "/src/spatial_SDM_RE")
        dyn.load(dynlib(version))
        run$obj$retape()
        run$obj$fn(run$opt$par)
      }
        
      ## Producing the marginal effect plots (done within the model for more flexibility & complex covaraite structure)
        ## Marginal SST effect 
          # SST_dat <- data.frame(expand.grid(SST_clim_0m_scl = seq(min(dat$SST_clim_0m_scl),max(dat$SST_clim_0m_scl),length.out=50),
          #                                   CHL_clim_0m_scl = seq(min(dat$CHL_clim_0m_scl),max(dat$CHL_clim_0m_scl),length.out=50),
          #                                   YEAR = 2020, X1000 = mean(datdat$data$X1000), Y1000=mean(datdat$data$Y1000), CPUE=1))
          SST_dat <- data.frame(expand.grid(SST_0m_scl = seq(min(dat$SST_0m_scl),max(dat$SST_0m_scl),length.out=50),
					                                  REPNRTCHL_0m_scl = 0, Area = 2,  
																						YEAR = 2020, X1000 = mean(datdat$data$X1000), Y1000=mean(datdat$data$Y1000), CPUE=1))
					Marginal_SST <- refit(run=run, data=datdat, pred_data=SST_dat, predict_type = 2, conf=conf) 
					SST_marg <- as.data.frame(Marginal_SST$dens_pred)
					colnames(SST_marg) <- paste0("Age",conf$ages)
					SST_dat <- cbind(SST_dat, SST_marg)
					SST_dat_df <- pivot_longer(SST_dat, cols =  starts_with("Age"), names_to="Age")
					fake = datdat$data
					fake$value = 1
					fake$Age = "Age1"
					p1 <- ggplot(SST_dat_df, aes(x=SST_0m_scl, y = value, col=Age)) + geom_line(size=1) + 
						scale_color_viridis_d() + theme_bw() + labs(y="Marginal effect (link scale)") + 
						geom_rug(data=fake, col="black")
          ggsave(p1, filename = paste0(sfolder, "/Marginal_SST.jpeg"), dpi ="retina", width = 6, height = 4, device = "jpeg")
 
				## Marginal CHL effect 
					# CHL_dat <- data.frame(expand.grid(SST_clim_0m_scl = 0,
					# 																	CHL_clim_0m_scl = seq(min(dat$CHL_clim_0m_scl),max(dat$CHL_clim_0m_scl),length.out=50),
					# 																	YEAR = 2020, X1000 = mean(datdat$data$X1000), Y1000=mean(datdat$data$Y1000), CPUE=1))
					CHL_dat <- data.frame(expand.grid(SST_0m_scl = 0,
																						REPNRTCHL_0m_scl = seq(min(dat$REPNRTCHL_0m_scl),max(dat$REPNRTCHL_0m_scl),length.out=50),
																						YEAR = 2020, Area = 2,
																						X1000 = mean(datdat$data$X1000), Y1000=mean(datdat$data$Y1000), CPUE=1))
					Marginal_CHL <- refit(run=run, data=datdat, pred_data=CHL_dat, predict_type = 2, conf=conf) 
					CHL_marg <- as.data.frame(Marginal_CHL$dens_pred)
					colnames(CHL_marg) <- paste0("Age",conf$ages)
					CHL_dat <- cbind(CHL_dat, CHL_marg)
					CHL_dat_df <- pivot_longer(CHL_dat, cols =  starts_with("Age"), names_to="Age")
					p1 <- ggplot(CHL_dat_df, aes(x=REPNRTCHL_0m_scl, y = value, col=Age)) + geom_line(size=1) + 
						scale_color_viridis_d() + theme_bw() + labs(y="Marginal effect (link scale)") + 
						geom_rug(data=fake, col="black")
          ggsave(p1, filename = paste0(sfolder, "/Marginal_CHL.jpeg"), dpi ="retina", width = 6, height = 4, device = "jpeg")
        
				## SST effect with CHL
					# SST_dat <- data.frame(expand.grid(SST_clim_0m_scl = seq(min(dat$SST_clim_0m_scl),max(dat$SST_clim_0m_scl),length.out=50),
																						# CHL_clim_0m_scl = seq(min(dat$CHL_clim_0m_scl),max(dat$CHL_clim_0m_scl),length.out=50),
																						# YEAR = 2020, X1000 = mean(datdat$data$X1000), Y1000=mean(datdat$data$Y1000), CPUE=1))
					# SST_dat <- data.frame(expand.grid(SST_0m_scl = seq(min(dat$SST_0m_scl),max(dat$SST_0m_scl),length.out=50),
																						# REPNRTCHL_0m_scl = seq(min(dat$REPNRTCHL_0m_scl),max(dat$REPNRTCHL_0m_scl),length.out=50),
																						# YEAR = 2020, X1000 = mean(datdat$data$X1000), Y1000=mean(datdat$data$Y1000), CPUE=1))
					# # data=datdat; pred_data=SST_dat; predict_type = 2 
					# Marginal_SST <- refit(run=run, data=datdat, pred_data=SST_dat, predict_type = 2, conf=conf) 
					# SST_marg <- as.data.frame(Marginal_SST$dens_pred)
					# colnames(SST_marg) <- paste0("Age", conf$ages)
					# SST_dat <- cbind(SST_dat, SST_marg)
					# SST_CHL_dat_df <- pivot_longer(SST_dat, cols = 7:(6+length(conf$ages)), names_to="Age")
					# SST_CHL_dat_df$Age <- factor(SST_CHL_dat_df$Age, levels=paste0("Age", conf$ages))
					
					# fake = datdat$data
					# fake$value = 1
					# p1 <- ggplot(SST_CHL_dat_df, aes(x=SST_0m_scl, y = REPNRTCHL_0m_scl, col=value, fill=value)) +
					# facet_wrap(. ~ Age) + geom_tile() + theme_bw() + 
						# scale_color_viridis_c(option="inferno", name="Effect") + scale_fill_viridis_c(option="inferno", name="Effect") + 
						# geom_point(data= fake, aes(x=SST_0m_scl, y = REPNRTCHL_0m_scl), col="black", fill="black", cex=0.5, alpha = 0.3)
					# ggsave(p1, filename = paste0(sfolder, "/Joint_SST_CHL_effects.jpeg"), dpi ="retina", width = 10, height = 8, device = "jpeg")
				


			

        
      ## Calculating the adreport variables: activating the prediction   
			# Mod_pred <- refit(run=run, pred_grid= NULL, predict_type=1, conf=conf) 
      # saveRDS(Mod_pred, file=paste0(sfolder, "/run_all.rds"))
      Mod_pred= run
			
      ## Index of abundance
      IA_IESSNS <- readRDS(paste0(getwd(), "/data/new_data/WGWIDE_IESSN_Biomass_billionTon.rds"))
      IA_IESSNS[1,9]<- 0.01
      ttt <- as.matrix(IA_IESSNS[conf$ages,conf$years-2007])
      qqq <- reshape2::melt(ttt)
      colnames(qqq) <- c("Age", "Year", "IA_IESSNS")
      qqq <- qqq %>% mutate(Year = as.numeric(Year), Age = as.numeric(as.character(factor(Age, labels=conf$ages))), IA_IESSNS = as.numeric(IA_IESSNS))
      # p1 <- ggplot(qqq, aes(x=Year, y=IA_IESSNS, col=Age)) + geom_line(size=2) + theme_bw() + gg_control
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/", model_number, "/IA_IESSN.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      
      IA2 <- as.data.frame(summary(Mod_pred$sdrep, "report"))
      if (ncol(IA2) ==2 ) colnames(IA2) <- c("IA", "IA_sd")
      if (ncol(IA2) ==4 ) colnames(IA2) <- c("IA_old", "IA_sd_old", "IA", "IA_sd")
      IA2$Year <- rep(1:conf$Nyear, conf$Nage)
      IA2$Age <- rep(conf$ages, each=conf$Nyear)
      IA2$IA_type <- rep(c("IA", "IA_corrected"), each=nrow(IA2)/2)
      write.table(IA2, file=paste0(sfolder, "/IA.txt"))
      
      ggplot(IA2 %>% filter(IA_type == "IA_corrected"), aes(x = Year, y = IA)) + theme_bw() + facet_wrap(.~Age, scale = "free_y") +
        geom_ribbon(aes(ymin = IA - IA_sd, ymax=IA + IA_sd), fill=gray(0.7, 0.5)) + geom_line(size=1.2) 
      
      
      plot_IA(mod = Mod_pred, conf = conf, IA_assessment = qqq, name = "IAnew", folder = sfolder)
      CG_plots(run = Mod_pred, name="IA", data=datdat, no_fish_zone=FALSE, folder = sfolder, conf=conf)
      CG_plots(run = Mod_pred, name="IA", data=datdat, no_fish_zone=TRUE, folder = sfolder, conf=conf)
      
      Plot_maps(run=Mod_pred, data=datdat, name="Best", conf=conf, no_fish_zone=FALSE, folder=sfolder)
 
      ## Plotting the env covariates
      Plot_env(run = Mod_pred, var = "SST_clim_0m_scl", truncate=TRUE, folder=sfolder)
      Plot_env(run = Mod_pred, var = "SST_0m", truncate=TRUE, folder=sfolder)
      Plot_env(run = Mod_pred, var = "SST_0m_scl", truncate=TRUE, folder=sfolder)
      # Plot_env(run = Mod_pred, var = "SST_clim_1m_scl", truncate=FALSE, folder=sfolder)
      Plot_env(run = Mod_pred, var = "CHL_clim_0m_scl", truncate=TRUE, folder=sfolder)
      Plot_env(run = Mod_pred, var = "REPNRTCHL_0m", truncate=TRUE, folder=sfolder)
      Plot_env(run = Mod_pred, var = "REPNRTCHL_0m_scl", truncate=TRUE, folder=sfolder)
      Plot_env(run = Mod_pred, var = "OMLT_0m_scl", truncate=TRUE, folder=sfolder)
      # Plot_env(run = Mod_pred, var = "CHL_clim_1m_scl", truncate=FALSE, folder=sfolder)
      # Plot_env(run = Mod_pred, var = "SSTfront_scl", truncate=FALSE, folder=sfolder)
      
      
 
      p1 <- ggplot(datdat$pred_grid) +
        geom_tile(aes(x=X1000, y=Y1000, fill = SSTfront_scl)) +
        geom_sf(data = Atlantic_proj$geometry, color ="grey27", size = .2)+
        xlab("") + ylab("") +
        coord_sf()+ theme_bw() + 
        facet_wrap(~YEAR)+
        # scale_fill_viridis_c() +
        scale_fill_gradient2(low="darkblue", mid="white", high="red2", midpoint=0) +
        ggtitle(paste0("SSTfront_scl", " distribution")) + 
        theme(axis.title = element_text(size=15),
              axis.text = element_text(size=14), 
              plot.title = element_text(hjust=0.5, size=18),
              strip.text = element_text(size=14))
      ggsave(p1, filename = paste0(sfolder, "/Env_map", "_", "SSTfront_scl", ".pdf"), dpi ="retina", width = 12, height = 10, device = "pdf")
      
           
      # asd <- hist(data$SST_scl, breaks=seq(floor(min(Mod_pred$test2$SST_scl, na.rm=T)), ceiling(max(Mod_pred$test2$SST_scl, na.rm=T)), by=0.5))
      # qwe <- hist(Mod_pred$test2$SST_scl,breaks=asd$breaks)
      # hist_val <- data.frame(breaks=c(asd$mids,asd$mids), source = c(rep("data", length(asd$mids)), rep("prediction", length(asd$mids))), val= c(asd$densiy, qwe$density))
      # ggplot(data=hist_val, aes(x=breaks, y=val)) + geom_line(size=1.5) + gg_control + theme_bw() +
      #   facet_grid(source~.)
      # asd <- hist(data$CHL_scl, breaks=seq(floor(min(Mod_pred$test2$CHL_scl, na.rm=T)), ceiling(max(Mod_pred$test2$CHL_scl, na.rm=T)), by=0.5))
      # qwe <- hist(Mod_pred$test2$CHL_scl,breaks=asd$breaks)
      # hist_val <- data.frame(breaks=c(asd$mids,asd$mids), source = c(rep("data", length(asd$mids)), rep("prediction", length(asd$mids))), val= c(asd$densiy, qwe$density))
      # ggplot(data=hist_val, aes(x=breaks, y=val)) + geom_line(size=1.5) + gg_control + theme_bw() +
      #   facet_grid(source~.) + coord_cartesian(xli=c(-5,10))
      # asd <- hist(data$mixed_scl, breaks=seq(floor(min(Mod_pred$test2$mixed_scl, na.rm=T)), ceiling(max(Mod_pred$test2$mixed_scl, na.rm=T)), by=0.5))
      # qwe <- hist(Mod_pred$test2$mixed_scl,breaks=asd$breaks)
      # hist_val <- data.frame(breaks=c(asd$mids,asd$mids), source = c(rep("data", length(asd$mids)), rep("prediction", length(asd$mids))), val= c(asd$densiy, qwe$density))
      # ggplot(data=hist_val, aes(x=breaks, y=val)) + geom_line(size=1.5) + gg_control + theme_bw() +
      #   facet_grid(source~.)
      # 
      
      
      
      
      
      
      
      
      
           
  # ## Case when we eliminate the extreme values of covariates (max = max(data))
  #     pred_grid_new <- pred_grid   
  #     pred_grid_new$CHL_scl <- ifelse(pred_grid_new$CHL_sc > max(data$CHL_sc), max(data$CHL_sc), pred_grid_new$CHL_sc)
  #     pred_grid_new$mixed_scl <- ifelse(pred_grid_new$mixed_scl > max(data$mixed_scl), max(data$mixed_scl), pred_grid_new$mixed_scl)
  #     pred_grid_new$SST_scl <- ifelse(pred_grid_new$SST_scl > max(data$SST_scl), max(data$SST_scl), pred_grid_new$SST_scl)
  #     
  #     Mod_extre <- refit(obj=tmb_obj_phase2, opt= opt_phase2, tmb_data=tmb_data, pred_grid=pred_grid_new) 
  #     
  #     plot_IA(report = Mod_extre$report, name = "rm_extreme", folder = sfolder)
  #     CG_plots(Mod_extre$test2, name="rm_extreme", no_fish_zone=FALSE, folder = sfolder)
  #     CG_plots(Mod_extre$test2, name="rm_extreme", no_fish_zone=TRUE, folder = sfolder)
  #     
  #     Plot_maps(pred_mu=Mod_extre$test2, pred_epsilon=Mod_extre$qwe2, name="rm_extreme", no_fish_zone=FALSE, folder=sfolder)
  #     
  #     

      
      
      
      
      
      
      
         
      
 ##### Final plots of the results
      
      # Create a continuous palette function
      Atlantic_proj2 <- st_transform(Atlantic, crs=projection_km)
      pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), test2$Pred,
                          na.color = "transparent")
      # Layers
      Years = levels(as.factor(test2$YEAR))
      Ages = levels(as.factor(test2$AGE))
      
      basemap <- leaflet() %>% 
        addProviderTiles(providers$Esri.OceanBasemap) %>%
        addSimpleGraticule(interval = 5)
      
      #Be patient the following takes time as it produces 16 maps per year
      library(leafsync)
      maplist_byYEAR <- list()
      
      for(q in Years){
        maplist_byYEAR[[q]] <- list()
        for(i in 1:13){
          sub_data <- subset(test2, subset=c(YEAR== q & AGE ==i))
          pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), sub_data$Pred,
                              na.color = "transparent")
          map <- basemap
          .df <- test2 %>% 
            filter(YEAR==q, AGE==i) %>% 
            mutate(x=X1000,y=Y1000, z=Pred) %>% 
            dplyr::select(x,y,z) 
          r <- rasterFromXYZ(.df, crs = projection_km)
          r_masked <- mask(r, Atlantic_proj2, inverse=T)
          bladat <- data %>%
            filter(YEAR==q) 
          maplist_byYEAR[[q]][[i]] <- map %>%
            addRasterImage(r_masked, colors = pal, opacity = 1) %>%
            addCircleMarkers(
              lat=bladat$LAT, lng=bladat$LON, radius=0.5,
              color= "red") %>% 
            addLegend(pal = pal, values = sub_data$Pred, opacity = 1,
                      title =  paste0("Predictions<br>", q, " Age: ", i))
        }
      }
      
      #You can visualize each Year plotting the following
      
      latticeview(maplist_byYEAR$`2007`)
      # latticeview(maplist_byYEAR$`2008`)
      latticeview(maplist_byYEAR$`2009`)
      latticeview(maplist_byYEAR$`2010`)
      latticeview(maplist_byYEAR$`2011`)
      latticeview(maplist_byYEAR$`2012`)
      latticeview(maplist_byYEAR$`2013`)
      latticeview(maplist_byYEAR$`2014`)
      latticeview(maplist_byYEAR$`2015`)
      latticeview(maplist_byYEAR$`2016`)
      latticeview(maplist_byYEAR$`2017`)
      latticeview(maplist_byYEAR$`2018`)
      latticeview(maplist_byYEAR$`2019`)
      latticeview(maplist_byYEAR$`2020`)
   
      
      
      
         
      
      maplist_byAGE <- list()
      
      for(i in Ages){
        maplist[[i]] <- list()
        for(q in Years){
          map <- basemap
          .df <- ALL_data %>% 
            filter(YEAR==q, AGE==i) %>% 
            mutate(x=X1000*1000,y=Y1000*1000, z=CPUE_pred) %>% 
            dplyr::select(x,y,z) 
          r <- rasterFromXYZ(.df, crs = CRS("+init=epsg:32630"))
          r_masked <- mask(r, Atlantic_proj, inverse=T)
          maplist_byAGE[[i]][[q]] <- map %>%
            addRasterImage(r_masked, colors = pal, opacity = 1) %>%
            addLegend(pal = pal, values = exp(predictions$data$est), opacity = 1,
                      title =  paste0("Predictions<br>", q, " Age: ", i))
        }
        
      }
      
      #You can visualize each Year plotting the following
      
      latticeview(maplist_byAGE$`1`)
      #latticeview(maplist_byAGE$`2`)
      #latticeview(maplist_byAGE$`3`)     
      
      
      
      #### More plots using leaflet 
      # # Create a continuous palette function
      # pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), test2$Pred,
      #                     na.color = "transparent")
      # # Layers
      # Years = levels(as.factor(test2$YEAR))
      # Ages = levels(as.factor(test2$AGE))
      # 
      # basemap <- leaflet() %>% 
      #   addProviderTiles(providers$Esri.OceanBasemap) %>%
      #   addSimpleGraticule(interval = 5)
      # 
      # #Be patient the following takes time as it produces 16 maps per year
      # library(leafsync)
      # maplist_byYEAR_PA <- list()
      # 
      # for(q in Years){
      #   maplist_byYEAR_PA[[q]] <- list()
      #   for(i in Ages){
      #     map <- basemap
      #     .df <- test2 %>% 
      #       filter(YEAR==q, AGE==i) %>% 
      #       mutate(x=X1000*1000,y=Y1000*1000, z=Pred) %>% 
      #       dplyr::select(x,y,z) 
      #     r <- rasterFromXYZ(.df, crs = CRS("+init=epsg:32630"))
      #     r_masked <- mask(r, Atlantic_proj, inverse=T)
      #     maplist_byYEAR_PA[[q]][[i]] <- map %>%
      #       addRasterImage(r_masked, colors = pal, opacity = 1) %>%
      #       addLegend(pal = pal, values = test2$Pred, opacity = 1,
      #                 title =  paste0("Predictions<br>", q, " P(age ", i, ")"))
      #   }
      # }
      # 
      # #You can visualize each Year plotting the following
      # 
      # latticeview(maplist_byYEAR_PA$`2012`)
      # latticeview(maplist_byYEAR_PA$`2013`)
      # latticeview(maplist_byYEAR_PA$`2014`)
      # latticeview(maplist_byYEAR_PA$`2015`)
      # latticeview(maplist_byYEAR_PA$`2016`)
      # latticeview(maplist_byYEAR_PA$`2017`)
      # latticeview(maplist_byYEAR_PA$`2018`)
      # latticeview(maplist_byYEAR_PA$`2019`)
      # latticeview(maplist_byYEAR_PA$`2020`)
      # 
      # # By age
      # maplist_byAGE_PA <- list()
      # 
      # for(i in Ages){
      #   maplist_byAGE_PA[[i]] <- list()
      #   for(q in Years){
      #     map <- basemap
      #     .df <- test2 %>% 
      #       filter(YEAR==q, AGE==i) %>% 
      #       mutate(x=X1000*1000,y=Y1000*1000, z=Pred) %>% 
      #       dplyr::select(x,y,z) 
      #     r <- rasterFromXYZ(.df, crs = CRS("+init=epsg:32630"))
      #     r_masked <- mask(r, Atlantic_proj, inverse=T)
      #     maplist_byAGE_PA[[i]][[q]] <- map %>%
      #       addRasterImage(r_masked, colors = pal, opacity = 1) %>%
      #       addLegend(pal = pal, values = test2$Pred, opacity = 1,
      #                 title =  paste0("Predictions<br>", "P(age ", i, ") Year: ", q))
      #   }
      #   
      # }
      # 
      # #You can visualize each Year plotting the following
      # 
      # latticeview(maplist_byAGE_PA$`1`)
      # latticeview(maplist_byAGE_PA$`2`)
      # latticeview(maplist_byAGE_PA$`3`)     
      # latticeview(maplist_byAGE_PA$`4`)     
      # latticeview(maplist_byAGE_PA$`5`)     
      # latticeview(maplist_byAGE_PA$`6`)     
      # latticeview(maplist_byAGE_PA$`7`)     
      # latticeview(maplist_byAGE_PA$`8`)     
      # latticeview(maplist_byAGE_PA$`9`)     
      # latticeview(maplist_byAGE_PA$`10`)     
      # latticeview(maplist_byAGE_PA$`11`)     
      # 
      
      # 
      # qres_tweedie <- function(TMBobject, y, mu, age, seed=NULL, to_remove) {
      #   if(!is.null(seed)) set.seed(seed)
      #   p <- stats::plogis(TMBobject$par[grep("thetaf", names(TMBobject$par))[age]]) + 1
      #   dispersion <- exp(TMBobject$par[grep("ln_phi", names(TMBobject$par))[age]])
      #   if (!is.null(to_remove)) u <- fishMod::pTweedie(q = y[-to_remove,age], p = p, mu = mu[-to_remove,age], phi = dispersion)
      #   if (is.null(to_remove)) u <- fishMod::pTweedie(q = y[,age], p = p, mu = mu[,age], phi = dispersion)
      #   if (p > 1 && p < 2)
      #     if (!is.null(to_remove)) u[y[-to_remove,age] == 0] <- stats::runif(sum(y[-to_remove,age] == 0), min = 0, max = u[y[-to_remove,age] == 0])
      #   if (is.null(to_remove)) u[y[,age] == 0] <- stats::runif(sum(y[,age] == 0), min = 0, max = u[y[,age] == 0])
      #   stats::qnorm(u)
      # }
      # 
      # do_qqplot <- function(age, opt=opt_phase3, obj = tmb_obj_phase3, seed=2,to_remove){
      #   res <- as.data.frame(yobs)
      #   res <- res[!is.na(yobs[,1]),]
      #   res$resid <- qres_tweedie(opt, y=as.matrix(yobs), mu=obj$report()$mu, age=age, seed=seed, to_remove)
      #   res$Year = data$YEAR[!is.na(yobs[,1])]
      #   res$X <- data$X1000[!is.na(yobs[,1])]
      #   res$Y <- data$Y1000[!is.na(yobs[,1])]
      #   ggplot(res, aes(sample=resid)) + facet_wrap(.~Year) + stat_qq() + stat_qq_line() + theme_bw() +
      #     coord_cartesian(xlim=c(-3,3),ylim=c(-3,3)) + ggtitle(paste0("Age ", age))
      # }
      # do_qqplot_grouped <- function(age, opt=opt_phase3, obj = tmb_obj_phase3, seed=2, to_remove){
      #   res <- as.data.frame(yobs)
      #   res <- res[!is.na(yobs[,1]),]
      #   res$resid <- qres_tweedie(opt, y=as.matrix(yobs), mu=obj$report()$mu, age=age, seed=seed, to_remove)
      #   res$Year = data$YEAR[!is.na(yobs[,1])]
      #   res$X <- data$X1000[!is.na(yobs[,1])]
      #   res$Y <- data$Y1000[!is.na(yobs[,1])]
      #   ggplot(res, aes(sample=resid)) + stat_qq() + stat_qq_line() + theme_bw() +
      #     coord_cartesian(xlim=c(-3,3),ylim=c(-3,3)) + ggtitle(paste0("Age ", age))
      # }
      # 
      # 
      # 
      # opt_phase_use <- opt_phase2
      # tmb_obj_phase_use <- tmb_obj_phase2
      # 
      # p1 <- do_qqplot(age=1, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) # OK
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age1.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot(age=2, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) # OK
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age2.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot(age=3, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) # OK-ish
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age3.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot(age=4, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) # hum...
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age4.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot(age=5, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) # hum...
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age5.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot(age=6, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) # hum...
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age6.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot(age=7, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) # hum...
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age7.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot(age=8, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) # OKish
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age8.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot(age=9, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) # OKish
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age9.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot(age=10, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) # OK
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age10.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # 
      # p1 <- do_qqplot_grouped(age=1, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) ; p1
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age1_aggregated.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot_grouped(age=2, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) ; p1
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age2_aggregated.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot_grouped(age=3, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) ; p1
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age3_aggregated.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot_grouped(age=4, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) ; p1
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age4_aggregated.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot_grouped(age=5, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) ; p1
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age5_aggregated.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot_grouped(age=6, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) ; p1
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age6_aggregated.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot_grouped(age=7, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) ; p1
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age7_aggregated.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot_grouped(age=8, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) ; p1
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age8_aggregated.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot_grouped(age=9, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove) ; p1
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age9_aggregated.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # p1 <- do_qqplot_grouped(age=10, opt=opt_phase_use, obj=tmb_obj_phase_use, to_remove=to_remove); p1
      # ggsave(p1, filename = paste0(getwd(), "/plots/VAST_type/", save_folder, "/res_age10_aggregated.pdf"), dpi ="retina", width = 10, height = 10, device = "pdf")
      # 
      # # ggplot(res) +
      # #   geom_point(aes_string("X", "Y", colour = "resid")) +
      # #   geom_sf(data = Atlantic_proj$geometry/1000, color ="grey27", size = .2)+
      # #   xlab("") + ylab("") +
      # #   facet_wrap(~Year) +
      # #   scale_color_gradient2() +
      # #   theme_minimal() + coord_sf(xlim=c(-2000,2000),ylim=c(6000,8500)) +
      # #   ggtitle("Residuals")
      # # 
      # 
      
      ## Residual maps
      # out$out$fixed_pred  <-  unlist(pivot_longer(data.frame(run$rep$fixed)[-out$to_remove,], cols=1:10, names_to = "AGE", values_to ="fixed_pred") %>% dplyr::select(fixed_pred))
      # out$out$RE_pred  <- unlist(pivot_longer(data.frame(run$rep$epsilon_st_A_mat)[-out$to_remove,], cols=1:10, names_to = "AGE", values_to ="RE_pred") %>% dplyr::select(RE_pred))
      res <- ggplot(Atlantic) + geom_sf() + geom_point(data = out$out %>% filter(age == 4), aes(x = LON,  y = LAT, col=residuals), size=1) + facet_wrap(.~YEAR) + theme_bw() + 
        scale_color_gradient2(low="blue",mid="white",high="red")
      ggsave(res, file=paste0(sfolder, "/res_map.pdf"), height=14, width=20, dpi=400)
      SST <- ggplot(Atlantic) + geom_sf() + geom_point(data = out$out %>% filter(age == 4), aes(x = LON,  y = LAT, col=SST_clim_0m_scl), size=1) + facet_wrap(.~YEAR) + theme_bw() + 
        scale_color_gradient2(low="blue",mid="white",high="red")
      ggsave(SST, file=paste0(sfolder, "/SST_map_july.pdf"), height=14, width=20, dpi=400)
      CHL <- ggplot(Atlantic) + geom_sf() + geom_point(data = out$out %>% filter(age == 4), aes(x = LON,  y = LAT, col=CHL_clim_0m_scl), size=1) + facet_wrap(.~YEAR) + theme_bw() + 
        scale_color_gradient2(low="blue",mid="white",high="red")
      ggsave(CHL, file=paste0(sfolder, "/CHL_map_july.pdf"), height=14, width=20, dpi=400)
      SST1 <- ggplot(Atlantic) + geom_sf() + geom_point(data = out$out %>% filter(age == 4), aes(x = LON,  y = LAT, col=SST_clim_1m_scl), size=1) + facet_wrap(.~YEAR) + theme_bw() + 
        scale_color_gradient2(low="blue",mid="white",high="red")
      ggsave(SST1, file=paste0(sfolder, "/SST_map_june.pdf"), height=14, width=20, dpi=400)
      CHL1 <- ggplot(Atlantic) + geom_sf() + geom_point(data = out$out %>% filter(age == 4), aes(x = LON,  y = LAT, col=CHL_clim_1m_scl), size=1) + facet_wrap(.~YEAR) + theme_bw() + 
        scale_color_gradient2(low="blue",mid="white",high="red")
      ggsave(CHL1, file=paste0(sfolder, "/CHL_map_june.pdf"), height=14, width=20, dpi=400)
      Front <- ggplot(Atlantic) + geom_sf() + geom_point(data = out$out %>% filter(age == 4), aes(x = LON,  y = LAT, col=SSTfront_scl), size=1) + facet_wrap(.~YEAR) + theme_bw() + 
        scale_color_gradient2(low="blue",mid="white",high="red")
      ggsave(Front, file=paste0(sfolder, "/Front_map_june.pdf"), height=14, width=20, dpi=400)
      fixed <- ggplot(Atlantic) + geom_sf() + geom_point(data = out$out %>% filter(age == 4), aes(x = LON,  y = LAT, col=fixed_pred), size=1) + facet_wrap(.~YEAR) + theme_bw() + 
        scale_color_gradient2(low="blue",mid="white",high="red")
      ggsave(fixed, file=paste0(sfolder, "/fixedeffect_map.pdf"), height=14, width=20, dpi=400)
      RE <- ggplot(Atlantic) + geom_sf() + geom_point(data = out$out %>% filter(age == 4), aes(x = LON,  y = LAT, col=RE_pred), size=1) + facet_wrap(.~YEAR) + theme_bw() + 
        scale_color_gradient2(low="blue",mid="white",high="red")
      ggsave(RE, file=paste0(sfolder, "/REeffect_map.pdf"), height=14, width=20, dpi=400)
      
      ## Resid vs covariate
      ggplot(data = out, aes(x = SSTfront_scl,  y = residuals)) + facet_wrap(.~age) + theme_bw() + geom_point() + geom_smooth()
      ggplot(data = out, aes(x = REPCHL_0m_scl,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth()+ 
        scale_x_continuous(limits = c(-2,4))
      ggplot(data = out, aes(x = SST_0m_scl,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth()
      ggplot(data = out, aes(x = vgpm_0m_scl,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth()
      ggplot(data = out, aes(x = SST_1m_scl,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth()
      ggplot(data = out, aes(x = REPCHL_1m_scl,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth()
      ggplot(data = out, aes(x = REPCHL_1m_scl,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth()
      ggplot(data = out, aes(x = SST_clim_0m_scl,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth()
      ggplot(data = out, aes(x = CHL_clim_0m_scl,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth() + 
        scale_x_continuous(limits = c(-2,4))
      ggplot(data = out, aes(x = BOTTDEPTH,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth()
      ggplot(data = out, aes(x = vgpm_0m,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth()
      ggplot(data = out, aes(x = CTDst,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth()
      ggplot(data = out, aes(x = WP2st,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth()
      ggplot(data = out, aes(x = VESSEL,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_boxplot() + geom_hline(yintercept=0, col="red", linetype=2)
      ggplot(data = out, aes(x = X1000,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth()
      ggplot(data = out, aes(x = Y1000,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_point() + geom_smooth()
      ggplot(data = out, aes(x = STTYPE,  y = residuals)) + facet_wrap(.~age) + theme_bw()+ geom_boxplot() + geom_hline(yintercept=0, col="red", linetype=2)
      
      