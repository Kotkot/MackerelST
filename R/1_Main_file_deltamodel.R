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
library(sdmTMBextra)						
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
                formula_pres = formula(CPUE ~ 1),
                # formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_0m_scl,k=3) + s(REPNRTCHL_0m_scl, k=3) + s(OMLT_0m_scl, k=3) + (1|VESSEL)),
                formula = formula(CPUE ~ -1 + as.factor(YEAR) + s(SST_0m_scl,k=3) + s(OMLT_0m_scl, k=3) + (1|VESSEL)),
                formula_mix = formula(CPUE ~ 1),
                mesh_type = "cutoff", # or "kmeans" "cutoff", or "own" (which modifies to the mackerel case)
                knots_mesh = 250,
								bias_correct = FALSE								
                )

## Set up data & Define parameters
datdat = Prepare_data(bio=bio, catch=dat, conf=conf, rerun_stack=FALSE)

saveRDS(datdat, file="delta_data.rds")
saveRDS(conf, file="delta_conf.rds")

## Fit model
startTime = Sys.time()
run = fit_mackerelST_delta(data=datdat, conf, conf_extra = NULL)
endTime = Sys.time()
endTime-startTime

betas <- matrix(run$opt$par[grep("beta", names(run$opt$par))], nrow=ncol(datdat$tmb_data$X), ncol= datdat$tmb_data$Nage, byrow=F)  

## Extract output
gg_control <- theme(axis.title = element_text(size=15),
                    axis.text = element_text(size=14), 
                    plot.title = element_text(hjust=0.5, size=18),
                    strip.text = element_text(size=14))  

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
                        
       