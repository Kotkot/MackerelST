##' Prepare/filter data for the TMB model run
##' @param bio is the biological data frame
##' @param catch is the catch data frame
##' @param conf  configurations
##' @param rerun_stack  force the rerun of the code to create the prediction grid 
##' @details Prepare and modifies the original data files to be able to run the spatio-temporal model 
##' @return Data to be provided to TMB
##' @export
##' 


Prepare_data = function(bio, catch, conf, rerun_stack, catchability_var=NULL) {
  
  ## Step 1: defining yobs (the CPUE at each station)
  ## Dividing the aggregated CPUE (kg/km2) to each age based on the sampled age comps weight (assuming random sample - only VALID outside North Sea)
  ## i.e. CPUE(a) at station = sum of the weight of individual in age a in the station   DIVIDED BY   the total sum of all aged individuals in the station
  ## If the aggregated CPUE = 0, then CPUE(a)=0 for all ages
  ## If no age samples but CPUE, flag it as NO_data to that the model handles it differently later on
    bio <- bio %>% mutate(AGEgroup = ifelse(AGE >= conf$plusgroup, conf$plusgroup, AGE)) %>% 
      filter(!is.na(AGEgroup), YEAR %in% conf$years)
    catch <- catch %>% filter(YEAR %in% conf$years) %>% drop_na(as.character(attributes(conf)$fixed_effect))
    catch$ID <- as.factor(catch$ID)
    
    bio_dat <- bio %>% mutate(AGEgroup = as.factor(AGEgroup),
                              ID = factor(ID, levels=levels(catch$ID))) %>% 
      group_by(ID, AGEgroup, .drop=FALSE) %>% summarize(sumweight=sum(WEIGHT,na.rm=T),n=n()) %>% mutate(prop_weight = prop.table(sumweight)) %>% 
      left_join(subset(catch, select=c("ID", "CPUE"))) %>% mutate(CPUE_age = CPUE * prop_weight) %>% dplyr::select(ID, AGEgroup, CPUE, CPUE_age, n, sumweight, prop_weight) %>% 
      mutate(CPUE_agenew = ifelse(is.na(prop_weight),NA, CPUE_age)) %>%
      mutate(CPUE_agenew = ifelse(CPUE == 0, 0, CPUE_agenew)) %>% 
      mutate(NO_data = ifelse(is.na(CPUE_agenew), 1, 0)) %>% filter(AGEgroup %in% conf$ages)
    
    # To check that I have done correctly the allocation of age info to ALL the survey stations that are included in the selected years
    # length(unique(bio_dat$ID)); length(levels(catch$ID))
    
    ydata <- bio_dat %>% dplyr::select(ID, AGEgroup, CPUE_agenew) %>% spread(key = AGEgroup, value = CPUE_agenew)
    data <- catch %>% distinct() %>% right_join(ydata) %>% drop_na(as.character(attributes(conf)$fixed_effect))

  ## Add UTM coordinate to the data 
    survey <- data %>% st_as_sf(crs = 4326, coords = c("LON", "LAT"))  %>%
      st_transform(projection_km)
    surv_utm_coords <- st_coordinates(survey)
    data$X1000 <- surv_utm_coords[,1] #/ 1000
    data$Y1000 <- surv_utm_coords[,2] #/ 1000

  ## Filtering the study area and year (i.e. excluding the NS) 
    data <- mutate(data, area = ifelse(LAT < 60 & LON > -10, "northsea", "other"))
    data <- data %>% filter(area != "northsea") %>% as.data.frame()
    
  ## Defining the yobs by age
  grep_val = paste0("^", conf$ages[1], "$")  
  yobs <- data[, grep(grep_val, colnames(data))+c(1:(conf$Nage))-1] # age 1 data starts on column 8
  Nobs = nrow(data)
  Nyear <- length(unique(data$YEAR))
  
  ## deriving the fixed effect design matrix for the presence-absence and the positive model
    mgcv_page <- NULL  # to initialize in case it is not used (because reported as output)
    mgcv_mod <- NULL  # to initialize in case it is not used (because reported as output)
    # presence-absence
    if (length(grep( "k =", attributes(conf)$fixed_effect_formula_pres)) == 0 )  X_pres <- model.matrix(attributes(conf)$fixed_effect_formula_pres, data)
    if (length(grep( "k =", attributes(conf)$fixed_effect_formula_pres)) > 0 )   {
      data$YEAR <- as.factor(data$YEAR) 
      mgcv_page <- mgcv::gam(attributes(conf)$fixed_effect_formula_pres, data = data)
      X_pres <- model.matrix(mgcv_page)
    }
    # positive model
    if (length(grep( "k =", attributes(conf)$fixed_effect_formula)) == 0 )  X <- model.matrix(attributes(conf)$fixed_effect_formula, data)
    if (length(grep( "k =", attributes(conf)$fixed_effect_formula)) > 0 )   {
      data$YEAR <- as.factor(data$YEAR) 
      mgcv_mod <- mgcv::gam(attributes(conf)$fixed_effect_formula, data = data)
      X <- model.matrix(mgcv_mod)
    }
  
  ## Now deriving the parameters needed for the random effect estimation
  RE_indexes <- as.matrix(rep(0L, nrow(data)))
  nobs_RE <- 0L
  ln_tau_G_index <- unlist(lapply(seq_along(nobs_RE), function(i) rep(i, each = nobs_RE[i]))) - 1L
  
  if (FALSE %in% is.na(attributes(conf)$RE_effects)) {
    RE_indexes <- sapply(attributes(conf)$RE_effects, function(x) as.numeric(factor(data[[x]], labels=1:length(unique(data[[x]]))))) - 1L
    nobs_RE <- unname(apply(RE_indexes, 2L, max)) + 1L
    if (length(nobs_RE) == 0L) nobs_RE <- 0L
    ln_tau_G_index <- unlist(lapply(seq_along(nobs_RE), function(i) rep(i, each = nobs_RE[i]))) - 1L
  }
  
  ## Doing similar thing for the mixture effect
  if (length(grep( "k =", attributes(conf)$mixture_formula)) == 0 )  X_mix <- model.matrix(attributes(conf)$mixture_formula, data)
  if (length(grep( "k =", attributes(conf)$mixture_formula)) > 0 )   {
    data$YEAR <- as.factor(data$YEAR) 
    mgcv_mod_mix <- mgcv::gam(attributes(conf)$mixture_formula, data = data)
    X_mix <- model.matrix(mgcv_mod_mix)
  }
  
  
  # Stuff needed for spatiotemporal A matrix. 
  # Construct our mesh:
  # simple mesh
  if (conf$mesh_type == "kmeans") spde <- sdmTMB::make_mesh(data, xy_cols = c("X1000", "Y1000"),
                            n_knots = conf$knots_mesh, type = "kmeans", refine = list(min.angle = 25, max.edge = Inf, max.n.strict = -1, max.n = 1000))
  
  if (conf$mesh_type == "cutoff") spde <- make_mesh(data, xy_cols = c("X1000", "Y1000"), n_knots = conf$knots_mesh, type = "cutoff_search")
  
	if (conf$mesh_type == "own") {
		spde <- make_mesh(data, xy_cols = c("X1000", "Y1000"), n_knots = 400, type = "cutoff_search")
		mesh_temp = INLA::inla.mesh.2d(rbind(spde$mesh$loc), max.edge = c(20000, 50000))
		plot(mesh_temp)
		spde <- make_mesh(data,  xy_cols = c("X1000", "Y1000"), mesh = mesh_temp)
		plot(spde$mesh)
		points(data[,c("X1000", "Y1000")])
  }
  
  # more refined mesh
  # boundary <- list(
  #   INLA::inla.nonconvex.hull(coordinates(all_strata_sf %>% sf_to_df() %>% dplyr::select(x,y)), 200),
  #   INLA::inla.nonconvex.hull(coordinates(all_strata_sf %>% sf_to_df() %>% dplyr::select(x,y)), 500))
  # 
  # mesh <- INLA::inla.mesh.2d(loc=data[,c("X1000", "Y1000")], boundary = boundary,  max.edge = c(150, 450), cutoff = 150, min.angle = 25)
  # plot(mesh)
  # points(data[,c("X1000", "Y1000")])
  # # 
  # spde <- make_mesh(data,  xy_cols = c("X1000", "Y1000"), mesh = mesh)
  # plot(spde$mesh)
  # points(data[,c("X1000", "Y1000")])
  
  
  # create the barrier mesh   
  Barrier_range_fraction = 0.1  
  
  bspde <- add_barrier_mesh(
    spde, Atlantic_proj, range_fraction = Barrier_range_fraction,
    proj_scaling = 1, plot = TRUE)


  # to modify the mesh structure to account for the barrier effect that is changing annually (both land and area with zero fish)
    bla <- sf_to_df(Atlantic_proj)
    no_fish_zone <- list()
    for (yr in 1:conf$Nyear) {
      bla1 <- st_transform(Exclude_by_year[[yr]], crs=projection_km) %>% sf_to_df() 
      bla1$multipolygon_id <- 2
      bla1$polygon_id <- bla1$polygon_id+1000
      bla1 <- bla1[, match(colnames(bla), colnames(bla1))]
      bla2 <- rbind(bla, bla1)
      bla2 <- bla2 %>% st_as_sf(coords = c("x", "y")) %>% st_set_crs(projection_km) %>% 
        group_by(polygon_id) %>% 
        summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
      if (conf$plotting == TRUE) ggplot(Atlantic_proj) + geom_sf() + theme_bw() + geom_sf(data = bla2, size = 0.2, fill="lightblue1")     
      no_fish_zone[[yr]] <- add_barrier_mesh(spde, bla2, range_fraction = Barrier_range_fraction,
                                             proj_scaling = 1000, plot = plotting)
    }
   
  
   
  data$sdm_orig_id <- seq(1, nrow(data))
  data$sdm_x <- spde$loc_xy[,1,drop=TRUE]
  data$sdm_y <- spde$loc_xy[,2,drop=TRUE]
  fake_data <- unique(data.frame(sdm_x = data$sdm_x, sdm_y = data$sdm_y))
  fake_data[["sdm_spatial_id"]] <- seq(1, nrow(fake_data))
  data <- base::merge(data, fake_data, by = c("sdm_x", "sdm_y"),
                      all.x = TRUE, all.y = FALSE)
  data <- data[order(data$sdm_orig_id),, drop=FALSE]

  
  # now calculating the empirical correlation 
    corr_list <- list()
    for (yr in seq_along(conf$years)) {
      corr_list[[yr]] <- cor(subset(data, YEAR == conf$years[yr])[,grep(grep_val, colnames(data))+c(1:(conf$Nage))-1], use = "complete.obs")
    }
    Corr_all <- do.call(rbind, corr_list)

  # If using only the same barrier across years
		A <- spde$A
		A_st <- INLA::inla.spde.make.A(bspde$mesh,
																	 loc = as.matrix(fake_data[, c("sdm_x", "sdm_y"), drop = FALSE]))
		
		n_s <- nrow(spde$mesh$loc)
		Nmesh <- n_s
  
  ### Now prepare the input data for the prediction (based on 10 x 10 km grid?) 
  ## Reading in the covariate information
  if (rerun_stack == TRUE | !exists("pred_grid_or")){    
    # env_pred <- readRDS(paste0(getwd(), "/data/new_data/stacks_per_year_with_clim_and_distance_from_SSTfronts for July at 3035.rds"))   # no mixed layer depth
    # env_pred <- readRDS(paste0(getwd(), "/data/new_data/stacks_per_year_with_clim_and_distance_from_SSTfronts for July at 3035_new.rds"))
    env_pred <- readRDS(paste0(getwd(), "/data/new_data/stacks_per_year_with_clim_and_distance_from_SSTfronts for July at 30352023-01-21.rds"))
    
    # for all env variable except front
    qqq <- extract(env_pred[["2019"]],c(1:ncell(env_pred[["2019"]])))
    ydf <- as.data.frame(qqq)
    head(ydf)
    coord <- coordinates(env_pred[["2019"]])
    ydf <- as.data.frame(cbind(ydf, coord))
    ydf_sf <- ydf %>% st_as_sf(coords = c("x","y"), crs=projection)  %>%
      st_transform(projection_km)
    ydf_sf$X1000 <- st_coordinates(ydf_sf)[,1]#/1000
    ydf_sf$Y1000 <- st_coordinates(ydf_sf)[,2]#/1000
    ydf_sf$sdm_spatial_id <- 1:nrow(ydf_sf)
    ydf_sf_cropped <- st_intersection(ydf_sf,all_strata_sf_union) # limit the projection area
    
    kept_ID <- ydf_sf_cropped$sdm_spatial_id  
    
    pred_grid <- c()
    for (yr in levels(data$YEAR)) {
      
      qqq <- extract(env_pred[[yr]],c(1:ncell(env_pred[[yr]])))
      ydf <- as.data.frame(qqq)
      head(ydf)
      coord <- coordinates(env_pred[[yr]])
      ydf <- as.data.frame(cbind(ydf, coord))
      ydf_sf <- ydf %>% st_as_sf(coords = c("x","y"), crs=projection)  %>%
        st_transform(projection_km)
      ydf_sf$X1000 <- st_coordinates(ydf_sf)[,1]#/1000
      ydf_sf$Y1000 <- st_coordinates(ydf_sf)[,2]#/1000
      ydf_sf <- ydf_sf %>% st_transform(crs=4326)
      ydf_sf$LON <- st_coordinates(ydf_sf)[,1]
      ydf_sf$LAT <- st_coordinates(ydf_sf)[,2]
      ydf_sf$sdm_spatial_id <- 1:nrow(ydf_sf)
      
      ### Now transforming all continuous variables
      ydf_sf1 <- ydf_sf %>% as_tibble() %>% mutate(across(where(is.double), ~ scale(.x))) 
      ydf_sf1 <- ydf_sf1 %>% dplyr::select(where(is.double)) 
      colnames(ydf_sf1) <- paste0(colnames(ydf_sf1), "_scl")
      ydf_sf <- cbind(ydf_sf,ydf_sf1)
      
      ydf_cropped <- ydf_sf[kept_ID,] %>% 
        st_drop_geometry() %>%
        mutate(YEAR = as.numeric(yr)) 
      
      pred_grid <- bind_rows(pred_grid, ydf_cropped)                                          
    }
    
    pred_grid$ID <- 1:nrow(pred_grid)
 
    ## Save the pred_grid to the global environment
    pred_grid_or <- pred_grid
    assign("pred_grid_or", pred_grid_or, envir = .GlobalEnv)
    
    ## remove large object
    rm(env_pred, pred_grid)
  }

  pred_grid_or$SSTfront_scl <- pred_grid_or$distance_SST_fronts_ostia_scl
  # pred_grid_or$YEAR_fct <- as.factor(pred_grid_or$YEAR)
  pred_grid_or$VESSEL <- "TFNA"
  
	# to add the YearArea variable
	if (length(grep("YearArea", colnames(pred_grid_or)))==0 & ("Area" %in% attr(conf, "fixed_effect"))){
		bla <- st_distance(polygons_df, st_as_sf(x=pred_grid_or, coords = c("LON", "LAT")) %>% st_set_crs(4326) %>% st_transform(projection_km) %>% st_cast("POINT"))
		if ("Area" %in% attr(conf, "fixed_effect")){
			pred_grid_or$Area = apply(bla, 2, which.min)
			pred_grid_or$YearArea = apply(pred_grid_or[,c('YEAR', 'Area')], 1, function(x) paste(x, collapse="_"))
		}  
  }
	
  ## Remove any locations (doing this for all years) with NAS --> ideally, there should be none...

  to_exclude_NA <- c()
  if (length(attr(conf, "fixed_effect"))>1){
    for (i in 2:length(attr(conf, "fixed_effect"))) {
      asd <- which(is.na(pred_grid_or[,as.character(attr(conf, "fixed_effect")[i])]) ==TRUE)
      to_exclude_NA <- union(to_exclude_NA, asd)
    }
  }
  to_exclude_NA <- sort(to_exclude_NA)
  pred_grid <- pred_grid_or %>% filter(! (sdm_spatial_id %in% unique(pred_grid_or[to_exclude_NA,'sdm_spatial_id'])))
  pred_grid$CPUE = 1
  
  nrow(pred_grid)
  
  ## Filter to years where we have data for (no projection)
  pred_grid <- pred_grid %>% filter(YEAR %in% conf$years)
  
  ## Now replace all the values of variables that are catchability by 0 when doing prediction 
  if (!is.null(catchability_var)) {
    pred_grid[,match(catchability_var, colnames(pred_grid))] <- 0
  }
  
  # design matrix for prediction for the presence-absence
  if (length(grep( "k =", attributes(conf)$fixed_effect_formula_pres)) == 0 )  X_proj_pres <- model.matrix(attributes(conf)$fixed_effect_formula_pres, pred_grid)
  if (length(grep( "k =", attributes(conf)$fixed_effect_formula_pres)) > 0 )   X_proj_pres = mgcv::predict.gam(mgcv_page, type = "lpmatrix", newdata = pred_grid)
  # design matrix for prediction for the positive model
  if (length(grep( "k =", attributes(conf)$fixed_effect_formula)) == 0 )  X_proj <- model.matrix(attributes(conf)$fixed_effect_formula, pred_grid)
  if (length(grep( "k =", attributes(conf)$fixed_effect_formula)) > 0 )   X_proj = mgcv::predict.gam(mgcv_mod, type = "lpmatrix", newdata = pred_grid)
  # design matrix for prediction for the mixture component
  if (length(grep( "k =", attributes(conf)$mixture_formula)) == 0 )  X_proj_mix <- model.matrix(attributes(conf)$mixture_formula, pred_grid)
  if (length(grep( "k =", attributes(conf)$mixture_formula)) > 0 )   X_proj_mix = mgcv::predict.gam(mgcv_mod_mix, type = "lpmatrix", newdata = pred_grid)
  
  proj_data <- subset(pred_grid, YEAR==2012)
  
  # keeping the random effect of interest (when involving year index or area)
  if ("YearArea" %in% attributes(conf)$RE_effects){
    where_index <- grep("YearArea", attributes(conf)$RE_effects)
    YearArea_index <- data.frame(index=RE_indexes[,where_index], val=data[["YearArea"]])
    YearArea_index <- YearArea_index[!duplicated(YearArea_index),]
    RE_indexes_proj <- data.frame(val = pred_grid[["YearArea"]])
    RE_indexes_proj <- RE_indexes_proj %>% left_join(YearArea_index)
    RE_indexes_proj <- RE_indexes_proj[,"index"] %>% replace_na(9999)
  }
  
  # create the A matrix: it is always the same because the underlying mesh does not change
  A_proj<- INLA::inla.spde.make.A(bspde$mesh,
                                  loc = as.matrix(proj_data[, c("X1000", "Y1000"), drop = FALSE]))
  
  # Npred
  Npred = nrow(A_proj)
  
  
  ### Calculating the location where you do not want to project (because of no-fish-zone)
  proj_data_NA_sf <- proj_data %>% mutate(ID = 1:nrow(proj_data)) %>% 
    st_as_sf(coords = c("X1000","Y1000"), crs=projection_km) 
  
  proj_NA <- matrix(1, nrow=nrow(proj_data), ncol=Nyear)
  for (ii in 1:Nyear) {
    where_NA <- st_intersection(proj_data_NA_sf,Exclude_by_year_simple[[ii]]) # limit the projection area
    NAs <- where_NA$ID
    proj_NA[NAs, ii] <- NA
  }  
  
  
  
  Nage  = conf$Nage 
  Nrowyobs = ifelse(conf$Nage == 1, length(yobs), nrow(yobs))
  Ncolyobs = ifelse(conf$Nage == 1, 1, ncol(yobs))
  temp_val= 0
  spatialmodel = 0 
  if (conf$keep_omega == TRUE | conf$keep_epsilon == TRUE) spatialmodel = 1
  
  if ("YearArea" %in% attributes(conf)$RE_effects) temp_val = RE_indexes_proj
  
  ### Starting model building 
  # tmb data specs
  if (conf$Model_type == "Normal" | conf$include_age_correlation == "with_epsilon") {
    tmb_data <- list(
      Nage             = Nage,
      X                = as.matrix(X),
      yobs             = as.matrix(yobs),
      Nobs             = Nobs,
      spde_barrier     = make_barrier_spde(bspde),
      barrier_scaling  = bspde$barrier_scaling,
      Aobs             = A,
      Ast              = A_st,
      A_spatial_index  = data$sdm_spatial_id - 1L,  
      year_i           = make_year_i(data$YEAR),
      Nyear            = Nyear,
      Nmesh            = Nmesh,
      Npred            = Npred,
      do_predict       = conf$Do_predict,
      calc_se					 = 0,
      X_proj					 = X_proj,
      A_proj					 = A_proj,
      family           = 0,            # does not matter right now because only the tweedie case is implemented
      link             = 1
    )
  }

  if (conf$Model_type == "Annual_barrier" & conf$mixture_model == 2) {
      tmb_data <- list(
        Nage             = Nage,
        X                = as.matrix(X),
        X_mix            = as.matrix(X_mix),
        yobs             = as.matrix(yobs),
        Nobs             = Nobs,
        spde_barrier2010     = make_barrier_spde(no_fish_zone[[1]]),   # between 2010-2020, 2011 has the smallest coverage
        spde_barrier2011     = make_barrier_spde(no_fish_zone[[2]]),   # between 2010-2020, 2011 has the smallest coverage
        spde_barrier2012     = make_barrier_spde(no_fish_zone[[3]]),   # between 2010-2020, 2011 has the smallest coverage
        spde_barrier2013     = make_barrier_spde(no_fish_zone[[4]]),   # between 2010-2020, 2011 has the smallest coverage
        spde_barrier2014     = make_barrier_spde(no_fish_zone[[5]]),   # between 2010-2020, 2011 has the smallest coverage
        spde_barrier2015     = make_barrier_spde(no_fish_zone[[6]]),   # between 2010-2020, 2011 has the smallest coverage
        spde_barrier2016     = make_barrier_spde(no_fish_zone[[7]]),   # between 2010-2020, 2011 has the smallest coverage
        spde_barrier2017     = make_barrier_spde(no_fish_zone[[8]]),   # between 2010-2020, 2011 has the smallest coverage
        spde_barrier2018     = make_barrier_spde(no_fish_zone[[9]]),   # between 2010-2020, 2011 has the smallest coverage
        spde_barrier2019     = make_barrier_spde(no_fish_zone[[10]]),  # between 2010-2020, 2011 has the smallest coverage
        spde_barrier2020     = make_barrier_spde(no_fish_zone[[11]]),  # between 2010-2020, 2011 has the smallest coverage
        barrier_scaling  = bspde$barrier_scaling,
        RE_indexes       = RE_indexes,
        nobs_RE          = nobs_RE,
        ln_tau_G_index   = ln_tau_G_index,
        Aobs             = A,
        Ast              = A_st,
        A_spatial_index  = data$sdm_spatial_id - 1L,  
        year_i           = make_year_i(data$YEAR),
        Nyear            = Nyear,
        Nmesh            = Nmesh,
        Npred            = Npred,
        do_predict       = conf$Do_predict,
        corr_str         = conf$corr_str,
        include_dd       = as.integer(conf$density_dependence),
        calc_se					 = 0,
        X_proj					 = X_proj,
        X_proj_mix			 = X_proj_mix,
        A_proj					 = A_proj,
        proj_NA					 = proj_NA,
        to_keep          = matrix(1, nrow=Nrowyobs, ncol=Ncolyobs), 
				family           = conf$family,            # does not matter right now because only the tweedie case is implemented
        link             = 1,
        sim              = 1
      )
    }
    
  if (conf$Model_type == "Annual_barrier" & conf$mixture_model == 3) {
    tmb_data <- list(
      Nage             = Nage,
      X_pres           = as.matrix(X_pres),
      X                = as.matrix(X),
      yobs             = as.matrix(yobs),
      Nobs             = Nobs,
      spde_barrier     = make_barrier_spde(bspde),  # between 2010-2020, 2011 has the smallest coverage
      barrier_scaling  = bspde$barrier_scaling,
      RE_indexes       = RE_indexes,
      nobs_RE          = nobs_RE,
      ln_tau_G_index   = ln_tau_G_index,
      COR              = Corr_all,
      Aobs             = A,
      Ast              = A_st,
      A_spatial_index  = data$sdm_spatial_id - 1L,  
      year_i           = make_year_i(data$YEAR),
      Nyear            = Nyear,
      Nmesh            = Nmesh,
      Npred            = Npred,
      do_predict       = conf$Do_predict,
      corr_str         = conf$corr_str,
      include_dd       = as.integer(conf$density_dependence),
      calc_se					 = 0,
      X_proj_pres      = X_proj_pres,
      X_proj					 = X_proj,
      A_proj					 = A_proj,
      proj_NA					 = proj_NA,
      to_keep          = matrix(1, nrow=Nrowyobs, ncol=Ncolyobs), 
      family           = conf$family,            # does not matter right now because only the tweedie case is implemented
      link             = 1,
      df               = 4,
      ARorIID          = conf$ARorIID,
      sim              = 1
    )
    }
    
  if (conf$Model_type == "Annual_barrier" & conf$mixture_model == 0) {
    tmb_data <- list(
      Nage             = Nage,
      X                = as.matrix(X),
      # yobs_vec         = as.numeric(unlist(yobs)),     # this is by row
      # yobs_index       = matrix(1:length(unlist(yobs)), nrow=Nobs, ncol=Nage)-1,     # this is by row
      yobs             = as.matrix(yobs),     # this is by row
      yobs_vec         = as.vector(as.matrix(yobs)),     # this is by row
      Nobs             = Nobs,
      spatialmodel     = spatialmodel,
      spde_barrier     = make_barrier_spde(bspde),  # between 2010-2020, 2011 has the smallest coverage
      barrier_scaling  = bspde$barrier_scaling,
      RE_indexes       = RE_indexes,
      nobs_RE          = nobs_RE,
      ln_tau_G_index   = ln_tau_G_index,
      COR              = Corr_all,
      Aobs             = A,
      Ast              = A_st,
      A_spatial_index  = data$sdm_spatial_id - 1L,  
      year_i           = make_year_i(data$YEAR),
      Nyear            = Nyear,
      Nmesh            = Nmesh,
      Npred            = Npred,
      do_predict       = conf$Do_predict,
      corr_str         = conf$corr_str,
      add_nugget       = conf$add_nugget,
      include_dd       = as.integer(conf$density_dependence),
      calc_se					 = 0,
      X_proj					 = X_proj,
      A_proj					 = A_proj,
      RE_indexes_proj  = temp_val,
      proj_NA					 = proj_NA,
      to_keep          = matrix(1, nrow=Nrowyobs, ncol=Ncolyobs), 
      family           = conf$family,            # does not matter right now because only the tweedie case is implemented
      link             = conf$link,     # 1 is log, 4 is boxcox
      df 							 = 3, 
      ARorIID          = conf$ARorIID,
      sim              = 1
    )
  }
  
  # tmb param specs    
  val <- 2
  if(conf$corr_str == 1) val <- rep(0, Nyear)
  if(conf$corr_str == 2) val <- rep(0, Nage)
  if(conf$corr_str == 3) val <- rep(0, Nage)
  tmb_phases <- list()
  
  if (conf$Model_type == "Normal") {
    tmb_params <- list(
      beta         = matrix(0, nrow=ncol(X), ncol=Nage),
      omega        = matrix(0, nrow=Nmesh, ncol=Nage),
      epsilon_st   = array(0, dim=c(Nmesh, Nyear, Nage)),
      epsilon_cov  = array(0, dim=c(Nmesh, Nyear, Nage)),
      transf_rho   = 0.0,
      logKappa     = rep(0, Nage),
      logKappa_age = 0,
      logTauO      = rep(0, Nage),
      logTauE      = rep(0.1, Nage),
      logTau_age   = 0,
      thetaf       = rep(0.1, Nage), 
      ln_phi       = rep(0.1, Nage)    
    )  
  }    
  if (conf$Model_type == "Annual_barrier") {
     if (conf$include_age_correlation == "allinone" & conf$mixture_model %in% c(0,2) & conf$density_dependence == FALSE){
       tmb_params <- list(
        beta          = matrix(1, nrow=ncol(X), ncol=Nage),
        omega         = matrix(0, nrow=Nmesh, ncol=Nage),
        epsilon_st    = array(0, dim=c(Nmesh, Nyear, Nage)),
        transf_rho    = val,
        transf_rho_age= 10,
        logKappa      = val,
        logTauO       = rep(0, Nage),
        logTauE       = val,
        logsds        = rep(0, Nage),
        thetaf        = rep(0.1, Nage), 
        ln_phi        = rep(3, Nage),
        ln_tau_G      = matrix(0.1, nrow=length(nobs_RE), ncol=Nage),
        RE            = matrix(0, nrow=sum(nobs_RE), ncol=Nage),
        s50           = 0,
        logslope      =-5
      )  
    }
     if (conf$include_age_correlation == "allinone" & conf$mixture_model ==0 & conf$density_dependence == TRUE){
       tmb_params <- list(
        beta          = matrix(1, nrow=ncol(X), ncol=Nage),
        teta_svc      = matrix(0, nrow=Nmesh, ncol=Nage),
        omega         = matrix(0, nrow=Nmesh, ncol=Nage),
        epsilon_st    = array(0, dim=c(Nmesh, Nyear, Nage)),
        transf_rho    = val,
        transf_rho_age= 10,
        logKappa      = val,
        logTauSVC     = 0,
        logTauO       = rep(0, Nage),
        logTauE       = val,
        logsds        = rep(0, Nage),
        thetaf        = rep(0.1, Nage), 
        ln_phi        = rep(0.1, Nage),
        ln_tau_G      = matrix(0.1, nrow=length(nobs_RE), ncol=Nage),
        RE            = matrix(0, nrow=sum(nobs_RE), ncol=Nage),
        s50           = 0,
        logslope      =-5
      )  
       tmb_phases <- list(
        beta          = 1,
        teta_svc      = 2,
        omega         = 100,
        epsilon_st    = 3,
        transf_rho    = 100,
        transf_rho_age= 100,
        logKappa      = 2,
        logTauSVC     = 2,
        logTauO       = 100,
        logTauE       = 3,
        logsds        = 3,
        thetaf        = 1, 
        ln_phi        = 1,
        ln_tau_G      = 1,
        RE            = 1,
        s50           = 100,
        logslope      = 100
      )  
    }
     if (conf$include_age_correlation == "allinone" & conf$mixture_model == 1){
      tmb_params <- list(
        beta          = matrix(0, nrow=ncol(X), ncol=Nage),
        omega         = matrix(0, nrow=Nmesh, ncol=Nage),
        epsilon_st    = array(0, dim=c(Nmesh, Nyear, Nage)),
        transf_rho    = val,
        transf_rho_age= 0,
        logKappa      = val,
        logTauO       = rep(0, Nage),
        logTauE       = val,
        ln_tau_G      = rep(0.1, length(nobs_RE)), #change to matrix(0.1, nrow=length(nobs_RE), ncol=Nage)
        RE            = matrix(0, nrow=sum(nobs_RE), ncol=Nage),
        s50           = 0,
        logslope      =-5,
        # beta1         = matrix(0, nrow=ncol(X), ncol=Nage),
        # omega1        = matrix(0, nrow=Nmesh, ncol=Nage),
        epsilon_st1   = array(0, dim=c(Nmesh, Nyear, Nage)),
        transf_rho1   = rep(0, Nyear),
        logKappa1     = rep(0, Nyear),
        # logTauO1      = rep(0, Nage),
        logTauE1      = rep(0.1, Nyear),
        # ln_tau_G1     = rep(0.1, length(nobs_RE)),
        # RE1           = matrix(0, nrow=sum(nobs_RE), ncol=Nage),
        # s501          = 0,
        # logslope1     = -5,
        beta_mix      = matrix(15, nrow=ncol(X_mix), ncol=Nage),
        thetaf        = rep(0.1, Nage), 
        ln_phi        = rep(0.1, Nage)
      )  
    }
     if (conf$include_age_correlation == "allinone" & conf$mixture_model == 3){
       tmb_params <- list(
         beta_absc          = matrix(0, nrow=ncol(X_pres), ncol=Nage),
         beta               = matrix(0, nrow=ncol(X), ncol=Nage),
         omega_absc         = matrix(0, nrow=Nmesh, ncol=Nage),
         omega              = matrix(0, nrow=Nmesh, ncol=Nage),
         epsilon_st_absc    = array(0, dim=c(Nmesh, Nyear, Nage)),
         epsilon_st         = array(0, dim=c(Nmesh, Nyear, Nage)),
         transf_rho_absc    = val,
         transf_rho         = val,
         transf_rho_age_absc= 10,
         transf_rho_age     = 10,
         logKappa_absc      = val,
         logKappa           = val,
         logTauO_absc       = rep(0, Nage),
         logTauO            = rep(0, Nage),
         logTauE_absc       = val,
         logTauE            = val,
         logsds             = rep(0, Nage),
         thetaf             = rep(0.1, Nage), 
         ln_phi             = rep(0.1, Nage),
         ln_tau_G_absc      = matrix(0.1, nrow=length(nobs_RE), ncol=Nage),
         ln_tau_G           = matrix(0.1, nrow=length(nobs_RE), ncol=Nage),
         RE_absc            = matrix(0, nrow=sum(nobs_RE), ncol=Nage),
         RE                 = matrix(0, nrow=sum(nobs_RE), ncol=Nage),
         s50_absc           = 0,
         s50                = 0,
         logslope_absc      =-5,
         logslope           =-5
       )  
    }
     if (!is.null(conf$cohort)){
      tmb_params <- list(
        init_year     = rep(0, Nyear), 
        init_age      = rep(0, Nage-1),
        RE_diff       = rep(0, (Nage-1)*(Nyear-1)), 
        logcohort_SD  = 0, 
        beta          = matrix(1, nrow=ncol(X), ncol=Nage),
        omega         = matrix(0, nrow=Nmesh, ncol=Nage),
        epsilon_st    = array(0, dim=c(Nmesh, Nyear, Nage)),
        transf_rho    = val,
        transf_rho_age= 10,
        logKappa      = val,
        logTauO       = rep(0, Nage),
        logTauE       = val,
        logsds        = rep(0, Nage),
        logTau_nugget = 0,
        transf_rho_nugget    = 0,
        nugget_effect = array(0, dim=c(Nobs, Nage)),
        thetaf        = rep(0.1, Nage), 
        ln_phi        = rep(3, Nage),
        ln_tau_G      = matrix(0.1, nrow=length(nobs_RE), ncol=Nage),
        RE            = matrix(0, nrow=sum(nobs_RE), ncol=Nage),
        s50           = 0,
        logslope      =-5
      )  
    }
  }
    output <- list()
    output$tmb_data <- tmb_data
    output$tmb_params <- tmb_params
    output$tmb_phases <- tmb_phases
    output$pred_grid <- pred_grid
    output$mgcv_page <- mgcv_page
    output$mgcv_mod <- mgcv_mod
    output$data <- data
    output$spde <- spde
    output$proj_NA <- proj_NA
    
    return(output)
    
}    



##' Prepare/filter data for the TMB model run for the model that combines 1. aggregated CPUE and 2. P(age)
##' @param bio is the biological data frame
##' @param catch is the catch data frame
##' @param conf  configurations
##' @param rerun_stack  force the rerun of the code to create the prediction grid 
##' @details Prepare and modifies the original data files to be able to run the spatio-temporal model 
##' @return Data to be provided to TMB
##' @export
##' 

Prepare_data_condlogit = function(bio, catch, conf, rerun_stack, catchability_var=NULL) {
  
  ## Step 1: defining yobs (the CPUE at each station)
  ## Dividing the aggregated CPUE (kg/km2) to each age based on the sampled age comps weight (assuming random sample - only VALID outside North Sea)
  ## i.e. CPUE(a) at station = sum of the weight of individual in age a in the station   DIVIDED BY   the total sum of all aged individuals in the station
  ## If the aggregated CPUE = 0, then CPUE(a)=0 for all ages
  ## If no age samples but CPUE, flag it as NO_data to that the model handles it differently later on
    bio <- bio %>% mutate(AGEgroup = ifelse(AGE >= conf$plusgroup, conf$plusgroup, AGE),
                          AGEgroup = ifelse(AGEgroup <= conf$minusgroup, conf$minusgroup, AGEgroup)) %>% 
      filter(!is.na(AGEgroup), YEAR %in% conf$years)
    catch <- catch %>% filter(YEAR %in% conf$years) %>% drop_na(as.character(attributes(conf)$fixed_effect))
    catch$ID <- as.factor(catch$ID)
    
    bio_dat <- bio %>% mutate(AGEgroup = as.factor(AGEgroup),
                              ID = factor(ID, levels=levels(catch$ID))) %>% 
      group_by(ID, AGEgroup, .drop=FALSE) %>% summarize(sumweight=sum(WEIGHT,na.rm=T),n=n()) %>% mutate(prop_weight = prop.table(sumweight)) %>% 
      left_join(subset(catch, select=c("ID", "CPUE"))) %>% mutate(CPUE_age = CPUE * prop_weight) %>% dplyr::select(ID, AGEgroup, CPUE, CPUE_age, n, sumweight, prop_weight) %>% 
      mutate(CPUE_agenew = ifelse(is.na(prop_weight),NA, CPUE_age)) %>%
      mutate(CPUE_agenew = ifelse(CPUE == 0, 0, CPUE_agenew)) %>% 
      mutate(NO_data = ifelse(is.na(CPUE_agenew), 1, 0)) %>% filter(AGEgroup %in% conf$ages)
    
    # To check that I have done correctly the allocation of age info to ALL the survey stations that are included in the selected years
    # length(unique(bio_dat$ID)); length(levels(catch$ID))
    
    ydata <- bio_dat %>% dplyr::select(ID, AGEgroup, prop_weight) %>% spread(key = AGEgroup, value = prop_weight)
    ydata1 <- bio_dat %>% dplyr::select(ID, AGEgroup, n) %>% group_by(ID) %>% summarize(n=sum(n))
    data <- catch %>% distinct() %>% right_join(ydata) %>% right_join(ydata1) %>% drop_na(as.character(attributes(conf)$fixed_effect))

  ## Add UTM coordinate to the data 
    survey <- data %>% st_as_sf(crs = 4326, coords = c("LON", "LAT"))  %>%
      st_transform(projection_km)
    surv_utm_coords <- st_coordinates(survey)
    data$X1000 <- surv_utm_coords[,1] #/ 1000
    data$Y1000 <- surv_utm_coords[,2] #/ 1000

  ## Filtering the study area and year (i.e. excluding the NS) 
    data <- mutate(data, area = ifelse(LAT < 60 & LON > -10, "northsea", "other"))
    data <- data %>% filter(area != "northsea") %>% as.data.frame()
    
  ## Defining the yobs by age
    grep_val = paste0("^", conf$ages[1], "$")  
    yobs <- data[, grep(grep_val, colnames(data))+c(1:(conf$Nage))-1] # age 1 data starts on column 8
    Nobs = nrow(data)
    Nyear <- length(unique(data$YEAR))
  
  ## deriving the fixed effect design matrix for the P(age) and the aggregated CPUE model
    mgcv_mod <- NULL  # to initialize in case it is not used (because reported as output)
    mgcv_mod_cpue <- NULL  # to initialize in case it is not used (because reported as output)
    # P(age)
    if (length(grep( "k =", attributes(conf)$fixed_effect_formula)) == 0 )  X <- model.matrix(attributes(conf)$fixed_effect_formula, data)
    if (length(grep( "k =", attributes(conf)$fixed_effect_formula)) > 0 )   {
      data$YEAR <- as.factor(data$YEAR) 
      mgcv_mod <- mgcv::gam(attributes(conf)$fixed_effect_formula, data = data)
      X <- model.matrix(mgcv_mod)
    }
    # aggregated CPUE model
    if (length(grep( "k =", attributes(conf)$fixed_effect_formula_cpue)) == 0 )  X_cpue <- model.matrix(attributes(conf)$fixed_effect_formula_cpue, data)
    if (length(grep( "k =", attributes(conf)$fixed_effect_formula_cpue)) > 0 )   {
      data$YEAR <- as.factor(data$YEAR) 
      mgcv_mod_cpue <- mgcv::gam(attributes(conf)$fixed_effect_formula_cpue, data = data)
      X_cpue <- model.matrix(mgcv_mod_cpue)
    }
  
  ## Now deriving the parameters needed for the random effect estimation
    RE_indexes <- as.matrix(rep(0L, nrow(data)))
    RE_indexes_cpue <- as.matrix(rep(0L, nrow(data)))
    nobs_RE <- 0L
    nobs_RE_cpue <- 0L
    ln_tau_G_index <- unlist(lapply(seq_along(nobs_RE), function(i) rep(i, each = nobs_RE[i]))) - 1L
    ln_tau_G_index_cpue <- unlist(lapply(seq_along(nobs_RE_cpue), function(i) rep(i, each = nobs_RE_cpue[i]))) - 1L
    
    if (FALSE %in% is.na(attributes(conf)$RE_effects)) {
      RE_indexes <- sapply(attributes(conf)$RE_effects, function(x) as.numeric(factor(data[[x]], labels=1:length(unique(data[[x]]))))) - 1L
      nobs_RE <- unname(apply(RE_indexes, 2L, max)) + 1L
      if (length(nobs_RE) == 0L) nobs_RE <- 0L
      ln_tau_G_index <- unlist(lapply(seq_along(nobs_RE), function(i) rep(i, each = nobs_RE[i]))) - 1L
    }
    if (FALSE %in% is.na(attributes(conf)$RE_effects_cpue)) {
      RE_indexes_cpue <- sapply(attributes(conf)$RE_effects_cpue, function(x) as.numeric(factor(data[[x]], labels=1:length(unique(data[[x]]))))) - 1L
      nobs_RE_cpue <- unname(apply(RE_indexes_cpue, 2L, max)) + 1L
      if (length(nobs_RE_cpue) == 0L) nobs_RE_cpue <- 0L
      ln_tau_G_index_cpue <- unlist(lapply(seq_along(nobs_RE_cpue), function(i) rep(i, each = nobs_RE_cpue[i]))) - 1L
    }
    
  
  # Stuff needed for spatiotemporal A matrix. 
  # Construct our mesh:
    if (conf$mesh_type == "kmeans") spde <- sdmTMB::make_mesh(data, xy_cols = c("X1000", "Y1000"),
                              n_knots = conf$knots_mesh, type = "kmeans", refine = list(min.angle = 25, max.edge = Inf, max.n.strict = -1, max.n = 1000))
    
    if (conf$mesh_type == "cutoff") spde <- make_mesh(data, xy_cols = c("X1000", "Y1000"), n_knots = conf$knots_mesh, type = "cutoff_search")
    
  	if (conf$mesh_type == "own") {
  		spde <- make_mesh(data, xy_cols = c("X1000", "Y1000"), n_knots = 400, type = "cutoff_search")
  		mesh_temp = INLA::inla.mesh.2d(rbind(spde$mesh$loc), max.edge = c(20000, 50000))
  		plot(mesh_temp)
  		spde <- make_mesh(data,  xy_cols = c("X1000", "Y1000"), mesh = mesh_temp)
  		plot(spde$mesh)
  		points(data[,c("X1000", "Y1000")])
    }
    
  # more refined mesh
    # boundary <- list(
    #   INLA::inla.nonconvex.hull(coordinates(all_strata_sf %>% sf_to_df() %>% dplyr::select(x,y)), 200),
    #   INLA::inla.nonconvex.hull(coordinates(all_strata_sf %>% sf_to_df() %>% dplyr::select(x,y)), 500))
    # 
    # mesh <- INLA::inla.mesh.2d(loc=data[,c("X1000", "Y1000")], boundary = boundary,  max.edge = c(150, 450), cutoff = 150, min.angle = 25)
    # plot(mesh)
    # points(data[,c("X1000", "Y1000")])
    # # 
    # spde <- make_mesh(data,  xy_cols = c("X1000", "Y1000"), mesh = mesh)
    # plot(spde$mesh)
    # points(data[,c("X1000", "Y1000")])
    
  
  # create the barrier mesh   
    Barrier_range_fraction = 0.1  
    
    bspde <- add_barrier_mesh(
      spde, Atlantic_proj, range_fraction = Barrier_range_fraction,
      proj_scaling = 1, plot = TRUE)
  

  # to modify the mesh structure to account for the barrier effect that is changing annually (both land and area with zero fish)
    bla <- sf_to_df(Atlantic_proj)
    no_fish_zone <- list()
    for (yr in 1:conf$Nyear) {
      bla1 <- st_transform(Exclude_by_year[[yr]], crs=projection_km) %>% sf_to_df() 
      bla1$multipolygon_id <- 2
      bla1$polygon_id <- bla1$polygon_id+1000
      bla1 <- bla1[, match(colnames(bla), colnames(bla1))]
      bla2 <- rbind(bla, bla1)
      bla2 <- bla2 %>% st_as_sf(coords = c("x", "y")) %>% st_set_crs(projection_km) %>% 
        group_by(polygon_id) %>% 
        summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
      if (conf$plotting == TRUE) ggplot(Atlantic_proj) + geom_sf() + theme_bw() + geom_sf(data = bla2, size = 0.2, fill="lightblue1")     
      no_fish_zone[[yr]] <- add_barrier_mesh(spde, bla2, range_fraction = Barrier_range_fraction,
                                             proj_scaling = 1000, plot = plotting)
    }
   
  
   
    data$sdm_orig_id <- seq(1, nrow(data))
    data$sdm_x <- spde$loc_xy[,1,drop=TRUE]
    data$sdm_y <- spde$loc_xy[,2,drop=TRUE]
    fake_data <- unique(data.frame(sdm_x = data$sdm_x, sdm_y = data$sdm_y))
    fake_data[["sdm_spatial_id"]] <- seq(1, nrow(fake_data))
    data <- base::merge(data, fake_data, by = c("sdm_x", "sdm_y"),
                        all.x = TRUE, all.y = FALSE)
    data <- data[order(data$sdm_orig_id),, drop=FALSE]

  
  # now calculating the empirical correlation 
    corr_list <- list()
    for (yr in seq_along(conf$years)) {
      corr_list[[yr]] <- cor(subset(data, YEAR == conf$years[yr])[,grep(grep_val, colnames(data))+c(1:(conf$Nage))-1], use = "complete.obs")
      # corr_list[[yr]] <- corr_list[[yr]][conf$ages, conf$ages]
    }
    Corr_all <- do.call(rbind, corr_list)
    # if (conf$Nage < 10) 
    # {
    #   Corr_all <- diag(conf$Nage)
    # }  

  # If using only the same barrier across years
  # if (conf$Model_type == "Normal") {
    A <- spde$A
    A_st <- INLA::inla.spde.make.A(bspde$mesh,
                                   loc = as.matrix(fake_data[, c("sdm_x", "sdm_y"), drop = FALSE]))
  
  # }
  # to modify the mesh structure to account for the barrier effect that is changing annually (both land and area with zero fish)
  # if (conf$Model_type == "Annual_barrier") {
  #   Make_spde_A <- function(yr) INLA::inla.spde.make.A(no_fish_zone[[yr]]$mesh, as.matrix(data[which(data$YEAR == (sort(unique(data$YEAR))[yr])), c("sdm_x", "sdm_y")]))
  #   A_st <- c()
  #   for (yr in seq_along(sort(unique(data$YEAR))))  A_st <- rbind(A_st, Make_spde_A(yr))
  # }
  # 
  n_s <- nrow(spde$mesh$loc)
  Nmesh <- n_s
  
  ### Now prepare the input data for the prediction (based on 10 x 10 km grid?) 
  ## Reading in the covariate information
  if (rerun_stack == TRUE | !exists("pred_grid_or")){    
    # env_pred <- readRDS(paste0(getwd(), "/data/new_data/stacks_per_year_with_clim_and_distance_from_SSTfronts for July at 3035.rds"))   # no mixed layer depth
    # env_pred <- readRDS(paste0(getwd(), "/data/new_data/stacks_per_year_with_clim_and_distance_from_SSTfronts for July at 3035_new.rds"))
    env_pred <- readRDS(paste0(getwd(), "/data/new_data/stacks_per_year_with_clim_and_distance_from_SSTfronts for July at 30352023-01-21.rds"))
    
    # for all env variable except front
    qqq <- extract(env_pred[["2019"]],c(1:ncell(env_pred[["2019"]])))
    ydf <- as.data.frame(qqq)
    head(ydf)
    coord <- coordinates(env_pred[["2019"]])
    ydf <- as.data.frame(cbind(ydf, coord))
    ydf_sf <- ydf %>% st_as_sf(coords = c("x","y"), crs=projection)  %>%
      st_transform(projection_km)
    ydf_sf$X1000 <- st_coordinates(ydf_sf)[,1]#/1000
    ydf_sf$Y1000 <- st_coordinates(ydf_sf)[,2]#/1000
    ydf_sf$sdm_spatial_id <- 1:nrow(ydf_sf)
    ydf_sf_cropped <- st_intersection(ydf_sf,all_strata_sf_union) # limit the projection area
    
    kept_ID <- ydf_sf_cropped$sdm_spatial_id  
    
    pred_grid <- c()
    for (yr in levels(data$YEAR)) {
      
      qqq <- extract(env_pred[[yr]],c(1:ncell(env_pred[[yr]])))
      ydf <- as.data.frame(qqq)
      head(ydf)
      coord <- coordinates(env_pred[[yr]])
      ydf <- as.data.frame(cbind(ydf, coord))
      ydf_sf <- ydf %>% st_as_sf(coords = c("x","y"), crs=projection)  %>%
        st_transform(projection_km)
      ydf_sf$X1000 <- st_coordinates(ydf_sf)[,1]#/1000
      ydf_sf$Y1000 <- st_coordinates(ydf_sf)[,2]#/1000
      ydf_sf <- ydf_sf %>% st_transform(crs=4326)
      ydf_sf$LON <- st_coordinates(ydf_sf)[,1]
      ydf_sf$LAT <- st_coordinates(ydf_sf)[,2]
      ydf_sf$sdm_spatial_id <- 1:nrow(ydf_sf)
      
      ### Now transforming all continuous variables
      ydf_sf1 <- ydf_sf %>% as_tibble() %>% mutate(across(where(is.double), ~ scale(.x))) 
      ydf_sf1 <- ydf_sf1 %>% dplyr::select(where(is.double)) 
      colnames(ydf_sf1) <- paste0(colnames(ydf_sf1), "_scl")
      ydf_sf <- cbind(ydf_sf,ydf_sf1)
      
      ydf_cropped <- ydf_sf[kept_ID,] %>% 
        st_drop_geometry() %>%
        mutate(YEAR = as.numeric(yr)) 
      
      pred_grid <- bind_rows(pred_grid, ydf_cropped)                                          
    }
    
    pred_grid$ID <- 1:nrow(pred_grid)
 
    ## Save the pred_grid to the global environment
    pred_grid_or <- pred_grid
    assign("pred_grid_or", pred_grid_or, envir = .GlobalEnv)
    
    ## remove large object
    rm(env_pred, pred_grid)
  }

  pred_grid_or$SSTfront_scl <- pred_grid_or$distance_SST_fronts_ostia_scl
  # pred_grid_or$YEAR_fct <- as.factor(pred_grid_or$YEAR)
  pred_grid_or$VESSEL <- "TFNA"
  
	# to add the YearArea variable
	if (length(grep("YearArea", colnames(pred_grid_or)))>0){
		bla <- st_distance(polygons_df, st_as_sf(x=pred_grid_or, coords = c("LON", "LAT")) %>% st_set_crs(4326) %>% st_transform(projection_km) %>% st_cast("POINT"))
		if ("Area" %in% attr(conf, "fixed_effect")){
			pred_grid_or$Area = apply(bla, 2, which.min)
			pred_grid_or$YearArea = apply(pred_grid_or[,c('YEAR', 'Area')], 1, function(x) paste(x, collapse="_"))
		}  
  }
	
  ## Remove any locations (doing this for all years) with NAS --> ideally, there should be none...

  to_exclude_NA <- c()
  if (length(attr(conf, "fixed_effect"))>1){
    for (i in 2:length(attr(conf, "fixed_effect"))) {
      asd <- which(is.na(pred_grid_or[,as.character(attr(conf, "fixed_effect")[i])]) ==TRUE)
      to_exclude_NA <- union(to_exclude_NA, asd)
    }
  }
  to_exclude_NA <- sort(to_exclude_NA)
  pred_grid <- pred_grid_or %>% filter(! (sdm_spatial_id %in% unique(pred_grid_or[to_exclude_NA,'sdm_spatial_id'])))
  pred_grid$CPUE = 1
  
  nrow(pred_grid)
  
  ## Filter to years where we have data for (no projection)
  pred_grid <- pred_grid %>% filter(YEAR %in% conf$years)
  
  ## Now replace all the values of variables that are catchability by 0 when doing prediction 
  if (!is.null(catchability_var)) {
    pred_grid[,match(catchability_var, colnames(pred_grid))] <- 0
  }
  
  # design matrix for prediction for the aggregated CPUE model
  if (length(grep( "k =", attributes(conf)$fixed_effect_formula_cpue)) == 0 )  X_proj_cpue <- model.matrix(attributes(conf)$fixed_effect_formula_cpue, pred_grid)
  if (length(grep( "k =", attributes(conf)$fixed_effect_formula_cpue)) > 0 )   X_proj_cpue = mgcv::predict.gam(mgcv_mod_cpue, type = "lpmatrix", newdata = pred_grid)
  # design matrix for prediction for the P(age) model
  if (length(grep( "k =", attributes(conf)$fixed_effect_formula)) == 0 )  X_proj <- model.matrix(attributes(conf)$fixed_effect_formula, pred_grid)
  if (length(grep( "k =", attributes(conf)$fixed_effect_formula)) > 0 )   X_proj = mgcv::predict.gam(mgcv_mod, type = "lpmatrix", newdata = pred_grid)

  proj_data <- subset(pred_grid, YEAR==2012)
  
  # keeping the random effect of interest (when involving year index or area)
  if ("YearArea" %in% attributes(conf)$RE_effects){
    where_index <- grep("YearArea", attributes(conf)$RE_effects)
    YearArea_index <- data.frame(index=RE_indexes[,where_index], val=data[["YearArea"]])
    YearArea_index <- YearArea_index[!duplicated(YearArea_index),]
    RE_indexes_proj <- data.frame(val = pred_grid[["YearArea"]])
    RE_indexes_proj <- RE_indexes_proj %>% left_join(YearArea_index)
    RE_indexes_proj <- RE_indexes_proj[,"index"] %>% replace_na(9999)
  }
  if ("YearArea" %in% attributes(conf)$RE_effects_cpue){
    where_index_cpue <- grep("YearArea", attributes(conf)$RE_effects_cpue)
    YearArea_index_cpue <- data.frame(index=RE_indexes_cpue[,where_index_cpue], val=data[["YearArea"]])
    YearArea_index_cpue <- YearArea_index[!duplicated(YearArea_index_cpue),]
    RE_indexes_proj_cpue <- data.frame(val = pred_grid[["YearArea"]])
    RE_indexes_proj_cpue <- RE_indexes_proj_cpue %>% left_join(YearArea_index_cpue)
    RE_indexes_proj_cpue <- RE_indexes_proj_cpue[,"index"] %>% replace_na(9999)
  }
  
  # create the A matrix: it is always the same because the underlying mesh does not change
  A_proj<- INLA::inla.spde.make.A(bspde$mesh,
                                  loc = as.matrix(proj_data[, c("X1000", "Y1000"), drop = FALSE]))
  
  # Npred
  Npred = nrow(A_proj)
  
  
  ### Calculating the location where you do not want to project (because of no-fish-zone)
  proj_data_NA_sf <- proj_data %>% mutate(ID = 1:nrow(proj_data)) %>% 
    st_as_sf(coords = c("X1000","Y1000"), crs=projection_km) 
  
  proj_NA <- matrix(1, nrow=nrow(proj_data), ncol=Nyear)
  for (ii in 1:Nyear) {
    where_NA <- st_intersection(proj_data_NA_sf,Exclude_by_year_simple[[ii]]) # limit the projection area
    NAs <- where_NA$ID
    proj_NA[NAs, ii] <- NA
  }  
  
  
  
  Nage  = conf$Nage 
  Nrowyobs = ifelse(conf$Nage == 1, length(yobs), nrow(yobs))
  Ncolyobs = ifelse(conf$Nage == 1, 1, ncol(yobs))
  temp_val= 0
  temp_val_cpue= 0
  spatialmodel = 0 
  if (conf$keep_omega == TRUE | conf$keep_epsilon == TRUE) spatialmodel = 1
  
  if ("YearArea" %in% attributes(conf)$RE_effects) temp_val = RE_indexes_proj
  if ("YearArea" %in% attributes(conf)$RE_effects_cpue) temp_val_cpue = RE_indexes_proj_cpue
  
  sumPobs = as.vector(apply(yobs, 1, sum))
  sumPobs[which(is.na(sumPobs)==TRUE)] = 0
  yobs_weight = rep(1, Nobs)
	yobs_weight[which(log(data$CPUE)>1)] <- log(data$CPUE)[which(log(data$CPUE)>1)]
	
  ### Starting model building 
  # tmb data specs

    tmb_data <- list(
      Nage                  = Nage,
      X                     = as.matrix(X),
      X_cpue                = as.matrix(X_cpue),
      yobs                  = as.matrix(yobs),     # this is by row
      sumPobs               = sumPobs,
      # yobs_weight           = as.numeric(as.vector(data$n)),
      yobs_weight           = yobs_weight,
			yobs_cpue             = data$CPUE,     # this is by row
      Nobs                  = Nobs,
      spatialmodel          = spatialmodel,
      spde_barrier          = make_barrier_spde(bspde),  # between 2010-2020, 2011 has the smallest coverage
      barrier_scaling       = bspde$barrier_scaling,
      RE_indexes            = RE_indexes,
      RE_indexes_cpue       = RE_indexes_cpue,
      nobs_RE               = nobs_RE,
      nobs_RE_cpue          = nobs_RE_cpue,
      ln_tau_G_index        = ln_tau_G_index,
      ln_tau_G_index_cpue   = ln_tau_G_index_cpue,
      COR                   = Corr_all,
      Aobs                  = A,
      Ast                   = A_st,
      A_spatial_index       = data$sdm_spatial_id - 1L,  
      year_i                = make_year_i(data$YEAR),
      Nyear                 = Nyear,
      Nmesh                 = Nmesh,
      Npred                 = Npred,
      do_predict            = conf$Do_predict,
      corr_str              = conf$corr_str,
      include_dd            = as.integer(conf$density_dependence),
      calc_se					      = 0,
      X_proj					      = X_proj,
      X_proj_cpue			     	= X_proj_cpue,
      A_proj					      = A_proj,
      RE_indexes_proj       = temp_val,
      RE_indexes_proj_cpue  = temp_val_cpue,
      proj_NA					      = proj_NA,
      to_keep               = rep(1, Nobs), 
      link                  = conf$link,     # 1 is log, 2 is logit, 4 is boxcox
      link_cpue             = conf$link_cpue,     # 1 is log, 4 is boxcox
      df 							      = 3, 
      ARorIID               = conf$ARorIID,
      ARorIID_cpue          = conf$ARorIID_cpue,
      sim                   = 1,
      boxcox_pow            = rep(0, Nage-1),                
      boxcox_pow_cpue       = 0            
    )

  # tmb param specs    
  val <- 2
  if(conf$corr_str == 1) val <- rep(-2, Nyear)
  if(conf$corr_str == 2) val <- rep(-2, Nage-1)
  if(conf$corr_str == 3) val <- rep(-2, Nage-1)
  tmb_phases <- list()
  

	if (conf$cohort %in% c(0,1)) {

		# tmb_params <- list(
		 # beta              = matrix(1, nrow=ncol(X), ncol=Nage-1),
		 # beta_cpue         = c(rep(5, Nyear), rep(0, ncol(X_cpue)-Nyear)),   # had to modify it manually because causing problems
		 # omega             = matrix(0, nrow=Nmesh, ncol=Nage-1),
		 # omega_cpue        = rep(0, Nmesh),
		 # epsilon_st        = array(0, dim=c(Nmesh, Nyear, Nage-1)),
		 # epsilon_st_cpue    = matrix(0, nrow=Nmesh, ncol=Nyear),
		 # transf_rho        = val,
		 # transf_rho_cpue    = 1,
		 # transf_rho_age    = 10,
		 # logKappa          = val,
		 # # logKappa_cpue     = 0,
		 # # logTauO           = rep(0, Nage),
		 # logTauE           = val,
		 # logTauE_cpue      = 0,
		 # logsds            = rep(0, Nage-1),
		 # thetaf_cpue       = 0.1, 
		 # ln_phi_cpue       = 5,
		 # ln_tau_G          = matrix(0.1, nrow=length(nobs_RE), ncol=Nage-1),
		 # ln_tau_G_cpue     = rep(0.1, nrow=length(nobs_RE_cpue)),
		 # RE                = matrix(0, nrow=max(sum(nobs_RE),1), ncol=Nage-1),
		 # RE_cpue           = rep(0, max(1,sum(nobs_RE_cpue))),
		 # s50               = 0,
		 # logslope          =-5
		# ) 
			
		tmb_params <- list(
		 # init_year         = runif(Nyear, -5, 5),
		 init_year         =  c(-1.5, -0.6, -2.8, -1.9, -0.6, -2.2, -3.2, -1.3, -6.4,-1.8, -2.4),
		 # init_year_mean    = 0.5, 
		 # init_year_SD      = 2, 
		 # init_age          = runif(Nage-2, -5, 5),
		 init_age          = c(0.6,  0.3, -0.4, -1.2, -1.5, -2.7 ,-7.0, -6.6),
		 # init_age_SD       = 0,
		 # RE_diff           = rep(0, (Nage-2)*(Nyear-1)), 
		 RE_diff           = runif((Nage-2)*(Nyear-1), -1, 1), 
		 logcohort_SD      = 0.2, 
		 beta              = matrix(5, nrow=ncol(X), ncol=Nage-1),
		 beta_cpue         = c(rnorm(Nyear,6,1), rep(0, ncol(X_cpue)-Nyear)),   # had to modify it manually because causing problems
		 omega             = matrix(0, nrow=Nmesh, ncol=Nage-1),
		 omega_cpue        = rep(0, Nmesh),
		 epsilon_st        = array(0, dim=c(Nmesh, Nyear, Nage-1)),
		 epsilon_st_cpue    = matrix(0, nrow=Nmesh, ncol=Nyear),
		 transf_rho        = rep(0, Nage-1),
		 transf_rho_cpue    = 1,
		 transf_rho_age    = 10,
		 logKappa          = rep(-4, Nage-1),
		 logKappa_cpue     = -4,
		 # logTauO           = rep(0, Nage),
		 logTauE           = rep(1, Nage-1),
		 logTauE_cpue      = 0,
		 logsds            = rep(0, Nage-1),
		 thetaf_cpue       = 0.7, 
		 ln_phi_cpue       = 3,
		 ln_tau_G          = matrix(0.1, nrow=length(nobs_RE), ncol=Nage-1),
		 ln_tau_G_cpue     = rep(0.1, nrow=length(nobs_RE_cpue)),
		 RE                = matrix(0, nrow=max(sum(nobs_RE),1), ncol=Nage-1),
		 RE_cpue           = rep(0, max(1,sum(nobs_RE_cpue))),
		 s50               = 0,
		 logslope          =-5
	 )  

	}

    output <- list()
    output$tmb_data <- tmb_data
    output$tmb_params <- tmb_params
    output$tmb_phases <- tmb_phases
    output$pred_grid <- pred_grid
    output$mgcv_mod_cpue <- mgcv_mod_cpue
    output$mgcv_mod <- mgcv_mod
    output$data <- data
    output$spde <- spde
    output$proj_NA <- proj_NA
    
    return(output)
    
}    

