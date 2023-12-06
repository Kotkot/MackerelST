### Useful code chunk to set up Rstudio
  # TMB:::setupRStudio()

### Loading some useful functions
  source("R/Strata.R")
  source("R/Exclude_strata.R")

### Some general configurations
  
  use_new_data <- TRUE
  use_new_env = TRUE
  explore_data <- FALSE

  
  ### Read in geographical extent of the study area and convert projection to UTM
  world <- st_read("Shapefiles/ne_10m_land.shp")
  world <- world %>% st_make_valid()
  Atlantic <- st_crop(world, c(xmin = -65, ymin = 46, xmax = 32, ymax = 80))
  Atlantic_proj <- st_transform(Atlantic, crs=projection_km)
  
  all_strata <- st_read("Shapefiles/all_strata.shp")
  all_strata_sf <- all_strata %>% st_transform(crs = projection_km)
  South_end <- data.frame(lon=c(-25,-5,-5,-25,-25), lat=c(67,67,62,62, 67))
  South_end_sf <- st_as_sf(x=South_end, coords = c("lon", "lat")) %>% st_set_crs(4326) %>% 
    summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
  South_end_sf <- South_end_sf %>% st_transform(crs = projection_km) 
  NorthSea <- st_as_sf(x=All_strata %>% filter(ID == 13), coords = c("lon", "lat")) %>%
    st_set_crs(4326) %>% summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON")
  NorthSea_sf <- NorthSea %>% st_transform(crs = projection_km)
  all_strata_sf <- st_difference(all_strata_sf, NorthSea_sf)
  all_strata_sf_union <- st_union(all_strata_sf, South_end_sf)
  all_strata_sf_union_LL <- all_strata_sf_union %>% st_transform(crs=4326)


### Load data
  
  # bla <- readRDS(paste0(getwd(), "/data/new_data/CatchData.rds"))  # no mixed layer depth
  # new_dat <- readRDS(paste0(getwd(), "/data/new_data/iessns with selected env 07022022.rds"))  # no mixed layer depth
  new_dat <- readRDS(paste0(getwd(), "/data/new_data/iessns with SST CHL OMLT fronts.rds"))
  new_dat1 <- readRDS(paste0(getwd(), "/data/new_data/iessns with OMLT.rds")) # the new OMLT values do not match the previous one
  new_dat$OMLT_0m <- new_dat1$OMLT_0m
  CTD <- readRDS(paste0(getwd(), "/data/new_data/CTD.rds"))
  Plankton <- readRDS(paste0(getwd(), "/data/new_data/WP2.rds"))
  new_bio <- readRDS(paste0(getwd(), "/data/new_data/Mackerels.rds"))

  if (use_new_data == TRUE) {
    dat <- new_dat 
    projection_use <- projection
  } else {
    projection_use <- projection
  }

  
  if (use_new_data == TRUE) bio <- new_bio
  glimpse(dat)
  glimpse(bio)
  
  if (explore_data == TRUE) {
    # Some exploration of the data
    bio_sf <- bio %>% mutate(ID = 1:nrow(bio)) %>% filter(!is.na(WEIGHT)) %>% st_as_sf(coords = c("LON", "LAT")) %>% st_set_crs(4326) %>%
      st_cast("POINT")
    bio_sf_age <- bio %>% mutate(ID = 1:nrow(bio)) %>% filter(!is.na(AGE)) %>% st_as_sf(coords = c("LON", "LAT")) %>% st_set_crs(4326) %>%
      st_cast("POINT")
    # ggplot(Atlantic) + geom_sf() + geom_sf(data = bio_sf) + geom_sf(data = bio_sf_age, col="red", pch=".") + theme_bw() + facet_wrap(~YEAR)
    
    bio %>% group_by(YEAR, STATION) %>% filter(!is.na(WEIGHT)) %>% summarize(n=n())
    bio %>% group_by(YEAR, STATION) %>% filter(!is.na(AGE)) %>% summarize(n=n())
    
    # ggplot(Atlantic) + geom_sf() + geom_sf(data = bio_sf, pch=20, size=2) + geom_sf(data = noage_sf, col="red", size=1, pch=20) + theme_bw() + facet_wrap(~YEAR)
    
  }
  
  # Creating an unique ID that is the common identifier between the catch & bio datasets
  bio$ID <- apply(cbind(as.character(bio$CRUISE), bio$COUNTRY, bio$VESSEL, bio$STATION, bio$YEAR), 1, function(x) paste(x,collapse="_"))
  dat$ID <- apply(cbind(as.character(dat$CRUISE), dat$COUNTRY, dat$VESSEL, dat$STATION, dat$YEAR), 1, function(x) paste(x,collapse="_"))
  
  # for which stations do we have no age data but mackerel are caught (i.e. length samples are taken)?
  bio <- bio %>% group_by(ID, LENGTH) %>% mutate(n=n()) %>% filter(n>0) %>% ungroup() %>% group_by(ID) %>% mutate(n=sum(!is.na(AGE))) %>% ungroup()
  check <- bio %>% group_by(ID, LENGTH) %>% mutate(n=n()) %>% filter(n>0) %>% ungroup() %>% group_by(ID) %>% mutate(n=sum(!is.na(AGE)))
  noage_sf <- check %>% filter(n==0) %>% st_as_sf(coords = c("LON", "LAT")) %>% st_set_crs(4326) %>%
    st_cast("POINT")
  
  # do we have catch info as well for these stations?
  bio <- bio %>% left_join(dat)
  aaa <- bio %>% filter(n==0) %>% dplyr::select(COUNTRY, VESSEL, CRUISE, STATION, YEAR, MONTH, DAY, HOUR, MIN, LAT, LON, SPECIES, LENGTH, WEIGHT, 
                                                SEX, CATCH, MATURATION, AGE) 
  
  range(aaa$CATCH)

  
  ## MANUALLY read in the area with "zeros" = the "zero" boundary line i.e. this is not the coastal zone but 
  ## areas based on expert knowledge that should not have mackerel in
  ## This is the list that is used for modeling purpose i.e. extended beyong the strata to take into account for the mesh 
  Exclude_by_year <- list()
  Strata_list_exclude <- list(#c(4L,5L,6L,10L,11L,12L,13L,"Triangle","Area9_74","South_end","West_11","West_10","West_Greenland"),  # 2007
    #c(10L,11L,13L,"Triangle","Area9_75","South_end","West_11","West_10","West_Greenland"),               # 2009
    c(10L,11L,13L,"South_end","West_11","West_10", "West_Greenland"),                                    # 2010
    c(10L,11L,13L,"Triangle","Area9_74","South_end","West_11","West_10", "West_Greenland"),              # 2011
    c(10L,11L,13L,"Area9_74","South_end","West_11","West_10", "West_Greenland"),                         # 2012
    c(11L,13L,"South_end","West_11", "West_Greenland"),                                                  # 2013
    c(13L, "West_Greenland"),                                                                            # 2014
    c(13L, "West_Greenland"),                                                                            # 2015
    c(13L, "West_Greenland"),                                                                            # 2016
    c(13L, "West_Greenland"),                                                                            # 2017
    c(13L, "West_Greenland"),                                                                            # 2018
    c(13L, "West_Greenland"),                                                                            # 2019
    c(13L, "West_Greenland")                                                                             # 2020
  )
  
  Exclude_by_year <- lapply(Strata_list_exclude, function(x) exclude_area(data=All_bounds, x))

  ### Now adding LAT and LON (with this case specific writing)
  dat_df <- st_drop_geometry(dat)
  dat$LAT <- dat_df[,grep("lat", colnames(dat_df), ignore.case=TRUE)]
  dat$LON <- dat_df[,grep("lon", colnames(dat_df), ignore.case=TRUE)]
  
  ### Now transforming all continuous variables
  dat1 <- dat %>% as_tibble() %>% mutate(across(where(is.double), ~ scale(.x))) 
  dat1 <- dat1 %>% dplyr::select(where(is.double)) 
  colnames(dat1) <- paste0(colnames(dat1), "_scl")
  dat <- cbind(dat, dat1)
  
  not_in_use <- function(plot=FALSE){  
    ### This is the strata list excluded by year (same as "Exclude_by_year") for plotting purpose
    Exclude_by_year_simple <- list()
    # # 2007
    # Exclude_by_year_simple[[1]] <- All_bounds %>% filter(ID %in% c(4L,5L,6L,10L,11L,12L,13L,"Triangle")) %>% 
    #   st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>% st_transform(crs = projection_km) %>% 
    #   group_by(ID) %>% summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
    # # 2009
    # Exclude_by_year_simple[[2]] <- All_bounds %>% filter(ID %in% c(10L,11L,13L,"Triangle")) %>% 
    #   st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>%  st_transform(crs = projection_km) %>% 
    #   group_by(ID) %>% summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
    # 2010
    Exclude_by_year_simple[[1]] <- All_bounds %>% filter(ID %in% c(10L,11L,13L)) %>% 
      st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>%  st_transform(crs = projection_km) %>% 
      group_by(ID) %>% summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
    # 2011
    Exclude_by_year_simple[[2]] <- All_bounds %>% filter(ID %in% c(10L,11L,13L,"Triangle")) %>% 
      st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>%  st_transform(crs = projection_km) %>% 
      group_by(ID) %>% summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
    # 2012
    Exclude_by_year_simple[[3]] <- Exclude_by_year_simple[[1]]
    # 2013
    Exclude_by_year_simple[[4]] <- All_bounds %>% filter(ID %in% c(11L,13L)) %>% 
      st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>%  st_transform(crs = projection_km) %>% 
      group_by(ID) %>% summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
    # 2014
    Exclude_by_year_simple[[5]] <- All_bounds %>% filter(ID %in% c(13L)) %>% 
      st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>%  st_transform(crs = projection_km) %>% 
      group_by(ID) %>% summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
    # 2015
    Exclude_by_year_simple[[6]] <- Exclude_by_year_simple[[5]]
    # 2016
    Exclude_by_year_simple[[7]] <- Exclude_by_year_simple[[5]]
    # 2017
    Exclude_by_year_simple[[8]] <- Exclude_by_year_simple[[5]]
    # 2018
    Exclude_by_year_simple[[9]] <- Exclude_by_year_simple[[5]]
    # 2019
    Exclude_by_year_simple[[10]] <- Exclude_by_year_simple[[5]]
    # 2020
    Exclude_by_year_simple[[11]] <- Exclude_by_year_simple[[5]] 
    
    if (plot == TRUE) {
      for (i in 1:Nyear){
        p1 <- ggplot(Atlantic_proj) + geom_sf() + theme_bw() + geom_sf(data = Exclude_by_year_simple[[i]], aes(fill=as.factor(ID)), size = 0.2) + 
          scale_fill_viridis_d(name="Strata ID")  # + geom_sf_label(data = Exclude_by_year[[1]], aes(label = ID)) 
        p1
      }
    }
    
    return(Exclude_by_year_simple)
  }  
  
  Exclude_by_year_simple <- not_in_use(plot=FALSE)
