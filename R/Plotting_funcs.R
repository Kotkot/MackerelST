### some code for plotting esthetics for later
gg_control <- theme(axis.title = element_text(size=15),
                    axis.text = element_text(size=14), 
                    plot.title = element_text(hjust=0.5, size=18),
                    strip.text = element_text(size=14))  

##' Plot indices of abundance
##' @param mod model run output obtained from fit_mackerelST
##' @param IA_assessment index of abundance data.frame with the following column names ("Age", "Year", "IA_IESSNS")
##' @param conf configuration file obtained from conf_TMB
##' @param name name to use for saving the output plots
##' @param folder defines save location
##' @details Plot indices of abundance
##' @return Plot indices of abundance
##' @export plot_IA
plot_IA <- function(mod, conf, IA_assessment, name="step2", folder = getwd(), ...){  
  IA2 <- as.data.frame(summary(Mod_pred$sdrep, "report"))
  if (ncol(IA2) ==2 ) colnames(IA2) <- c("IA", "IA_sd")
  if (ncol(IA2) ==4 ) colnames(IA2) <- c("IA_old", "IA_sd_old", "IA", "IA_sd")
  IA2$Year <- rep(1:conf$Nyear, conf$Nage)
  IA2$Age <- rep(conf$ages, each=conf$Nyear)
  IA2$IA_type <- rep(c("IA", "IA_log", "IA_corrected", "IA_log_corrected"), each=nrow(IA2)/4)
  IA2$Year <- rep(1:conf$Nyear, conf$Nage)
  IA2$Age <- rep(conf$ages, each=conf$Nyear)
  IA2 <- IA2 %>% mutate(Index = as.numeric(IA/1000000), Year = as.numeric(as.character(factor(as.factor(Year), labels=conf$years))))
  IA2 <- IA2 %>% mutate(YC = Year - Age)   # adding the cohort 
  
  p1 <- ggplot(IA2 %>% filter(IA_type == "IA"), aes(x=Year, y=Index, col=factor(Age))) + geom_line(size=2) + theme_bw() + gg_control
  ggsave(p1, filename = paste0(folder, "/IAspatial", name, "_raw.png"), dpi ="retina", width = 10, height = 7, device = "png")
  p2 <- ggplot(IA2 %>% filter(IA_type == "IA_corrected"), aes(x=Year, y=Index, col=factor(Age))) + geom_line(size=2) + theme_bw() + gg_control
  ggsave(p2, filename = paste0(folder, "/IAspatial", name, "_corrected.png"), dpi ="retina", width = 10, height = 7, device = "png")
  
  IA2 <- IA2 %>% left_join(IA_assessment)
  
  
  # by year
  p1 <- ggplot(IA2 %>% filter(IA_type == "IA"), aes(x=IA_IESSNS, y=Index)) + facet_wrap(~Age, scales="free") + geom_text(aes(label=Year, col=factor(Year))) + 
    geom_point() + theme_bw() + geom_smooth(method="lm", se=FALSE) + gg_control + 
    labs(x="IESSN estimate (billion ton)", y="Spatial model estimate (in CPUE unit)") +
    stat_poly_eq(formula = y ~x, aes(label = paste(..rr.label.., sep = "~~~")), parse = TRUE) 
  ggsave(p1, filename = paste0(folder, "/IA_IESSN_vs_spatial", name, ".png"), dpi ="retina", width = 15, height = 10, device = "png")
  p1 <- ggplot(IA2 %>% filter(IA_type == "IA_corrected"), aes(x=IA_IESSNS, y=Index)) + facet_wrap(~Age, scales="free") + geom_text(aes(label=Year, col=factor(Year))) + 
    geom_point() + theme_bw() + geom_smooth(method="lm", se=FALSE) + gg_control + 
    labs(x="IESSN estimate (billion ton)", y="Spatial model estimate (in CPUE unit)") +
    stat_poly_eq(formula = y ~x, aes(label = paste(..rr.label.., sep = "~~~")), parse = TRUE) 
  ggsave(p1, filename = paste0(folder, "/IA_IESSN_vs_spatial", name, "_corrected.png"), dpi ="retina", width = 15, height = 10, device = "png")
  
  # Internal consistency
  p1 <- ggplot(IA2 %>% filter(IA_type == "IA"), aes(x=YC, y=Index, col=factor(Age))) + 
    geom_line(size=2) + theme_bw() + gg_control + 
    labs(x="Cohort year", y="Spatial model estimate (in CPUE unit)") 
  ggsave(p1, filename = paste0(folder, "/IA_cohort", name, ".png"), dpi ="retina", width = 10, height = 7, device = "png")
  p1 <- ggplot(IA2 %>% filter(IA_type == "IA_corrected"), aes(x=YC, y=Index, col=factor(Age))) + 
    geom_line(size=2) + theme_bw() + gg_control + 
    labs(x="Cohort year", y="Spatial model estimate (in CPUE unit)") 
  ggsave(p1, filename = paste0(folder, "/IA_cohort_corrected", name, ".png"), dpi ="retina", width = 10, height = 7, device = "png")
  p1 <- ggplot(IA2 %>% filter(IA_type == "IA"), aes(x=Age, y=Index, col=factor(YC))) + 
    geom_line(size=2) + theme_bw() + gg_control + 
    labs(x="Age", y="Spatial model estimate (in CPUE unit)") 
  ggsave(p1, filename = paste0(folder, "/IA_cohortplot", name, ".png"), dpi ="retina", width = 10, height = 7, device = "png")
  p1 <- ggplot(IA2 %>% filter(IA_type == "IA_corrected"), aes(x=Age, y=Index, col=factor(YC))) + 
    geom_line(size=2) + theme_bw() + gg_control + 
    labs(x="Age", y="Spatial model estimate (in CPUE unit)") 
  ggsave(p1, filename = paste0(folder, "/IA_cohortplot_corrected", name, ".png"), dpi ="retina", width = 10, height = 7, device = "png")
  p1 <- ggcorr(IA2 %>% filter(IA_type == "IA_corrected") %>% dplyr::select(Year, Age, Index) %>%
                 pivot_wider(names_from = Age, values_from = Index) %>% dplyr::select(! Year), label=TRUE)
  ggsave(p1, filename = paste0(folder, "/IA_cohort_corr", name, ".png"), dpi ="retina", width = 10, height = 7, device = "png")
    
}



##' Plot various map of predicted density, epsilon, etc on prediction grid 
##' @param run model run output obtained from fit_mackerelST function
##' @param data data output from prepare_data function
##' @param conf configuration file obtained from conf_TMB function
##' @param name name to use for saving the output plots
##' @param folder defines save location
##' @details Plots various maps of outputs
##' @return Plots various maps of outputs
##' @export Plot_maps
# save map outputs 
Plot_maps <- function(run=Mod_pred, data=datdat, name="step2", no_fish_zone=FALSE, plot_data=FALSE, folder = getwd(), conf, ...)  {
  
  ## Prediction maps
  if (no_fish_zone == TRUE){
    run$dens_pred[which(is.na(data$data$proj_NAisNA))] <- NA
    run$epsilon_pred[which(is.na(pred_mu$isNA))] <- NA
  }
  
  dens_pred = run$dens_pred %>% mutate(AGE = AGE + min(conf$ages) - min(AGE))
  if (length(run$epsilon_pred) != 0 ) epsilon_pred = run$epsilon_pred %>% mutate(AGE = AGE + min(conf$ages) - min(AGE))

  # by year
  for (yr in seq_along(sort(unique(data$data$YEAR))))
  {
    yy = sort(unique(data$data$YEAR))[yr]
    p1 <- ggplot(subset(dens_pred, YEAR == yy) %>% mutate(AGE_fake = as.factor(paste0("Age ", AGE)), 
                                                        AGE_fake2 = factor(AGE_fake, levels=paste0("Age ", conf$ages)),
                                                        Pred2 = ifelse(Pred>quantile(Pred,0.995, na.rm=T), quantile(Pred,0.995, na.rm=T), Pred))) +
      geom_raster(aes_string("X1000", "Y1000", fill = "logPred")) +
      geom_sf(data = Atlantic_proj$geometry, color ="grey27", size = .2)+
      geom_sf(data = Exclude_by_year_simple[[yr]]$geometry, col="red", fill=NA, size = 0.2) +
      xlab("") + ylab("") +
      facet_wrap(~ AGE_fake2) +
      coord_sf(xlim = c(1000,5000), ylim = c(4000, 6200),
               datum = projection_km, expand = FALSE) + theme_bw() + 
      # scale_fill_viridis_c() +
      scale_fill_gradient2(low="darkblue", mid="white", high="red2", midpoint=0) +
      ggtitle(paste0("CPUE distribution in ", yy)) + 
      theme(axis.title = element_text(size=15),
            axis.text = element_text(size=14), 
            plot.title = element_text(hjust=0.5, size=18),
            strip.text = element_text(size=14))
    ggsave(p1, filename = paste0(folder, "/Prediction", yy, "_", name, "nofishzone", no_fish_zone, ".png"), dpi ="retina", width = 12, height = 10, device = "png")
  }
  # by age
  for (age in conf$ages)
  {
    p1 <- ggplot(dens_pred %>% filter(AGE == age)) +
      geom_raster(aes_string("X1000", "Y1000", fill = "logPred")) +
      geom_sf(data = Atlantic_proj$geometry, color ="grey27", size = .2)+
      # geom_sf(data = Exclude_by_year_simple[[yr]]$geometry, col="red", fill=NA, size = 0.2) +
      xlab("") + ylab("") +
      facet_wrap(~ YEAR) +
      coord_sf(xlim = c(1000,5000), ylim = c(4000, 6200),
               datum = projection_km, expand = FALSE) + theme_bw() + 
      # scale_fill_viridis_c(option = "plasma") +
      scale_fill_gradient2(low="darkblue", mid="white", high="red2", midpoint=1) +
      ggtitle(paste0("CPUE distribution for age ", age)) + 
      theme(axis.title = element_text(size=15),
            axis.text = element_text(size=14), 
            plot.title = element_text(hjust=0.5, size=18),
            strip.text = element_text(size=14))
    ggsave(p1, filename = paste0(folder, "/Prediction_age", age,"_", name,  "nofishzone", no_fish_zone, ".png"), dpi ="retina", width = 12, height = 10, device = "png")
  }
  # by age (truncated at values below 1)
  for (age in conf$ages)
  {
    p1 <- ggplot(dens_pred %>% filter(AGE == age) %>% mutate(logPred_trunc = ifelse(logPred<1, 1, logPred))) +
      geom_raster(aes_string("X1000", "Y1000", fill = "logPred_trunc")) +
      geom_sf(data = Atlantic_proj$geometry, color ="grey27", size = .2)+
      # geom_sf(data = Exclude_by_year_simple[[yr]]$geometry, col="red", fill=NA, size = 0.2) +
      xlab("") + ylab("") +
      facet_wrap(~ YEAR) +
      coord_sf(xlim = c(1000,5000), ylim = c(4000, 6200),
               datum = projection_km, expand = FALSE) + theme_bw() + 
      scale_fill_viridis_c(option = "plasma") +
      # scale_fill_gradient2(low="darkblue", mid="white", high="red2", midpoint=1) +
      ggtitle(paste0("CPUE distribution for age ", age)) + 
      theme(axis.title = element_text(size=15),
            axis.text = element_text(size=14), 
            plot.title = element_text(hjust=0.5, size=18),
            strip.text = element_text(size=14))
    ggsave(p1, filename = paste0(folder, "/Prediction_truncated_age", age,"_", name,  "nofishzone", no_fish_zone, ".png"), dpi ="retina", width = 12, height = 10, device = "png")
  }
  
  if (plot_data  == TRUE){
    ## For the actual data itself
    for (age in conf$ages)
    {
      asd <- data
      asd$logobs <- log(asd[,grep("ID", colnames(data$data))+age]+0.01)
      p1 <- ggplot(asd) +
        geom_point(aes_string("X1000", "Y1000", fill = "logobs", col="logobs"), size=0.2) +
        geom_sf(data = Atlantic_proj$geometry, color ="grey27", size = .2)+
        # geom_sf(data = Exclude_by_year_simple[[yr]]$geometry, col="red", fill=NA, size = 0.2) +
        xlab("") + ylab("") +
        facet_wrap(~ YEAR) +
        coord_sf(xlim = c(1000,5000), ylim = c(4000, 6200),
                 datum = projection_km, expand = FALSE) + theme_bw() + 
        scale_fill_viridis_c(option = "plasma") +
        scale_color_viridis_c(option = "plasma") +
        ggtitle(paste0("CPUE distribution for age ", age)) + 
        theme(axis.title = element_text(size=15),
              axis.text = element_text(size=14), 
              plot.title = element_text(hjust=0.5, size=18),
              strip.text = element_text(size=14))
      ggsave(p1, filename = paste0(folder, "/Obs_age", age, "nofishzone", no_fish_zone, ".png"), dpi ="retina", width = 12, height = 10, device = "png")
    }
    ## observation trucated at 1
    for (age in conf$ages)
    {
      asd <- data
      asd$logobs <- ifelse(asd[,grep("ID", colnames(data$data))+age]>1, log(asd[,grep("ID", colnames(data$data))+age]), 0)
      p1 <- ggplot(asd) +
        geom_point(aes_string("X1000", "Y1000", fill = "logobs", col="logobs"), size=0.2) +
        geom_sf(data = Atlantic_proj$geometry, color ="grey27", size = .2)+
        # geom_sf(data = Exclude_by_year_simple[[yr]]$geometry, col="red", fill=NA, size = 0.2) +
        xlab("") + ylab("") +
        facet_wrap(~ YEAR) +
        coord_sf(xlim = c(1000,5000), ylim = c(4000, 6200),
                 datum = projection_km, expand = FALSE) + theme_bw() + 
        scale_fill_viridis_c(option = "plasma") +
        scale_color_viridis_c(option = "plasma") +
        ggtitle(paste0("CPUE distribution for age ", age)) + 
        theme(axis.title = element_text(size=15),
              axis.text = element_text(size=14), 
              plot.title = element_text(hjust=0.5, size=18),
              strip.text = element_text(size=14))
      ggsave(p1, filename = paste0(folder, "/Obs_trunc_age", age, "nofishzone", no_fish_zone, ".png"), dpi ="retina", width = 12, height = 10, device = "png")
    }
  }
  
  ## For the estimated spatio-temporal field
  if (length(run$epsilon_pred) != 0) {
    for (yr in seq_along(sort(unique(data$data$YEAR))))
    {
      yy = sort(unique(data$data$YEAR))[yr]
      p1 <- ggplot(subset(epsilon_pred, YEAR == yy) %>% mutate(AGE_fake = as.factor(paste0("Age ", AGE)), AGE_fake2 = factor(AGE_fake, levels=paste0("Age ", 1:15)))) +
        geom_raster(aes_string("X1000", "Y1000", fill = "Pred")) +
        geom_sf(data = Atlantic_proj$geometry, color ="grey27", size = .2)+
        geom_sf(data = Exclude_by_year_simple[[yr]]$geometry, col="red", fill=NA, size = 0.2) +
        xlab("") + ylab("") +
        facet_wrap(~ AGE_fake2) +
        coord_sf(xlim = c(1000,5000), ylim = c(4000, 6200),
                 datum = projection_km, expand = FALSE) + theme_bw() + 
        scale_fill_gradient2(low="darkblue", mid="white", high="red") +
        ggtitle(paste0("CPUE distribution in ", yy)) + 
        theme(axis.title = element_text(size=15),
              axis.text = element_text(size=14), 
              plot.title = element_text(hjust=0.5, size=18),
              strip.text = element_text(size=14))
      ggsave(p1, filename = paste0(folder, "/Epsilon_st", yy, "_", name, "nofishzone", no_fish_zone, ".png"), dpi ="retina", width = 12, height = 10, device = "png")
    }
    # by age
    for (age in conf$ages)
    {
      p1 <- ggplot(subset(epsilon_pred, AGE == age)) +
        geom_raster(aes_string("X1000", "Y1000", fill = "Pred")) +
        geom_sf(data = Atlantic_proj$geometry, color ="grey27", size = .2)+
        geom_sf(data = Exclude_by_year_simple[[yr]]$geometry, col="red", fill=NA, size = 0.2) +
        xlab("") + ylab("") +
        facet_wrap(~ YEAR) +
        coord_sf(xlim = c(1000,5000), ylim = c(4000, 6200),
                 datum = projection_km, expand = FALSE) + theme_bw() + 
        scale_fill_gradient2(low="darkblue", mid="white", high="red") +
        ggtitle(paste0("CPUE distribution age ", age)) + 
        theme(axis.title = element_text(size=15),
              axis.text = element_text(size=14), 
              plot.title = element_text(hjust=0.5, size=18),
              strip.text = element_text(size=14))
      ggsave(p1, filename = paste0(folder, "/Epsilon_st_age", age, "_", name, "nofishzone", no_fish_zone, ".png"), dpi ="retina", width = 12, height = 10, device = "png")
    }
  }

}        



##' Model diagnostics: QQplot of residuals based on PIT residuals
##' @param run model run output obtained from fit_mackerelST function
##' @param data data output from prepare_data function
##' @param conf configuration file obtained from conf_TMB function
##' @param name name to use for saving the output plots
##' @param folder defines save location
##' @details residual qqplot
##' @return residual qqplot
##' @export Plot_residuals
#### Model diagnostics
Plot_residuals <- function(run, data, parval = NULL, conf, name="step2", folder = getwd(), do_all=TRUE, ...){
  out <- c()
  resids <- list()
  
  if (is.null(parval) == TRUE) parval = run$obj$env$last.par.best
  
  to_remove <- as.numeric(which(is.na(data$tmb_data$yobs[,1]) == TRUE))
  if(length(to_remove)==0) to_remove=NULL
  
  for (age in seq_along(conf$ages)){
    # age = 4
    model_output <- list()
    model_output$response <- data$tmb_data$yobs[-to_remove,age]
    model_output$mu <- run$obj$report(parval)$mu[-to_remove,age]
    model_output$family <- "tweedie"
    if (is.null(parval) == TRUE) {
      if (!is.null(conf$conf_model$ln_phi)) model_output$ln_phi <- run$opt$par[grep("ln_phi", names(run$opt$par))][as.numeric(conf$conf_model$ln_phi)[age]]
      if (is.null(conf$conf_model$ln_phi)) model_output$ln_phi <- run$opt$par[grep("ln_phi", names(run$opt$par))][age]
      if (!is.null(conf$conf_model$thetaf)) model_output$thetaf <- run$opt$par[grep("thetaf", names(run$opt$par))][as.numeric(conf$conf_model$thetaf)[age]]
      if (is.null(conf$conf_model$thetaf)) model_output$thetaf <- run$opt$par[grep("thetaf", names(run$opt$par))][age]
    }
    if (is.null(parval) == FALSE) {
      if (!is.null(conf$conf_model$ln_phi)) model_output$ln_phi <- parval[grep("ln_phi", names(parval))][as.numeric(conf$conf_model$ln_phi)[age]]
      if (is.null(conf$conf_model$ln_phi)) model_output$ln_phi <- parval[grep("ln_phi", names(parval))][age]
      if (!is.null(conf$conf_model$thetaf)) model_output$thetaf <- parval[grep("thetaf", names(parval))][as.numeric(conf$conf_model$thetaf)[age]]
      if (is.null(conf$conf_model$thetaf)) model_output$thetaf <- parval[grep("thetaf", names(parval))][age]
    }
    model_output$obj <- run$obj
    model_output$tmb_data <- data$tmb_data
    
    d <- data$data[-to_remove,]
    
    d$residuals <- get.residuals(model_output)
    d <- d[order(d$residuals),]
    d$qq <- (qnorm(ppoints(nrow(d))))
    p1 <- ggplot(d, aes(x= qq, y=residuals)) + geom_point(size=1) + geom_abline(slope=1) + theme_bw() + 
      scale_color_viridis_c(option = "plasma") + labs(x="Theoretical", y="Sample") + ggtitle(paste0("Age:", conf$ages[age]))
    resids[[age]] <- p1
    if (do_all == TRUE) {
      ggsave(p1, filename = paste0(folder, "/res_age", conf$ages[age], name, "_aggregated.png"), dpi ="retina", width = 7, height = 7, device = "png")
      p2 <- ggplot(d, aes(x= qq, y=residuals)) + geom_point(size=1) + geom_abline(slope=1) + theme_bw() + facet_wrap(~YEAR) +
        scale_color_viridis_c(option = "plasma") + labs(x="Theoretical", y="Sample") + ggtitle(paste0("Age:", conf$ages[age]))
      p2
      ggsave(p2, filename = paste0(folder, "/res_age", conf$ages[age], name, "_byyear.png"), dpi ="retina", width = 7, height = 7, device = "png")
    }
    
    d$age = conf$ages[age]
    out <- rbind(out, d)
  }
  
  if (conf$Nage <= 9) pp <- do.call("grid.arrange", c(resids, ncol=3, nrow=3))
  if (conf$Nage > 9) pp <- do.call("grid.arrange", c(resids, ncol=4, nrow=3))
  
  ggsave(pp, filename = paste0(folder, "/res_age_all", name, "_aggregated.png"), dpi ="retina", width = 8, height = 6, device = "png")
  
  
  return(list(out=out, to_remove=to_remove))
}



##' Model diagnostics: QQplot of residuals as in DHARMa (based on simulation)
##' @param simulated_response this is a matrix of simulated response from the model: row for observation, columns for simulation index
##' @param pred model best prediction for the observation vector
##' @param obs observation vector
##' @param plot name to use for saving the output plots
##' @details residual qqplot as in DHARMa
##' @return residual qqplot as in DHARMa
##' @export Plot_residuals
#### Model diagnostics
dharma_residuals <- function(simulated_response, pred, obs, plot = TRUE, ...) {
  
  res <- DHARMa::createDHARMa(
    simulatedResponse = simulated_response,
    observedResponse = obs,
    fittedPredictedResponse = pred,
    ...
  )
  u <- res$scaledResiduals
  n <- length(u)
  m <- seq_len(n) / (n + 1)
  z <- stats::qqplot(m, u, plot.it = FALSE)
  if (plot) {
    DHARMa::plotQQunif(
      res,
      testUniformity = FALSE,
      testOutliers = FALSE, testDispersion = FALSE
    )
  }
  invisible(data.frame(observed = z$y, expected = z$x))
}



##' Plot center of gravity of predicted distribution over the grid
##' @param run model run output obtained from fit_mackerelST function
##' @param data data output from prepare_data function
##' @param conf configuration file obtained from conf_TMB function
##' @param name name to use for saving the output plots
##' @param folder defines save location
##' @details CG plot
##' @return CG plot
##' @export CG_plots
#### Plotting the center of gravity
CG_plots <- function(run=Mod_pred, data = datdat, name="step2", no_fish_zone=FALSE, folder = getwd(), conf, ...){
  
  loc_NAs <- reshape2::melt(data$proj_NA)
  colnames(loc_NAs) <- c("LOC", "YEAR", "isNA")
  loc_NAs$YEAR <- as.factor(loc_NAs$YEAR)
  loc_NAs$YEAR <- factor(loc_NAs$YEAR, labels=conf$years)
  bla <- run$dens_pred 
  bla$AGE = bla$AGE + min(conf$ages) - min(bla$AGE)
  if (no_fish_zone == TRUE){
    bla <- run$dens_pred %>% left_join(loc_NAs) %>% mutate(Pred = ifelse(is.na(isNA) == TRUE, NA, Pred),
                                                   logPred  = ifelse(is.na(isNA) == TRUE, NA, logPred))
  }
  bla <- bla %>% mutate(YC = (as.numeric(YEAR) - as.numeric(AGE)+2009))
  CG <- bla %>% group_by(YEAR, AGE) %>% summarize(meanX = sum(LON*Pred/sum(Pred, na.rm=T), na.rm=T), meanY = sum(LAT*Pred/sum(Pred, na.rm=T), na.rm=T)) 
  CG <- bla %>% group_by(YEAR, AGE) %>% summarize(meanX = sum(LON*Pred/sum(Pred, na.rm=T), na.rm=T), meanY = sum(LAT*Pred/sum(Pred, na.rm=T), na.rm=T)) 
  CG <- CG %>% mutate(YC = (as.numeric(YEAR) - as.numeric(AGE)+2009)) 

  ## Each panel is age, and we track the center of gravity by year
  p1 <- ggplot(Atlantic) + geom_sf() + coord_sf(xlim = c(-30, 15), ylim = c(60, 71), expand = FALSE) +
    geom_path(data=CG, aes(x=meanX, y=meanY), col=grey(0.8), linewidth=1.2) + 
    geom_point(data=CG, aes(x=meanX, y=meanY, col=factor(YEAR)), size=3) + 
    theme_bw() + gg_control + facet_wrap(~factor(AGE)) + geom_text_repel(data=CG, aes(x=meanX, y=meanY, label=YEAR), size=4, max.overlaps=20, segment.color = NA) + 
    labs(x="Longitude", y="Latitude") + scale_color_viridis_d(name = "Year") + 
    theme(legend.title = element_text(size = 15),legend.text = element_text(size = 14))
  p1
  ggsave(p1, filename = paste0(folder, "/CG_byage_nofishzone", no_fish_zone, ".png"), dpi ="retina", width = 12, height = 8, device = "png")
  
  ## Each panel is Year, and we track the center of gravity by age
  p2 <- ggplot(Atlantic) + geom_sf() + coord_sf(xlim = c(-30, 15), ylim = c(60, 71), expand = FALSE) +
    geom_path(data=CG, aes(x=meanX, y=meanY), col=grey(0.8), linewidth=1.2) + 
    geom_point(data=CG, aes(x=meanX, y=meanY, col=AGE), size=3) +  
    theme_bw() + gg_control + facet_wrap(~YEAR, ncol=3) + geom_text_repel(data=CG, aes(x=meanX, y=meanY, label=AGE), size=4, max.overlaps=20, segment.color = NA) + 
    labs(x="Longitude", y="Latitude") + scale_color_viridis_c(name = "Age") + scale_fill_viridis_c(name = "Age") + 
    theme(legend.title = element_text(size = 15),legend.text = element_text(size = 14))
  
  p2
  ggsave(p2, filename = paste0(folder, "/CG_byyear_nofishzone", no_fish_zone, ".png"), dpi ="retina", width = 10, height = 8, device = "png")
  
  ## Each panel is Year, and we track the center of gravity by age
  p3 <-  ggplot(Atlantic) + geom_sf() + coord_sf(xlim = c(-30, 15), ylim = c(60, 71), expand = FALSE) +
    geom_path(data=CG, aes(x=meanX, y=meanY), col=grey(0.2), linewidth=1.5) + facet_wrap(.~YC) + 
    geom_point(data=CG, aes(x=meanX, y=meanY, col=factor(AGE)), size=2) +  
    theme_bw() + gg_control + 
    geom_text_repel(data=CG, aes(x=meanX, y=meanY, label=AGE), size=3, max.overlaps=20, segment.color = NA) + 
    labs(x="Longitude", y="Latitude") + scale_color_manual(values=c(scales::viridis_pal()(8)), name="Age") + 
    scale_fill_manual(values=c(scales::viridis_pal()(8)), name="Age") 
  
  p3
  ggsave(p3, filename = paste0(folder, "/CG_YC_all_nofishzone", no_fish_zone, ".png"), dpi ="retina", width = 12, height = 8, device = "png")
  
  set.seed(1)
  pre_2012 = bla %>% filter(YC <= 2012) %>% group_by(AGE, YC) %>% 
    summarize(meanX = sum(LON*Pred/sum(Pred, na.rm=T), na.rm=T), meanY = sum(LAT*Pred/sum(Pred, na.rm=T), na.rm=T)) %>% 
    ungroup() %>% group_by(AGE) %>% summarize(meanX = mean(meanX), meanY = mean(meanY))
  pre_2012$YC = "Pre2012"
  post_2013 = bla %>% filter(YC > 2012) %>% group_by(AGE, YC) %>% 
    summarize(meanX = sum(LON*Pred/sum(Pred, na.rm=T), na.rm=T), meanY = sum(LAT*Pred/sum(Pred, na.rm=T), na.rm=T)) %>% 
    ungroup() %>% group_by(AGE) %>% summarize(meanX = mean(meanX), meanY = mean(meanY))
  post_2013$YC = "Post2013"
  p3 <- ggplot(Atlantic) + geom_sf() + coord_sf(xlim = c(-30, 15), ylim = c(60, 71), expand = FALSE) +
    geom_path(data=pre_2012, aes(x=meanX, y=meanY), col=grey(0.7), linewidth=2) + 
    geom_point(data=pre_2012, aes(x=meanX, y=meanY, col=factor(AGE)), size=3) + 
    # geom_path(data=CG %>% filter(YC <= 2012), aes(x=meanX, y=meanY, col=YC), col=grey(0.8, alpha=0.5), linewidth=1) +
    # geom_point(data=CG %>% filter(YC <= 2012), aes(x=meanX, y=meanY), size=2, col=grey(0.8, alpha=0.5)) +
    # geom_point(data=CG %>% filter(YC <= 2012, AGE ==3), aes(x=meanX, y=meanY), size=3, col="#440154FF") +
    # geom_point(data=CG %>% filter(YC <= 2012, AGE %in% c(3,7,10)), aes(x=meanX, y=meanY, col=factor(AGE)), size=2) +
    geom_convexhull(data=CG %>% filter(YC <= 2012, AGE ==3), aes(x=meanX, y=meanY), alpha = 0.7, fill=NA, col="#440154FF", linewidth=1.1) +
    geom_convexhull(data=CG %>% filter(YC <= 2012, AGE ==5), aes(x=meanX, y=meanY), alpha = 0.7, fill=NA, col="#365C8DFF", linewidth=1.1) +
    geom_convexhull(data=CG %>% filter(YC <= 2012, AGE ==10), aes(x=meanX, y=meanY), alpha = 0.7, fill=NA, col="#FDE725FF", linewidth=1.1) +
    theme_bw() + gg_control + 
    geom_text_repel(data=pre_2012, aes(x=meanX, y=meanY, label=AGE, col=factor(AGE)), size=4, max.overlaps=20, segment.color = NA, nudge_x=c(0.4,0.1,-0.1,-0.1,-0.5,-0.5,0.2,0), nudge_y=c(0,0.4,0,0,0.4,0.4,0.5,0)) + 
    labs(x="Longitude", y="Latitude") + scale_color_manual(values=c(scales::viridis_pal()(8)), name="Age") + scale_fill_manual(values=c(scales::viridis_pal()(8)), name="Age")+ 
    ggtitle("(a) Average pre-2012")# + labs(subtitle = "(a)")
  p30 <- p3 + theme(legend.position="none")
  p31 <- ggplot(Atlantic) + geom_sf() + coord_sf(xlim = c(-30, 15), ylim = c(60, 71), expand = FALSE) +
    geom_path(data=post_2013, aes(x=meanX, y=meanY), col=grey(0.7), linewidth=2) +
    geom_point(data=post_2013, aes(x=meanX, y=meanY, col=factor(AGE)), size=3) +  
    # geom_point(data=CG %>% filter(YC > 2012, AGE %in% c(3,7,10)), aes(x=meanX, y=meanY, col=factor(AGE)), size=2) +
    geom_convexhull(data=CG %>% filter(YC > 2012, AGE ==3), aes(x=meanX, y=meanY), alpha = 0.7, fill=NA, col="#440154FF", linewidth=1.1) +
    geom_convexhull(data=CG %>% filter(YC > 2012, AGE ==5), aes(x=meanX, y=meanY), alpha = 0.7, fill=NA, col="#365C8DFF", linewidth=1.1) +
    geom_convexhull(data=CG %>% filter(YC > 2012, AGE ==10), aes(x=meanX, y=meanY), alpha = 0.7, fill=NA, col="#FDE725FF", linewidth=1.1) +
    # geom_path(data=CG %>% filter(YC > 2012), aes(x=meanX, y=meanY, col=YC), col=grey(0.8, alpha=0.5), linewidth=1) +
    # geom_point(data=CG %>% filter(YC > 2012), aes(x=meanX, y=meanY), size=2, col=grey(0.8, alpha=0.5)) +
    # geom_point(data=CG %>% filter(YC > 2012, AGE ==3), aes(x=meanX, y=meanY), size=3, col="#440154FF") +
    # geom_point(data=CG %>% filter(YC > 2012, AGE ==10), aes(x=meanX, y=meanY), size=3, col="#FDE725FF") +
    theme_bw() + gg_control + 
    geom_text_repel(data=post_2013, aes(x=meanX, y=meanY, label=AGE, col=factor(AGE)), size=4, max.overlaps=20, segment.color = NA, nudge_x=c(-0.5,-0.2,-0.5,-0.5,-0.5)) + 
    labs(x="Longitude", y="Latitude") + scale_color_manual(values=c(scales::viridis_pal()(8))[1:5], name="Age") + scale_fill_manual(values=c(scales::viridis_pal()(8))[1:5], name="Age")  + 
    ggtitle("(b) Average post-2013")+ theme(legend.position="none")# + labs(subtitle = "(b)")
  p32 <- ggplot(Atlantic) + geom_sf() + coord_sf(xlim = c(-30, 15), ylim = c(60, 71), expand = FALSE) +
    geom_path(data=CG %>% filter(YC == 2007), aes(x=meanX, y=meanY), col=grey(0.7), linewidth=2) +
    geom_point(data=CG %>% filter(YC == 2007), aes(x=meanX, y=meanY, col=factor(AGE)), size=3) +  
    theme_bw() + gg_control + 
    geom_text_repel(data=CG %>% filter(YC == 2007), aes(x=meanX, y=meanY, label=AGE, col=factor(AGE)), size=4, max.overlaps=20, segment.color = NA, nudge_y=c(0.5,0,0.4, 0.5,0.5,0.2,0.7,-0.5), nudge_x=c(0,-0.5,0.4,-0.1,-0.4,0.6,0,0)) + 
    labs(x="Longitude", y="Latitude") + scale_color_manual(values=c(scales::viridis_pal()(8)), name="Age") + scale_fill_manual(values=c(scales::viridis_pal()(8)), name="Age") + 
    ggtitle("(c) 2007") + theme(legend.position="none")# + labs(subtitle = "(c)")
  p33 <- ggplot(Atlantic) + geom_sf() + coord_sf(xlim = c(-30, 15), ylim = c(60, 71), expand = FALSE) +
    geom_path(data=CG %>% filter(YC == 2009), aes(x=meanX, y=meanY), col=grey(0.7), linewidth=2) +
    geom_point(data=CG %>% filter(YC == 2009), aes(x=meanX, y=meanY, col=factor(AGE)), size=3) +  
    theme_bw() + gg_control + 
    geom_text_repel(data=CG %>% filter(YC == 2009), aes(x=meanX, y=meanY, label=AGE, col=factor(AGE)), size=4, max.overlaps=20, segment.color = NA, nudge_y=c(0,0.2,0.1,-0.2,0.45,0.5,0.5,0.5), nudge_x=c(0,0,0.3,0,0.2,0,0,0)) + 
    labs(x="Longitude", y="Latitude") + scale_color_manual(values=c(scales::viridis_pal()(8)), name="Age") + scale_fill_manual(values=c(scales::viridis_pal()(8)), name="Age")  + 
    ggtitle("(d) 2009")+ theme(legend.position="none")# + labs(subtitle = "(d)")
  p34 <- ggplot(Atlantic) + geom_sf() + coord_sf(xlim = c(-30, 15), ylim = c(60, 71), expand = FALSE) +
    geom_path(data=CG %>% filter(YC == 2011), aes(x=meanX, y=meanY), col=grey(0.7), linewidth=2) +
    geom_point(data=CG %>% filter(YC == 2011), aes(x=meanX, y=meanY, col=factor(AGE)), size=3) +  
    theme_bw() + gg_control + 
    geom_text_repel(data=CG %>% filter(YC == 2011), aes(x=meanX, y=meanY, label=AGE, col=factor(AGE)), size=4, max.overlaps=20, segment.color = NA, nudge_y=0.5) + 
    labs(x="Longitude", y="Latitude") + scale_color_manual(values=c(scales::viridis_pal()(8))[1:7], name="Age") + scale_fill_manual(values=c(scales::viridis_pal()(8))[1:7], name="Age")  + 
    ggtitle("(e) 2011")+ theme(legend.position="none")# + labs(subtitle = "(e)")
  p35 <- ggplot(Atlantic) + geom_sf() + coord_sf(xlim = c(-30, 15), ylim = c(60, 71), expand = FALSE) +
    geom_path(data=CG %>% filter(YC == 2013), aes(x=meanX, y=meanY), col=grey(0.7), linewidth=2) +
    geom_point(data=CG %>% filter(YC == 2013), aes(x=meanX, y=meanY, col=factor(AGE)), size=3) +  
    theme_bw() + gg_control + 
    geom_text_repel(data=CG %>% filter(YC == 2013), aes(x=meanX, y=meanY, label=AGE, col=factor(AGE)), size=4, max.overlaps=20, segment.color = NA, nudge_y=c(0.5,0.7,0.5,0.7,0.5)) + 
    labs(x="Longitude", y="Latitude") + scale_color_manual(values=c(scales::viridis_pal()(8)), name="Age") + scale_fill_manual(values=c(scales::viridis_pal()(8)), name="Age") + 
    ggtitle("(f) 2013")+ theme(legend.position="none")# + labs(subtitle = "(f)")
  
  p_legend <- get_legend(p3 + theme(legend.title = element_text(size = 15),legend.text = element_text(size = 14)), position = NULL) %>% as_ggplot()
  
  
  pp <- grid.arrange(grobs=list(p30,p31,p32,p33,p34,p35,p_legend), layout_matrix=matrix(c(1,1,1,2,2,2,7,
                                                                                          3,3,3,4,4,4,7,
                                                                                          5,5,5,6,6,6,7),nrow=3,byrow=T))
  ggsave(pp, filename = "D:\\Dropbox\\IMR_projects\\Mackerel_distribution\\MS\\Figs\\Fig5new.png", dpi ="retina", width = 13, height = 11, device = "png")
  # ggsave(pp, filename = paste0(folder, "/CG_YC_select_nnofishzone", no_fish_zone, "new.png"), dpi ="retina", width = 11, height = 10, device = "png")
} 



##' Plot indices of abundance
##' @param run model run output obtained from fit_mackerelST function
##' @param data data output from prepare_data function
##' @param conf configuration file obtained from conf_TMB function
##' @param name name to use for saving the output plots
##' @param folder defines save location
##' @details Plot indices of abundance
##' @return Plot indices of abundance
##' @export Plot_env
#### Plotting the distribution of the environmental covariates to explain their distribution shift
Plot_env <- function(run = Mod_pred, var = "SST_0m_scl", truncate=FALSE, folder = getwd()) {
  
  # disp_win_wgs84 <- st_sfc(st_point(c(-25, 32)), st_point(c(70, 80)),
                         # crs = 4326)
	# disp_win_trans <- st_transform(disp_win_wgs84, crs = projection_km)
	# disp_win_coord <- st_coordinates(disp_win_trans)
											 
	dat_use = run$dens_pred %>% filter(AGE == 1)
  dat_use$resp <- dat_use[,eval(parse(text=var))]
  if (var == "SST_0m") dat_use$resp <- dat_use$resp-273.15
  dat_use <- dat_use %>% mutate(resp_trunc = ifelse(resp > quantile(resp, 0.99, na.rm=T), quantile(resp, 0.99, na.rm=T), ifelse(resp < quantile(resp, 0.01, na.rm=T), quantile(resp, 0.01, na.rm=T), resp)))
  if (truncate==TRUE) dat_use$resp_val = dat_use$resp_trunc
  if (truncate==FALSE) dat_use$resp_val = dat_use$resp
  p1 <- ggplot(dat_use) +
    geom_raster(aes(x=X1000, y=Y1000, fill = resp_val)) +
    geom_sf(data = Atlantic_proj$geometry, color ="grey27", size = .2)+
    xlab("Easting (km)") + ylab("Northing (km)") +
    coord_sf(xlim = c(1000,5000), ylim = c(4000, 6200),
             datum = projection_km, expand = FALSE) + theme_bw() + 
    facet_wrap(~YEAR, nrow=4)+
    # scale_fill_viridis_c() +
    scale_fill_gradient2(low="darkblue", mid="white", high="red2", midpoint=0) +
    ggtitle(paste0(var, " distribution")) + 
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=14), 
          plot.title = element_text(hjust=0.5, size=18),
          strip.text = element_text(size=14)) + 
    scale_x_continuous(breaks = c(1000, 3000, 5000))
  
  ggsave(p1, filename = paste0(folder, "/Env_map", "_", var, ".png"), dpi ="retina", width = 15, height = 12, device = "png")
}






##' Plot marginal effect of different variables 
##' @param run model run output obtained from fit_mackerelST function
##' @param data data output from prepare_data function
##' @param conf configuration file obtained from conf_TMB function
##' @param variable name of the variable to plot the marginal effect plot for
##' @param ref_year reference year to choose for the marginal effect calculation
##' @param folder defines save location
##' @details Plot indices of abundance
##' @return Plot indices of abundance
##' @export Marginal_plots
# DEPRACATED (better predict directly through the model)
Marginal_plots <- function(run, data_list, conf, variable="SST_0m", ref_year = 2018, folder = getwd(), ...) {
  
  # selected variable 
  variables <- as.character(attributes(conf)$fixed_effect)
  which.var <- which(variables == variable)
  whichnot_var_year <- which(! variables %in% c("YEAR", variable))
  
  dat <- data.frame(YEAR = factor(ref_year), t(colMeans(as.matrix(data_list$data[,variables[whichnot_var_year]]))))
  colnames(dat) <- c("YEAR", variables[whichnot_var_year])
  dat1 <- cbind(dat, X= seq(min(data_list$data[,variable]), max(data_list$data[,variable]), length.out=100))
  colnames(dat1) <- c(colnames(dat), variable)
  
  X_var = mgcv::predict.gam(data_list$mgcv_mod, type = "lpmatrix", newdata = dat1)
  
  Marginal_df <- c()
  for (age in 1:conf$Nage){
    betas <- run$opt$par[(1:ncol(X_var)+ncol(X_var)*(age-1))] 
    cov = run$sdrep$cov.fixed[(1:ncol(X_var)+ncol(X_var)*(age-1)),(1:ncol(X_var)+ncol(X_var)*(age-1))] 
    vars <- rep(0, nrow(X_var))
    for (obs in 1:nrow(X_var)){
      xx <- X_var[obs,]
      for (i in 1:nrow(cov)){
        for (j in 1:ncol(cov)){
          vars[obs] = vars[obs] + xx[i]*xx[j]*exp(betas[i]*xx[i])*exp(betas[j]*xx[j])*cov[i,j]
        }
      }
    }
    dat1$pred <- as.numeric(exp(X_var %*% betas))
    dat1$pred_se <- sqrt(vars)
    dat1$CI_lower <- dat1$pred - 1.96*dat1$pred_se
    dat1$CI_upper <- dat1$pred + 1.96*dat1$pred_se
    dat1$age <- conf$ages[age]
    Marginal_df <- rbind(Marginal_df, dat1)
  }
  Marginal_df$age <- as.factor(Marginal_df$age)
  
  psst <- ggplot(Marginal_df, aes_string(x=variable, y="pred", col="age")) + geom_line(size=1.5) + 
    theme_bw() + scale_color_viridis_d() + labs(x=variable, y = "Marginal effect") + coord_cartesian(ylim=c(0, 300))
  ggsave(psst, filename = paste0(folder, "/Marginal_", variable, ".png"), 
         dpi ="retina", width = 6, height = 4, device = "png")
  
}





