##' Extract outputs from model run 
##' @param Model_type defines the study years to keep in the analysis
##' @param keep_omega defines plus group cut-off
##' @param include_age_correlation defines the age groups to focus the analysis
##' @param years the study years to keep in the analysis
##' @param plusgroup defines plus group cut-off
##' @param ages defines the age groups to focus the analysis
##' @param formula formula as in lme4 to look at covariate effects. 'Smooth' effect as in gam is possible
##' @param Do_predict Predict outputs or not
##' @param plotting do the model results
##' @details this is the main function that configures the type of TMB model to run, the data to include (years, age), and
##' the general configuration about plotting, prediction, etc
##' @return Configurations to set up the model and the data
##' @export extract_output

extract_output <- function(run, data){

    pred_grid <- data$pred_grid
  # extract the estimates predicted values
    mu_proj2 <- run$rep$mu_proj
    test2 <- data.table::as.data.table(mu_proj2)
    colnames(test2) <- c("LOC", "YEAR", "AGE", "Pred")
    test2$YEAR <- as.factor(test2$YEAR)
    test2$YEAR <- as.factor(as.numeric(as.character(factor(test2$YEAR, labels=sort(unique(data$tmb_data$YEAR))))))
    pred_grid$YEAR <- as.factor(as.numeric(as.character(factor(pred_grid$YEAR, labels=sort(unique(data$tmb_data$YEAR))))))
    pred_grid$LOC <- rep(1:(nrow(pred_grid)/conf$Nyear), conf$Nyear)
    pred <- test2 %>% left_join(pred_grid)
  
  # extract the estimated spatio-temporal field
    epsilon2 <- run$rep$epsilon_st_A_proj
    qwe2 <- data.table::as.data.table(epsilon2)
    colnames(qwe2) <- c("LOC", "YEAR", "AGE", "Pred")
    qwe2$YEAR <- as.factor(qwe2$YEAR)
    qwe2$YEAR <- as.factor(as.numeric(as.character(factor(qwe2$YEAR, labels=sort(unique(data$tmb_data$YEAR))))))
    RE_st <- qwe2 %>% left_join(pred_grid)
  
  out <- list(predicted=pred, RE_st = RE_st)
    
  return(out)
}

