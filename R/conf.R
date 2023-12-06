##' Configure data and model for the TMB run
##' @param Model_type defines the study years to keep in the analysis
##' @param conf_model TMB model run configuration i.e. fixing parameters, etc
##' @param keep_omega defines plus group cut-off
##' @param include_age_correlation defines the age groups to focus the analysis
##' @param density_dependence determines whether we include the RE density-dependent effect on fixed effects
##' @param spatio_temporal defines the type of spatio-temporal correlation to use: "AR1", "IID"
##' @param years the study years to keep in the analysis
##' @param plusgroup defines plus group cut-off
##' @param ages defines the age groups to focus the analysis
##' @param formula formula as in lme4 to look at covariate effects. 'Smooth' effect as in gam is possible
##' @param formula_mix formula as in lme4 to look at covariate effects on the mixture probability
##' @param cohort whether to include the cohort effect
##' @param Do_predict Predict outputs or not
##' @param plotting do the model results
##' @details this is the main function that configures the type of TMB model to run, the data to include (years, age), and
##' the general configuration about plotting, prediction, etc
##' @return Configurations to set up the model and the data
##' @export


conf_TMB = function(Model_type, keep_omega,  keep_epsilon=TRUE, include_age_correlation, 
                    years=2010:2020, plusgroup=16, minusgroup = 0,
                    ages=1:10, 
                    formula_pres, formula, formula_mix, cohort = 0L,
                    conf_model = NULL, 
                    mixture_model = TRUE, 
                    family = 2, 
                    link = 1, 
                    ARorIID = c(0,0), 
                    add_nugget = 1,
                    density_dependence = TRUE, 
                    corr_str = 1, 
                    knots_mesh = 250, 
										bias_correct = FALSE,
                    mesh_type = "cutoff", 
                    Do_predict, plotting){
  conf = list()
  conf$Model_type = Model_type
  conf$keep_omega = keep_omega
  conf$keep_epsilon = keep_epsilon
  conf$include_age_correlation = include_age_correlation
  conf$years = years
  conf$plusgroup = plusgroup
  conf$minusgroup = minusgroup
  conf$ages = ages
  conf$formula_pres = formula_pres
  conf$formula = formula
  conf$formula_mix = formula_mix
  conf$family = family
  conf$ARorIID = ARorIID
  conf$add_nugget = add_nugget
  conf$link = link
  conf$mesh_type =  mesh_type
  conf$knots_mesh = knots_mesh
  conf$Do_predict = Do_predict 
  conf$mixture_model = mixture_model 
  conf$density_dependence = density_dependence 
  conf$cohort = cohort 
  conf$corr_str = corr_str 
  conf$plotting = plotting
  conf$save_folder = "New_run"
  conf$years = years
  conf$Nyear = length(years)
  conf$Nage = length(ages)
  conf$conf_model =  conf_model
  conf$bias_correct =  bias_correct

  if (min(ages) < minusgroup) cli_warn(c("the minimum age in ages should be bigger than the minusgroup"))
  if (max(ages) > plusgroup) cli_warn(c("the maximum age in ages should be smaller than the plusgroup"))
  
  # Now checking if the formula contains any random effect terms and split to fixed & random effects  
  attributes(conf)$fixed_effect_formula_pres = nobars(formula_pres)
  attributes(conf)$fixed_effect_formula = nobars(formula)
  if (nobars(formula) != "CPUE ~ -1" & nobars(formula) != "CPUE ~ 1")  attributes(conf)$fixed_effect =  extract_variables(nobars(formula))
  if (nobars(formula) == "CPUE ~ -1" & nobars(formula) == "CPUE ~ 1")  attributes(conf)$fixed_effect =  NULL
  attributes(conf)$RE_effects = NA 
  if (!is.null(findbars(formula))) attributes(conf)$RE_effects = barnames(findbars(formula))
  if (length(grep("VESSEL", attributes(conf)$RE_effects))>0){
    if (grep("VESSEL", attributes(conf)$RE_effects) != 1) cli_warn(c("the VESSEL effect should always be put as the firt random effect (for some internal modeling reason)"))
  }  
  
  attributes(conf)$mixture_formula = nobars(formula_mix)
  return(conf)
}



##' Configure data and model for the TMB run for the 2-step model that models P(age) and aggregated CPUE
##' @param Model_type defines the study years to keep in the analysis
##' @param conf_model TMB model run configuration i.e. fixing parameters, etc
##' @param keep_omega defines plus group cut-off
##' @param include_age_correlation defines the age groups to focus the analysis
##' @param density_dependence determines whether we include the RE density-dependent effect on fixed effects
##' @param spatio_temporal defines the type of spatio-temporal correlation to use: "AR1", "IID"
##' @param years the study years to keep in the analysis
##' @param plusgroup defines plus group cut-off
##' @param ages defines the age groups to focus the analysis
##' @param formula formula as in lme4 to look at covariate effects. 'Smooth' effect as in gam is possible
##' @param formula_mix formula as in lme4 to look at covariate effects on the mixture probability
##' @param cohort whether to include the cohort effect
##' @param Do_predict Predict outputs or not
##' @param plotting do the model results
##' @details this is the main function that configures the type of TMB model to run, the data to include (years, age), and
##' the general configuration about plotting, prediction, etc
##' @return Configurations to set up the model and the data
##' @export


conf_TMB_condlogit = function(Model_type, keep_omega,  keep_epsilon=TRUE, include_age_correlation, years=2010:2020, 
                              plusgroup=16, ages=1:10, minusgroup=0, 
                    formula, formula_cpue, cohort = 0L,
                    conf_model = NULL, 
                    mixture_model = TRUE, 
                    family = 2, 
                    link = 2, 
                    link_cpue = 1, 
                    ARorIID = c(0,0), 
                    ARorIID_cpue = 0, 
                    density_dependence = TRUE, 
                    corr_str = 1, 
                    knots_mesh = 250, 
										bias_correct = FALSE,
                    mesh_type = "cutoff", 
                    Do_predict, plotting){
  conf = list()
  conf$Model_type = Model_type
  conf$keep_omega = keep_omega
  conf$keep_epsilon = keep_epsilon
  conf$include_age_correlation = include_age_correlation
  conf$years = years
  conf$plusgroup = plusgroup
  conf$minusgroup = minusgroup
  conf$ages = ages
  conf$formula = formula
  conf$formula_cpue = formula_cpue
  conf$family = family
  conf$ARorIID = ARorIID
  conf$ARorIID_cpue = ARorIID_cpue
  conf$link = link
  conf$link_cpue = link_cpue
  conf$mesh_type =  mesh_type
  conf$knots_mesh = knots_mesh
  conf$Do_predict = Do_predict 
  conf$mixture_model = mixture_model 
  conf$density_dependence = density_dependence 
  conf$cohort = cohort 
  conf$corr_str = corr_str 
  conf$plotting = plotting
  conf$save_folder = "New_run"
  conf$years = years
  conf$Nyear = length(years)
  conf$Nage = length(ages)
  conf$conf_model =  conf_model
  conf$bias_correct =  bias_correct

  # Now checking if the formula contains any random effect terms and split to fixed & random effects  
  attributes(conf)$fixed_effect_formula_cpue = nobars(formula_cpue)
  attributes(conf)$fixed_effect_formula = nobars(formula)
  attributes(conf)$fixed_effect_cpue =  extract_variables(nobars(formula_cpue))
  # attributes(conf)$fixed_effect =  extract_variables(nobars(formula))
  if (nobars(formula) != "CPUE ~ -1" & nobars(formula) != "CPUE ~ 1") attributes(conf)$fixed_effect =  extract_variables(nobars(formula))
  if (nobars(formula) == "CPUE ~ -1" | nobars(formula) == "CPUE ~ 1") attributes(conf)$fixed_effect =  NULL
  attributes(conf)$RE_effects_cpue = NA 
  attributes(conf)$RE_effects = NA 
  if (!is.null(findbars(formula_cpue))) attributes(conf)$RE_effects_cpue = barnames(findbars(formula_cpue))
  if (!is.null(findbars(formula))) attributes(conf)$RE_effects = barnames(findbars(formula))
  
  if (min(ages) < minusgroup) cli_warn(c("the minimum age in ages should be bigger than the minusgroup"))
  if (max(ages) > plusgroup) cli_warn(c("the maximum age in ages should be smaller than the plusgroup"))
  return(conf)
}

