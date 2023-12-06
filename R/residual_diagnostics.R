qres_tweedie <- function(object, y, mu) {
  p <- stats::plogis(object$thetaf) + 1
  dispersion <- exp(object$ln_phi)
  u <- fishMod::pTweedie(q = y, p = p, mu = mu, phi = dispersion)
  if (p > 1 && p < 2)
    u[y == 0] <- stats::runif(sum(y == 0), min = 0, max = u[y == 0])
  stats::qnorm(u)
}


qres_gamma <- function(object, y, mu) {
  phi <- exp(object$ln_phi)
  s1 <- phi
  s2 <- mu / s1
  u <- stats::pgamma(q = y, shape = s1, scale = s2)
  stats::qnorm(u)
}

qres_gaussian <- function(object, y, mu) {
  dispersion <- exp(object$ln_phi)
  u <- stats::pnorm(q = y, mean = mu, sd = dispersion)
  stats::qnorm(u)
}

qres_lognormal <- function(object, y, mu) {
  dispersion <- exp(object$ln_phi)
  u <- stats::plnorm(q = y, mean = log(mu) - dispersion^2/2, sd = dispersion)
  stats::qnorm(u)
}

# https://en.wikipedia.org/wiki/Location%E2%80%93scale_family
pt_ls <- function(q, df, mu, sigma) stats::pt((q - mu)/sigma, df)

qres_student <- function(object, y, mu) {
  dispersion <- exp(object$ln_phi)
  u <- pt_ls(q = y, df = object$df, mu = mu, sigma = dispersion)
  stats::qnorm(u)
}


#' @export
#' @importFrom stats predict
get.residuals <- function(object, ...) {
  res_func <- switch(object$family,
                     gaussian = qres_gaussian,
                     binomial = qres_binomial,
                     tweedie  = qres_tweedie,
                     Beta     = qres_beta,
                     gamma    = qres_gamma,
                     nbinom2  = qres_nbinom2,
                     poisson  = qres_pois,
                     student  = qres_student,
                     lognormal= qres_lognormal
  )
  y <- object$response
  mu <- object$mu
  res_func(object, y, mu, ...)
}






