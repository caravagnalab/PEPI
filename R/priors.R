#' Fit a Gamma distribution to two quantiles.
#'
#' Finds the shape and rate of the Gamma distribution whose `p1` and `p2`
#' quantiles are `q1` and `q2`, by least squares on the CDF. Useful to turn
#' a plausible range of a rate (e.g. a switch rate) into Gamma prior
#' hyperparameters.
#'
#' @param q1,q2 Lower and upper quantile values.
#' @param p1,p2 Probabilities of the two quantiles.
#' @return A list with elements `shape` and `rate`.
#' @examples
#' g = fit_gamma_quantiles(0.01, 0.5)
#' stats::qgamma(c(0.05, 0.95), shape = g$shape, rate = g$rate)
#' @export

fit_gamma_quantiles = function(q1, q2, p1 = 0.05, p2 = 0.95){

  # moment-matching starting point
  mu0 = (q1 + q2)/2
  sd0 = (q2 - q1)/3.29

  shape0 = mu0^2 / sd0^2
  rate0 = mu0 / sd0^2

  objective = function(par){

    shape = exp(par[1])
    rate = exp(par[2])

    err1 = stats::pgamma(q1, shape = shape, rate = rate) - p1
    err2 = stats::pgamma(q2, shape = shape, rate = rate) - p2

    err1^2 + err2^2
  }

  fit = stats::optim(par = log(c(shape0, rate0)), fn = objective)

  list(shape = exp(fit$par[1]), rate = exp(fit$par[2]))

}
