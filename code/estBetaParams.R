estBetaParams <- function(mu, var) {
  #  shape parameters a > 0; b > 0
  #  y ~ dbeta(a, b)
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}