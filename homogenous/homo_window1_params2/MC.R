# Minimum-contrast estimation for the bivariate LGCP defined in the Appendix.
# est_BMC_pcf: Estimates parameters using pair correlation functions (pcf).
# est_BMC_K: Estimates parameters using K-functions.

est_BMC_pcf = function(r, par0, g11_hat, g22_hat, g12_hat, weight1, weight2, weight3, q=0.25){
  
  g11_hat[!is.finite(g11_hat) | g11_hat <= 0] <- 1e-6
  g22_hat[!is.finite(g22_hat) | g22_hat <= 0] <- 1e-6
  g12_hat[!is.finite(g12_hat) | g12_hat <= 0] <- 1e-6
  
  cov_exp <- function(r, sigma2, xi) {sigma2 * exp(-r / xi)}
  
  g11_model <- function(r, par) {
    sigma2_y  <- par[1]
    xi_y   <- par[2]
    sigma2_u1 <- par[3]
    xi_u1  <- par[4]
    exp(cov_exp(r, sigma2_y, xi_y) + cov_exp(r, sigma2_u1, xi_u1))
  }
  
  g22_model <- function(r, par) {
    sigma2_y  <- par[1]
    xi_y   <- par[2]
    sigma2_u2 <- par[5]
    xi_u2  <- par[6]
    exp(cov_exp(r, sigma2_y, xi_y) + cov_exp(r, sigma2_u2, xi_u2))
  }
  
  g12_model <- function(r, par) {
    sigma2_y <- par[1]
    xi_y  <- par[2]
    exp(cov_exp(r, sigma2_y, xi_y))
  }
  
  BMC_L <- function(par, r, g11_hat, g22_hat, g12_hat, weight1, weight2, weight3, q = 0.25) {
    term11 <- sum(weight1 *(g11_model(r, par)^q - g11_hat^q)^2)
    term22 <- sum(weight2 *(g22_model(r, par)^q - g22_hat^q)^2)
    term12 <- sum(weight3 *(g12_model(r, par)^q - g12_hat^q)^2)
    
    term11 + term22 + term12
  }
  
  fit <- optim(par = par0, fn = BMC_L, r = r,
               g11_hat = g11_hat, g22_hat = g22_hat, g12_hat = g12_hat,
               weight1 = weight1, weight2 = weight2, weight3 = weight3, q = q, method = "L-BFGS-B",
               lower = c(0.1, 0.001, 0.1, 0.001, 0.1, 0.001),
               upper = c(6, 0.3, 6, 0.3, 6, 0.3)
  )
  
  fit$par
}


est_BMC_K = function(r, par0, K11_hat, K22_hat, K12_hat, weight1, weight2, weight3, q = 0.25) {
  
  K11_hat[!is.finite(K11_hat) | K11_hat <= 0] <- 1e-8
  K22_hat[!is.finite(K22_hat) | K22_hat <= 0] <- 1e-8
  K12_hat[!is.finite(K12_hat) | K12_hat <= 0] <- 1e-8
  
  cov_exp <- function(r, sigma2, xi) {
    sigma2 * exp(-r / xi)
  }
  
  g11_model <- function(r, par) {
    sigma2_y  <- par[1]
    xi_y      <- par[2]
    sigma2_u1 <- par[3]
    xi_u1     <- par[4]
    exp(cov_exp(r, sigma2_y, xi_y) + cov_exp(r, sigma2_u1, xi_u1))
  }
  
  g22_model <- function(r, par) {
    sigma2_y  <- par[1]
    xi_y      <- par[2]
    sigma2_u2 <- par[5]
    xi_u2     <- par[6]
    exp(cov_exp(r, sigma2_y, xi_y) + cov_exp(r, sigma2_u2, xi_u2))
  }
  
  g12_model <- function(r, par) {
    sigma2_y <- par[1]
    xi_y     <- par[2]
    exp(cov_exp(r, sigma2_y, xi_y))
  }
  
  cumtrapz0 <- function(x, y) {
    n <- length(x)
    out <- numeric(n)
    if (n >= 2) {
      dx <- diff(x)
      out[-1] <- cumsum(dx * (head(y, -1) + tail(y, -1)) / 2)
    }
    out
  }
  
  K_from_g <- function(r, g) {
    2 * pi * cumtrapz0(r, r * g)
  }
  
  K11_model <- function(r, par) K_from_g(r, g11_model(r, par))
  K22_model <- function(r, par) K_from_g(r, g22_model(r, par))
  K12_model <- function(r, par) K_from_g(r, g12_model(r, par))
  
  BMC_L <- function(par, r, K11_hat, K22_hat, K12_hat, weight1, weight2, weight3, q = 0.25) {
    term11 <- sum(weight1 * (K11_model(r, par)^q - K11_hat^q)^2)
    term22 <- sum(weight2 * (K22_model(r, par)^q - K22_hat^q)^2)
    term12 <- sum(weight3 * (K12_model(r, par)^q - K12_hat^q)^2)
    term11 + term22 + term12
  }
  
  fit <- optim(par = par0, fn = BMC_L, r = r,
               K11_hat = K11_hat, K22_hat = K22_hat, K12_hat = K12_hat,
               weight1 = weight1, weight2 = weight2, weight3 = weight3,
               q = q,
               method = "L-BFGS-B",
               lower = c(0.1, 0.001, 0.1, 0.001, 0.1, 0.001),
               upper = c(6, 0.3, 6, 0.3, 6, 0.3))
  
  fit$par
}


