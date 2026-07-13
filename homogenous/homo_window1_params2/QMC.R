
Intg_diag = function(s, sg1, t1, sg2, t2){
  if (any(s < 0)){
    stop("Distance should be non-negative.")
  }
  C1 = sg1^2 * exp( - s / t1)
  C2 = sg2^2 * exp( - s / t2)
  result = s * (exp(C1 + C2))
  
  return(result)
}

Intg_diag_s1 = function(s, sg1, t1, sg2, t2){
  if (any(s < 0)){
    stop("Distance should be non-negative.")
  }
  
  C1 = sg1^2 * exp( - s / t1)
  C2 = sg2^2 * exp( - s / t2)
  result = s * (exp(C1 + C2)) * 2 * sg1 * exp(- s / t1)
  
  return(result)
}

Intg_diag_t1 = function(s, sg1, t1, sg2, t2){
  if (any(s < 0)){
    stop("Distance should be non-negative.")
  }
  
  C1 = sg1^2 * exp( - s / t1)
  C2 = sg2^2 * exp( - s / t2)
  result = s^2 * (exp(C1 + C2)) * C1 / (t1^2)
  
  return(result)
}

Intg_diag_s2 = function(s, sg1, t1, sg2, t2){
  if (any(s < 0)){
    stop("Distance should be non-negative.")
  }
  
  C1 = sg1^2 * exp( - s / t1)
  C2 = sg2^2 * exp( - s / t2)
  result = s * (exp(C1 + C2)) * 2 * sg2 * exp(- s / t2)
  
  return(result)
}


Intg_diag_t2 = function(s, sg1, t1, sg2, t2){
  if (any(s < 0)){
    stop("Distance should be non-negative.")
  }
  
  C1 = sg1^2 * exp( - s / t1)
  C2 = sg2^2 * exp( - s / t2)
  result = s^2 * (exp(C1 + C2)) * C2 / (t2^2)
  
  return(result)
}


Intg_off = function(s, sigma, theta, Corr){
  if (any(s < 0)){
    stop("Distance should be non-negative.")
  }

  if (Corr == "Neg"){
    C = - sigma^2 * exp( - s / theta)   # cov(log(Lambda1), log(Lambda2)) < 0
  }else if (Corr == "Pos"){
    C = sigma^2 * exp( - s / theta)   # cov(log(Lambda1), log(Lambda2)) > 0
  }else{
    stop("Corr can only be Neg or Pos.")
  }
  result = s * (exp(C))
  
  return(result)
}


Intg_off_s = function(s, sigma, theta, Corr){
  if (any(s < 0)){
    stop("Distance should be non-negative.")
  }
  
  if (Corr == "Neg"){
    C = - sigma^2 * exp( - s / theta)   # cov(log(Lambda1), log(Lambda2)) < 0
    result = - s * (exp(C)) * 2 * sigma * exp( - s / theta)
  }else if (Corr == "Pos"){
    C = sigma^2 * exp( - s / theta)   # cov(log(Lambda1), log(Lambda2)) > 0
    result = s * (exp(C)) * 2 * sigma * exp( - s / theta)
  }else{
    stop("Corr can only be Neg or Pos.")
  }
  
  return(result)
}


Intg_off_t = function(s, sigma, theta, Corr){
  if (any(s < 0)){
    stop("Distance should be non-negative.")
  }
  
  if (Corr == "Neg"){
    C = - sigma^2 * exp( - s / theta)   # cov(log(Lambda1), log(Lambda2)) < 0
  }else if (Corr == "Pos"){
    C = sigma^2 * exp( - s / theta)   # cov(log(Lambda1), log(Lambda2)) > 0
  }else{
    stop("Corr can only be Neg or Pos.")
  }
  result = s^2 * (exp(C)) * C / (theta^2)
  
  return(result)
}



#' @title Theoretical Ripley's K function matrix
#' @description
#' calculate the integrand at a sequence of distances
#' the first 2 columns: diagonal terms of Ripley's K
#' the 3rd column: off-diagonal terms of Ripley's K
#' dimension: (length(R) - 5) by 3 because we don't calculate the values near distance 0
#'
#' @param sigma1 sigma1^2: variance parameter of Y1
#' @param theta1 spatial scale of Y1
#' @param sigma2 sigma2^2: variance parameter of Y2
#' @param theta2 spatial scale of Y2
#' @param sigma3 sigma3^2: variance parameter of Y3
#' @param theta3 spatial scale of Y3
#' @param R vector of a sequence of distances
#' @param Corr character of correlation of log intensities of two point patterns (negative (default) or positive)
#'
#' @return (length(R) - 5) by 3 matrix, theoretical Ripley's K matrix
#' @export

Ktheo = function(sigma1, theta1, sigma2, theta2, sigma3, theta3, R, Corr){
  ## compatibility check
  # if (sigma1 <=0 | theta1 <=0 | sigma2 <= 0 |theta2 <= 0 | sigma3 <= 0 | theta3 <= 0){
  #   stop(paste("Parameters sigma_i and theta_i should be positive."))
  # }
  
  if (any(duplicated(R))){
    stop("The vector R contains repetitive elements.")
  }
  
  if (min(R) < 0){
    stop("Components of R should be non-negative.")
  }
  
  temp = matrix(NA, length(R), 3)
  delta = c(R[1],diff(R))
  temp[, 1] = 2 * pi * cumsum(Intg_diag(s = R, sg1 = sigma1, t1 = theta1,
                                        sg2 = sigma3, t2 = theta3) * delta)
  temp[, 2] = 2 * pi * cumsum(Intg_diag(s = R, sg1 = sigma2, t1 = theta2,
                                        sg2 = sigma3, t2 = theta3) * delta)
  temp[, 3] = 2 * pi * cumsum(Intg_off(s = R, sigma = sigma3, theta = theta3, Corr) * delta)
  
  return(temp[-c(1:5), ])
}



#' @title Estimated Ripley's K function matrix
#' @description
#' Calculate the estimated Ripley's K and cross-K values at a sequence of distances
#' the first 2 columns: diagonal terms of estimated Ripley's K
#' the 3rd column: off-diagonal terms of estimated Ripley's K
#' dimension: (length(R) - 5) by 3 because we don't calculate the values near distance 0
#'
#' @param data bivariate point pattern, ppp object with two marks
#' @param rmax Maximum desired value of the distances sequence r
#'
#' @return R: the distances sequence; Kh: (length(R) - 5) by 3 matrix, estimated Ripley's K matrix
#' @export
#' @importFrom spatstat Kest Kcross is.ppp

Khat = function(data, rmax, correction = "isotropic", Inhom = F){
  ## compatibility check
  if (!is.ppp(data) | is.null(data$marks) | length(levels(data$marks)) != 2){
    stop("data should be an object of ppp with two marks.")
    #  should correspond to a bivariate point pattern
  }
  
  if (rmax <= 0){
    stop("rmax should be positive.")
  }
  
  r = seq(0, rmax, length = 513)
  
  ## diagonal estimated K
  if (Inhom == F){
    K1 = Kest(split(data)[[1]], r = r, 
              nlarge = Inf, correction = correction)
    K2 = Kest(split(data)[[2]], r = r, 
              nlarge = Inf, correction = correction)
    ## off-diagonal estimated cross-K
    Kc1 = Kcross(data, r = r, nlarge = Inf,
                 i = levels(data$marks)[1],
                 j = levels(data$marks)[2],correction = correction)
    Kc2 = Kcross(data, r = r, nlarge = Inf,
                 i = levels(data$marks)[2],
                 j = levels(data$marks)[1],correction = correction)
  }else{
    K1 = Kinhom(split(data)[[1]], r = r, lambda = density(split(data)[[1]]),
                nlarge = Inf, correction = correction)
    K2 = Kinhom(split(data)[[2]], r = r, lambda = density(split(data)[[2]]),
                nlarge = Inf, correction = correction)
    ## off-diagonal estimated cross-K
    Kc1 = Kcross.inhom(data, r = r, 
                       i = levels(data$marks)[1],
                       j = levels(data$marks)[2],correction = correction,
                       lambdaI=density(split(data)[[1]]),lambdaJ=density(split(data)[[2]]))
    Kc2 = Kcross.inhom(data, r = r, 
                       i = levels(data$marks)[2],
                       j = levels(data$marks)[1],correction = correction,
                       lambdaI=density(split(data)[[2]]),lambdaJ=density(split(data)[[1]]))
  }
  
  n1 = length(which(data$marks == levels(data$marks)[1]))   # number of the points in the 1st point pattern
  n2 = length(which(data$marks == levels(data$marks)[2]))   # number of the points in the 2nd point pattern
  
  ## estimated K matrix
  r = K1$r                            # K1, K2, Kc1, Kc2 should have the same distances sequence r
  temp = matrix(NA, length(r), 3)
  ## Consider two common edge correction methods "isotropic" and "border"
  if (correction == "isotropic"){
    temp[ , 1] = as.vector(K1$iso)
    temp[ , 2] = as.vector(K2$iso)
    temp[ , 3] = as.vector(Kc1$iso * n2 + Kc2$iso * n1) / (n1 + n2)
  }
  if (correction == "border"){
    temp[ , 1] = as.vector(K1$border)
    temp[ , 2] = as.vector(K2$border)
    temp[ , 3] = as.vector(Kc1$border * n2 + Kc2$border * n1) / (n1 + n2)
  }
  if (correction == "bord.modif"){
    temp[ , 1] = as.vector(K1$bord.modif)
    temp[ , 2] = as.vector(K2$bord.modif)
    temp[ , 3] = as.vector(Kc1$bord.modif * n2 + Kc2$bord.modif * n1) / (n1 + n2)
  }
  
  return(list(Kh = temp[-c(1:5), ], R = r))
}



#' @title Derivative of each entry of theoretical Ripley's K function matrix
#' with respect to model parameters (sigma1, theta1, sigma2, theta2, sigma3, theta3)
#' It's a necessary term for estimating S(theta) in equation (40)
#' @description
#' calculate the integrand at a sequence of distances
#' K11_deriv: derivative of the 1st diagonal term of Ripley's K matrix
#' K22_deriv: derivative of the 2nd diagonal term of Ripley's K matrix
#' K12_deriv: derivative of the diagonal term of Ripley's K matrix
#' dimension: 3 matrices of (length(R) - 5) by 3 because we don't calculate the values near distance 0
#'
#' @param sigma1 sigma1^2: variance parameter of Y1
#' @param theta1 spatial scale of Y1
#' @param sigma2 sigma2^2: variance parameter of Y2
#' @param theta2 spatial scale of Y2
#' @param sigma3 sigma3^2: variance parameter of Y3
#' @param theta3 spatial scale of Y3
#' @param R vector of a sequence of distances
#' @param Corr character of correlation of log intensities of two point patterns (negative (default) or positive)
#'
#' @return 3 matrices of (length(R) - 5) by 3, derivative of theoretical Ripley's K matrix
#' @export

Ktheo_deriv = function(sigma1, theta1, sigma2, theta2, sigma3, theta3, R, Corr){
  K11_deriv = K12_deriv = K22_deriv = matrix(0, length(R), 6)
  delta = c(R[1],diff(R))
  K11_deriv[, 1] = 2 * pi * cumsum(Intg_diag_s1(s = R, sg1 = sigma1, t1 = theta1,
                                                sg2 = sigma3, t2 = theta3) * delta)
  K11_deriv[, 2] = 2 * pi * cumsum(Intg_diag_t1(s = R, sg1 = sigma1, t1 = theta1,
                                                sg2 = sigma3, t2 = theta3) * delta)
  K11_deriv[, 5] = 2 * pi * cumsum(Intg_diag_s2(s = R, sg1 = sigma1, t1 = theta1,
                                                sg2 = sigma3, t2 = theta3) * delta)
  K11_deriv[, 6] = 2 * pi * cumsum(Intg_diag_t2(s = R, sg1 = sigma1, t1 = theta1,
                                                sg2 = sigma3, t2 = theta3) * delta)
  
  K22_deriv[, 3] = 2 * pi * cumsum(Intg_diag_s1(s = R, sg1 = sigma2, t1 = theta2,
                                                sg2 = sigma3, t2 = theta3) * delta)
  K22_deriv[, 4] = 2 * pi * cumsum(Intg_diag_t1(s = R, sg1 = sigma2, t1 = theta2,
                                                sg2 = sigma3, t2 = theta3) * delta)
  K22_deriv[, 5] = 2 * pi * cumsum(Intg_diag_s2(s = R, sg1 = sigma2, t1 = theta2,
                                                sg2 = sigma3, t2 = theta3) * delta)
  K22_deriv[, 6] = 2 * pi * cumsum(Intg_diag_t2(s = R, sg1 = sigma2, t1 = theta2,
                                                sg2 = sigma3, t2 = theta3) * delta)
  
  K12_deriv[, 5] = 2 * pi * cumsum(Intg_off_s(s = R, sigma = sigma3, theta = theta3, Corr) * delta)
  K12_deriv[, 6] = 2 * pi * cumsum(Intg_off_t(s = R, sigma = sigma3, theta = theta3, Corr) * delta)
  return(list(K11_deriv = K11_deriv[-c(1:5), ], K22_deriv = K22_deriv[-c(1:5), ], K12_deriv = K12_deriv[-c(1:5), ]))
}


#' @title Derivative of each entry of theoretical Q function matrix
#' with respect to model parameters (sigma1, theta1, sigma2, theta2, sigma3, theta3, beta1, beta2)
#' It's a necessary term for estimating S(theta) in equation (40)
#' @description
#' calculate the integrand at a sequence of distances
#' Q11_deriv: derivative of the 1st diagonal term of Q matrix
#' Q22_deriv: derivative of the 2nd diagonal term of Q matrix
#' Q12_deriv: derivative of the diagonal term of Q matrix
#' dimension: 3 matrices of (length(R) - 5) by 3 because we don't calculate the values near distance 0
#'
#' @param sigma1 sigma1^2: variance parameter of Y1
#' @param theta1 spatial scale of Y1
#' @param sigma2 sigma2^2: variance parameter of Y2
#' @param theta2 spatial scale of Y2
#' @param sigma3 sigma3^2: variance parameter of Y3
#' @param theta3 spatial scale of Y3
#' @param beta1 lambda1 = exp(beta1)
#' @param beta2 lambda2 = exp(beta2)
#' @param R vector of a sequence of distances
#' @param Corr character of correlation of log intensities of two point patterns (negative (default) or positive)
#'
#' @return 3 matrices of (length(R) - 5) by 3, derivative of theoretical Ripley's K matrix
#' @export

Qtheo_deriv = function(sigma1, theta1, sigma2, theta2, sigma3, theta3, beta1, beta2, R, Corr){
  Q11_deriv = Q12_deriv = Q22_deriv = matrix(0, length(R), 8)
  delta = c(R[1],diff(R))
  Q11_deriv[, 1] = 2 * pi * cumsum(Intg_diag_s1(s = R, sg1 = sigma1, t1 = theta1,
                                                sg2 = sigma3, t2 = theta3) * delta) * exp(2 * beta1)
  Q11_deriv[, 2] = 2 * pi * cumsum(Intg_diag_t1(s = R, sg1 = sigma1, t1 = theta1,
                                                sg2 = sigma3, t2 = theta3) * delta) * exp(2 * beta1)
  Q11_deriv[, 5] = 2 * pi * cumsum(Intg_diag_s2(s = R, sg1 = sigma1, t1 = theta1,
                                                sg2 = sigma3, t2 = theta3) * delta) * exp(2 * beta1)
  Q11_deriv[, 6] = 2 * pi * cumsum(Intg_diag_t2(s = R, sg1 = sigma1, t1 = theta1,
                                                sg2 = sigma3, t2 = theta3) * delta) * exp(2 * beta1)
  Q11_deriv[, 7] = 2 * pi * cumsum(Intg_diag(s = R, sg1 = sigma1, t1 = theta1,
                                             sg2 = sigma3, t2 = theta3) * delta) * 2 * exp(2 * beta1)
  
  
  Q22_deriv[, 3] = 2 * pi * cumsum(Intg_diag_s1(s = R, sg1 = sigma2, t1 = theta2,
                                                sg2 = sigma3, t2 = theta3) * delta) * exp(2 * beta2)
  Q22_deriv[, 4] = 2 * pi * cumsum(Intg_diag_t1(s = R, sg1 = sigma2, t1 = theta2,
                                                sg2 = sigma3, t2 = theta3) * delta) * exp(2 * beta2)
  Q22_deriv[, 5] = 2 * pi * cumsum(Intg_diag_s2(s = R, sg1 = sigma2, t1 = theta2,
                                                sg2 = sigma3, t2 = theta3) * delta) * exp(2 * beta2)
  Q22_deriv[, 6] = 2 * pi * cumsum(Intg_diag_t2(s = R, sg1 = sigma2, t1 = theta2,
                                                sg2 = sigma3, t2 = theta3) * delta) * exp(2 * beta2)
  Q22_deriv[, 8] = 2 * pi * cumsum(Intg_diag(s = R, sg1 = sigma2, t1 = theta2,
                                             sg2 = sigma3, t2 = theta3) * delta) * 2 * exp(2 * beta2)
  
  Q12_deriv[, 5] = 2 * pi * cumsum(Intg_off_s(s = R, sigma = sigma3, theta = theta3, Corr) * delta) * exp(beta1 + beta2)
  Q12_deriv[, 6] = 2 * pi * cumsum(Intg_off_t(s = R, sigma = sigma3, theta = theta3, Corr) * delta) * exp(beta1 + beta2)
  Q12_deriv[, 7] = 2 * pi * cumsum(Intg_off(s = R, sigma = sigma3, theta = theta3, Corr) * delta) * exp(beta1 + beta2)
  Q12_deriv[, 8] = 2 * pi * cumsum(Intg_off(s = R, sigma = sigma3, theta = theta3, Corr) * delta) * exp(beta1 + beta2)
  
  return(list(Q11_deriv = Q11_deriv[-c(1:5), ], Q22_deriv = Q22_deriv[-c(1:5), ], Q12_deriv = Q12_deriv[-c(1:5), ]))
}


#' @title L2 type discrepancy between theoretical and estimated K
#'
#' @param data bivariate point pattern
#' @param par vector of log transfrom of sigma1,theta1,sigma2,theta2,sigma3,theta3
#' @param power transform applied to theoretical or estimated K
#' @param rmax Maximum desired value of the distances sequence r
#' @param Corr character of correlation of log intensities of two point patterns (negative (default) or positive)
#'
#' @return a scalar of L2 type discrepancy between theoretical and estimated K
#' @export

L = function(par, data, power, rmax, Corr, correction = "isotropic", Inhom = F, 
             power_equal = T, power1 = NULL, power2 = NULL, power3 = NULL){
  
  sigma1 = exp(par[1])
  theta1 = exp(par[2])
  sigma2 = exp(par[3])
  theta2 = exp(par[4])
  sigma3 = exp(par[5])
  theta3 = exp(par[6])
  
  Ke = Khat(data, rmax, correction = correction, Inhom = Inhom)
  R = Ke$R
  
  Kt = Ktheo(sigma1, theta1, sigma2, theta2, sigma3, theta3, R, Corr)
  Kh = Ke$Kh
  if (power_equal == T){
    diff2 = ((Kt)^{power} - (Kh)^{power})^2
    diff2[ , 3] = diff2[ , 3] * 2
    Emp = (Kh)^(power * 2)
    Emp[ , 3] = Emp[ , 3] * 2
  }else{
    diff2 = Emp = matrix(NA, nrow(Kt), ncol(Kt))
    diff2[, 1] = ((Kt[, 1])^{power1} - (Kh[, 1])^{power1})^2
    diff2[, 2] = ((Kt[, 2])^{power2} - (Kh[, 2])^{power2})^2
    diff2[, 3] = ((Kt[, 3])^{power3} - (Kh[, 3])^{power3})^2
    diff2[ , 3] = diff2[ , 3] * 2
    Emp[, 1] = (Kh[, 1])^(power1 * 2)
    Emp[, 2] = (Kh[, 2])^(power2 * 2)
    Emp[, 3] = (Kh[, 3])^(power3 * 2) 
    Emp[ , 3] = Emp[ , 3] * 2
  }
  result = log(sum(diff2) * R[2])
  
  return((result))
}


#' @title L2 type discrepancy between unscaled theoretical and estimated K
#'
#' @param lm1 intensity of the first point pattern
#' @param lm2 intensity of the second point pattern
#' @param data bivariate point pattern
#' @param par vector of log transfrom of sigma1,theta1,sigma2,theta2,sigma3,theta3
#' @param power transform applied to theoretical or estimated K
#' @param rmax Maximum desired value of the distances sequence r
#' @param Corr character of correlation of log intensities of two point patterns (negative (default) or positive)
#'
#' @return a scalar of L2 type discrepancy between theoretical and estimated K
#' @export

Q = function(par, factor, data, power, rmax, Corr, correction = "isotropic", Inhom = F, 
             power_equal = T, power1 = NULL, power2 = NULL, power3 = NULL){
  
  sigma1 = exp(par[1])
  theta1 = exp(par[2])
  sigma2 = exp(par[3])
  theta2 = exp(par[4])
  sigma3 = exp(par[5])
  theta3 = exp(par[6])
  beta1 = par[7]
  beta2 = par[8]
  
  Ke = Khat(data, rmax, correction = correction, Inhom = Inhom)
  R = Ke$R
  Kh = Ke$Kh
  
  Kt = Ktheo(sigma1, theta1, sigma2, theta2, sigma3, theta3, R, Corr)
  # in case of large values, so we scale by factor
  Qt = Kt %*% diag(c(exp(beta1)^2, exp(beta2)^2, exp(beta1 + beta2)) / factor)
  
  X1 = split(data)[[1]]
  npts1 = npoints(X1)
  lmh1 = npts1 * (npts1 - 1)/((area.owin(Window(X1)))^2)
  X2 = split(data)[[2]]
  npts2 = npoints(X2)
  lmh2 = npts2 * (npts2 - 1)/((area.owin(Window(X2)))^2)
  lmh12 = (npts1 * npts2) / ( area.owin(Window(X1)) * area.owin(Window(X2)) )
  Qh = Kh %*% diag(c(lmh1, lmh2, lmh12) / factor) # scale by factor in case of large value
  
  if (power_equal == T){
    diff2 = ((Qt)^{power} - (Qh)^{power})^2
    diff2[ , 3] = diff2[ , 3] * 2
    Emp = (Qh)^(power * 2)
    Emp[ , 3] = Emp[ , 3] * 2
  }else{
    diff2 = Emp = matrix(NA, nrow(Qt), ncol(Qt))
    diff2[, 1] = ((Qt[, 1])^{power1} - (Qh[, 1])^{power1})^2
    diff2[, 2] = ((Qt[, 2])^{power2} - (Qh[, 2])^{power2})^2
    diff2[, 3] = ((Qt[, 3])^{power3} - (Qh[, 3])^{power3})^2
    diff2[ , 3] = diff2[ , 3] * 2
    Emp[, 1] = (Qh[, 1])^(power1 * 2)
    Emp[, 2] = (Qh[, 2])^(power2 * 2)
    Emp[, 3] = (Qh[, 3])^(power3 * 2) 
    Emp[ , 3] = Emp[ , 3] * 2
  }
  result = log(sum(diff2) * R[2]) 
  
  
  return((result))
}

## Calculate covariance matrix of model parameters (use scaled K function)
cal_Stheta_useK = function(par, data, power, rmax, Corr, correction = "isotropic", Inhom = F){
  
  sigma1 = exp(par[1])
  theta1 = exp(par[2])
  sigma2 = exp(par[3])
  theta2 = exp(par[4])
  sigma3 = exp(par[5])
  theta3 = exp(par[6])
  
  Ke = Khat(data, rmax, correction = correction, Inhom = Inhom)
  R = Ke$R
  Kt = Ktheo(sigma1, theta1, sigma2, theta2, sigma3, theta3, R, Corr)
  Kh = Ke$Kh
  Kt_deriv = Ktheo_deriv(sigma1, theta1, sigma2, theta2, sigma3, theta3, R, Corr)
  
  Jt = Kt 
  Jh = Kh
  J_deriv11 = Kt_deriv$K11_deriv
  J_deriv22 = Kt_deriv$K22_deriv
  J_deriv12 = J_deriv21 = Kt_deriv$K12_deriv
  
  Cb11 = (Jh[,1] - Jt[,1]) * (Jt[,1])^(2 * power - 2) * J_deriv11
  Cb22 = (Jh[,2] - Jt[,2]) * (Jt[,2])^(2 * power - 2) * J_deriv22
  Cb12 = Cb21 = (Jh[,3] - Jt[,3]) * (Jt[,3])^(2 * power - 2) * J_deriv12
  
  Cb11 = colSums(Cb11) * diff(R)[1]; Cb22 = colSums(Cb22) * diff(R)[1]; Cb12 = Cb21 = colSums(Cb12) * diff(R)[1]
  return(list(Cb11 = Cb11, Cb22 = Cb22, Cb12 = Cb12, Cb21 = Cb21))
}

cal_Btheta_useK = function(par, power, rmax, Corr, correction = "isotropic", Inhom = F){
  
  sigma1 = exp(par[1])
  theta1 = exp(par[2])
  sigma2 = exp(par[3])
  theta2 = exp(par[4])
  sigma3 = exp(par[5])
  theta3 = exp(par[6])
  
  R = seq(0, rmax, length = 513)
  Kt = Ktheo(sigma1, theta1, sigma2, theta2, sigma3, theta3, R, Corr)
  Kt_deriv = Ktheo_deriv(sigma1, theta1, sigma2, theta2, sigma3, theta3, R, Corr)
  
  Jt = Kt
  J_deriv11 = Kt_deriv$K11_deriv
  J_deriv22 = Kt_deriv$K22_deriv
  J_deriv12 = J_deriv21 = Kt_deriv$K12_deriv
  
  Cb11 = t(J_deriv11) %*% diag((Jt[,1])^(2 * power - 2)) %*% J_deriv11 * diff(R)[1] 
  Cb22 = t(J_deriv22) %*% diag((Jt[,2])^(2 * power - 2)) %*% J_deriv22 * diff(R)[1]
  Cb12 = t(J_deriv12) %*% diag((Jt[,3])^(2 * power - 2)) %*% J_deriv12 * diff(R)[1]
  
  result = Cb11 + Cb22 + 2 * Cb12
  
  return(result)
}

## Calculate covariance matrix of model parameters (use unscaled K function)
cal_Stheta_useQ = function(par, factor, data, power, rmax, Corr, correction = "isotropic", Inhom = F){
  
  sigma1 = exp(par[1])
  theta1 = exp(par[2])
  sigma2 = exp(par[3])
  theta2 = exp(par[4])
  sigma3 = exp(par[5])
  theta3 = exp(par[6])
  beta1 = par[7]
  beta2 = par[8]
  
  Ke = Khat(data, rmax, correction = correction, Inhom = Inhom)
  R = Ke$R
  Kt = Ktheo(sigma1, theta1, sigma2, theta2, sigma3, theta3, R, Corr)
  Qt_deriv = Qtheo_deriv(sigma1, theta1, sigma2, theta2, sigma3, theta3, beta1, beta2, R, Corr)
  Jt = Kt %*% diag(c(exp(beta1)^2, exp(beta2)^2, exp(beta1 + beta2)) / factor)
  J_deriv11 = Qt_deriv$Q11_deriv / factor
  J_deriv22 = Qt_deriv$Q22_deriv / factor
  J_deriv12 = J_deriv21 = Qt_deriv$Q12_deriv / factor
  
  Kh = Ke$Kh
  X1 = split(data)[[1]]
  npts1 = npoints(X1)
  lmh1 = npts1 * (npts1 - 1)/((area.owin(Window(X1)))^2)
  X2 = split(data)[[2]]
  npts2 = npoints(X2)
  lmh2 = npts2 * (npts2 - 1)/((area.owin(Window(X2)))^2)
  lmh12 = (npts1 * npts2) / ( area.owin(Window(X1)) * area.owin(Window(X2)) )
  Jh = Kh %*% diag(c(lmh1, lmh2, lmh12) / factor) # in case of large values
  
  Cb11 = (Jh[,1] - Jt[,1]) * (Jt[,1])^(2 * power - 2) * J_deriv11
  Cb22 = (Jh[,2] - Jt[,2]) * (Jt[,2])^(2 * power - 2) * J_deriv22
  Cb12 = Cb21 = (Jh[,3] - Jt[,3]) * (Jt[,3])^(2 * power - 2) * J_deriv12
  
  Cb11 = colSums(Cb11) * diff(R)[1]; Cb22 = colSums(Cb22) * diff(R)[1]; Cb12 = Cb21 = colSums(Cb12) * diff(R)[1]
  return(list(Cb11 = Cb11, Cb22 = Cb22, Cb12 = Cb12, Cb21 = Cb21))
}

cal_Btheta_useQ = function(par, factor, power, rmax, Corr, correction = "isotropic", Inhom = F){
  
  sigma1 = exp(par[1])
  theta1 = exp(par[2])
  sigma2 = exp(par[3])
  theta2 = exp(par[4])
  sigma3 = exp(par[5])
  theta3 = exp(par[6])
  beta1 = par[7]
  beta2 = par[8]
  
  R = seq(0, rmax, length = 513)
  Kt = Ktheo(sigma1, theta1, sigma2, theta2, sigma3, theta3, R, Corr)
  Qt_deriv = Qtheo_deriv(sigma1, theta1, sigma2, theta2, sigma3, theta3, beta1, beta2, R, Corr)
  Jt = Kt %*% diag(c(exp(beta1)^2, exp(beta2)^2, exp(beta1 + beta2)) / factor)
  J_deriv11 = Qt_deriv$Q11_deriv / factor
  J_deriv22 = Qt_deriv$Q22_deriv / factor
  J_deriv12 = J_deriv21 = Qt_deriv$Q12_deriv / factor
  
  Cb11 = t(J_deriv11) %*% diag((Jt[,1])^(2 * power - 2)) %*% J_deriv11 * diff(R)[1] 
  Cb22 = t(J_deriv22) %*% diag((Jt[,2])^(2 * power - 2)) %*% J_deriv22 * diff(R)[1]
  Cb12 = t(J_deriv12) %*% diag((Jt[,3])^(2 * power - 2)) %*% J_deriv12 * diff(R)[1]
  
  result = Cb11 + Cb22 + 2 * Cb12
  
  return(result)
}

