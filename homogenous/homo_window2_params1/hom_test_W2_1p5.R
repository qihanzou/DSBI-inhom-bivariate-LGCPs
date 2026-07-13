set.seed(2026)
t0 = Sys.time()
library(spatstat.geom)
library(spatstat.explore)
library(spatstat.random)
library(spatstat.model)
library(reticulate)
library(future)
library(future.apply)
library(ggplot2)
library(terra)
library(abind)
library(progressr)
library(INLA)
library(inlabru)
library(sf)
library(fmesher)
source("C:/Users/qihan/Desktop/hom_W2/MC.R")
source("C:/Users/qihan/Desktop/hom_W2/QMC.R")
use_python("C:/Users/qihan/anaconda3/envs/py39env/python.exe", required = TRUE)
plan(multisession, workers = parallelly::availableCores() - 1)

simulate_grf = function(nx, ny, L, range, var) {
  mu0 = as.im(0, owin(xrange = c(0, L), yrange = c(0, L)), dimyx = c(ny, nx))
  sim_field = rLGCP(model = "exponential", mu = mu0, var = var, scale = range, win = as.owin(mu0), saveLambda = TRUE)
  field_im = attr(sim_field, "Lambda")
  field_im = eval.im(log(field_im))
  as.matrix(field_im)
}

simulate_MLGCP_point_pattern = function(params) {
  Y = simulate_grf(params$nx, params$ny, params$L, range = params$Y_scale, var = 1)
  U_list = vector("list", params$p)
  for (j in seq_len(params$p)) {
    U_list[[j]] = simulate_grf(params$nx, params$ny, params$L, range = params$U_scale[j], var = 1)
  }
  Lambda_list = vector("list", params$p)
  X_list = vector("list", params$p)
  for (i in seq_len(params$p)) {
    shared_term = params$Y_var * Y - (params$Y_var^2) / 2
    individual_term = params$U_var[i] * U_list[[i]] - (params$U_var[i]^2) / 2
    Lambda_list[[i]] = exp(params$beta0[i]) * exp(shared_term + individual_term)
    Lambda_im = im(mat = Lambda_list[[i]], xrange = c(0, params$L), yrange = c(0, params$L))
    X_list[[i]] = rpoispp(Lambda_im)
  }
  list(X_list = X_list, Lambda_list = Lambda_list, fields = list(Y = Y, U_list = U_list), params = params)
}

has_points = function(sim) {
  all(sapply(sim$X_list, function(x) x$n > 0))
}

fit_kppm_safe = function(Xi) {
  fit = tryCatch(
    kppm(unmark(Xi) ~ 1, clusters = "LGCP"),
    error = function(e) NULL
  )
  fit
}

make_counts = function(P, xbreaks, ybreaks, nxy) {
  cx = cut(P$x, xbreaks, include.lowest = TRUE, labels = FALSE)
  cy = cut(P$y, ybreaks, include.lowest = TRUE, labels = FALSE)
  idx = (cy - 1) * nxy + cx
  matrix(tabulate(idx, nbins = nxy * nxy), nrow = nxy, byrow = FALSE)
}

make_expected_counts = function(lambda_im, W, nxy) {
  df = as.data.frame(lambda_im)
  names(df)[1:3] = c("x", "y", "lambda")
  xbreaks = seq(W$xrange[1], W$xrange[2], length.out = nxy + 1)
  ybreaks = seq(W$yrange[1], W$yrange[2], length.out = nxy + 1)
  cx = cut(df$x, xbreaks, include.lowest = TRUE, labels = FALSE)
  cy = cut(df$y, ybreaks, include.lowest = TRUE, labels = FALSE)
  pixel_area = mean(diff(lambda_im$xcol)) * mean(diff(lambda_im$yrow))
  idx = (cy - 1) * nxy + cx
  E_sum = rowsum(df$lambda * pixel_area, idx, reorder = FALSE)
  out = numeric(nxy * nxy)
  out[as.integer(rownames(E_sum))] = E_sum[, 1]
  matrix(out, nrow = nxy, byrow = FALSE)
}

pearson_residual_image = function(N, E, eps = 1e-6) {
  (N - E) / sqrt(pmax(E, eps))
}

resize_count_to_base = function(count_mat, base_nxy) {
  nxy = nrow(count_mat)
  
  if (nxy == base_nxy)
    return(count_mat)
  
  if (base_nxy %% nxy != 0)
    stop("base_nxy must be divisible by nxy for multi-scale count images")
  
  fact = base_nxy / nxy
  kronecker(count_mat, matrix(1, nrow = fact, ncol = fact))
}

make_count_image = function(X, W, nxy_set, base_nxy = max(nxy_set)) {
  X1 = X[marks(X) == "Type1"]
  X2 = X[marks(X) == "Type2"]
  img_list = vector("list", length(nxy_set))
  for (k in seq_along(nxy_set)) {
    nxy = nxy_set[k]
    xbreaks = seq(W$xrange[1], W$xrange[2], length.out = nxy + 1)
    ybreaks = seq(W$yrange[1], W$yrange[2], length.out = nxy + 1)
    counts_type1 = make_counts(X1, xbreaks, ybreaks, nxy)
    counts_type2 = make_counts(X2, xbreaks, ybreaks, nxy)
    counts_all = make_counts(X, xbreaks, ybreaks, nxy)
    counts_type1 = resize_count_to_base(counts_type1, base_nxy)
    counts_type2 = resize_count_to_base(counts_type2, base_nxy)
    counts_all = resize_count_to_base(counts_all, base_nxy)
    img_list[[k]] = abind(counts_type1, counts_type2, counts_all, along = 3)
  }
  do.call(abind, c(img_list, along = 3))
}

make_residual_count_image = function(X, W, Lambda_est, nxy_set = c(50, 10), base_nxy = max(nxy_set), eps = 1e-6) {
  X1 = X[marks(X) == "Type1"]
  X2 = X[marks(X) == "Type2"]
  img_list = vector("list", length(nxy_set))
  for (k in seq_along(nxy_set)) {
    nxy = nxy_set[k]
    xbreaks = seq(W$xrange[1], W$xrange[2], length.out = nxy + 1)
    ybreaks = seq(W$yrange[1], W$yrange[2], length.out = nxy + 1)
    N1 = make_counts(X1, xbreaks, ybreaks, nxy)
    N2 = make_counts(X2, xbreaks, ybreaks, nxy)
    Nall = N1 + N2
    E1 = make_expected_counts(Lambda_est$Type1, W, nxy)
    E2 = make_expected_counts(Lambda_est$Type2, W, nxy)
    Eall = E1 + E2
    R1 = resize_count_to_base(pearson_residual_image(N1, E1, eps = eps), base_nxy)
    R2 = resize_count_to_base(pearson_residual_image(N2, E2, eps = eps), base_nxy)
    Rall = resize_count_to_base(pearson_residual_image(Nall, Eall, eps = eps), base_nxy)
    img_list[[k]] = abind(R1, R2, Rall, along = 3)
  }
  do.call(abind, c(img_list, along = 3))
}

build_param_list = function(idx, beta0_1, beta0_2, sigma_y, sigma_u_1, sigma_u_2, scale_y, scale_u_1, scale_u_2) {
  list(p = 2, L = 1.5, nx = 50, ny = 50, beta0 = c(beta0_1[idx], beta0_2[idx]), Y_var = sigma_y[idx], U_var = c(sigma_u_1[idx], sigma_u_2[idx]), Y_scale = scale_y[idx], U_scale = c(scale_u_1[idx], scale_u_2[idx]))
}

make_MC_data = function(sim, params, types = c("Type1", "Type2")) {
  W = owin(c(0, params$L), c(0, params$L))
  X = superimpose(Type1 = sim$X_list[[1]], Type2 = sim$X_list[[2]], W = W)
  marks(X) = factor(marks(X), levels = types)
  X
}

make_MC_init = function(par_ini_QBMC, beta0_init) {
  log(c(par_ini_QBMC[3], par_ini_QBMC[4], par_ini_QBMC[5], par_ini_QBMC[6], par_ini_QBMC[1], par_ini_QBMC[2], exp(beta0_init[1]), exp(beta0_init[2])))
}

MC_to_QBMC_est = function(par_MC) {
  c(exp(par_MC[5]), exp(par_MC[6]), exp(par_MC[1]), exp(par_MC[3]), exp(par_MC[2]), exp(par_MC[4]))
}

make_MC_reference_case = function(params, types = c("Type1", "Type2"), max_try = 100) {
  for (a in 1:max_try) {
    sim = simulate_MLGCP_point_pattern(params)
    if (!has_points(sim)) next
    X = make_MC_data(sim, params, types)
    return(list(X = X, sim = sim))
  }
  stop("Error")
}



fit_MC_Q_once = function(X, par0_MC, power, rmax, Corr, factor) {
  
  lower_MC = c(
    log(0.5),  # sigma_U1
    log(0.001),  # scale_U1
    log(0.5),  # sigma_U2
    log(0.001),  # scale_U2
    log(0.5),  # sigma_Y
    log(0.001),  # scale_Y
    5,        # beta0_1
    5         # beta0_2
  )
  
  upper_MC = c(
    log(2),      # sigma_U1
    log(0.3),    # scale_U1
    log(2),      # sigma_U2
    log(0.3),    # scale_U2
    log(2),      # sigma_Y
    log(0.3),    # scale_Y
    7,         # beta0_1
    7          # beta0_2
  )
  
  optim(par = par0_MC, fn = Q, factor = factor, data = X, power = power, rmax = rmax, Corr = Corr, method = "L-BFGS-B", lower = lower_MC, upper = upper_MC, hessian = TRUE, control = list(maxit = 5000))
}



extract_one_case_combined = function(params, r,par_ini_BMC = NULL, use_BMC = FALSE, image_nxy_set = NULL, par_ini_QBMC, QBMC_power, QBMC_rmax, types = c("Type1", "Type2"), Corr = "Pos") {
  W = owin(c(0, 1.5), c(0, 1.5))
  if (is.null(image_nxy_set))
    image_nxy_set = params$nx
  
  repeat {
    sim = simulate_MLGCP_point_pattern(params)
    if (has_points(sim)) break
  }
  X = superimpose(Type1 = sim$X_list[[1]], Type2 = sim$X_list[[2]], W = W)
  marks(X) = factor(marks(X), levels = types)
  beta0_init = params$beta0
  factor_Q = exp(2 * mean(beta0_init))
  par0_MC = make_MC_init(par_ini_QBMC, beta0_init)
  QBMC_fit = tryCatch(fit_MC_Q_once(X = X, par0_MC = par0_MC, power = QBMC_power, rmax = QBMC_rmax, Corr = Corr, factor = factor_Q),
                      error = function(e) NULL)
  QBMC_ok = !is.null(QBMC_fit)
  QBMC_est = NULL
  Est_MC = NULL
  Hes = NULL
  QBMC_value = NULL
  if (QBMC_ok) {
    QBMC_est = MC_to_QBMC_est(QBMC_fit$par)
    Est_MC = QBMC_fit$par
    Hes = eigen(QBMC_fit$hessian)$values
    QBMC_value = QBMC_fit$value
  }
  
  fit1 = fit_kppm_safe(X[marks(X) == "Type1"])
  fit2 = fit_kppm_safe(X[marks(X) == "Type2"])
  if (is.null(fit1) || is.null(fit2)) {
    return(list(ok = FALSE, X = X, ok_QBMC = QBMC_ok, QBMC_est = QBMC_est, Est_MC = Est_MC, Hes = Hes, value = QBMC_value, Np_type1 = sim$X_list[[1]]$n, Np_type2 = sim$X_list[[2]]$n, Np_all = sim$X_list[[1]]$n + sim$X_list[[2]]$n))
  }
  beta_kppm1 = coef(fit1)
  beta_kppm2 = coef(fit2)
  Lambda_est = list(Type1 = predict(fit1, type = "trend", dimyx = c(params$ny, params$nx)), Type2 = predict(fit2, type = "trend", dimyx = c(params$ny, params$nx)))
  
  pcf_cross_list = list()
  pcfc_list = list()
  for (a in 1:(length(types) - 1)) {
    for (b in (a + 1):length(types)) {
      i = types[a]
      j = types[b]
      name = paste0(i, "_", j)
      gc = tryCatch(
        Lcross.inhom(X, i = i, j = j, lambdaI = Lambda_est[[i]], lambdaJ = Lambda_est[[j]], r = r, correction = "border"),
        error = function(e) NULL
      )
      if (is.null(gc)) return(list(ok = FALSE, X = X, ok_QBMC = QBMC_ok, QBMC_est = QBMC_est, Est_MC = Est_MC, Hes = Hes, value = QBMC_value, Np_type1 = sim$X_list[[1]]$n, Np_type2 = sim$X_list[[2]]$n, Np_all = sim$X_list[[1]]$n + sim$X_list[[2]]$n))
      pcf_cross_list[[name]] = gc$border - gc$r
      
      if (use_BMC) {
        pcfc = tryCatch(
          pcfcross.inhom(X, i = i, j = j, lambdaI = Lambda_est[[i]], lambdaJ = Lambda_est[[j]], r = r, correction = "translate"),
          error = function(e) NULL
        )
        if (is.null(pcfc)) return(list(ok = FALSE, X = X, ok_QBMC = QBMC_ok, QBMC_est = QBMC_est, Est_MC = Est_MC, Hes = Hes, value = QBMC_value, Np_type1 = sim$X_list[[1]]$n, Np_type2 = sim$X_list[[2]]$n, Np_all = sim$X_list[[1]]$n + sim$X_list[[2]]$n))
        pcfc_list[[name]] = pcfc$trans
      }
    }
  }
  
  pcf_within_list = list()
  pcf_list = list()
  for (i in types) {
    Xi = X[marks(X) == i]
    gi = tryCatch(
      Linhom(Xi, lambda = Lambda_est[[i]], r = r, correction = "border"),
      error = function(e) NULL
    )
    if (is.null(gi)) return(list(ok = FALSE, X = X, ok_QBMC = QBMC_ok, QBMC_est = QBMC_est, Est_MC = Est_MC, Hes = Hes, value = QBMC_value, Np_type1 = sim$X_list[[1]]$n, Np_type2 = sim$X_list[[2]]$n, Np_all = sim$X_list[[1]]$n + sim$X_list[[2]]$n))
    pcf_within_list[[i]] = gi$border - gi$r
    
    if (use_BMC) {
      pcfi = tryCatch(
        pcfinhom(Xi, lambda = Lambda_est[[i]], r = r, correction = "translate"),
        error = function(e) NULL
      )
      if (is.null(pcfi)) return(list(ok = FALSE, X = X, ok_QBMC = QBMC_ok, QBMC_est = QBMC_est, Est_MC = Est_MC, Hes = Hes, value = QBMC_value, Np_type1 = sim$X_list[[1]]$n, Np_type2 = sim$X_list[[2]]$n, Np_all = sim$X_list[[1]]$n + sim$X_list[[2]]$n))
      pcf_list[[i]] = pcfi$trans
    }
  }
  
  count_img = make_residual_count_image(X, W, Lambda_est, image_nxy_set, base_nxy = max(image_nxy_set))
  beta_vec = c(beta_kppm1["(Intercept)"], beta_kppm2["(Intercept)"])
  
  BMC_est = NULL
  if (use_BMC) {
    BMC_est = tryCatch(
      est_BMC_pcf(r, par_ini_BMC, pcf_list$Type1, pcf_list$Type2, pcfc_list$Type1_Type2, weight1 = rep(1, length(r)), weight2 = rep(1, length(r)), weight3 = rep(1, length(r)), q = 0.25),
      error = function(e) NULL
    )
  }
  
  list(
    ok = TRUE,
    X = X,
    pcf = pcf_within_list,
    pcfcross = pcf_cross_list,
    beta = beta_vec,
    BMC_est = BMC_est,
    count_img = count_img,
    Np_type1 = sim$X_list[[1]]$n,
    Np_type2 = sim$X_list[[2]]$n,
    Np_all = sim$X_list[[1]]$n + sim$X_list[[2]]$n,
    ok_QBMC = QBMC_ok,
    QBMC_est = QBMC_est,
    Est_MC = Est_MC,
    Hes = Hes,
    value = QBMC_value
  )
}


stack_feature_results = function(res_list, response_par = NULL) {
  keep_true = sapply(res_list, function(x) isTRUE(x$ok))
  res_ok = res_list[keep_true]
  X_ok = lapply(res_ok, function(x) x$X)
  pcf_within = lapply(res_ok, function(x) x$pcf)
  pcf_cross = lapply(res_ok, function(x) x$pcfcross)
  pcf_list1 = do.call(cbind, lapply(pcf_within, function(u) u[["Type1"]]))
  pcf_list2 = do.call(cbind, lapply(pcf_within, function(u) u[["Type2"]]))
  pcf_cross_list12 = do.call(cbind, lapply(pcf_cross, function(u) u[["Type1_Type2"]]))
  keep_pcf = (colSums(is.na(pcf_list1)) == 0) & (colSums(is.na(pcf_list2)) == 0) & (colSums(is.na(pcf_cross_list12)) == 0)
  pcf_list1 = pcf_list1[, keep_pcf, drop = FALSE]
  pcf_list2 = pcf_list2[, keep_pcf, drop = FALSE]
  pcf_cross_list12 = pcf_cross_list12[, keep_pcf, drop = FALSE]
  pcf1 = array(as.numeric(t(pcf_list1)), dim = c(ncol(pcf_list1), nrow(pcf_list1), 1))
  pcf2 = array(as.numeric(t(pcf_list2)), dim = c(ncol(pcf_list2), nrow(pcf_list2), 1))
  pcf12 = array(as.numeric(t(pcf_cross_list12)), dim = c(ncol(pcf_cross_list12), nrow(pcf_cross_list12), 1))
  pcf_Y = abind(pcf12, along = 3)
  d1_12 = pcf1 - pcf12
  d2_21 = pcf2 - pcf12
  U_pcf = abind(d1_12, d2_21, along = 3)
  pcf_input = abind(pcf_Y, U_pcf, along = 3)
  
  count_list = lapply(res_ok, function(x) x$count_img)
  count_keep = count_list[keep_pcf]
  count_img = aperm(do.call(abind, c(count_keep, along = 4)), c(4, 1, 2, 3))
  beta_mat = do.call(rbind, lapply(res_ok, function(x) x$beta))
  beta_mat = beta_mat[keep_pcf, , drop = FALSE]
  Np_type1 = as.numeric(lapply(res_ok, function(x) x$Np_type1))[keep_pcf]
  Np_type2 = as.numeric(lapply(res_ok, function(x) x$Np_type2))[keep_pcf]
  Np_all = as.numeric(lapply(res_ok, function(x) x$Np_all))[keep_pcf]
  Np = cbind(Np_all, Np_type1, Np_type2)
  aux_input = beta_mat
  out = list(X = X_ok[keep_pcf], keep_true = keep_true, keep_pcf = keep_pcf, pcf_input = pcf_input, count_img = count_img, aux_input = aux_input, beta_input = beta_mat, Np = Np)
  #out$X = X_ok[keep_pcf]
  
  if (!is.null(response_par)) {
    response_par_ok = response_par[keep_true, , drop = FALSE]
    out$Y = sqrt(response_par_ok[keep_pcf, , drop = FALSE])
    out$true_par = response_par_ok[keep_pcf, , drop = FALSE]
  }
  BMC_list = lapply(res_ok, function(x) x$BMC_est)
  has_BMC = sapply(BMC_list, function(x) !is.null(x))
  if (length(BMC_list) > 0 && any(has_BMC)) {
    BMC_keep = BMC_list[keep_pcf]
    has_BMC_keep = sapply(BMC_keep, function(x) !is.null(x))
    pred_par_BMC = matrix(NA, nrow = length(BMC_keep), ncol = 6)
    colnames(pred_par_BMC) = c("sigma_Y", "scale_Y", "sigma_U1", "sigma_U2", "scale_U1", "scale_U2")
    if (any(has_BMC_keep)) {
      pred_par_BMC[has_BMC_keep, ] = do.call(rbind, lapply(BMC_keep[has_BMC_keep], function(x) {
        c(x[1], x[2], x[3], x[5], x[4], x[6])
      }))
    }
    out$pred_par_BMC = pred_par_BMC
    out$keep_BMC = has_BMC_keep
  }
  out
}

scale_by_train = function(train_x, test_x) {
  scale_mean = colMeans(train_x, na.rm = TRUE)
  scale_sd = apply(train_x, 2, sd, na.rm = TRUE)
  scale_sd[scale_sd == 0] = 1
  train_scaled = sweep(train_x, 2, scale_mean, "-")
  train_scaled = sweep(train_scaled, 2, scale_sd, "/")
  test_scaled = sweep(test_x, 2, scale_mean, "-")
  test_scaled = sweep(test_scaled, 2, scale_sd, "/")
  list(train = train_scaled, test = test_scaled, mean = scale_mean, sd = scale_sd)
}

scale_pcf_by_train = function(train_x, test_x) {
  n_pos = dim(train_x)[2]
  n_channels = dim(train_x)[3]
  train_scaled = train_x
  test_scaled = test_x
  scale_mean = matrix(NA, nrow = n_pos, ncol = n_channels)
  scale_sd = matrix(NA, nrow = n_pos, ncol = n_channels)
  for (ch in seq_len(n_channels)) {
    scale_mean[, ch] = colMeans(train_x[, , ch], na.rm = TRUE)
    scale_sd[, ch] = apply(train_x[, , ch], 2, sd, na.rm = TRUE)
    scale_sd[!is.finite(scale_sd[, ch]), ch] = 1
    scale_sd[scale_sd[, ch] == 0, ch] = 1
    train_scaled[, , ch] = sweep(train_x[, , ch], 2, scale_mean[, ch], "-")
    train_scaled[, , ch] = sweep(train_scaled[, , ch], 2, scale_sd[, ch], "/")
    test_scaled[, , ch] = sweep(test_x[, , ch], 2, scale_mean[, ch], "-")
    test_scaled[, , ch] = sweep(test_scaled[, , ch], 2, scale_sd[, ch], "/")
  }
  list(train = train_scaled, test = test_scaled, mean = scale_mean, sd = scale_sd)
}

compare_par = function(true_par, pred_par, par_names = NULL) {
  true_par = as.matrix(true_par)
  pred_par = as.matrix(pred_par)
  true_mean = colMeans(true_par, na.rm = TRUE)
  est_mean  = colMeans(pred_par, na.rm = TRUE)
  bias = est_mean - true_mean
  data.frame(
    parameter = par_names,
    true      = true_mean,
    estimate  = est_mean,
    bias      = bias,
    rmse      = sqrt(colMeans((pred_par - true_par)^2, na.rm = TRUE)),
    mape      = colMeans(abs((pred_par - true_par) / pmax(true_par, 1e-8)), na.rm = TRUE) * 100
  )
}

# ------------------------------------------------------------------------------
set.seed(2026)
W = owin(c(0, 1.5), c(0, 1.5))
rmax = rmax.rule("K", W)
r = seq(0, rmax, length.out = 513)
par_ini_BMC = c(sigma2_y = 1, xi_y = 0.1, sigma2_u1 = 1, xi_u1 = 0.1, sigma2_u2 = 1, xi_u2 = 0.1)
par_ini_QBMC = c(1.5, 0.2, 1, 0.15, 1, 0.15)
types = c("Type1", "Type2")
image_nxy_set = 50
# ------------------------------------------------------------------------------
# Test
ntest = 500
set.seed(2026)
beta0_1_test = runif(ntest, 6, 6)
beta0_2_test = runif(ntest, 6, 6)
sigma_y_test = runif(ntest, 1.5, 1.5)
scale_y_test = runif(ntest, 0.2, 0.2)
sigma_u_1_test = runif(ntest, 1, 1)
sigma_u_2_test = runif(ntest, 1, 1)
scale_u_1_test = runif(ntest, 0.15, 0.15)
scale_u_2_test = runif(ntest, 0.15, 0.15)
test_par_full = cbind(sigma_y_test, scale_y_test, sigma_u_1_test, sigma_u_2_test, scale_u_1_test, scale_u_2_test)

# ------------------------------------------------------------------------------
# QMC Monte Carlo tuning 
QBMC_MC_Nsim = 300
QBMC_power_grid = c(0.1, 0.2, 0.3, 0.4, 0.5)
QBMC_rmax_grid = seq(0.15, 0.35, by = 0.025)
tun_par = data.frame(power = rep(QBMC_power_grid, each = length(QBMC_rmax_grid)), rmax = rep(QBMC_rmax_grid, length(QBMC_power_grid)))
Corr = "Pos"
Corr_sign = 1

params_QBMC_tune = list(
  p = 2, L = 1.5, nx = 50, ny = 50,
  beta0 = c(mean(range(beta0_1_test)), mean(range(beta0_2_test))),
  Y_var = mean(range(sigma_y_test)),
  U_var = c(mean(range(sigma_u_1_test)), mean(range(sigma_u_2_test))),
  Y_scale = mean(range(scale_y_test)),
  U_scale = c(mean(range(scale_u_1_test)), mean(range(scale_u_2_test)))
)

plan(sequential)
gc()
plan(multisession, workers = max(1, parallelly::availableCores() - 1))
handlers(global = TRUE)
handlers("txtprogressbar")
QBMC_tune_out = with_progress({
  p_QBMC = progressor(steps = nrow(tun_par))
  ref_QBMC_tune = make_MC_reference_case(params_QBMC_tune, types = types)
  Data_ref = ref_QBMC_tune$X
  bh1 = log(subset(Data_ref, marks == "Type1")$n / area.owin(Data_ref$window))
  bh2 = log(subset(Data_ref, marks == "Type2")$n / area.owin(Data_ref$window))
  factor_Q = exp(2 * mean(params_QBMC_tune$beta0))
  WL = params_QBMC_tune$L
  simu_est = vector("list", nrow(tun_par))
  T1 = Sys.time()
  for (i in seq_len(nrow(tun_par))) {
    p_QBMC(sprintf("tuning %d/%d", i, nrow(tun_par)))
    power = tun_par$power[i]
    rmax = tun_par$rmax[i]
    init = make_MC_init(par_ini_QBMC, params_QBMC_tune$beta0)
    t1 = Sys.time()
    m2 = tryCatch(fit_MC_Q_once(X = Data_ref, par0_MC = init, power = power, rmax = rmax, Corr = Corr, factor = factor_Q),
                  error = function(e) NULL)
    t2 = Sys.time()
    if (is.null(m2)) {
      simu_est[[i]] = list(Est = rep(NA, 8), init = init, power = power, rmax = rmax, Hes = rep(NA, 8), det = NA, sd = rep(NA, 9), est_time = t2 - t1, est_units = units(t2 - t1), cal_cov_time = NA, cal_cov_units = NA)
      next
    }
    Btheta = cal_Btheta_useQ(m2$par, factor = factor_Q, power, rmax, Corr, correction = "isotropic", Inhom = FALSE)
    Par_est = c(exp(m2$par[1:6]), m2$par[7:8])
    result = tryCatch({
      t3 = Sys.time()
      cal_cov = vector("list", QBMC_MC_Nsim)
      for (j in seq_len(QBMC_MC_Nsim)) {
        set.seed(j)
        params_cov = params_QBMC_tune
        params_cov$beta0 = c(bh1, bh2)
        params_cov$Y_var = Par_est[5]
        params_cov$Y_scale = Par_est[6]
        params_cov$U_var = c(Par_est[1], Par_est[3])
        params_cov$U_scale = c(Par_est[2], Par_est[4])
        sim_cov = simulate_MLGCP_point_pattern(params_cov)
        Data_cov = make_MC_data(sim_cov, params_cov, types)
        Stheta = cal_Stheta_useQ(par = m2$par, factor = factor_Q, data = Data_cov, power, rmax, Corr, correction = "isotropic", Inhom = FALSE)
        cal_cov[[j]] = list(Cb11 = Stheta$Cb11, Cb22 = Stheta$Cb22, Cb12 = Stheta$Cb12, Cb21 = Stheta$Cb21)
      }
      Output = list(Cb11 = sapply(cal_cov, function(x) x$Cb11), Cb22 = sapply(cal_cov, function(x) x$Cb22), Cb12 = sapply(cal_cov, function(x) x$Cb12), Cb21 = sapply(cal_cov, function(x) x$Cb21))
      Stheta_matrix = cov(t(Output$Cb11)) + cov(t(Output$Cb22)) + 2 * cov(t(Output$Cb12)) +
        cov(t(Output$Cb11), t(Output$Cb12)) + cov(t(Output$Cb11), t(Output$Cb21)) + cov(t(Output$Cb11), t(Output$Cb22)) +
        cov(t(Output$Cb12), t(Output$Cb11)) + cov(t(Output$Cb12), t(Output$Cb21)) + cov(t(Output$Cb12), t(Output$Cb22)) +
        cov(t(Output$Cb21), t(Output$Cb11)) + cov(t(Output$Cb21), t(Output$Cb12)) + cov(t(Output$Cb21), t(Output$Cb22)) +
        cov(t(Output$Cb22), t(Output$Cb11)) + cov(t(Output$Cb22), t(Output$Cb12)) + cov(t(Output$Cb22), t(Output$Cb21))
      Cov_matrix = solve(Btheta) %*% Stheta_matrix %*% solve(Btheta) / (WL^2)
      t4 = Sys.time()
      Sgm1 = Par_est[1]
      Sgm2 = Par_est[3]
      Sgm3 = Par_est[5]
      temp = (Sgm1^2 + Sgm3^2) * (Sgm2^2 + Sgm3^2)
      Corr_deriv_theta = c(-Sgm1 * Sgm3^2 * temp^(-1.5) * (Sgm2^2 + Sgm3^2), 0, -Sgm2 * Sgm3^2 * temp^(-1.5) * (Sgm1^2 + Sgm3^2), 0, 2 * Sgm3 * temp^(-0.5) - Sgm3^3 * temp^(-1.5) * (Sgm1^2 + Sgm2^2 + 2 * Sgm3^2), 0, 0, 0)
      var_Corr = t(Corr_deriv_theta) %*% Cov_matrix %*% Corr_deriv_theta
      det = det(Cov_matrix[1:6, 1:6])
      sd = sqrt(c(diag(Cov_matrix), var_Corr))
      list(Est = m2$par, init = init, power = power, rmax = rmax, Hes = eigen(m2$hessian)$values, det = det, sd = sd, est_time = t2 - t1, est_units = units(t2 - t1), cal_cov_time = t4 - t3, cal_cov_units = units(t4 - t3))
    }, error = function(e) {
      list(Est = m2$par, init = init, power = power, rmax = rmax, Hes = eigen(m2$hessian)$values, det = NA, sd = rep(NA, 9), est_time = t2 - t1, est_units = units(t2 - t1), cal_cov_time = NA, cal_cov_units = NA)
    })
    simu_est[[i]] = result
  }
  T2 = Sys.time()
  Output2_pos = list(Corr_sign = Corr_sign, T_total = T2 - T1,
                     power = sapply(simu_est, function(x) x$power),
                     rmax = sapply(simu_est, function(x) x$rmax),
                     Est = sapply(simu_est, function(x) x$Est),
                     sd = sapply(simu_est, function(x) x$sd),
                     Hes = sapply(simu_est, function(x) x$Hes),
                     init = sapply(simu_est, function(x) x$init),
                     det = sapply(simu_est, function(x) x$det),
                     est_units = sapply(simu_est, function(x) x$est_units),
                     est_time = sapply(simu_est, function(x) x$est_time),
                     cal_cov_units = sapply(simu_est, function(x) x$cal_cov_units),
                     cal_cov_time = sapply(simu_est, function(x) x$cal_cov_time))
  ok_det = which(is.finite(Output2_pos$det))
  best_id = ok_det[which.min(Output2_pos$det[ok_det])]
  Output2_pos$best = data.frame(power = Output2_pos$power[best_id], rmax = Output2_pos$rmax[best_id], det = Output2_pos$det[best_id])
  Output2_pos
})
QBMC_power_opt = QBMC_tune_out$best$power[1]
QBMC_rmax_opt = QBMC_tune_out$best$rmax[1]
print(QBMC_tune_out$best)

#power rmax          det
#1   0.5  0.2 8.937144e-16

# ------------------------------------------------------------------------------
plan(sequential)
gc()
plan(multisession, workers = max(1, parallelly::availableCores() - 1))
handlers(global = TRUE)
handlers("txtprogressbar")
with_progress({
  p = progressor(along = seq_len(ntest))
  res_test = future_lapply(
    seq_len(ntest),
    function(j) {
      p(sprintf("combined iteration %d/%d", j, ntest))
      params = build_param_list(
        j,
        beta0_1_test, beta0_2_test,
        sigma_y_test, sigma_u_1_test, sigma_u_2_test,
        scale_y_test, scale_u_1_test, scale_u_2_test
      )
      extract_one_case_combined(params, r, par_ini_BMC = par_ini_BMC, use_BMC = TRUE, image_nxy_set = image_nxy_set, par_ini_QBMC = par_ini_QBMC, QBMC_power = QBMC_power_opt, QBMC_rmax = QBMC_rmax_opt, types = types, Corr = Corr)
    },
    future.seed = 2025,
    future.packages = c("spatstat.geom", "spatstat.explore", "spatstat.random", "spatstat.model", "terra", "abind", "future.apply")
  )
})
test_data = stack_feature_results(res_test, response_par = test_par_full)
pcf_test = test_data$pcf_input
img_test = test_data$count_img
aux_test = test_data$aux_input
test_par = test_data$true_par
pred_par_BMC = test_data$pred_par_BMC
keep_BMC = test_data$keep_BMC
scaled_aux = scale_by_train(aux_train, aux_test)
aux_train_scaled = scaled_aux$train
aux_test_scaled = scaled_aux$test

scaled_pcf = scale_pcf_by_train(pcf_train, pcf_test)
pcf_train_scaled = scaled_pcf$train
pcf_test_scaled = scaled_pcf$test




###

compare_par = function(true_par, pred_par, par_names = NULL) {
  true_par = as.matrix(true_par)
  pred_par = as.matrix(pred_par)
  true_mean = colMeans(true_par, na.rm = TRUE)
  est_mean  = colMeans(pred_par, na.rm = TRUE)
  bias = est_mean - true_mean
  data.frame(
    parameter = par_names,
    true      = true_mean,
    estimate  = est_mean,
    bias      = bias,
    rmse      = sqrt(colMeans((pred_par - true_par)^2, na.rm = TRUE)),
    mape      = colMeans(abs((pred_par - true_par) / pmax(true_par, 1e-8)), na.rm = TRUE) * 100
  )
}

#QMC
keep_QBMC = sapply(res_test, function(x) isTRUE(x$ok_QBMC))
pred_par_QBMC = matrix(NA, nrow = ntest, ncol = 6)
colnames(pred_par_QBMC) = c("sigma_Y", "scale_Y", "sigma_U1", "sigma_U2", "scale_U1", "scale_U2")
pred_par_QBMC[keep_QBMC, ] = do.call(rbind, lapply(res_test[keep_QBMC], function(x) x$QBMC_est))
idx = apply(pred_par_QBMC[keep_QBMC, , drop = FALSE] < 10e20, 1, all)
compare_QBMC = compare_par(test_par_full[keep_QBMC, , drop = FALSE][idx, , drop = FALSE], pred_par_QBMC[keep_QBMC, , drop = FALSE][idx, , drop = FALSE], c("sigma_Y", "scale_Y", "sigma_U1", "sigma_U2", "scale_U1", "scale_U2"))

compare_QBMC[] <- lapply(compare_QBMC, function(x) if(is.numeric(x)) sprintf("%.3f", x) else x)
print(compare_QBMC, row.names = FALSE)


# ------------------------------------------------------------------------------
# DSBI 
source_python("BLGCP_NN.py")
img_train_50 = img_train[, , , 1:3, drop = FALSE]
img_test_50  = img_test[, , , 1:3, drop = FALSE]
n_runs = 10
res_all_50 = vector("list", n_runs)
predict_batch_size = 2000
for (i in 1:n_runs) {
  seed_j = 2027 + i
  set.seed(seed_j)
  seed_i = sample.int(10000, 1)
  set_all_seeds(seed_i)
  model_path_i = sprintf("C:/Users/qihan/Desktop/hom_W2/W2_sim1_model_50_run_%d.pth", i)
  model_hybrid = train_and_save_hybrid_NN_model(pcf_train_scaled, img_train_50, Y_train, aux_train_scaled, model_path = model_path_i, batch_size = 100, epochs = 30, lr = 1e-3, requested_device = "auto", seed = seed_i, verbose = TRUE, standardize_img = TRUE)
  pred_par = load_hybrid_NN_model_predict(pcf_test_scaled, img_test_50, aux_test_scaled, model_path_i, requested_device = "auto", predict_batch_size = predict_batch_size)^2
  res_all_50[[i]] = pred_par
  rm(model_hybrid, pred_par)
  clear_all()
  gc()
}
pred_par_50 = get_mean_sd(res_all_50)$mean

compare_hybrid_50 = compare_par(test_par, pred_par_50, c("sigma_Y", "scale_Y", "sigma_U1", "sigma_U2", "scale_U1", "scale_U2"))
compare_hybrid_50[] = lapply(compare_hybrid_50, function(x) if (is.numeric(x)) sprintf("%.3f", x) else x)
print(compare_hybrid_50, row.names = FALSE)




save.image("C:/Users/qihan/Desktop/hom_W2/sim1_W2.RData")





# INLA
param_names = c("sigma_Y", "scale_Y", "sigma_U1", "sigma_U2", "scale_U1", "scale_U2")
X_inla_list = test_data$X
true_par_INLA = test_data$true_par
colnames(true_par_INLA) = param_names



L = 1.5
boundary_poly = matrix(c(0, 0, L, 0, L, L, 0, L, 0, 0), ncol = 2, byrow = TRUE)
boundary_INLA = sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(boundary_poly))), crs = NA)
mesh_INLA = fmesher::fm_mesh_2d(boundary = boundary_INLA, max.edge = c(0.04, 0.12), cutoff = 0.005)


sigma_y0  = as.numeric(quantile(sigma_y_train, 0.90))
sigma_u10 = as.numeric(quantile(sigma_u_1_train, 0.90))
sigma_u20 = as.numeric(quantile(sigma_u_2_train, 0.90))

scale_y0  = as.numeric(quantile(scale_y_train, 0.90))
scale_u10 = as.numeric(quantile(scale_u_1_train, 0.90))
scale_u20 = as.numeric(quantile(scale_u_2_train, 0.90))

prior_INLA = list(
  shared = list(
    range = c(2 * scale_y0, 0.90),   # P(scale_Y < scale_y0) = 0.90
    sigma = c(sigma_y0, 0.10)        # P(sigma_Y > sigma_y0) = 0.10
  ),
  u1 = list(
    range = c(2 * scale_u10, 0.90),  # P(scale_U1 < scale_u10) = 0.90
    sigma = c(sigma_u10, 0.10)
  ),
  u2 = list(
    range = c(2 * scale_u20, 0.90),  # P(scale_U2 < scale_u20) = 0.90
    sigma = c(sigma_u20, 0.10)
  )
)


print(prior_INLA)




plan(sequential)
gc()
plan(multisession, workers = max(1, 1))
handlers(global = TRUE)
handlers("txtprogressbar")

INLA_res = with_progress({
  p_INLA = progressor(along = seq_along(X_inla_list))
  future_lapply(
    seq_along(X_inla_list),
    function(j) {
      p_INLA(sprintf("INLA iteration %d/%d", j, length(X_inla_list)))
      
      tryCatch({
        
        Xj = X_inla_list[[j]]
        X_type1 = Xj[marks(Xj) == "Type1"]
        X_type2 = Xj[marks(Xj) == "Type2"]
        dat_type1 = sf::st_as_sf(data.frame(x = X_type1$x, y = X_type1$y), coords = c("x", "y"), crs = NA)
        dat_type2 = sf::st_as_sf(data.frame(x = X_type2$x, y = X_type2$y), coords = c("x", "y"), crs = NA)
        
        spde_shared = INLA::inla.spde2.pcmatern(
          mesh = mesh_INLA,
          alpha = 1.5,
          prior.range = prior_INLA$shared$range,
          prior.sigma = prior_INLA$shared$sigma
        )
        
        spde_u1 = INLA::inla.spde2.pcmatern(
          mesh = mesh_INLA,
          alpha = 1.5,
          prior.range = prior_INLA$u1$range,
          prior.sigma = prior_INLA$u1$sigma
        )
        
        spde_u2 = INLA::inla.spde2.pcmatern(
          mesh = mesh_INLA,
          alpha = 1.5,
          prior.range = prior_INLA$u2$range,
          prior.sigma = prior_INLA$u2$sigma
        )
        
        cmp = ~
          Shared(geometry, model = spde_shared) +
          Type1_field(geometry, model = spde_u1) +
          Type2_field(geometry, model = spde_u2) +
          Intercept_type1(1) +
          Intercept_type2(1)
        
        fml_type1 = geometry ~ Intercept_type1 + Shared + Type1_field
        fml_type2 = geometry ~ Intercept_type2 + Shared + Type2_field
        
        lik_type1 = inlabru::bru_obs(
          "cp",
          formula = fml_type1,
          data = dat_type1,
          samplers = boundary_INLA,
          domain = list(geometry = mesh_INLA),
          tag = "type1"
        )
        
        lik_type2 = inlabru::bru_obs(
          "cp",
          formula = fml_type2,
          data = dat_type2,
          samplers = boundary_INLA,
          domain = list(geometry = mesh_INLA),
          tag = "type2"
        )
        
        fit = inlabru::bru(
          cmp,
          lik_type1,
          lik_type2,
          options = list(
            control.inla = list(
              int.strategy = "eb"
            ),
            bru_max_iter = 1
          )
        )
        
        fit = tryCatch(
          inlabru::bru_rerun(fit),
          error = function(e) fit
        )
        
        hyper = fit$summary.hyperpar
        
        est = c(
          sigma_Y  = hyper$mean[2],
          scale_Y  = hyper$mean[1] / 2,
          sigma_U1 = hyper$mean[4],
          sigma_U2 = hyper$mean[6],
          scale_U1 = hyper$mean[3] / 2,
          scale_U2 = hyper$mean[5] / 2
        )
        
        list(
          ok = all(is.finite(est)),
          est = est,
          fit = fit,
          hyper = hyper,
          error = NULL
        )
        
      }, error = function(e) {
        list(
          ok = FALSE,
          est = rep(NA, 6),
          fit = NULL,
          error = conditionMessage(e)
        )
      })
    },
    future.seed = 2025,
    future.packages = c("sf", "spatstat.geom", "INLA", "inlabru")
  )
})



keep_INLA = sapply(INLA_res, function(x) isTRUE(x$ok))
pred_par_INLA = matrix(NA, nrow = ntest, ncol = 6)
colnames(pred_par_INLA) = c("sigma_Y", "scale_Y", "sigma_U1", "sigma_U2", "scale_U1", "scale_U2")
pred_par_INLA[keep_INLA, ] = do.call(rbind, lapply(INLA_res[keep_INLA], function(x) x$est))
idx_INLA = apply(pred_par_INLA[keep_INLA, , drop = FALSE] < 10, 1, all)
compare_INLA = compare_par(test_par_full[keep_INLA, , drop = FALSE][idx_INLA, , drop = FALSE], pred_par_INLA[keep_INLA, , drop = FALSE][idx_INLA, , drop = FALSE], c("sigma_Y", "scale_Y", "sigma_U1", "sigma_U2", "scale_U1", "scale_U2"))
compare_INLA[] = lapply(compare_INLA, function(x) if (is.numeric(x)) sprintf("%.3f", x) else x)
print(compare_INLA, row.names = FALSE)






useful_output <- list(
  settings = list(
    ntest = ntest,
    types = types,
    image_nxy_set = image_nxy_set,
    QBMC_power_opt = QBMC_power_opt,
    QBMC_rmax_opt = QBMC_rmax_opt
  ),
  
  tuning = list(
    QBMC_best = QBMC_tune_out$best,
    QBMC_grid = data.frame(
      power = QBMC_tune_out$power,
      rmax  = QBMC_tune_out$rmax,
      det   = QBMC_tune_out$det
    )
  ),
  
  true_parameters = list(
    test_par_full = test_par_full,
    test_par = test_par
  ),
  
  predictions = list(
    DSBI_50 = pred_par_50,
    QBMC    = pred_par_QBMC,
    INLA    = pred_par_INLA
  ),
  
  valid_indices = list(
    keep_QBMC = keep_QBMC,
    idx_QBMC  = idx,
    keep_INLA = keep_INLA,
    idx_INLA  = idx_INLA
  ),
  
  comparison_tables = list(
    DSBI_50 = compare_hybrid_50,
    QBMC    = compare_QBMC,
    INLA    = compare_INLA
  )
)

save(
  useful_output,
  file = "C:/Users/qihan/Desktop/hom_W2/sim1_useful_output_W2.RData"
)



load("C:/Users/qihan/Desktop/hom_W2/sim1_useful_output_W2.RData")

# ------------------------------------------------------------
# Extract results
# ------------------------------------------------------------

pred <- useful_output$predictions
truth <- useful_output$true_parameters
valid <- useful_output$valid_indices

par_names <- c(
  "sigma_Y", "scale_Y", "sigma_U1", "sigma_U2", "scale_U1", "scale_U2"
)

add_par_names <- function(x) {
  colnames(x) <- par_names
  x
}

pred_par_50   <- add_par_names(pred$DSBI_50)
pred_par_QBMC <- add_par_names(pred$QBMC)
pred_par_INLA <- add_par_names(pred$INLA)

test_par_full <- add_par_names(truth$test_par_full)
test_par      <- add_par_names(truth$test_par)

keep_QBMC <- valid$keep_QBMC
idx_QBMC  <- valid$idx_QBMC

keep_INLA <- valid$keep_INLA


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

print_compare <- function(compare_obj) {
  compare_print <- compare_obj
  
  compare_print[] <- lapply(
    compare_print,
    function(x) if (is.numeric(x)) sprintf("%.3f", x) else x
  )
  
  print(compare_print, row.names = FALSE)
}

compare_methods <- function(true, pred) {
  compare_par(true, pred, par_names)
}


# ------------------------------------------------------------
# Parameter bounds for INLA post-processing
# ------------------------------------------------------------

lower_train <- c(
  sigma_Y  = 0.5,
  scale_Y  = 0.001,
  sigma_U1 = 0.5,
  sigma_U2 = 0.5,
  scale_U1 = 0.001,
  scale_U2 = 0.001
)

upper_train <- c(
  sigma_Y  = 2,
  scale_Y  = 0.3,
  sigma_U1 = 2,
  sigma_U2 = 2,
  scale_U1 = 0.3,
  scale_U2 = 0.3
)

in_bounds <- function(x, par_k) {
  is.finite(x) &
    x >= lower_train[par_k] &
    x <= upper_train[par_k]
}


# ------------------------------------------------------------
# INLA comparison
# Keep each parameter only when its INLA estimate is valid
# ------------------------------------------------------------

pred_par_INLA_pp <- matrix(
  NA,
  nrow = nrow(pred_par_INLA),
  ncol = length(par_names),
  dimnames = list(NULL, par_names)
)

true_par_INLA_pp <- matrix(
  NA,
  nrow = nrow(test_par_full),
  ncol = length(par_names),
  dimnames = list(NULL, par_names)
)

for (par_k in par_names) {
  idx_k <- keep_INLA & in_bounds(pred_par_INLA[, par_k], par_k)
  
  pred_par_INLA_pp[idx_k, par_k] <- pred_par_INLA[idx_k, par_k]
  true_par_INLA_pp[idx_k, par_k] <- test_par_full[idx_k, par_k]
}

compare_INLA <- compare_methods(true_par_INLA_pp, pred_par_INLA_pp)
print_compare(compare_INLA)


# ------------------------------------------------------------
# DSBI comparison
# ------------------------------------------------------------

compare_hybrid_50 <- compare_methods(test_par, pred_par_50)
print_compare(compare_hybrid_50)


# ------------------------------------------------------------
# QMC comparison
# ------------------------------------------------------------

pred_QBMC_valid <- pred_par_QBMC[keep_QBMC, , drop = FALSE][idx_QBMC, , drop = FALSE]
true_QBMC_valid <- test_par_full[keep_QBMC, , drop = FALSE][idx_QBMC, , drop = FALSE]

compare_QBMC <- compare_methods(true_QBMC_valid, pred_QBMC_valid)
print_compare(compare_QBMC)







plot_order <- c(
  "sigma_Y", "sigma_U1", "sigma_U2",
  "scale_Y", "scale_U1", "scale_U2"
)

ylim_list <- list(
  sigma_Y  = c(0.5, 2),
  sigma_U1 = c(0.5, 2),
  sigma_U2 = c(0.5, 2),
  scale_Y  = c(0, 0.3),
  scale_U1 = c(0, 0.3),
  scale_U2 = c(0, 0.3)
)

true_values <- c(
  sigma_Y  = 1.5,
  sigma_U1 = 1,
  sigma_U2 = 1,
  scale_Y  = 0.2,
  scale_U1 = 0.15,
  scale_U2 = 0.15
)

png(
  filename = "sim1_boxplot_param2_W2_final.png",
  width = 2800,
  height = 1800,
  res = 300
)

par(mfrow = c(2, 3))

for (par_k in plot_order) {
  idx_INLA_k <- keep_INLA & in_bounds(pred_par_INLA[, par_k], par_k)
  
  boxplot(
    pred_par_50[, par_k],
    pred_QBMC_valid[, par_k],
    pred_par_INLA[idx_INLA_k, par_k],
    names = c("DSBI", "MC", "INLA"),
    main = par_k,
    ylim = ylim_list[[par_k]]
  )
  
  abline(
    h = true_values[par_k],
    lty = 2,
    lwd = 1.5,
    col = "red"
  )
}

dev.off()











outside_INLA <- sapply(par_names, function(par_k) {
  keep_INLA &
    is.finite(pred_par_INLA[, par_k]) &
    (
      pred_par_INLA[, par_k] < lower_train[par_k] |
        pred_par_INLA[, par_k] > upper_train[par_k]
    )
})

colnames(outside_INLA) <- par_names

outside_counts <- colSums(outside_INLA, na.rm = TRUE)

outside_counts





INLA_range_summary <- data.frame(
  parameter = par_names,
  n_successful_INLA = sum(keep_INLA),
  n_outside_range = outside_counts,
  n_NA_after_filtering = colSums(is.na(pred_par_INLA_pp[keep_INLA, , drop = FALSE])),
  prop_outside_range = outside_counts / sum(keep_INLA)
)

print(INLA_range_summary)





































# ------------------------------------------------------------
# Load W1 and W2 without overwriting useful_output
# ------------------------------------------------------------

load_useful <- function(file) {
  env <- new.env()
  load(file, envir = env)
  env$useful_output
}

useful_W1 <- load_useful("C:/Users/qihan/Desktop/hom_W1/sim1_useful_output_W1.RData")
useful_W2 <- load_useful("C:/Users/qihan/Desktop/hom_W2/sim1_useful_output_W2.RData")


# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

par_names <- c(
  "sigma_Y", "scale_Y", "sigma_U1", "sigma_U2", "scale_U1", "scale_U2"
)

plot_order <- c(
  "sigma_Y", "sigma_U1", "sigma_U2",
  "scale_Y", "scale_U1", "scale_U2"
)

ylim_list <- list(
  sigma_Y  = c(0.5, 2),
  sigma_U1 = c(0.5, 2),
  sigma_U2 = c(0.5, 2),
  scale_Y  = c(0, 0.3),
  scale_U1 = c(0, 0.3),
  scale_U2 = c(0, 0.3)
)

true_values <- c(
  sigma_Y  = 1.5,
  sigma_U1 = 1,
  sigma_U2 = 1,
  scale_Y  = 0.2,
  scale_U1 = 0.15,
  scale_U2 = 0.15
)


# ------------------------------------------------------------
# INLA bounds
# ------------------------------------------------------------

lower_train <- c(
  sigma_Y  = 0.5,
  scale_Y  = 0.001,
  sigma_U1 = 0.5,
  sigma_U2 = 0.5,
  scale_U1 = 0.001,
  scale_U2 = 0.001
)

upper_train <- c(
  sigma_Y  = 2,
  scale_Y  = 0.3,
  sigma_U1 = 2,
  sigma_U2 = 2,
  scale_U1 = 0.3,
  scale_U2 = 0.3
)

in_bounds <- function(x, par_k) {
  is.finite(x) &
    x >= lower_train[par_k] &
    x <= upper_train[par_k]
}


# ------------------------------------------------------------
# Extract one window
# ------------------------------------------------------------

add_par_names <- function(x) {
  colnames(x) <- par_names
  x
}

extract_window_predictions <- function(useful_output) {
  
  pred  <- useful_output$predictions
  valid <- useful_output$valid_indices
  
  pred_par_50   <- add_par_names(pred$DSBI_50)
  pred_par_QBMC <- add_par_names(pred$QBMC)
  pred_par_INLA <- add_par_names(pred$INLA)
  
  keep_QBMC <- valid$keep_QBMC
  idx_QBMC  <- valid$idx_QBMC
  
  keep_INLA <- valid$keep_INLA
  
  pred_QBMC_valid <- pred_par_QBMC[keep_QBMC, , drop = FALSE][idx_QBMC, , drop = FALSE]
  
  pred_par_INLA_pp <- matrix(
    NA_real_,
    nrow = nrow(pred_par_INLA),
    ncol = length(par_names),
    dimnames = list(NULL, par_names)
  )
  
  for (par_k in par_names) {
    idx_INLA_k <- keep_INLA & in_bounds(pred_par_INLA[, par_k], par_k)
    pred_par_INLA_pp[idx_INLA_k, par_k] <- pred_par_INLA[idx_INLA_k, par_k]
  }
  
  list(
    DSBI = pred_par_50,
    MC   = pred_QBMC_valid,
    INLA = pred_par_INLA_pp
  )
}

res_W1 <- extract_window_predictions(useful_W1)
res_W2 <- extract_window_predictions(useful_W2)


# ------------------------------------------------------------
# Convert to long data frame for ggplot
# ------------------------------------------------------------

mat_to_long <- function(mat, window, method) {
  
  df <- data.frame(
    Window = window,
    Method = method,
    mat,
    check.names = FALSE
  )
  
  long <- reshape(
    df,
    varying = par_names,
    v.names = "Estimate",
    timevar = "Parameter",
    times = par_names,
    direction = "long"
  )
  
  rownames(long) <- NULL
  long
}

plot_df <- rbind(
  mat_to_long(res_W1$DSBI, "W1", "DSBI"),
  mat_to_long(res_W2$DSBI, "W2", "DSBI"),
  mat_to_long(res_W1$MC,   "W1", "MC"),
  mat_to_long(res_W2$MC,   "W2", "MC"),
  mat_to_long(res_W1$INLA, "W1", "INLA"),
  mat_to_long(res_W2$INLA, "W2", "INLA")
)

plot_df <- plot_df[is.finite(plot_df$Estimate), ]

plot_df$Parameter <- factor(plot_df$Parameter, levels = plot_order)
plot_df$Method <- factor(plot_df$Method, levels = c("DSBI", "MC", "INLA"))
plot_df$Window <- factor(plot_df$Window, levels = c("W1", "W2"))


# ------------------------------------------------------------
# Data for true-value horizontal lines
# ------------------------------------------------------------

true_df <- data.frame(
  Parameter = factor(names(true_values), levels = plot_order),
  true = as.numeric(true_values)
)


# ------------------------------------------------------------
# Data to force parameter-specific y-limits
# ------------------------------------------------------------

ylim_df <- do.call(
  rbind,
  lapply(plot_order, function(par_k) {
    data.frame(
      Parameter = factor(par_k, levels = plot_order),
      Method = "DSBI",
      y = ylim_list[[par_k]]
    )
  })
)

ylim_df$Method <- factor(ylim_df$Method, levels = c("DSBI", "MC", "INLA"))


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

library(ggplot2)

p_param_compare <- ggplot(
  plot_df,
  aes(x = Method, y = Estimate, fill = Window)
) +
  geom_boxplot(
    outlier.size = 0.6,
    colour = "black",
    position = position_dodge(width = 0.75)
  ) +
  geom_hline(
    data = true_df,
    aes(yintercept = true),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "red"
  ) +
  geom_blank(
    data = ylim_df,
    aes(x = Method, y = y),
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 2) +
  scale_fill_manual(
    values = c("W1" = "#9ecae1", "W2" = "#fdae6b")
  ) +
  theme_bw() +
  labs(
    title = bquote("Parameter estimate boxplots for two windows: " ~
                     W[1] == "[0,1]"^2 * " and " ~
                     W[2] == "[0,1.5]"^2),
    x = "Method",
    y = "Estimate",
    fill = "Window"
  ) +
  theme(
    strip.text = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    legend.position = "bottom"
  )

png(
  filename = "sim1_boxplot_param_compare_W1_W2_colour.png",
  width = 3000,
  height = 1800,
  res = 300
)

print(p_param_compare)

dev.off()