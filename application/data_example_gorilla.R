rm(list = ls())
set.seed(2026)
library(spatstat.geom)
library(spatstat.explore)
library(spatstat.random)
library(spatstat.model)
library(spatstat.data)
library(reticulate)
library(future)
library(future.apply)
library(progressr)
library(abind)
library(ggplot2)
use_python("C:/Users/qihan/anaconda3/envs/py39env/python.exe", required = TRUE)
plan(multisession, workers = max(1, parallelly::availableCores() - 1))
handlers(global = TRUE)
handlers("txtprogressbar")

data(gorillas, package = "spatstat.data")
X = gorillas
X_major = unmark(X[X$marks$group == "major"])
X_minor = unmark(X[X$marks$group == "minor"])
types = c("Type1", "Type2")
scale_im = function(z) {
  v = as.vector(z$v)
  eval.im((z - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
}
elevation = scale_im(gorillas.extra$elevation)
waterdist = scale_im(gorillas.extra$waterdist)
slopeangle = scale_im(gorillas.extra$slopeangle)
heat = gorillas.extra$heat
slopetype = gorillas.extra$slopetype
vegetation = gorillas.extra$vegetation

cov_list = list(
  elevation = elevation,
  waterdist = waterdist,
  slopeangle = slopeangle,
  heat = heat,
  slopetype = slopetype,
  vegetation = vegetation)

fit_kppm_safe = function(Xi, cov_list, group = "major") {
  fit = tryCatch({
    if (group == "major") {
      kppm(unmark(Xi) ~ elevation + waterdist + slopeangle + heat + vegetation, clusters = "LGCP", data = cov_list)
    } else {
      kppm(unmark(Xi) ~ elevation + waterdist + slopeangle + heat + slopetype + vegetation, clusters = "LGCP", data = cov_list)
    }
  }, error = function(e) NULL)
  fit
}

make_counts = function(P, xbreaks, ybreaks, nxy) {
  if (P$n == 0)
    return(matrix(0, nrow = nxy, ncol = nxy))
  
  cx = cut(P$x, xbreaks, include.lowest = TRUE, labels = FALSE)
  cy = cut(P$y, ybreaks, include.lowest = TRUE, labels = FALSE)
  keep = !is.na(cx) & !is.na(cy)
  idx = (cy[keep] - 1) * nxy + cx[keep]
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
  keep = !is.na(idx) & is.finite(df$lambda)
  E_sum = rowsum(df$lambda[keep] * pixel_area, idx[keep], reorder = FALSE)
  out = numeric(nxy * nxy)
  out[as.integer(rownames(E_sum))] = E_sum[, 1]
  matrix(out, nrow = nxy, byrow = FALSE)
}

pearson_residual_image = function(N, E, eps = 1e-6) {
  R = (N - E) / sqrt(pmax(E, eps))
  R
}

resize_count_to_base = function(count_mat, base_nxy) {
  nxy = nrow(count_mat)
  if (nxy == base_nxy) {
    return(count_mat)
  }
  if (base_nxy %% nxy != 0) {
    stop("base_nxy must be divisible by nxy (nrow(count_mat))")
  }
  fact = base_nxy / nxy
  kronecker(count_mat, matrix(1, nrow = fact, ncol = fact))
}

make_residual_count_image = function(X, W, Lambda_est, nxy_set = c(50), base_nxy = max(nxy_set), eps = 1e-6) {
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
    R1 = pearson_residual_image(N1, E1, eps = eps)
    R2 = pearson_residual_image(N2, E2, eps = eps)
    Rall = pearson_residual_image(Nall, Eall, eps = eps)
    R1 = resize_count_to_base(R1, base_nxy)
    R2 = resize_count_to_base(R2, base_nxy)
    Rall = resize_count_to_base(Rall, base_nxy)
    img_list[[k]] = abind(R1, R2, Rall, along = 3)
  }
  do.call(abind, c(img_list, along = 3))
}

get_lambda_trend = function(fit, beta = NULL, dimyx = c(50, 50)) {
  fit_ppm = as.ppm(fit)
  if (is.null(beta)) {
    out = predict(fit_ppm, type = "trend", dimyx = dimyx)
  } else {
    names(beta) = names(coef(fit_ppm))
    out = predict(fit_ppm, type = "trend", new.coef = beta, dimyx = dimyx)
  }
  out
}

get_beta_vec = function(fit1, fit2) {
  b1 = coef(fit1)
  b2 = coef(fit2)
  c(
    b1["(Intercept)"], b2["(Intercept)"],
    b1["elevation"], b2["elevation"],
    b1["waterdist"], b2["waterdist"],
    b1["slopeangle"], b2["slopeangle"],
    b1["heatModerate"], b2["heatModerate"],
    b1["heatCoolest"], b2["heatCoolest"],
    b2["slopetypeToe"],
    b2["slopetypeFlat"],
    b2["slopetypeMidslope"],
    b2["slopetypeUpper"],
    b2["slopetypeRidge"],
    b1["vegetationColonising"], b2["vegetationColonising"],
    b1["vegetationGrassland"], b2["vegetationGrassland"],
    b1["vegetationPrimary"], b2["vegetationPrimary"],
    b1["vegetationSecondary"], b2["vegetationSecondary"],
    b1["vegetationTransition"], b2["vegetationTransition"]
  )
}

get_features_from_pattern = function(X, r, cov_list, image_nxy = 50) {
  fit1 = fit_kppm_safe(X[marks(X) == "Type1"], cov_list, "major")
  fit2 = fit_kppm_safe(X[marks(X) == "Type2"], cov_list, "minor")
  if (is.null(fit1) || is.null(fit2))
    return(list(ok = FALSE))
  Lambda_est = list(
    Type1 = predict(fit1, type = "trend", dimyx = c(image_nxy, image_nxy)),
    Type2 = predict(fit2, type = "trend", dimyx = c(image_nxy, image_nxy))
  )
  
  Lc = tryCatch(
    Lcross.inhom(X, i = "Type1", j = "Type2", lambdaI = Lambda_est$Type1, lambdaJ = Lambda_est$Type2, r = r, correction = "border"),
    error = function(e) NULL
  )
  if (is.null(Lc)) return(list(ok = FALSE))
  Lc = as.numeric(Lc$border - Lc$r)
  
  L1 = tryCatch(Linhom(unmark(X[marks(X) == "Type1"]), lambda = Lambda_est$Type1, r = r, correction = "border"),
                error = function(e) NULL
  )
  if (is.null(L1)) return(list(ok = FALSE))
  L1 = as.numeric(L1$border - L1$r)
  
  L2 = tryCatch(Linhom(unmark(X[marks(X) == "Type2"]), lambda = Lambda_est$Type2, r = r, correction = "border"),
                error = function(e) NULL
  )
  if (is.null(L2)) return(list(ok = FALSE))
  L2 = as.numeric(L2$border - L2$r)
  
  count_img = make_residual_count_image(X, Window(X), Lambda_est, nxy_set = c(image_nxy), base_nxy = image_nxy)
  beta_vec = get_beta_vec(fit1, fit2)
  
  list(ok = TRUE, L1 = L1, L2 = L2, Lc = Lc, beta = beta_vec, count_img = count_img, Np_type1 = sum(marks(X) == "Type1"), Np_type2 = sum(marks(X) == "Type2"), Np_all = X$n)
}

stack_feature_results = function(res_list, response_par = NULL) {
  keep_true = sapply(res_list, function(x) isTRUE(x$ok))
  res_ok = res_list[keep_true]
  L1 = do.call(cbind, lapply(res_ok, function(x) x$L1))
  L2 = do.call(cbind, lapply(res_ok, function(x) x$L2))
  Lc = do.call(cbind, lapply(res_ok, function(x) x$Lc))
  keep_pcf = (colSums(is.na(L1)) == 0) & (colSums(is.na(L2)) == 0) & (colSums(is.na(Lc)) == 0)
  L1 = L1[, keep_pcf, drop = FALSE]
  L2 = L2[, keep_pcf, drop = FALSE]
  Lc = Lc[, keep_pcf, drop = FALSE]
  L1_arr = array(as.numeric(t(L1)), dim = c(ncol(L1), nrow(L1), 1))
  L2_arr = array(as.numeric(t(L2)), dim = c(ncol(L2), nrow(L2), 1))
  Lc_arr = array(as.numeric(t(Lc)), dim = c(ncol(Lc), nrow(Lc), 1))
  L_Y = abind(Lc_arr, along = 3)
  d1_12 = L1_arr - Lc_arr
  d2_21 = L2_arr - Lc_arr
  U_L = abind(d1_12, d2_21, along = 3)
  pcf_input = abind(L_Y, U_L, along = 3)
  count_list = lapply(res_ok, function(x) x$count_img)
  count_img = aperm(do.call(abind, c(count_list[keep_pcf], along = 4)), c(4, 1, 2, 3))
  beta_mat = do.call(rbind, lapply(res_ok, function(x) x$beta))[keep_pcf, , drop = FALSE]
  Np_type1 = as.numeric(lapply(res_ok, function(x) x$Np_type1))[keep_pcf]
  Np_type2 = as.numeric(lapply(res_ok, function(x) x$Np_type2))[keep_pcf]
  Np_all = as.numeric(lapply(res_ok, function(x) x$Np_all))[keep_pcf]
  aux_input = beta_mat
  out = list(keep_true = keep_true, keep_pcf = keep_pcf, pcf_input = pcf_input, count_img = count_img, aux_input = aux_input, beta_input = beta_mat, Np = cbind(Np_all, Np_type1, Np_type2))
  
  if (!is.null(response_par)) {
    response_par_ok = response_par[keep_true, , drop = FALSE]
    out$Y = sqrt(response_par_ok[keep_pcf, , drop = FALSE])
    out$true_par = response_par_ok[keep_pcf, , drop = FALSE]
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

get_mean_sd = function(res_list) {
  arr = simplify2array(res_list)
  list(mean = apply(arr, c(1, 2), mean), sd = apply(arr, c(1, 2), sd))
}

fit_major = fit_kppm_safe(X_major, cov_list, "major")
fit_minor = fit_kppm_safe(X_minor, cov_list, "minor")
ppm_major = as.ppm(fit_major)
ppm_minor = as.ppm(fit_minor)
beta_major = coef(ppm_major)
beta_minor = coef(ppm_minor)
image_nxy = 50

lambda_mean_trend_major = get_lambda_trend(fit_major, dimyx = c(image_nxy, image_nxy))
lambda_mean_trend_minor = get_lambda_trend(fit_minor, dimyx = c(image_nxy, image_nxy))
mean_trend_major = eval.im(log(lambda_mean_trend_major))
mean_trend_minor = eval.im(log(lambda_mean_trend_minor))
W_major = Window(mean_trend_major)
W_minor = Window(mean_trend_minor)
W = union.owin(W_major, W_minor)
zero_mu = eval.im(0 * mean_trend_major)
X_obs = superimpose(Type1 = X_major, Type2 = X_minor, W = W)
marks(X_obs) = factor(marks(X_obs), levels = types)
rmax = rmax.rule("K", W)
r = seq(0, rmax, length.out = 513)

obs_one = get_features_from_pattern(X_obs, r, cov_list, image_nxy)
obs_data = stack_feature_results(list(obs_one))
pcf_obs = obs_data$pcf_input
img_obs = obs_data$count_img
aux_obs = obs_data$aux_input
img_obs_1 = img_obs[1, , , 1]
img_obs_2 = img_obs[1, , , 2]
img_obs_3 = img_obs[1, , , 3]


ntrain = 10000
summary_major = summary(ppm_major)
summary_minor = summary(ppm_minor)
CIs_major = cbind(summary_major$coefs.SE.CI$CI95.lo, summary_major$coefs.SE.CI$CI95.hi)
CIs_minor = cbind(summary_minor$coefs.SE.CI$CI95.lo, summary_minor$coefs.SE.CI$CI95.hi)
beta0_1_train = runif(ntrain, CIs_major[1, 1], CIs_major[1, 2])
beta_elevation_1_train = runif(ntrain, CIs_major[2, 1], CIs_major[2, 2])
beta_waterdist_1_train = runif(ntrain, CIs_major[3, 1], CIs_major[3, 2])
beta_slopeangle_1_train = runif(ntrain, CIs_major[4, 1], CIs_major[4, 2])
beta_heatModerate_1_train = runif(ntrain, CIs_major[5, 1], CIs_major[5, 2])
beta_heatCoolest_1_train = runif(ntrain, CIs_major[6, 1], CIs_major[6, 2])
beta_vegetationColonising_1_train = runif(ntrain, CIs_major[7, 1], CIs_major[7, 2])
beta_vegetationGrassland_1_train = runif(ntrain, CIs_major[8, 1], CIs_major[8, 2])
beta_vegetationPrimary_1_train = runif(ntrain, CIs_major[9, 1], CIs_major[9, 2])
beta_vegetationSecondary_1_train = runif(ntrain, CIs_major[10, 1], CIs_major[10, 2])
beta_vegetationTransition_1_train = runif(ntrain, CIs_major[11, 1], CIs_major[11, 2])

beta0_2_train = runif(ntrain, CIs_minor[1, 1], CIs_minor[1, 2])
beta_elevation_2_train = runif(ntrain, CIs_minor[2, 1], CIs_minor[2, 2])
beta_waterdist_2_train = runif(ntrain, CIs_minor[3, 1], CIs_minor[3, 2])
beta_slopeangle_2_train = runif(ntrain, CIs_minor[4, 1], CIs_minor[4, 2])
beta_heatModerate_2_train = runif(ntrain, CIs_minor[5, 1], CIs_minor[5, 2])
beta_heatCoolest_2_train = runif(ntrain, CIs_minor[6, 1], CIs_minor[6, 2])
beta_slopetypeToe_2_train = runif(ntrain, CIs_minor[7, 1], CIs_minor[7, 2])
beta_slopetypeFlat_2_train = runif(ntrain, CIs_minor[8, 1], CIs_minor[8, 2])
beta_slopetypeMidslope_2_train = runif(ntrain, CIs_minor[9, 1], CIs_minor[9, 2])
beta_slopetypeUpper_2_train = runif(ntrain, CIs_minor[10, 1], CIs_minor[10, 2])
beta_slopetypeRidge_2_train = runif(ntrain, CIs_minor[11, 1], CIs_minor[11, 2])
beta_vegetationColonising_2_train = runif(ntrain, CIs_minor[12, 1], CIs_minor[12, 2])
beta_vegetationGrassland_2_train = runif(ntrain, CIs_minor[13, 1], CIs_minor[13, 2])
beta_vegetationPrimary_2_train = runif(ntrain, CIs_minor[14, 1], CIs_minor[14, 2])
beta_vegetationSecondary_2_train = runif(ntrain, CIs_minor[15, 1], CIs_minor[15, 2])
beta_vegetationTransition_2_train = runif(ntrain, CIs_minor[16, 1], CIs_minor[16, 2])

all_dists = nndist(X)
all_dists_major = nndist(X_major)
all_dists_minor = nndist(X_minor)
scale_interval = range(all_dists[all_dists != 0])
scale_interval_major = range(all_dists_major[all_dists_major != 0])
scale_interval_minor = range(all_dists_minor[all_dists_minor != 0])

sigma_y_train = runif(ntrain, 0.1, 1.5)
sigma_u_1_train = runif(ntrain, 0.1, 1.5)
sigma_u_2_train = runif(ntrain, 0.1, 1.5)
scale_y_train = runif(ntrain, scale_interval[1], scale_interval[2])
scale_u_1_train = runif(ntrain, scale_interval_major[1], scale_interval_major[2])
scale_u_2_train = runif(ntrain, scale_interval_minor[1], scale_interval_minor[2])

train_par = cbind(sigma_y_train, scale_y_train, sigma_u_1_train, sigma_u_2_train, scale_u_1_train, scale_u_2_train)
colnames(train_par) = c("sigma_Y", "scale_Y", "sigma_U1", "sigma_U2", "scale_U1", "scale_U2")

simulate_one_train = function(j) {
  beta1 = c(beta0_1_train[j], beta_elevation_1_train[j], beta_waterdist_1_train[j], beta_slopeangle_1_train[j], beta_heatModerate_1_train[j], beta_heatCoolest_1_train[j], beta_vegetationColonising_1_train[j], beta_vegetationGrassland_1_train[j], beta_vegetationPrimary_1_train[j], beta_vegetationSecondary_1_train[j], beta_vegetationTransition_1_train[j])
  beta2 = c(beta0_2_train[j], beta_elevation_2_train[j], beta_waterdist_2_train[j], beta_slopeangle_2_train[j], beta_heatModerate_2_train[j], beta_heatCoolest_2_train[j], beta_slopetypeToe_2_train[j], beta_slopetypeFlat_2_train[j], beta_slopetypeMidslope_2_train[j], beta_slopetypeUpper_2_train[j], beta_slopetypeRidge_2_train[j], beta_vegetationColonising_2_train[j], beta_vegetationGrassland_2_train[j], beta_vegetationPrimary_2_train[j], beta_vegetationSecondary_2_train[j], beta_vegetationTransition_2_train[j])
  lambda1 = get_lambda_trend(fit_major, beta = beta1, dimyx = c(image_nxy, image_nxy))
  lambda2 = get_lambda_trend(fit_minor, beta = beta2, dimyx = c(image_nxy, image_nxy))
  mean1 = eval.im(log(lambda1))
  mean2 = eval.im(log(lambda2))
  repeat {
    zero_mu_j = eval.im(0 * mean1)
    Y_shared = log(attr(rLGCP(model = "exponential", mu = zero_mu_j, var = sigma_y_train[j]^2, scale = scale_y_train[j], win = W, saveLambda = TRUE), "Lambda"))
    mu1 = eval.im(mean1 - 0.5 * sigma_y_train[j]^2 - 0.5 * sigma_u_1_train[j]^2 + Y_shared)
    mu2 = eval.im(mean2 - 0.5 * sigma_y_train[j]^2 - 0.5 * sigma_u_2_train[j]^2 + Y_shared)
    X1 = rLGCP(model = "exponential", mu = mu1, var = sigma_u_1_train[j]^2, scale = scale_u_1_train[j], win = W, saveLambda = FALSE)
    X2 = rLGCP(model = "exponential", mu = mu2, var = sigma_u_2_train[j]^2, scale = scale_u_2_train[j], win = W, saveLambda = FALSE)
    if (X1$n > 0 && X2$n > 0) break
  }
  Xj = superimpose(Type1 = X1, Type2 = X2, W = W)
  marks(Xj) = factor(marks(Xj), levels = types)
  get_features_from_pattern(Xj, r, cov_list, image_nxy)
}

plan(sequential)
gc()
plan(multisession, workers = max(1, 6))
with_progress({
  p = progressor(along = seq_len(ntrain))
  res_train = future_lapply(
    seq_len(ntrain),
    function(j) {
      p(sprintf("iteration %d/%d", j, ntrain))
      simulate_one_train(j)
    },
    future.seed = 2026,
    future.packages = c("spatstat.geom", "spatstat.explore", "spatstat.random", "spatstat.model", "abind")
  )
})

train_data = stack_feature_results(res_train, response_par = train_par)
pcf_train = train_data$pcf_input
img_train = train_data$count_img
Y_train = train_data$Y
aux_train = train_data$aux_input
scaled_aux = scale_by_train(aux_train, aux_obs)
aux_train_scaled = scaled_aux$train
aux_obs_scaled = scaled_aux$test

R_address = "C:/Users/qihan/Desktop/data_example"
setwd(R_address)
source_python("BLGCP_NN.py")
print(get_device_info("auto"))

n_runs = 10
results = vector("list", n_runs)
for (i in 1:n_runs) {
  seed_j = 3053 + i
  set.seed(seed_j)
  seed_i = sample.int(10000, 1)
  set_all_seeds(seed_i)
  model_hybrid = train_and_save_hybrid_NN_model(pcf_train, img_train, Y_train, aux_train_scaled, model_path = sprintf("C:/Users/qihan/Desktop/data_example/gorilla_model_run_%d.pth", i), batch_size = 100, epochs = 30, lr = 1e-3, requested_device = "auto", seed = seed_i, verbose = TRUE, standardize_img = TRUE)
  pred_i = load_hybrid_NN_model_predict(pcf_obs, img_obs, aux_obs_scaled, sprintf("C:/Users/qihan/Desktop/data_example/gorilla_model_run_%d.pth", i), requested_device = "auto", predict_batch_size = 2000)^2
  results[[i]] = pred_i
}

pred_out = get_mean_sd(results)
pred_mean = pred_out$mean
pred_sd = pred_out$sd
pred_min = apply(simplify2array(results), c(1, 2), min)
pred_max = apply(simplify2array(results), c(1, 2), max)
colnames(pred_mean) = colnames(train_par)
colnames(pred_sd) = colnames(train_par)
colnames(pred_min) = colnames(train_par)
colnames(pred_max) = colnames(train_par)
pred_table = data.frame(parameter = colnames(train_par), estimate = as.numeric(pred_mean[1, ]), sd = as.numeric(pred_sd[1, ]), min = as.numeric(pred_min[1, ]), max = as.numeric(pred_max[1, ]))
print(pred_table)

env_mlgcp_cross_fast = function(theta, method_name = "DSBI", nsim = 99, r_use = r) {
  sigma_y = theta[1]
  scale_y = theta[2]
  sigma_u1 = theta[3]
  sigma_u2 = theta[4]
  scale_u1 = theta[5]
  scale_u2 = theta[6]
  env_cross = envelope(X_obs,
                       fun = function(Y, r) {
                         Kcross.inhom(Y, i = "Type1", j = "Type2", r = r, lambdaI = lambda_mean_trend_major, lambdaJ = lambda_mean_trend_minor, correction = "border")
                       },
                       r = r_use,
                       simulate = expression({
                         Y_shared = log(attr(rLGCP(model = "exponential", mu = zero_mu, var = sigma_y^2, scale = scale_y, win = W, saveLambda = TRUE), "Lambda"))
                         mu_major = eval.im(mean_trend_major - 0.5 * sigma_y^2 - 0.5 * sigma_u1^2 + Y_shared)
                         mu_minor = eval.im(mean_trend_minor - 0.5 * sigma_y^2 - 0.5 * sigma_u2^2 + Y_shared)
                         X_major_sim = rLGCP(model = "exponential", mu = mu_major, win = W, var = sigma_u1^2, scale = scale_u1, saveLambda = FALSE)
                         X_minor_sim = rLGCP(model = "exponential", mu = mu_minor, win = W, var = sigma_u2^2, scale = scale_u2, saveLambda = FALSE)
                         Xsim = superimpose(Type1 = X_major_sim, Type2 = X_minor_sim, W = W)
                         marks(Xsim) = factor(marks(Xsim), levels = c("Type1", "Type2"))
                         Xsim
                       }),
                       nsim = nsim, savefuns = TRUE, global = FALSE, verbose = FALSE)
  test_cross = dclf.test(env_cross)
  plot(env_cross, main = paste0("Cross-type envelope (", method_name, ")"))
  list(env = env_cross, test = test_cross)
}

env_mlgcp_major_fast = function(theta, method_name = "DSBI", nsim = 99, r_use = r) {
  sigma_y = theta[1]
  scale_y = theta[2]
  sigma_u1 = theta[3]
  scale_u1 = theta[5]
  X_obs_major = unmark(X_obs[marks(X_obs) == "Type1"])
  env_major = envelope(X_obs_major,
                       fun = function(Y, r) {
                         Kinhom(Y, r = r, lambda = lambda_mean_trend_major, correction = "border")
                       }, r = r_use,
                       simulate = expression({Y_shared = log(attr(rLGCP(model = "exponential", mu = zero_mu, var = sigma_y^2, scale = scale_y, win = W, saveLambda = TRUE), "Lambda"))
                         mu_major = eval.im(mean_trend_major - 0.5 * sigma_y^2 - 0.5 * sigma_u1^2 + Y_shared)
                         rLGCP(model = "exponential", mu = mu_major, win = W, var = sigma_u1^2, scale = scale_u1, saveLambda = FALSE)
                       }),
                       nsim = nsim, savefuns = TRUE, global = FALSE, verbose = FALSE)
  test_major = dclf.test(env_major)
  plot(env_major, main = paste0("Major group envelope (", method_name, ")"))
  list(env = env_major, test = test_major)
}

env_mlgcp_minor_fast = function(theta, method_name = "DSBI", nsim = 99, r_use = r) {
  sigma_y = theta[1]
  scale_y = theta[2]
  sigma_u2 = theta[4]
  scale_u2 = theta[6]
  X_obs_minor = unmark(X_obs[marks(X_obs) == "Type2"])
  env_minor = envelope(X_obs_minor,
    fun = function(Y, r) {
      Kinhom(Y, r = r, lambda = lambda_mean_trend_minor, correction = "border")
    }, r = r_use,
    simulate = expression({Y_shared = log(attr(rLGCP(model = "exponential", mu = zero_mu, var = sigma_y^2, scale = scale_y, win = W, saveLambda = TRUE), "Lambda"))
    mu_minor = eval.im(mean_trend_minor - 0.5 * sigma_y^2 - 0.5 * sigma_u2^2 + Y_shared)
    rLGCP(model = "exponential", mu = mu_minor, win = W, var = sigma_u2^2, scale = scale_u2, saveLambda = FALSE)
    }), nsim = nsim, savefuns = TRUE, global = FALSE, verbose = FALSE)
  test_minor = dclf.test(env_minor)
  plot(env_minor, main = paste0("Minor group envelope (", method_name, ")"))
  list(env = env_minor, test = test_minor)
}

theta_hat = as.numeric(pred_mean[1, ])
set.seed(2025)
plot.res = 600
png("gorilla_env_cross_DSBI.png", width = 6.2 * plot.res, height = 4.5 * plot.res, res = plot.res)
res_cross = env_mlgcp_cross_fast(theta_hat, "DSBI")
dev.off()
set.seed(2025)
png("gorilla_env_major_DSBI.png", width = 6.2 * plot.res, height = 4.5 * plot.res, res = plot.res)
res_major = env_mlgcp_major_fast(theta_hat, "DSBI")
dev.off()
set.seed(2025)
png("gorilla_env_minor_DSBI.png", width = 6.2 * plot.res, height = 4.5 * plot.res, res = plot.res)
res_minor = env_mlgcp_minor_fast(theta_hat, "DSBI")
dev.off()

save.image("gorilla_app.RData")



res_cross$env
res_major$env
res_minor$env

library(GET)
res_get1 <- global_envelope_test(curve_sets = list(Cross = res_cross$env, Major = res_major$env, Minor = res_minor$env), type = "erl")
res_get1
plot(res_get1)


plot.res <- 600
png("G_env.png", width = 6.2*plot.res, height = 4.5*plot.res, res = plot.res)
plot(res_get1)
dev.off()


data(gorillas, package = "spatstat.data")
X1 = gorillas
X1_major = unmark(X1[X1$marks$group == "major"])
X1_minor = unmark(X1[X1$marks$group == "minor"])
win <- as.polygonal(gorillas$window)
win_df <- do.call(rbind, lapply(seq_along(win$bdry), function(i) {
  data.frame(
    x = win$bdry[[i]]$x,
    y = win$bdry[[i]]$y,
    id = i
  )
}))
major_df <- data.frame(x = X1_major$x, y = X1_major$y, group = "Major")
minor_df <- data.frame(x = X1_minor$x, y = X1_minor$y, group = "Minor")
points_df <- rbind(major_df, minor_df)
plot.res = 600
png("gorilla_groups.png", width = 6.2 * plot.res, height = 4.5 * plot.res, res = 600)
ggplot() +
  geom_polygon(data = win_df, aes(x = x, y = y, group = id), fill = NA, color = "black") +
  geom_point(data = points_df, aes(x = x, y = y, color = group), size = 1, alpha = 0.5) +
  labs(
    title = "Gorillas nests with groups",
    x = "X coordinate",
    y = "Y coordinate",
    color = "Group"
  ) +
  scale_color_manual(values = c("Major" = "black", "Minor" = "red")) +
  coord_equal() +
  theme_minimal()
dev.off()


save.image("gorilla_application_full.RData")



