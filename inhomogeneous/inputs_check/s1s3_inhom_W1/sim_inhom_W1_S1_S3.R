set.seed(2026)
library(reticulate)
use_python("C:/Users/qihan/anaconda3/envs/py39env/python.exe", required = TRUE)
t0 = Sys.time()
library(spatstat.geom)
library(spatstat.random)
library(spatstat.model)
library(future)
library(future.apply)
library(ggplot2)
library(terra)
library(abind)
library(progressr)
plan(multisession, workers = parallelly::availableCores() - 1)
R_address = "C:/Users/qihan/Desktop/sim_inhom_S1_S3"
setwd(R_address)

simulate_grf = function(nx, ny, L, range, var) {
  mu0 = as.im(0, owin(xrange = c(0, L), yrange = c(0, L)), dimyx = c(ny, nx))
  sim_field = rLGCP(model = "exponential", mu = mu0, var = var, scale = range, win = as.owin(mu0), saveLambda = TRUE)
  field_im = attr(sim_field, "Lambda")
  field_im = eval.im(log(field_im))
  as.matrix(field_im)
}

simulate_MLGCP_point_pattern = function(params) {
  Y = simulate_grf(params$nx, params$ny, params$L, range = params$Y_scale, var = 1)
  Z1 = params$Z_field1
  Z2 = params$Z_field2
  U_list = vector("list", params$p)
  for (j in seq_len(params$p)) {
    U_list[[j]] = simulate_grf(params$nx, params$ny, params$L, range = params$U_scale[j], var = 1)
  }
  Lambda_list = vector("list", params$p)
  X_list = vector("list", params$p)
  for (i in seq_len(params$p)) {
    shared_term = params$Y_var * Y - (params$Y_var^2) / 2
    individual_term = params$U_var[i] * U_list[[i]] - (params$U_var[i]^2) / 2
    Lambda_list[[i]] = exp(params$beta0[i] + params$beta1[i] * Z1 + params$beta2[i] * Z2) *exp(shared_term + individual_term)
    Lambda_im = im(mat = Lambda_list[[i]], xrange = c(0, params$L), yrange = c(0, params$L))
    X_list[[i]] = rpoispp(Lambda_im)
  }
  list(X_list = X_list, Lambda_list = Lambda_list, fields = list(Y = Y, Z1 = Z1, Z2 = Z2, U_list = U_list), params = params)
}

has_points = function(sim) {
  all(sapply(sim$X_list, function(x) x$n > 0))
}

fit_kppm_safe = function(Xi, Z_im1, Z_im2) {
  fit = tryCatch(
    kppm(unmark(Xi) ~ Z_im1 + Z_im2, clusters = "LGCP"),
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

extract_one_case_features = function(params, image_nxy_set = NULL) {
  W = owin(c(0, 1), c(0, 1))
  types = c("Type1", "Type2")
  
  if (is.null(image_nxy_set))
    image_nxy_set = params$nx
  
  repeat {
    sim = simulate_MLGCP_point_pattern(params)
    if (has_points(sim)) break
  }
  
  X = superimpose(Type1 = sim$X_list[[1]], Type2 = sim$X_list[[2]], W = W)
  marks(X) = factor(marks(X), levels = types)
  Z_im1 = as.im(sim$fields$Z1, owin(xrange = c(0, params$L), yrange = c(0, params$L)), dimyx = c(params$ny, params$nx))
  Z_im2 = as.im(sim$fields$Z2, owin(xrange = c(0, params$L), yrange = c(0, params$L)), dimyx = c(params$ny, params$nx))
  fit1 = fit_kppm_safe(X[marks(X) == "Type1"], Z_im1, Z_im2)
  fit2 = fit_kppm_safe(X[marks(X) == "Type2"], Z_im1, Z_im2)
  
  if (is.null(fit1) || is.null(fit2)) {
    return(list(ok = FALSE))
  }
  
  beta_kppm1 = coef(fit1)
  beta_kppm2 = coef(fit2)
  Lambda_est = list(Type1 = predict(fit1, type = "trend", dimyx = c(params$ny, params$nx)), Type2 = predict(fit2, type = "trend", dimyx = c(params$ny, params$nx)))
  
  count_img = make_residual_count_image(X, W, Lambda_est, image_nxy_set, base_nxy = max(image_nxy_set))
  beta_vec = c(beta_kppm1["(Intercept)"], beta_kppm2["(Intercept)"], beta_kppm1["Z_im1"], beta_kppm2["Z_im1"], beta_kppm1["Z_im2"], beta_kppm2["Z_im2"])
  
  list(ok = TRUE, beta = beta_vec, count_img = count_img)
}


build_param_list = function(idx, beta0_1, beta0_2, beta1_1, beta1_2, beta2_1, beta2_2, sigma_y, sigma_u_1, sigma_u_2, scale_y, scale_u_1, scale_u_2, Z1, Z2) {
  list(p = 2, L = 1, nx = 50, ny = 50, beta0 = c(beta0_1[idx], beta0_2[idx]), beta1 = c(beta1_1[idx], beta1_2[idx]), beta2 = c(beta2_1[idx], beta2_2[idx]), Y_var = sigma_y[idx], U_var = c(sigma_u_1[idx], sigma_u_2[idx]), Y_scale = scale_y[idx], U_scale = c(scale_u_1[idx], scale_u_2[idx]), Z_field1 = Z1, Z_field2 = Z2)
}

stack_feature_results = function(res_list, response_par = NULL) {
  keep_true = sapply(res_list, function(x) isTRUE(x$ok))
  res_ok = res_list[keep_true]
  count_list = lapply(res_ok, function(x) x$count_img)
  count_img = aperm(do.call(abind, c(count_list, along = 4)), c(4, 1, 2, 3))
  beta_mat = do.call(rbind, lapply(res_ok, function(x) x$beta))
  aux_input = beta_mat
  out = list(keep_true = keep_true, count_img = count_img, aux_input = aux_input, beta_input = beta_mat)
  
  if (!is.null(response_par)) {
    response_par_ok = response_par[keep_true, , drop = FALSE]
    out$Y = sqrt(response_par_ok)
    out$true_par = response_par_ok
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

compare_par = function(true_par, pred_par, par_names = NULL) {
  true_par = as.matrix(true_par)
  pred_par = as.matrix(pred_par)
  data.frame(
    parameter = par_names,
    bias = apply(pred_par - true_par, 2, mean),
    rmse = sqrt(apply((pred_par - true_par)^2, 2, mean)),
    mape = apply(abs((pred_par - true_par) / pmax(true_par, 1e-8)), 2, mean) * 100,
    cor = sapply(1:ncol(true_par), function(j) cor(true_par[, j], pred_par[, j]))
  )
}

plot_all_par = function(true_par, pred_par, main_prefix = "NN") {
  true_par = as.matrix(true_par)
  pred_par = as.matrix(pred_par)
  sigma_id = c(1, 3, 4)
  scale_id = c(2, 5, 6)
  plot_id = c(sigma_id, scale_id)
  par_names = c("sigma Y", "scale Y", "sigma U1", "sigma U2", "scale U1", "scale U2")
  old_par = par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  for (k in seq_along(plot_id)) {
    j = plot_id[k]
    x = true_par[, j]
    y = pred_par[, j]
    if (j %in% sigma_id) {
      x_lim = c(1, 1.5)
      y_lim = c(0.5, 2)
    } else {
      x_lim = c(0.07, 0.13)
      y_lim = c(0.04, 0.18)
    }
    plot(x, y, xlab = "True", ylab = "Predicted", main = paste0(par_names[j], " (", main_prefix, ")"), pch = 19, cex = 0.6, col = rgb(0, 0, 0, 0.4), xlim = x_lim, ylim = y_lim)
    abline(0, 1, col = 2, lwd = 2)
    rmse_j = sqrt(mean((y - x)^2, na.rm = TRUE))
    legend("topleft", legend = paste0("RMSE = ", round(rmse_j, 3)), bty = "n", text.col = "red", cex = 1.1)
  }
}

plot_beta_ppm = function(true_beta, beta_hat) {
  true_beta = as.matrix(true_beta)
  beta_hat = as.matrix(beta_hat)
  beta_names = c("beta0", "beta1", "beta2")
  beta_id = rbind(c(1, 3, 5), c(2, 4, 6))
  old_par = par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  for (type_i in 1:2) {
    for (k in 1:3) {
      id = beta_id[type_i, k]
      x = as.numeric(true_beta[, id])
      y = as.numeric(beta_hat[, id])
      if (k == 1) {
        x_lim = c(5.5, 6.5)
        y_lim = c(4, 8)
      } else {
        x_lim = c(-0.5, 0.5)
        y_lim = c(-1, 1)
      }
      plot(x, y, xlab = "True", ylab = "Estimated", main = paste0(beta_names[k], " (type ", type_i, ")"), pch = 19, cex = 0.6, col = rgb(0, 0, 0, 0.4), xlim = x_lim, ylim = y_lim)
      abline(0, 1, col = 2, lwd = 2)
      rmse_j = sqrt(mean((y - x)^2, na.rm = TRUE))
      bias_j = mean(y - x, na.rm = TRUE)
      legend("topleft", inset = c(0, 0), legend = c(paste0("RMSE = ", round(rmse_j, 3)), paste0("Bias = ", round(bias_j, 3))), bty = "n", text.col = "red", cex = 1.1)
    }
  }
}

get_mean_sd = function(res_list) {
  arr = simplify2array(res_list)
  list(mean = apply(arr, c(1, 2), mean), sd = apply(arr, c(1, 2), sd))
}
# ------------------------------------------------------------------------------
# Covariates
set.seed(123)
Z1 = simulate_grf(50, 50, 1, range = 0.01, var = 1)
Z1 = Z1 - mean(Z1)
set.seed(456)
Z2 = simulate_grf(50, 50, 1, range = 0.01, var = 1)
Z2 = Z2 - mean(Z2)

set.seed(2026)
image_nxy_set = c(50) 

# ------------------------------------------------------------------------------
# Training
ntrain = 100000
beta0_1_train = runif(ntrain, 5, 7)
beta0_2_train = runif(ntrain, 5, 7)
beta1_1_train = runif(ntrain, -1, 1)
beta1_2_train = runif(ntrain, -1, 1)
beta2_1_train = runif(ntrain, -1, 1)
beta2_2_train = runif(ntrain, -1, 1)
sigma_y_train = runif(ntrain, 0.5, 2)
scale_y_train = runif(ntrain, 0.001, 0.2)
sigma_u_1_train = runif(ntrain, 0.5, 2)
sigma_u_2_train = runif(ntrain, 0.5, 2)
scale_u_1_train = runif(ntrain, 0.001, 0.2)
scale_u_2_train = runif(ntrain, 0.001, 0.2)
train_par = cbind(sigma_y_train, scale_y_train, sigma_u_1_train, sigma_u_2_train, scale_u_1_train, scale_u_2_train)

plan(sequential)
gc()
plan(multisession, workers = max(1, parallelly::availableCores() - 1))
handlers(global = TRUE)
handlers("txtprogressbar")
with_progress({
  p = progressor(along = seq_len(ntrain))
  res_train = future_lapply(
    seq_len(ntrain),
    function(j) {
      p(sprintf("iteration %d/%d", j, ntrain))
      params = build_param_list(j, beta0_1_train, beta0_2_train, beta1_1_train, beta1_2_train, beta2_1_train, beta2_2_train, sigma_y_train, sigma_u_1_train, sigma_u_2_train, scale_y_train, scale_u_1_train, scale_u_2_train, Z1, Z2)
      extract_one_case_features(params, image_nxy_set = image_nxy_set)
    },
    future.seed = 2026,
    future.packages = c("spatstat.geom", "spatstat.random", "spatstat.model", "terra", "abind")
  )
})

train_data = stack_feature_results(res_train, response_par = train_par)
img_train = train_data$count_img
Y_train = train_data$Y
aux_train = train_data$aux_input

save.image("C:/Users/qihan/Desktop/sim_inhom_S1_S3/sim2_training_W1_S1_S3.RData")
# ------------------------------------------------------------------------------
# Test
ntest = 5000
set.seed(2027)
beta0_1_test = runif(ntest, 5.5, 6.5)
beta0_2_test = runif(ntest, 5.5, 6.5)
beta1_1_test = runif(ntest, -0.5, 0.5)
beta1_2_test = runif(ntest, -0.5, 0.5)
beta2_1_test = runif(ntest, -0.5, 0.5)
beta2_2_test = runif(ntest, -0.5, 0.5)
sigma_y_test   = runif(ntest, 1, 1.5)
scale_y_test   = runif(ntest, 0.07, 0.13)
sigma_u_1_test = runif(ntest, 1, 1.5)
sigma_u_2_test = runif(ntest, 1, 1.5)
scale_u_1_test = runif(ntest, 0.07, 0.13)
scale_u_2_test = runif(ntest, 0.07, 0.13)
test_par_full = cbind(sigma_y_test, scale_y_test, sigma_u_1_test, sigma_u_2_test, scale_u_1_test, scale_u_2_test)

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
      p(sprintf("iteration %d/%d", j, ntest))
      params = build_param_list(j, beta0_1_test, beta0_2_test, beta1_1_test, beta1_2_test, beta2_1_test, beta2_2_test, sigma_y_test, sigma_u_1_test, sigma_u_2_test, scale_y_test, scale_u_1_test, scale_u_2_test, Z1, Z2)
      extract_one_case_features(params, image_nxy_set = image_nxy_set)
    },
    future.seed = 2027, 
    future.packages = c("spatstat.geom", "spatstat.random", "spatstat.model", "terra", "abind")
  )
})
test_data = stack_feature_results(res_test, response_par = test_par_full)
img_test = test_data$count_img
aux_test = test_data$aux_input
test_par = test_data$true_par
scaled_aux = scale_by_train(aux_train, aux_test)
aux_train_scaled = scaled_aux$train
aux_test_scaled = scaled_aux$test

# ------------------------------------------------------------------------------
# DSBI 
source_python("BLGCP_NN_S1_S3_no_Ls.py")
print(get_device_info("auto"))

# ------------------------------------------------------------------------------
img_train_50 = img_train[, , , 1:3, drop = FALSE]
img_test_50  = img_test[, , , 1:3, drop = FALSE]
n_runs = 10
res_all_50 = vector("list", n_runs)
predict_batch_size = 2000
for (i in 1:n_runs) {
  seed_j = 2025 + i
  set.seed(seed_j)
  seed_i = sample.int(10000, 1)
  set_all_seeds(seed_i)
  model_path_i = sprintf("C:/Users/qihan/Desktop/sim_inhom_S1_S3/W1_model_S1_S3_run_%d.pth", i)
  model_hybrid = train_and_save_hybrid_NN_model(img_train_50, Y_train, aux_train_scaled, model_path = model_path_i, batch_size = 100, epochs = 30, lr = 1e-3, requested_device = "auto", seed = seed_i, verbose = TRUE, standardize_img = TRUE)
  pred_par = load_hybrid_NN_model_predict(img_test_50, aux_test_scaled, model_path_i, requested_device = "auto", predict_batch_size = predict_batch_size)^2
  res_all_50[[i]] = pred_par
  rm(model_hybrid, pred_par)
  clear_all()
  gc()
}
pred_par_50 = get_mean_sd(res_all_50)$mean
png(filename = "dsbi_50.png", width = 2400, height = 1800, res = 300)
plot_all_par(test_par, pred_par_50, main_prefix = "DSBI")
dev.off()
compare_hybrid_50 = compare_par(test_par, pred_par_50, c("sigma_Y", "scale_Y", "sigma_U1", "sigma_U2", "scale_U1", "scale_U2"))
compare_hybrid_50[] = lapply(compare_hybrid_50, function(x) if (is.numeric(x)) sprintf("%.3f", x) else x)
print(compare_hybrid_50, row.names = FALSE)




save.image("C:/Users/qihan/Desktop/sim_inhom_S1_S3/sim2_W1_S1_S3.RData")

diff_50 = pred_par_50 - test_par

par_names = c("sigma_Y", "scale_Y", "sigma_U1", "sigma_U2", "scale_U1", "scale_U2")
colnames(diff_50) = par_names

diff_df = data.frame(Method = "50x50 only", diff_50)

diff_long = reshape(
  diff_df,
  varying = par_names,
  v.names = "Error",
  timevar = "Parameter",
  times = par_names,
  direction = "long"
)

rownames(diff_long) = NULL

diff_long$Group = ifelse(grepl("^sigma", diff_long$Parameter), "sigma", "scale")

diff_long$Parameter = factor(
  diff_long$Parameter,
  levels = par_names
)

diff_long$Group = factor(
  diff_long$Group,
  levels = c("sigma", "scale")
)

library(ggplot2)

p_diff_50 = ggplot(diff_long, aes(x = Parameter, y = Error)) +
  geom_boxplot(
    outlier.size = 0.6,
    colour = "black",
    fill = "grey80"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Group, scales = "free") +
  theme_bw() +
  labs(
    title = "Prediction error boxplot (Window = [0, 1] X [0, 1])",
    x = "Parameter",
    y = "Prediction error"
  )


png(filename = "inhom_boxplot.png", width = 2400, height = 1800, res = 300)
print(p_diff_50)
dev.off()






