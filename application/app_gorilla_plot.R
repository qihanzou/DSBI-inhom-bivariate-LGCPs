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
  aux_input = cbind(beta_mat, Np_all, Np_type1, Np_type2)
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



library(ggplot2)

obs_one = get_features_from_pattern(X_obs, r, cov_list, image_nxy)
obs_data = stack_feature_results(list(obs_one))

pcf_obs = obs_data$pcf_input
img_obs = obs_data$count_img
aux_obs = obs_data$aux_input

imgs = list(
  "Major Type" = img_obs[1, , , 1],
  "Minor Type" = img_obs[1, , , 2],
  "Pooled"     = img_obs[1, , , 3]
)

# Original data-based limits
zlim_range = range(unlist(imgs), na.rm = TRUE)

make_df = function(z, name) {
  df = as.data.frame(as.table(z))
  names(df) = c("x", "y", "value")
  df$x = as.numeric(df$x)
  df$y = as.numeric(df$y)
  df$panel = name
  df
}

plot_df = do.call(rbind, Map(make_df, imgs, names(imgs)))

res_col = colorRampPalette(c(
  "#241255",
  "#005f73",
  "#0a8f6a",
  "#e8e8e8",
  "#7a9a01",
  "#c2a500",
  "#f0a97b",
  "#e8b4d3"
))(256)

p_resid = ggplot(plot_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  facet_wrap(~ panel, nrow = 1) +
  coord_equal(expand = FALSE) +
  scale_fill_gradientn(
    colours = res_col,
    limits = zlim_range,
    name = "Residual"
  ) +
  theme_void() +
  theme(
    strip.text = element_text(
      size = 14,
      face = "bold",
      margin = margin(b = 12)
    ),
    strip.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    panel.spacing = grid::unit(1.0, "lines"),
    plot.margin = margin(12, 12, 12, 12)
  )

print(p_resid)

ggsave(
  "gorilla_residual_count_images.png",
  plot = p_resid,
  width = 9.5,
  height = 3.5,
  dpi = 300
)

