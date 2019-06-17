x <- wine[, -1]

rG <- 1:3
rq <- 1:5

fit <- mcgfa(x, rG, rq, init_method = 'emEM', parallel = TRUE, cores = 4)
fit_km <- mcgfa(x, rG, rq, init_method = 'kmeans', parallel = TRUE, cores = 4)

fit$bic
fit_km$bic