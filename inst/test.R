
x <- wine[, -1]

mcgfa(x, c(4, 1, 2), 5)
N <- nrow(x)
rG <- 1:3
rq <- 1:5

fit <- mcgfa(x, rG, rq, init_method = 'soft', init_class = inits)
fit_km <- mcgfa(x, rG, rq, init_method = 'kmeans')
fit_pgmm <- mcgfa(x, rG, rq, init_method = 'pgmm')

fit$bic
