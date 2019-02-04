# This function:
#
# (1) Performs all of the set up for one run
# of the AECM algorithm given one triple {G, q, model}.
#
# (2) Calls the c wrapper function.
#
# (3) Packages all of the c output in a list.
#
# It is hidden from the user.

mcgfa_EM <- function(
    X, # matrix of data
    G, # number of groups
    q, # number of latent factors
    model, # covariance structure of model to fit
    z, # initial z matrix
    known = NULL, # known classes vector
    tol = 1.0e-03, # stopping rule in the Aitken rule
    eta_max = 1000, # maximum contamination factor
    alpha_min = 0.5, # minimum proportion of 'good' points / inliers
    max_it = 400 # maximum number of EM iterations
) {
    #--------#
    # SET UP #
    #--------#

    cov_model <- substr(model, 1, 3)
    fixed_alpha <- ifelse(substr(model, 4, 4) == 'C', 1L, 0L)
    fixed_eta <- ifelse(substr(model, 5, 5) == 'C', 1L, 0L)
    
    cov_model_list <- c("CCC", "CCU", "CUC", "CUU", "UCC", "UCU", "UUC", "UUU")
    N <- nrow(X)
    p <- ncol(X) 
    iterations <- 0

    class_ind <- if (!is.null(known) && (!all(known == 0))) 1 else 0

    # parameter initialization/allocation
    v <- array(0.99, c(N, G))
    mu <- array(0, c(p, G))
    alpha <- numeric(G)
    eta <- rep(1.01, G)

    bic <- 0
    cov_model_number <- match(cov_model, cov_model_list)

    # Initialize lambda & psi
    inits <- init_load(X, z, G, N, p, q, cov_model)
    Lambda <- inits$lambda
    Psi <- inits$psi

    #---------------------#
    # CALL AECM ALGORITHM #
    # --------------------#

    EM.out <- run_mcgfa(X, z, v, bic, known, q, p, G, N, mu, cov_model_number, fixed_alpha, fixed_eta,
                        class_ind, Lambda, Psi, eta, alpha, eta_max, alpha_min, tol, max_it, iterations)

    # ---------------#
    # PREPARE OUTPUT #
    # ---------------#

    init_z <- z
    z <- matrix(EM.out$z, N, G, byrow = TRUE)
    v <- matrix(EM.out$v, N, G, byrow = TRUE)
    mu <- matrix(EM.out$mu, p, G)
    group <- apply(z, 1, which.max)
    isBad.soft <- numeric(N)
    for(i in 1:N) isBad.soft[i] <- v[i, group[i]]
    isBad <- ifelse(isBad.soft < 0.5, 1, 0)

    # reshape psi into array of matrices
    psi <- array(0, c(p, p, G))
    if (cov_model %in% c("CCC", "UCC")) { psi_diag <- rep(EM.out$psi, p*p*G) }
    if (cov_model %in% c("CCU", "UCU")) { psi_diag <- rep(EM.out$psi, G) }
    if (cov_model %in% c("CUC", "UUC")) { psi_diag <- rep(EM.out$psi, each = p) }
    if (cov_model %in% c("CUU", "UUU")) { psi_diag <- EM.out$psi }
    for (g in 0:(G-1)) psi[,,g+1] <- diag(psi_diag[(g*p+1):(g*p+p)])

    # reshape lambda into array of matrices
    lambda <- array(0, c(p, q, G))
    if (substr(cov_model, 1, 1) == "C") { lambda_temp <- rep(EM.out$lambda, G) }
    if (substr(cov_model, 1, 1) == "U") { lambda_temp <- EM.out$lambda }

    for (g in 0:(G-1)) lambda[,,g+1] <- matrix(lambda_temp[(g*p*q+1):(g*p*q+p*q)], p, q, byrow = TRUE)

    # covariance matrix: sigma <- lambda * lambda' + psi
    sigma <- array(0, c(p, p, G))
    for (g in 1:G) sigma[,,g] <- lambda[,,g] %*% t(lambda[,,g]) + psi[,,g]

    # number of model parameters:
    npar <- switch(cov_model,
                  CCC = (G-1) + p*G + ((p*q - q*(q-1)/2) + 1),
                  CCU = (G-1) + p*G + ((p*q - q*(q-1)/2) + p),
                  CUC = (G-1) + p*G + ((p*q - q*(q-1)/2) + G),
                  CUU = (G-1) + p*G + ((p*q - q*(q-1)/2) + G*p),
                  UCC = (G-1) + p*G + (G*(p*q - q*(q-1)/2) + 1),
                  UCU = (G-1) + p*G + (G*(p*q - q*(q-1)/2) + p),
                  UUC = (G-1) + p*G + (G*(p*q - q*(q-1)/2) + G),
                  UUU = (G-1) + p*G + (G*(p*q - q*(q-1)/2) + G*p))
    
    npar <- npar + ifelse(fixed_alpha, 1, G) + ifelse(fixed_eta, 1, G)

    # final list
    result <- list(
        model = model,
        G = G,
        q = q,
        z = z,
        group = group,
        bic = EM.out$bic,
        isBad = isBad,
        isBad.soft = isBad.soft,
        mu = mu,
        alpha = EM.out$alpha,
        eta = EM.out$eta,
        lambda = lambda,
        psi = psi,
        sigma = sigma,
        npar = npar,
        iterations = EM.out$iterations,
        init_z = init_z
    )

    return(result)

}
