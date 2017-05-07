predict.mcgfa = function(object, y, ...) {
    fit = object
    if (ncol(y) != ncol(fit$X)) y = t(y)
    if (ncol(y) != ncol(fit$X)) stop("Input vector has incorrect number of dimensions.")

    m = nrow(y)

    if (!is.null(attr(fit$X,"scaled:scale")))
        y = scale(y, center = attr(fit$X,"scaled:center"), scale = attr(fit$X,"scaled:scale"))

    dens = matrix(0,m,fit$G)
    for (g in 1:fit$G)
        dens[,g] = cn_dens(y, fit$mu[,g], fit$sigma[,,g], fit$alpha[g], fit$eta[g])
    soft = dens/matrix(rowSums(dens),m,fit$G)
    hard = apply(soft,1,which.max)
    list(soft = soft, hard = hard)
}

plot.mcgfa = function(x, vars, hideBadPoints=F, ...) {
    fit = x
    if (hideBadPoints)
        pairs(fit$X[fit$isBad==0,vars], col=fit$group)
    else
        pairs(fit$X[,vars], col=fit$group, pch = fit$isBad+1)
}

cn_dens = function(x,mu,sigma,alpha,eta) {
    p = length(mu)
    good_dens = (2*pi)^(-p/2) * det(sigma)^(-1/2) * exp(-0.5*mahalanobis(x,mu,sigma))
    bad_dens = (2*pi)^(-p/2) * det(sigma*eta)^(-1/2) * exp(-0.5*mahalanobis(x,mu,sigma*eta))
    alpha*good_dens + (1-alpha)*bad_dens
}

