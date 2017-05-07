init_load = function(x, z, G, N, p, q, model) {

    # this function produces initial values for Lambda and Phi given values for z and a model type

    mu = matrix(0,G,p)
    n = colSums(z)
    pi = colSums(z) / N
    for (g in 1:G) mu[g,] = colSums(x*z[,g]) / n[g]
    sampcov = array(0,c(p,p,G))
    for (g in 1:G)
        sampcov[,,g] = cov.wt(x, wt=(z[,g]), center=mu[g,])$cov * (sum(z[,g]^2)-1)/sum(z[,g])

    if (substr(model,1,1) == "U") {
        lambda = array(0,c(p,q,G))
        for(g in 1:G) {
            eigvec = eigen(sampcov[,,g])$vector
            eigval = eigen(sampcov[,,g])$values
            for (i in 1:p)
                for (j in 1:q)
                    lambda[i,j,g] = sqrt(eigval[j]) * eigvec[i]
        }

        if (model %in% c("UCU","UCC")) {
            psi = c(rep(0,p))
            for(g in 1:G){
                psi = psi + pi[g] * abs(diag(sampcov[,,g] - lambda[,,g] %*% t(lambda[,,g])))
            }

            if (model=="UCC") psi = sum(psi)/p
        }

        if (model %in% c("UUU","UUC")) {
            psi_tmp = matrix(0,G,p)
            for (g in 1:G){
                psi_tmp[g,]<-abs(diag(sampcov[,,g] - lambda[,,g] %*% t(lambda[,,g])))
            }
            if (model=="UUU") psi = as.vector(t(psi_tmp))
            if (model=="UUC") psi = rowMeans(psi_tmp)
        }
    } else if (substr(model,1,1) == "C") {
        stilde = matrix(0,p,p)
        for(g in 1:G){
            stilde = stilde + pi[g] * sampcov[,,g];
        }
        eigvec = eigen(stilde)$vector
        eigval = eigen(stilde)$values
        lambda = matrix(0,p,q)

        for (i in 1:p)
            for (j in 1:q)
                lambda[i,j] = sqrt(eigval[j]) * eigvec[i]

        if (model == "CCU") psi = abs(diag(stilde - lambda %*% t(lambda)))
        if (model == "CCC") psi = sum(abs(diag(stilde - lambda %*% t(lambda))))/p

        if (model %in% c("CUC","CUU")) {
            psi_tmp = matrix(0,G,p)
            for (g in 1:G){
                psi_tmp[g,] = abs(diag(sampcov[,,g] - lambda%*%t(lambda)))
            }

            if (model == "CUC") psi = rowMeans(psi_tmp)
            if (model == "CUU") psi = as.vector(t(psi_tmp))
        }
    }

    # now store lambda in row major order
    if (is.na(dim(lambda)[3]))
        lambda = t(lambda)
    else
        for (i in 1:dim(lambda)[3])
            lambda[,,i] = t(lambda[,,i])

    list(lambda=lambda,psi=psi)

}
