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
    known=NULL, # known classes vector
    init_method="kmeans", # initialization method
    init_class, # intial class labels (only for init_method = 'hard' or 'soft')
    tol=1.0e-03, # stopping rule in the Aitken rule
    eta_max=1000, # maximum contamination factor
    alpha_min=0.5, # minimum proportion of 'good' points / inliers
    max_it=400 # maximum number of EM iterations
)
{
    #--------#
    # SET UP #
    #--------#

    model_list = c("CCC","CCU","CUC","CUU","UCC","UCU","UUC","UUU")
    N = nrow(X)    # sample size
    p = ncol(X)    # number of variables
    iterations = 0

    class_ind = if (!is.null(known) && (!all(known==0))) 1 else 0

    # parameter initialization/allocation
    v = array(1,c(N,G))
    mu = array(0,c(p,G))
    alpha = numeric(G)
    eta = rep(1.01,G)

    bic = 0;
    model_number = which(model_list==model)

    #------------------#
    # Z INITIALIZATION #
    #------------------#

    if (init_method == "pgmm") {
        # Class labels initialized using the PGMM package.

        dummy = capture.output(
            suppressWarnings(McNi <- pgmm::pgmmEM(x=X,zstart=2,rq=q,rG=G,modelSubset=model,class=known))
        )

        z = McNi$zhat

        # Now process Lambda and Psi from the pgmm format:

        # LAMBDA
        if(substr(model,1,1)=="C") {
            # lambda constrained to be identical across groups
            Lambda = t(McNi$load)
        }
        if(substr(model,1,1)=="U") {
            # lambda free to vary across groups
            Lambda = array(0,c(p,q,G))
            for(j in 1:G) Lambda[,,j] <- t(McNi$load[[j]])
        }

        # PSI
        if(model %in% c("CCC", "UCC")) {
            # psi constrained & isotropic
            Psi = McNi$noisev
        }
        if(model %in% c("CCU", "UCU")) {
            # psi constrained & non-isotropic
            Psi = McNi$noisev
        }
        if(model %in% c("CUU","UUU")) {
            # psi unconstrained & non-isotropic
            Psi = array(0,c(p,G))
            for(j in 1:G) Psi[,j] <- diag(McNi$noisev[[j]])
        }
        if(model %in% c("CUC","UUC")) {
            # psi unconstrained & isotropic
            Psi = McNi$noisev
        }

    } else if (init_method == "kmeans") {
        # Class labels initialized via kmeans clustering

        set.seed(123456) # for deterministic output
        z_ind = kmeans(X,G,nstart=10)$cluster # ten random starts
        z = matrix(0,N,G)
        for (i in 1:N) z[i,z_ind[i]] = 1

    } else if (init_method=="hard") {
        # Class labels initialized

        z_ind = init_class
        z = matrix(0,N,G)
        for (i in 1:N) z[i,init_class[i]] = 1

    } else if (init_method=="soft") {

        z = init_class

    } else if (init_method=="supervised") {
        # Some observations have known labels;
        # Others assigned equal probability of arising from any group

        z = matrix(0,N,G)

        for (i in 1:N){
            if (known[i])
                z[i,known[i]] = 1
            else z[i,] = 1/G
        }
    }

    # Except when initializing with pgmm, the covariance parameters
    # Lambda and Psi must now be estimated from the initial z.

    if (init_method != "pgmm") {
        inits = init_load(X, z, G, N, p, q, model)
        Lambda = inits$lambda
        Psi = inits$psi
    }

    #---------------------#
    # CALL AECM ALGORITHM #
    # --------------------#

    EM.out = run_mcgfa(X,z,v,bic,known,q,p,G,N,mu,model_number,class_ind,Lambda,Psi,eta,alpha,eta_max,alpha_min,tol,max_it,iterations)

    # ---------------#
    # PREPARE OUTPUT #
    # ---------------#

    z = matrix(EM.out$z,N,G,byrow=T)
    v = matrix(EM.out$v,N,G,byrow=T)
    mu = matrix(EM.out$mu,p,G)
    group = apply(z,1,which.max)
    isBad.soft = numeric(N)
    for(i in 1:N) isBad.soft[i] = v[i,group[i]]
    isBad = ifelse(isBad.soft<0.5,1,0)

    # reshape psi into array of matrices
    psi = array(0,c(p,p,G))
    if (model %in% c("CCC","UCC")) psi_diag = rep(EM.out$psi,p*p*G)
    if (model %in% c("CCU","UCU")) psi_diag = rep(EM.out$psi,G)
    if (model %in% c("CUC","UUC")) psi_diag = rep(EM.out$psi,each=p)
    if (model %in% c("CUU","UUU")) psi_diag = EM.out$psi
    for (g in 0:(G-1)) psi[,,g+1] = diag(psi_diag[(g*p+1):(g*p+p)])

    # reshape lambda into array of matrices
    lambda = array(0,c(p,q,G))
    if (substr(model,1,1)=="C") lambda_temp = rep(EM.out$lambda,G)
    if (substr(model,1,1)=="U") lambda_temp = EM.out$lambda

    for (g in 0:(G-1)) lambda[,,g+1] = matrix(lambda_temp[(g*p*q+1):(g*p*q+p*q)],p,q,byrow=T)

    # covariance matrix: sigma = lambda * lambda' + psi
    sigma = array(0,c(p,p,G))
    for (g in 1:G) sigma[,,g] = lambda[,,g] %*% t(lambda[,,g]) + psi[,,g]

    # number of model parameters:
    npar = switch(model,
                  CCC = (G-1) + p*G + ((p*q - q*(q-1)/2) + 1) + 2*G,
                  CCU = (G-1) + p*G + ((p*q - q*(q-1)/2) + p) + 2*G,
                  CUC = (G-1) + p*G + ((p*q - q*(q-1)/2) + G) + 2*G,
                  CUU = (G-1) + p*G + ((p*q - q*(q-1)/2) + G*p) + 2*G,
                  UCC = (G-1) + p*G + (G*(p*q - q*(q-1)/2) + 1) + 2*G,
                  UCU = (G-1) + p*G + (G*(p*q - q*(q-1)/2) + p) + 2*G,
                  UUC = (G-1) + p*G + (G*(p*q - q*(q-1)/2) + G) + 2*G,
                  UUU = (G-1) + p*G + (G*(p*q - q*(q-1)/2) + G*p) + 2*G)

    # final list
    result = list(
        model=model,
        G=G,
        q=q,
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
        iterations = EM.out$iterations
    )

    return(result)

}
