# This function takes input data in the R format,
# makes the foreign function call to C, and then
# packages the C output in a list.

run_mcgfa = function(X,z,v,bic,known,q,p,G,N,mu,model,class_ind,Lambda,Psi,eta,alpha,eta_max,alpha_min,tol,max_it,iterations) {

    # R uses column-major order => transpose matrices for C
    X = as.vector(t(X))
    z = as.vector(t(z))
    mu = as.vector(t(mu))

    c_out = .C("mcgfa_c", as.double(X), as.double(z), as.double(v), as.double(bic),
               as.integer(known), as.integer(q), as.integer(p), as.integer(G), as.integer(N),
               as.double(mu),as.integer(model), as.integer(class_ind), as.double(Lambda), as.double(Psi),
               as.double(eta), as.double(alpha), as.double(eta_max), as.double(alpha_min), as.double(tol),
               as.integer(max_it), as.integer(iterations),
               PACKAGE="mcgfa")
    list(z = c_out[[2]],
         v = c_out[[3]],
         bic = c_out[[4]],
         mu = c_out[[10]],
         lambda = c_out[[13]],
         psi = c_out[[14]],
         eta = c_out[[15]],
         alpha = c_out[[16]],
         iterations = c_out[[21]]
    )
}
