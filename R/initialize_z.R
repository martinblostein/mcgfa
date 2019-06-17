initialize_z_list <- function(init_method, N, x, rG, init_z, ememargs, known) {
    
    if (init_method == "given") {
        return(init_z)
    } 
    
    lapply(seq_along(rG), function(i) {
        initialize_z(init_method, N, x, rG[[i]], init_z[[i]], ememargs, known)
    })
}

initialize_z <- function(init_method, N, x, G, init_z, ememargs, known) {
    
    if (G == 1) {
        return(matrix(rep(1, N), N, 1))
    }
    
    if (init_method == "emEM") {
        
        ii <- seq_len(ememargs$numstart)
        
        z_candidates <- lapply(ii, function(i) {
            generate_random_z(N, G)
        })
        
        bic <- sapply(ii, function(i) {
            mcgfa_EM(x, G = G, q = ememargs$q, model = ememargs$model, z = z_candidates[[i]], max_it = 5)$bic  
        })
        
        z <- z_candidates[[which.max(bic)]]
        
    } else if (init_method == "kmeans") {
        # Class labels initialized via kmeans clustering
        
        set.seed(123456) # for deterministic output
        z_ind <- kmeans(x, G, nstart = 10)$cluster # ten random starts
        z <- matrix(0, N, G)
        for (i in 1:N) z[i, z_ind[i]] <- 1
        
    } else if (init_method == "given") {
        
        z <- init_z
        
    } else if (init_method == "supervised") {
        # Some observations have known labels;
        # Others assigned equal probability of arising from any group
        
        z <- matrix(0, N, G)
        
        for (i in 1:N){
            if (known[i]) {
                z[i, known[i]] <- 1
            } else {
                z[i, ] <- 1/G
            }
        }
    }
    
    z
    
}

generate_random_z <- function(N, G) {
    z <- matrix(runif(N * G), N, G)
    z / rowSums(z) 
}