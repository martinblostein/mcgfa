# This function checks to ensure that the user input to init_z
# is of the appropriate form before computation begins.

verify_init_z <- function(init_method, init_z, N, rG) {
    
    if (init_method != 'given') {
        if (!is.null(init_z)) {
            warning("init_z argument ignored unless init_method = 'given'.")
        }
        return()
    }
    
    if (is.null(init_z)) {
        stop("If using 'given' initialization method, must provide init_z.")
    }
    if (!is_list_of_matrices(init_z)) {
        stop("init_z must be a list of matrices.")
    }
    if (length(init_z) != length(rG)) {
        stop('init_z must be a list of length equal to the length of rG')
    }
    if (!all(sapply(init_z,nrow) == N)) {
        stop(paste("Elements of init_z must have number of rows equal",
                   "to number of observations (for init_method = 'soft')."))
    }
    if (!all(sapply(init_z,ncol) == rG)) {
        stop(paste("Elements of init_z must have number of columns corresponding",
                   "to the values of rG (for init_method = 'soft')."))
    }
    if (any(unlist(init_z) < 0)) {
        stop(paste("Elements of init_z must have only positive entries."))
    }
    if (any(abs(sapply(init_z,rowSums)-1) > 1e-5)) {
        stop(paste("Elements of init_z must have rows that sum to 1.",
                   "(for init_method = 'soft')."))
    }
}

is_list_of_matrices <- function(x, G) {
    is.list(x) && all(sapply(x, is.matrix))
}

