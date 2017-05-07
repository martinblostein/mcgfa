# This function checks to ensure that the user input to init_class
# is of the appropriate form before computation begins.

check_init_class = function(init_class, init_method, X, rG, num_G, original_rG_order) {

    if (is.null(init_class))
        stop("If using manual initialization ('hard' or 'soft'), must provide init_class.")

    if (!is.vector(init_class) && !is.matrix(init_class) && !is.list(init_class))
        stop("init_class must be a list of vectors or matrices.")

    if (!is.list(init_class))
        init_class = list(init_class)

    init_class = init_class[original_rG_order]

    if (!is.list(init_class) && length(init_class) != num_G)
        stop("init_class must be a list with length equal to that of rG.")

    if (init_method == "hard") {
        if (!is.list(init_class)) init_class = list(init_class)
        if (!all(sapply(init_class,is.vector)))
            stop("init_class must be a list of vectors (for init_method = 'hard').")
        if (!all(sapply(init_class,length) == nrow(X)))
            stop(paste("Length of init_class elements must match number of observations",
                       "(for init_method = 'hard')."))
        if (any(sapply(sapply(init_class,unique),length) > rG))
            stop(paste("Cannot manually initialize with more classes than corresponding",
                       "value in in rG."))
        if (min(unlist(sapply(init_class,table))) < 2)
            stop(paste("Initialization vectors must assign at least 2 observations to",
                       "each class."))
    }

    if (init_method == "soft") {
        if (!all(sapply(init_class,is.matrix)))
            stop("init_class must be a list of matrices (for init_method = 'hard').")
        if (!all(sapply(init_class,nrow) == nrow(X)))
            stop(paste("Elements of init_class must have number of rows equal",
                       "to number of observations (for init_method = 'soft')."))
        if (!all(sapply(init_class,ncol) == rG))
            stop(paste("Elements of init_class must have number of columns corresponding",
                       "to the values of rG (for init_method = 'soft')."))
        if (any(unlist(init_class) < 0))
            stop(paste("Elements of init_class must have only positive entries."))
        if (any(abs(sapply(init_class,rowSums)-1) > 1e-5))
            stop(paste("Elements of init_class must have rows that sum to 1.",
                       "(for init_method = 'soft')."))
    }
}
