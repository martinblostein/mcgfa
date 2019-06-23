# mcgfa()
#
# This function:
#
# (1) Performs all of the set up for a set of different
# runs of the AECM algorithm, given the sets {rG, rq, models}.
#
# (2) Calls the mcgfa_EM function for each individual model.
#
# (3) Selects best model by BIC.
#
# (4) Displays results and formats output.
#
# This function is the sole function callable by the user.
# It thus has full documentation available by typing ?mcgfa.
#
# As much type checking as possible takes place in this
# function so that any issues are found before computation
# begins.

mcgfa <- function(
    X, # data
    rG = 1:3, # number of groups
    rq = 1:3, # number of latent factors
    models = "all", # the subset of models to be applied
    known = NULL, # known labelled data for semi
    init_method = "emEM", # initialization method for AECM
    init_z = NULL, # intial z matrix (only for init_method = "given")
    max_it = 400, # maximum number of AECM iterations
    tol = 1e-3, # tolerance for Aitken criterion
    alpha_min = 0.5, # minimum proportion of 'good' points
    eta_max = 1000, # maximum contamination factor
    scale = T, # if true, scale data before beginning model fit
    parallel = F, # run code in parallel?
    cores = NULL, # number of parallel processes to run
    silent = F, # print info when done?
    ememargs = list(numstart = 25, iterations = 5, model = "UUUUU", q = max(rq)) # options for the emEM init method 
) {

    #---------------------#
    # CHECK & CLEAN INPUT #
    #---------------------#

    # check X
    tryCatch(X <- as.matrix(X), error = function(e) stop("X must be matrix-like."))

    if (any(is.na(X)))  {
        stop("No NAs allowed.")
    }
    if (!is.numeric(X)) {
        stop("X must be numeric.")
    }
    if (ncol(X) == 1) {
        stop("The model cannot be applied to univariate data.")
    }

    # check rG and rq
    if (!all((rG %% 1 == 0) & (rG > 0))) {
        stop ("Only positive integer numbers of components are possible (rG).")
    }
    if (!all((rq %% 1 == 0) & (rq > 0))) {
        stop ("Only postive integer numbers of latent factors are possible (rq).")
    }
    if (any(duplicated(rG))) {
        rG <- unique(rG)
        message("Removed duplicated group size choices.")
    }
    if (any(duplicated(rq))) {
        rG <- unique(rG)
        message("Removed duplicated choices of number of factors.")
    }
  
    

    if (with(list(p = ncol(X), q = max(rq)), (p - q) ^ 2 <= (p + q))) {
        stop("Cannot fit models with (p - q)^2 <= (p + q).")
    }

    N <- nrow(X)
    init_z <- init_z[order(rG)]
    rG <- sort(rG)
    rq <- sort(rq)
    
    num_G <- length(rG)
    num_q <- length(rq)
    
    models <- parse_models_argument(models)

    num_model <- length(models)

    # check known class labels
    if (!is.null(known)) {
        if (length(known) != nrow(X)) {
            stop("Length of known classes vector must equal the number of observations.")
        }
        known_class_numbers <- unique(known[known != 0])
        if (length(known_class_numbers) > min(rG)) {
            stop("Cannot have more known classes than model with fewest components allows.")
        }
        if (max(known_class_numbers) > length(known_class_numbers)) {
            stop("Known class numbers must be sequential beginning with 1.")
        }

        init_method <- "supervised"
    }

    # check initialization method
    if (!(init_method %in% c("emEM", "kmeans", "given", "supervised"))) {
        stop("Invalid initialization type. Choose one of 'emEM', 'kmeans', 'given', 'supervised'.")
    }

    # check manual initialization
    
    verify_init_z(init_method, init_z, N, rG)

    # other arguments
    if (max_it < 1) {
        stop ("Must allow at least one AECM iteration (max_it).")
    }
    if (tol < 0) {
        stop ("Stopping criterion tolerance must be positive (tol).")
    }
    if (alpha_min <= 0 || alpha_min >= 1) {
        stop ("alpha_min must lie in (0, 1).")
    }
    if (eta_max <= 0) {
        stop ("eta_max must be positive.")
    }

    # parallelization warning
    if (parallel && is.null(cores)) {
        cores <- parallel::detectCores()
        cat(paste("Number of parallel cores not specified.", cores,
                  "cores detected automatically. Is this number appropriate",
                  "for your hardware? (y/n)"))
        line <- NULL
        repeat {
            line <- tolower(readline(">>> "))
            if (line %in% c("n", "no")) {
                cat("Exiting, please selected appropriate number of cores using ``cores'' argument, or perform serial computation.")
                return()
            } else if (line %in% c("y", "yes")) {
                cat("Proceeding with automatically detected number of cores...")
                break()
            } else {
                cat("(y/n) format required. Please try again.")
            }
        }
    }

    # scaling data by mean & standard deviation if reqested
    # X will hold the centering and scale amounts as attributes
    if (scale) {
        X <- scale(X)
    }

    BIC <- array(NA, c(num_G, num_q, num_model), dimnames = list(
        paste("# of groups=", rG, sep=""),
        paste("# of factors=", rq, sep=""),
        paste("parsimonious model=", models, sep=""))
    )
    
    # --------------------------- #
    # ----- INITIALIZATION ------ #
    # --------------------------- #
    
    z_list <- initialize_z_list(init_method, N, X, rG, init_z, ememargs, known)

    # --------------------------- #
    # CALL MODEL FITTING FUNCTION #
    # --------------------------- #

    if (parallel) {

        # PARALLEL COMPUTATION #
        GQM <- expand.grid(model = models, q = rq, G = rG, equiv = NA_character_, stringsAsFactors = FALSE)
        GQM$i <- 1:nrow(GQM)
       
        equiv_CCCCC <- (GQM$G == 1) & substr(GQM$model, 3, 3) == 'C'
        equiv_CCUCC <- (GQM$G == 1) & substr(GQM$model, 3, 3) == 'U' 
        
        GQM$equiv[equiv_CCCCC] <- 'CCCCC'
        GQM$equiv[equiv_CCUCC] <- 'CCUCC'
        
        GQM_run <- GQM[is.na(GQM$equiv),]
        
        if (nrow(GQM_run) > 0) {
            GQM_run$equiv <- NA_character_
        }
        
        if (any(equiv_CCCCC)) {
            GQM_run <- rbind(GQM_run,
                             expand.grid(model = 'CCCCC', q = rq, G = 1, equiv = 'CCCCC',
                                         i = NA_integer_,  stringsAsFactors = FALSE))
        }
        if (any(equiv_CCUCC)) {
            GQM_run <- rbind(GQM_run,
                             expand.grid(model = 'CCUCC', q = rq, G = 1, equiv = 'CCUCC',
                                         i = NA_integer_, stringsAsFactors = FALSE))
        }
        
        parallel_output <- parallel::mclapply(1:nrow(GQM_run), function(i) {
            row <- GQM_run[i,]
            iG <- which(row$G == rG)
            tryCatch({
                mcgfa_EM(X, row$G, row$q, row$model, z_list[[iG]], known, tol, eta_max, alpha_min, max_it)
            }, error = function(e) cat(paste("\n Model not estimated.\n", e)))
        }, mc.cores = cores, mc.preschedule = FALSE)

        ccccc_runs <- parallel_output[sapply(GQM_run$equiv, identical, y = 'CCCCC')]
        ccucc_runs <- parallel_output[sapply(GQM_run$equiv, identical, y = 'CCUCC')]
        
        expanded_output <- lapply(1:nrow(GQM), function(i) {
            row <- GQM[i,]
            
            switch(row$equiv,
                   CCCCC = ccccc_runs[[match(row$q, rq)]],
                   CCUCC = ccucc_runs[[match(row$q, rq)]],
                   parallel_output[[match(i, GQM_run$i)]])
        })
        
        fits <- list()
        ii <- 1

        for (i in 1:num_G) {
            fits[[i]] <- list()
            for (j in 1:num_q) {
                fits[[i]][[j]] <- list()
                for (k in 1:num_model) {
                    fits[[i]][[j]][[k]] <- list()

                    # populate this list with the results from above
                    fits[[i]][[j]][[k]] <- expanded_output[[ii]]
                    BIC[i, j, k] <- expanded_output[[ii]]$bic
                    ii <- ii + 1
                }
            }
        }
    } else {

        # SERIAL COMPUTATION #
        
        if (!silent) {
            # set up progress bar
            pb_counter <- 0
            pb <- create_progress_bar(rG,num_G,num_q,num_model,models)
        }

        # For a one-component model, all of the models are equivalent to either CCCCC or CCUCC.
        # Run one from each equivalence class & store them in equiv_fit list.
        if (1 %in% rG) {
            equiv_fit <- list(CCCCC = list(), CCUCC = list())
            
            if (any(substr(models, 3, 3) == "C")) {
                for (j in 1:num_q) {
                    equiv_fit$CCCCC[[j]] <- mcgfa_EM(X, 1, rq[j], "CCCCC", z_list[[match(1, rG)]],
                                                     known, tol, eta_max, alpha_min, max_it)
                    if (!silent) {
                        pb_counter <- pb_counter + 1
                        setTxtProgressBar(pb, pb_counter)
                    }
                }
                
            }
            if (any(substr(models, 3, 3) == "U")) {
                for (j in 1:num_q) {
                    equiv_fit$CCUCC[[j]] <- mcgfa_EM(X, 1, rq[j], "CCUCC", z_list[[match(1, rG)]],
                                                     known, tol, eta_max, alpha_min, max_it)
                    if (!silent) {
                        pb_counter <- pb_counter + 1
                        setTxtProgressBar(pb, pb_counter)
                    }
                }
            }
        }

        # Now fit the rest of the models & populate the nested list of model fits:

        fits <- list()
        for (i in 1:num_G) {
            fits[[i]] <- list()
            for (j in 1:num_q) {
                fits[[i]][[j]] <- list()
                for (k in 1:num_model) {
                    fits[[i]][[j]][[k]] <- list()

                    tryCatch({
                        if (rG[i] == 1) {
                            # in one component case, just take model fit that has already been stored
                            equiv_class <- ifelse(substr(models[k], 3, 3) == "C", "CCCCC", "CCUCC")
                            fits[[i]][[j]][[k]] <- equiv_fit[[equiv_class]][[j]]
                        } else {
                            # otherwise, fit the multi-component model
                            fits[[i]][[j]][[k]] <- mcgfa_EM(X, rG[i], rq[j], models[k], z_list[[i]],
                                                            known, tol, eta_max, alpha_min, max_it)
                            if (!silent) {
                                pb_counter <- pb_counter + 1
                                setTxtProgressBar(pb, pb_counter)
                            }
                        }
                        # store BIC values corresponding to model fits
                        BIC[i, j, k] <- fits[[i]][[j]][[k]]$bic
                    }, error = function(e) cat(paste("\n Model not estimated.\n", e)))
                }
            }
        }
        
        if (!silent) {
            close(pb)
        }
    }
    
    successful_fits <- any(!is.na(BIC) & !is.infinite(BIC))

    if (successful_fits) {
        
        # --------------------------- #
        # DETERMINE BEST MODEL BY BIC #
        # --------------------------- #

        best.ind <- which(BIC==max(BIC[!is.infinite(BIC) & !is.na(BIC)]),arr.ind=T)
        best.ind <- best.ind[1,]
        best.fit <- fits[[best.ind]]
        
        # --------------- #
        # DISPLAY RESULTS #
        # --------------- #
        
        if (!silent) {
            if (num_model == 1) {
                cat(paste0("The MCGFA-", best.fit$model, " model with G = ",best.fit$G,
                           " components and q = ",best.fit$q,
                           " latent factors gives a BIC value of ", round(best.fit$bic, 3), ".", "\n"))
            } else {
                cat("\n")
                cat("  # ------------------------------- #", "\n")
                cat("  #  MCGFA Model Selection Results  #", "\n")
                cat("  # ------------------------------- #", "\n\n")
                cat("  Best BIC value of",best.fit$bic,
                    "obtained for the",best.fit$model,
                    "model with G =",best.fit$G,
                    "components and q =",best.fit$q,
                    "latent factors.", "\n")
            }
        }
    } else {
        if (!silent) {
            cat('No successful model fits.\n')
        }
        
        best.fit <- NA
        
    }

    return(
        structure(
            c(
                list(
                    X = X,
                    all.bic = BIC),
                best.fit
            ),
            class = "mcgfa"
        )
    )
}

all_models <- apply(do.call(expand.grid, rep(list(c('C', 'U')), 5)), 1, paste, collapse = '')

parse_models_argument <- function(models) {

  if (is.null(models) || length(models) == 0 || 'all' %in% models) {
    return(all_models)
  }
  
  models <- toupper(models)
  
  invalid_model_strings <- grep('^[UCX]{5}$', models, invert = TRUE, value = TRUE)
  
  if (length(invalid_model_strings) > 0) {
    stop(paste("Invalid model specification(s):", paste(invalid_model_strings, collapse = ', ')))
  }
  
  expanded_models <- lapply(models, expand_model_string)
  expanded_models <- do.call(c, expanded_models)
  
  unique_models <- unique(expanded_models)
  
  unique_models
  
}

expand_model_string <- function(model_string) {
  chars <- strsplit(model_string, '')[[1]]
  
  chars <- lapply(chars, function(chr) {
    if (chr == 'X') {
      c('C', 'U')
    } else {
      chr
    }
  })
  
  apply(expand.grid(chars), 1, paste, collapse = '')
}
