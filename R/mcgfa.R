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
    rG=1:3, # number of groups
    rq=1:3, # number of latent factors
    models="all", # the subset of models to be applied
    known=NULL, # known labelled data for semi
    init_method="kmeans", # initialization method for AECM
    init_class=NULL, # intial class labels (only for init_method = 'hard' or 'soft')
    max_it=400, # maximum number of AECM iterations
    tol=1e-3, # tolerance for Aitken criterion
    alpha_min=0.5, # minimum proportion of 'good' points
    eta_max=1000, # maximum contamination factor
    scale=T, # if true, scale data before beginning model fit
    parallel=F, # run code in parallel?
    cores=NULL, # number of parallel processes to run
    silent=F # print info when done?
) {

    #---------------------#
    # CHECK & CLEAN INPUT #
    #---------------------#

    # check X
    tryCatch(X <- as.matrix(X),
             error = function(e) cat(paste("X must be matrix-like.")))
    if(any(is.na(X)))  stop("No NAs allowed.")
    if (!is.numeric(X)) stop("X must be numeric.")
    if (ncol(X) == 1) stop("The model cannot be applied to univariate data.")

    # check rG and rq
    if (!all(rG%%1==0 && rG>0)) stop ("Only positive integer numbers of components are possible (rG).")
    if (!all(rq%%1==0 && rq>0)) stop ("Only postive integer numbers of latent factors are possible (rq).")
    if (any(duplicated(rG))) { rG = unique(rG); print("Removed duplicated group size choices.") }
    if (any(duplicated(rq))) { rG = unique(rG); print("Removed duplicated choices of number of factors.") }
    if (any(duplicated(models))) { models = unique(rG); print("Removed duplicated model choices.") }

    if (max(rq) > ncol(X)/2)
        stop("Cannot fit models with q > p/2.")

    original_rG_order = order(rG)
    rG = sort(rG)
    rq = sort(rq)

    num_G = length(rG)
    num_q = length(rq)

    # check models
    all_models = c("CCC","CCU","CUC","CUU","UCC","UCU","UUC","UUU")
    if (models == "all" || is.null(models)) models = all_models
    if (!all(models %in% all_models))
        stop("Invalid model name(s).")
    num_model = length(models)

    # check known class labels
    if (!is.null(known)) {
        if (length(known) != nrow(X))
            stop("Length of known classes vector must equal the number of observations.")
        known_class_numbers = unique(known[known != 0])
        if (length(known_class_numbers) > min(rG))
            stop("Cannot have more known classes than model with fewest components allows.")
        if (max(known_class_numbers) > length(known_class_numbers))
            stop("Known class numbers must be sequential beginning with 1.")

        init_method = "supervised"
    }

    # check initialization method
    if (!(init_method %in% c("kmeans","pgmm","hard","soft","supervised")))
        stop("Invalid initialization type. Choose one of 'kmeans','pgmm', 'hard', 'soft' or 'supervised'.")

    # check manual initialization
    if (init_method %in% c("hard", "soft"))
        check_init_class(init_class, init_method, X, rG, num_G, original_rG_order)

    # other arguments
    if (max_it < 1) stop ("Must allow at least one AECM iteration (max_it).")
    if (tol < 0) stop ("Stopping criterion tolerance must be positive (tol).")
    if (alpha_min <= 0 || alpha_min >= 1) stop ("alpha_min must lie in (0,1).")
    if (eta_max <= 0) stop ("eta_max must be positive.")

    # parallelization warning
    if (parallel && is.null(cores)) {
        cores = parallel::detectCores()
        cat(paste("Number of parallel cores not specified.", cores,
                  "cores detected automatically. Is this number appropriate",
                  "for your hardware? (y/n)"))
        line = NULL
        repeat {
            line = tolower(readline(">>> "))
            if (line %in% c("n","no")) {
                cat("Exiting, please selected appropriate number of cores using ``cores'' argument, or perform serial computation.")
                return()
            }
            else if (line %in% c("y","yes")) {
                cat("Proceeding with automatically detected number of cores...")
                break()
            }
            else
                cat("(y/n) format required. Please try again.")
        }


    }

    # scaling data by mean & standard deviation if reqested
    # X will hold the centering and scale amounts as attributes
    if (scale) X = scale(X)

    BIC = array(NA,c(num_G,num_q,num_model),dimnames=list(
        paste("# of groups=",rG,sep=""),
        paste("# of factors=",rq,sep=""),
        paste("parsimonious model=",models,sep=""))
    )

    # --------------------------- #
    # CALL MODEL FITTING FUNCTION #
    # --------------------------- #

    parallel_wrapper <- function(GQM) {
        # GQM is a list of the form: {G, q, modelname}
        i = which(GQM[[1]] == rG)
        tryCatch({
            mcgfa_EM(X, rG[GQM[1]], rq[GQM[2]], models[GQM[3]], known, init_method,
                     init_class[[i]], tol, eta_max, alpha_min, max_it)
        }, error = function(e) cat(paste("\n Model not estimated.\n", e)))
    }

    if (parallel) {

        # PARALLEL COMPUTATION #

        GQM = list()
        for (g in 1:num_G) {
            for (q in 1:num_q) {
                for (m in 1:num_model) {
                    GQM[[(g-1)*(num_q*num_model) + (q-1)*num_model + m]] = c(g,q,m)
                }
            }
        }
        #if (is.null(cores)) cores = parallel::detectCores()

        parallel_output = parallel::mclapply(GQM, parallel_wrapper, mc.cores=cores, mc.preschedule=FALSE)

        fits = list()
        ii = 1

        for(i in 1:num_G) {
            fits[[i]] = list()
            for(j in 1:num_q) {
                fits[[i]][[j]] = list()
                for(k in 1:num_model) {
                    fits[[i]][[j]][[k]] = list()

                    # populate this list with the results from above
                    fits[[i]][[j]][[k]] = parallel_output[[ii]]
                    BIC[i,j,k] = parallel_output[[ii]]$bic
                    ii = ii+1
                }
            }
        }
    }

    else {

        # SERIAL COMPUTATION #

        # set up progress bar
        pb_counter = 0
        pb = create_progress_bar(rG,num_G,num_q,num_model,models)

        # For a one-compnent model, all of the models are equivalent to either CCC or CCU:
        CCC_equiv = c("CCC", "CUC", "UCC", "UUC")
        CCU_equiv = c("CCU", "CUU", "UCU", "UUU")

        # Run one from each equivalence class & store them in equiv_fit list.
        if (1 %in% rG) {
            equiv_fit = list(CCC=list(),CCU=list())
            for (equiv_class in list(CCC_equiv,CCU_equiv)) {
                if (any(models %in% equiv_class)) {
                    for (j in 1:num_q) {
                        equiv_fit[[equiv_class[1]]][[j]] = mcgfa_EM(X, 1, rq[j], equiv_class[1], known, init_method,
                                                                    init_class[[1]], tol, eta_max, alpha_min, max_it)
                        pb_counter = pb_counter + 1
                        setTxtProgressBar(pb, pb_counter)
                    }
                }
            }
        }

        # Now fit the rest of the models & populate the nested list of model fits:

        fits = list()
        for (i in 1:num_G) {
            fits[[i]] = list()
            for(j in 1:num_q) {
                fits[[i]][[j]] = list()
                for(k in 1:num_model) {
                    fits[[i]][[j]][[k]] = list()

                    tryCatch({
                        if (rG[i] == 1) {
                            # in one component case, just take model fit that has already been stored
                            if (models[k] %in% CCC_equiv) equiv_class = "CCC"
                            else                          equiv_class = "CCU"
                            fits[[i]][[j]][[k]] = equiv_fit[[equiv_class]][[j]]
                        }
                        else {
                            # otherwise, fit the multi-component model
                            fits[[i]][[j]][[k]] = mcgfa_EM(X, rG[i], rq[j], models[k], known, init_method,
                                                           init_class[[i]], tol, eta_max, alpha_min, max_it)
                            pb_counter = pb_counter + 1
                            setTxtProgressBar(pb, pb_counter)
                        }
                        # store BIC values corresponding to model fits
                        BIC[i,j,k] = fits[[i]][[j]][[k]]$bic
                    }
                    ,error = function(e) cat(paste("\n Model not estimated.\n",e)))
                }
            }
        }
        close(pb)
    }

    # --------------------------- #
    # DETERMINE BEST MODEL BY BIC #
    # --------------------------- #

    best.ind = which(BIC==max(BIC,na.rm=T),arr.ind=T)
    best.ind = best.ind[1,]
    best.fit = fits[[best.ind]]

    # --------------- #
    # DISPLAY RESULTS #
    # --------------- #

    if(!silent) {
        if (num_model==1) {
            cat(paste0("The MCGFA-", best.fit$model, " model with G = ",best.fit$G,
                       " components and q = ",best.fit$q,
                       " latent factors gives a BIC value of ", round(best.fit$bic,3), ".", "\n"))
        }
        else {
            cat("\n")
            cat("  # ------------------------------- #","\n")
            cat("  #  MCGFA Model Selection Results  #","\n")
            cat("  # ------------------------------- #","\n\n")
            cat("  Best BIC value of",best.fit$bic,
                "obtained for the",best.fit$model,
                "model with G =",best.fit$G,
                "components and q =",best.fit$q,
                "latent factors.","\n")
        }
    }

    return(
        structure(
            c(
                list(
                    call = match.call(),
                    X = X,
                    all.bic = BIC),
                best.fit
            ),
            class = "mcgfa"
        )
    )
}
