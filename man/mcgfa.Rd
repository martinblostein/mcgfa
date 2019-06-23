\name{mcgfa}
\alias{mcgfa}
\title{Model Fitting for Mixtures of Contaminated Gaussian Factor Analyzers}

\description{Performs clustering and classification using the Mixtures of Contaminated Gaussian Factor Analyzers model. This model allows for automatic detection of outliers and noise, and automatically detects outliers. It is appropriate for high-dimensional numerical data.}

\usage{
mcgfa(X, rG=1:3, rq=1:3, models="all", known=NULL, init_method="emEM",
      init_z, max_it=400, tol=1e-3, alpha_min=0.5, eta_max=1000,
      scale=T, parallel=F, cores=NULL, silent=F, ememargs = 
      list(numstart = 25, iterations = 5, model = "UUUUU", q = max(rq)))
}

\arguments{
    \item{X}{The data matrix. Rows correspond to observations and columns to variables/features; thus \eqn{X} is an \eqn{N}-by-\eqn{p} matrix. \eqn{X} must be numerical.}
    \item{rG}{The set of values used for the number of components.}
    \item{rq}{The set of values used for the number of latent factors.}
    \item{models}{The set of parsimonious models used. Either the string "all", or a character vector specifying a subset of models. See Details.}
    \item{known}{If NULL, clustering or ``unsupervised learning'' is performed. If a vector of length equal to the number of observations, semi-supervised learning is performed. In this case, the \code{i}-th entry of the argument \code{class} is either \eqn{0}, or some number in \eqn{{1, 2, \ldots , G}}, where \eqn{G} is the number of components. A value of \eqn{0} indicates that the label for observation \code{i} is unknown, while a nonzero value gives the known label. Note that if all values are nonzero, the model parameters are simply fit and fully supervised classification of new observations may be performed using the \code{predict} method of the \code{mcgfa} class.}
    \item{init_method}{Determines how starting values for the AECM algorithm are generated.  Valid options are "emEM", "kmeans", "given" or "supervised". If \code{known} parameter is provided, this parameter is automatically set to "supervised".
    
        \itemize{
            \item{\code{emEM}: initialization determined using the emEM method, where initialization candidates are randomly generated, and the AECM algorithm is applied for a small number of iterations. Whichever random initialization produces the best model in the short AECM runs is selected as the initialization for the full AECM runs.}
            \item{\code{kmeans}: begin by clustering observations using built-in \code{kmeans} R function.}
            \item{\code{given}: initialization manually specified in \code{init_z} argument}
            \item{\code{supervised}: Observations with known class labels are assigned to their known component with probability 1, while unlabelled components are initialized with equal prior probability of membership for each component}
        }
        See Details for more information.
    }
    \item{init_z}{Used when \code{init_method} is set to \code{given} to manually prescribe the initialization used by the AECM algorithm. A list of row-stochastic matrices of same length as \code{rG}, where \eqn{z[i,g]} represents the prior probability that observation \eqn{i} is a member of group \eqn{g}.}
    \item{max_it}{The maximum number of iterations for the AECM algorithm.}
    \item{tol}{The tolerance for the Aitken acceleration stopping criterion.}
    \item{alpha_min}{The minimum allowable value for alpha, which represents the proportion of ``good'' points in a given group.}
    \item{eta_max}{The maximum allowable value for eta, which is the covariance inflation factor for outlying points.}
    \item{scale}{If TRUE, the data is scaled before the algorithm is begun. Recommended.}
    \item{parallel}{If TRUE, computation takes place in parallel on several processors.}
    \item{cores}{Only relevent if \code{parallel=TRUE}. Determines the number of cores used in parallel computation. If left undefined, number of available cores is determined automatically.}
    \item{silent}{If TRUE, function will not print any output at completion.}
    \item{ememargs}{A list used to set options for the emEM initialization method:
      \itemize{
        \item{\code{numstart}: The number of random starting values}
        \item{\code{iterations}: The number of AECM iterations applied to each starting value}
        \item{\code{model}: The covariance model used in the emEM iterations}
        \item{\code{q}: The number of latent factors used in the emEM iterations}
      }
    }
}

\details{

    This function implements the Mixtures of Contaminated Gaussian Factor Analyzers (MCGFA) model, for model-based clustering and classification. The approach is meant to be applied on on high-dimensional, and noisy, data. A description of the model can be found in Punzo, Blostein and McNicholas (2017). Parameter estimation is performed using an Alternating Expectation-Conditional Maximization (AECM) algorithm (Meng and Van Dyk, 1997).

    Besides clustering into components, this algorithm also automatically detects outliers through the use of the contaminated Gaussian distribution. So, in the end each observation is classified in a nested fashion: first into components, then into as inliers/outliers.

    To specify an individual MCGFA model, one must determine: the number of components (\eqn{G}), the number of latent factors (\eqn{q}), and the model name. The model name is one of the following thirty-two options:
    \itemize{
    \item{CCCCC, CCUCC, CUCCC, CUUCC, UCCCC, UCUCC, UUCCC, UUUCC \cr
CCCCU, CCUCU, CUCCU, CUUCU, UCCCU, UCUCU, UUCCU, UUUCU \cr
CCCUC, CCUUC, CUCUC, CUUUC, UCCUC, UCUUC, UUCUC, UUUUC \cr
CCCUU, CCUUU, CUCUU, CUUUU, UCCUU, UCUUU, UUCUU, UUUUU}
    }
    
    Explanation of the model naming scheme can be found in the next subsection, as well as instructions on how to conveniently generate of subset of the full set of models.

    The user may choose to fit a single model. However, usually many models are fit to the same data, and a model selection criterion is used to determine the best one. By multiple values for the \code{rG}, \code{rq} and \code{models} parameters, many models can be fit at once. The \code{mcgfa} function then selects the best model according to the Bayesian Information Criterion (BIC).

    When fitting many models to large data, parallel computation may be employed by \code{mcgfa}, using the \code{parallel} parameter. This parallelization is provided by the \code{mcapply} function from the \code{parallel} package for R.

    \subsection{Model Names & Constraints}{
    The model name indicates which constraints are to be imposed on to the covariance structure of the factor analysis model, and as well as on to the parameters \eqn{\eta} and \eqn{\alpha}.
    
    Because the MCGFA is a mixture of factor analyzers model, the covariance matrix of the \eqn{g}-th group, \eqn{\Sigma_g}, can be decomposed as follows:
    \deqn{\Sigma_g = \Lambda_g\Lambda_g' + \Psi_g}
    where \eqn{\Lambda_g} is a \eqn{p} by \eqn{q} factor loading matrix, and \eqn{\Psi_g} is a diagonal matrix that determines the noise variance for each of the \eqn{p} variables. The family of eight models is formed by introducing different sets of constraints on \eqn{\Lambda} and \eqn{\Psi}.

    The five-letter model names are interpreted as follows. ``C'' indicates that the constraint is imposed, and ``U'' that it is not.
    
    \enumerate{
      \item{\eqn{\Lambda} constrained to be equal across groups (\eqn{\Lambda_g = \Lambda})}
      \item{\eqn{\Psi} constrained to be equal across groups (\eqn{\Psi_g = \Psi})}
      \item{Error variances constrained to be to be equal \emph{within} groups (\eqn{\Psi_g = I\psi})}
      \item{\eqn{\alpha} constrained to be equal across groups (\eqn{\alpha_g = \alpha})}
      \item{\eqn{\eta} constrained to be equal across groups (\eqn{\eta_g = \eta})}
    }
    
    The subset of models to fit is specified in the \code{models} argument as a character vector. For convenience, the user may also use the character 'X' in a model name to indicate a wildcard. For example, \code{'UUUXX'} is equivalent to the set \code{c('UUUCC','UUUCU','UUUUC','UUUUU')}. These wildcards may be combined, so for example \code{c('UUUUU','CCCCX')} is equivalent to \code{c('UUUUU','CCCCU','CCCCC')}. Any duplicate models generated by wildcards will be removed.
}

    \subsection{Initialization Methods}{
        Because the AECM algorithm cannot guarantee convergence to a global maximum, success of the algorithm in reaching a good fit depends heavily on starting values. Starting values in this case refer to the prior probability of class membership for each observation. This is represented in an \eqn{n}-by-\eqn{G} stochastic matrix \eqn{z}, where \eqn{n} is the number of observations and G is the number of components. The element \eqn{z[i,j]} indicates the prior probability that observation \eqn{i} lies in group \eqn{G}. Several different options are provided for the generation of this initialization matrix, through the \code{init_method} argument.
        
        The default initialization method is \code{emEM}. In this case, several candidates for initial classification matrices are generated, and then the AECM algorithm is applied to each for a small number of iterations. Whichever candidate achieves the best BIC value is selected as the initialization for the full AECM runs. This process occurs separately for each value of \eqn{G}, but only using one parsimonious model and one value of \eqn{q}. The number of candidates, number of iterations, parsimonious model and \eqn{q} used for the emEM initialization can be specified in the \code{ememargs} argument. Options for the emEM initialization can be provided in the \code{ememargs} argument.
        
        The second option is \code{kmeans}, which uses the \eqn{k}-means clustering method to generate a ``hard'' initial classification. \eqn{k}-means is a fast, simple clustering algorithm that returns a hard classification.
        
        Finally, using the \code{given} option, the user may provide a specific initialization matrix for each value of \code{G} in \code{rG}. The initialization is provided in the \code{init_z} argument. That argument must be a list of the same length as \code{rG}, where each element is a matrix with N rows and number of columns corresponding the matching value of \code{rG}. The \eqn{[i,j]-th} element of each matrix indicates the initial probability that observation \eqn{i} lies in class \eqn{j}.
        
    }

    \subsection{Classification}{
    This function provides two methods for model-based classification. The first is usual fully-supervised classification. For this method, the model is fit to fully labelled data (by providing a \code{known} vector argument with no zeros). Then, after model fitting takes place, the \code{predict} generic function can be applied to the \code{mcgfa} object returned by this function, to predict the class labels of a matrix of unlabelled observations. This method has the advantage that as new observations arise, predictions can be made very quickly, without refitting the model.

    The other form of classification is a form of semi-supervised learning. Semi-supervised classification makes use of both unlabelled and labelled data for training. The information that the unlabelled data provides about the structure of the dataset can improve classification. For this method, the \code{known} vector argument is simply provided with elements equal to zero, indicated that the corresponding observation has no known label. The MCGFA model is then fit and the unknown elements are classified as usual.
    The second method requires that the model be fitted again each time new data arrives. However, it is also completely possible to apply the \code{predict} method a model fit using semi-supervised classification.


}

}

\value{
    \item{X}{The data the model was fit to. If the data was scaled, then this matrix will be as well, with the centering and scale factors applied to it stored as attributes.}
    \item{all.bic}{An array of all of the BIC values for every model fitted, indexed by model name, number of groups and number of latent factors.}
    \item{model}{The name of the best model: determines the constraints on the covariance structure. See Details for more information.}
    \item{G}{The number of groups in the best model.}
    \item{q}{The number of latent factors for the best model.}
    \item{z}{The ``soft'' clustering matrix. The element in the \code{i}-th row and \code{g}-th column gives the posterior probability that observation \code{i} is in group \code{g}.}
    \item{group}{The ``hard'' classifications into groups. A vector of integers in the range \code{rG}, giving the a posteriori classification of each observation for the best model.}
    \item{isBad}{Similar, to \code{group}, the maximum a posteriori classification of each observation as ``good'' or ``bad''. If the \code{i}-th value is 1, the \code{i}-th observation is labelled as an outlier or noise in the best model.}
    \item{isBad.soft}{The ``soft'' classification vector of each observation as ``good'' or ``bad''. The closer the \code{i}-th value is to 1, the more likely the \code{i}-th observation is an outlier or noise.}
    \item{mu}{The matrix of group means (each column corresponds to a group), for the best model.}
    \item{alpha}{The proportion of ``good'' points for each group, for the best model.}
    \item{eta}{The contamination inflation factors for each group, for the best model.}
    \item{lambda}{The factor loading matrices of each group, for the best model.}
    \item{psi}{The error variance matrices of each group, for the best model.}
    \item{sigma}{The covariance matrices of each group, for the best model.}
    \item{npar}{The total number of free parameters in the best model.}
    \item{iterations}{The number of iterations of the AECM algorithm performed, for the fitting of the best model.}
    \item{init_z}{The initial z matrix used for initialization of the AECM algorithm, for each value of G in rG.}
}

\references{
Meng, X.-L. and Van Dyk, D. (1997). The EM Algorithm - An old folksong sung to a fast new tune. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{59(3)}, 511-567.
Punzo, A., Blostein, M. and McNicholas, P. D. (2017). High-dimensional clustering with the contaminated Gaussian distribution \emph{arXiv preprint arXiv:1408.2128v2.}
}

\examples{

# TOY EXAMPLE: 2 GAUSSIAN COMPONENTS WITH NOISE
# small example for easy visualization
data(gaussnoise)
toy_classes <- gaussnoise[,1]
is_noise <- gaussnoise[,2]
X <- gaussnoise[,-c(1,2)]

set.seed(1)
toy_fit <- mcgfa(X, rG=1:3, rq=1, models="XXXCC")
# plot: outliers => triangle
plot(toy_fit)
# check clustering performance: 2 errors
table(toy_classes[1:175],toy_fit$group[1:175])
# check noise detection performance: 2 false positive, 3 false negative
table(is_noise,toy_fit$isBad)

\donttest{
# REAL DATA EXAMPLE: WINE DATA
data(wine)
# simulate clustering by completely hiding known classes from algorithm
X <- wine[,-1]
wine_classes <- wine[,1]
set.seed(1)
wine_fit <- mcgfa(X, rG=1:3, rq=1:4)

# check performance:
#  - correctly selected 3 groups
#  - 2 errors
table(wine_classes,wine_fit$group)

# CLASSIFICATION EXAMPLE: OLIVE DATA

# load data
data(olive)
X <- olive[,-c(1,2)]
# classes correspond to regions of Italy
regions <- olive[,1]

# take eighth of observations to form training data
set.seed(1)
train_ind <- sample.int(nrow(X),nrow(X)/8)

known <- rep(0,nrow(X))
known[train_ind] <- regions[train_ind]

# FULL SUPERVISION (only labelled data available to classifier)
# In this case, the model is formed using only the 71 labelled
# observations. Then this model is used to predict the class
# membership of the remaining 501 observations.
fit <- mcgfa(X[train_ind,], rq = 3:4, rG = 3)
pred <- predict(fit,X[-train_ind,])
# check classification performance: poor separation
table(regions[-train_ind],pred$hard)

# PARTIAL SUPERVISION (unlabelled data available to classifier)
# In this case, the model is formed using the 71 labelled 
# observations, as well the the 501 unlabelled observations.
# As a part of this initial model fit, class predictions for the
# unlabelled observations are automatically generated.
fit2 <- mcgfa(X,rq=3:4,rG=3,known=known)
# check classification performance: extremely good, only 1 error
table(regions[-train_ind],fit2$group[-train_ind])
# This shows that the use of the  unlabelled observations during
# model fitting can significantly improve classification performance.
}
}

\keyword{multivariate}
\keyword{cluster}
\keyword{classif}