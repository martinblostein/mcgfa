\name{mcgfa}
\alias{mcgfa}
\title{Model Fitting for Mixtures of Contaminated Gaussian Factor Analyzers}

\description{Performs clustering and classification using the Mixtures of Contaminated Gaussian Factor Analyzers model. This model allows for automatic detection of outliers and noise, and automatically detects outliers. It is appropriate for high-dimensional numerical data.}

\usage{
mcgfa(X, rG=1:3, rq=1:3, models="all", known=NULL, init_method="kmeans",
      init_class=NULL, max_it=400, tol=1e-3, alpha_min=0.5, eta_max=1000,
      scale=T, parallel=F, cores=NULL, silent=F)
}

\arguments{
    \item{X}{The data matrix. Rows correspond to observations and columns to variables/features; thus \eqn{X} is an \eqn{N}-by-\eqn{p} matrix. \eqn{X} must be numerical.}
    \item{rG}{The set of values used for the number of components.}
    \item{rq}{The set of values used for the number of latent factors.}
    \item{models}{The set of models used. Some subset of \{CCC, CCU, CUC, CUU, UCC, UCU, UUC, UUU\}.}
    \item{known}{If NULL, clustering or ``unsupervised learning'' is performed. If a vector of length equal to the number of observations, semi-supervised learning is performed. In this case, the \code{i}-th entry of the argument \code{class} is either \eqn{0}, or some number in \eqn{{1, 2, \ldots , G}}, where \eqn{G} is the number of components. A value of \eqn{0} indicates that the label for observation \code{i} is unknown, while a nonzero value gives the known label. Note that if all values are nonzero, the model parameters are simply fit and fully supervised classification of new observations may be performed using the \code{predict} method of the \code{mcgfa} class.}
    \item{init_method}{Selects how starting values for the AECM algorithm are determined. Valid options are "kmeans", "pgmm", "hard", and "soft".

        \itemize{
            \item{\code{kmeans}: begin by clustering observations using built-in \code{kmeans} R function.}
            \item{\code{pgmm}: apply the corresponding PGMM model of the \code{pgmm} R package and use the resulting model fit as a starting point.}
            \item{\code{hard}: begin with observations hard-clustered according to the \code{init_class} parameter.}
            \item{\code{soft}: begin with observations soft-clustered according to the \code{init_class} parameter.}
        }
        See Details for more information.
    }
    \item{init_class}{Used when \code{init_method} is set to \code{hard} or \code{soft} to manually prescribe initial classification for the AECM algorithm. If \code{init_method=hard}, a list of vectors of the same length as \code{rG}. Each vector should have length equal to the number of observations. In this case, the elements of the vector must be numbers in \code{{1,2,...,G}}, where \eqn{G} is the number of groups. If \code{init_method=hard}, a list of row-stochastic matrices of same length as \code{rG}. These matrices have number of rows equal to the number of observations, and number of columns correspoding to the values in \code{rG}.}
    \item{max_it}{The maximum number of iterations for the AECM algorithm.}
    \item{tol}{The tolerance for the Aitken acceleration stopping criterion.}
    \item{alpha_min}{The minimum allowable value for alpha, which represents the proportion of ``good'' points in a given group.}
    \item{eta_max}{The maximum allowable value for eta, which is the covariance inflation factor for outlying points.}
    \item{scale}{If TRUE, the data is scaled before the algorithm is begun. Reccomended.}
    \item{parallel}{If TRUE, computation takes place in parallel on several processors.}
    \item{cores}{Only relevent if \code{parallel=TRUE}. Determines the number of cores used in parallel computation. If left undefined, number of available cores is determined automatically.}
    \item{silent}{If TRUE, function will not print any output at completion.}
}

\details{

    This function implements the Mixtures of Contaminated Gaussian Factor Analyzers (MCGFA) model, for model-based clustering and classification. The algorithm is meant to be applied on on high-dimensional, and noisy, data. A description of the model can be found in Punzo and McNicholas (2014). Parameter estimation is performed using an Alternating Expectation-Conditional Maximization (AECM) algorithm (Meng and Van Dyk, 1997).

    Besides clustering into components, this algorithm also automatically detects outliers through the use of the contaminated Gaussian distribution. So, in the end each observation is classified in a nested fashion: first into components, then into as inliers/outliers.

    To specify an individual MCGFA model, one must determine: the number of components (\eqn{G}), the number of latent factors (\eqn{q}), and the model name. The model name is one of the following eight:
    \itemize{
    \item{\{CCC, CCU, CUC, CUU, UCC, UCU, UUC, UUU\}}
    }
    Explaination of the model names can be found in the next subsection.

    The user may choose to fit a single model. However, usually many models are fit to the same data, and a model selection criterion is used to determine the best one. By entering ranges for the \code{rG}, \code{rq} and \code{models} parameters, many models can be fit at once. The \code{mcgfa} function then selects the best model according to the Bayesian Information Criterion (BIC).

    When fitting many models to large data, parallel computation may be employed by \code{mcgfa}, using the \code{parallel} parameter. This parallelization is provided by the \code{mcapply} function from the \code{parallel} package for R.

    \subsection{Model Names & Covariance Constraints}{
    The model name indicates which covariance constraints are to be imposed on to the covariance structure of the factor analysis model. Because the MCGFA is a mixture of factor analyzers model, the covariance matrix of the \eqn{g}-th group, \eqn{\Sigma_g}, can be decomposed as follows:
    \deqn{\Sigma_g = \Lambda_g\Lambda_g' + \Psi_g}
    where \eqn{\Lambda_g} is a \eqn{p} by \eqn{q} matrix of factor loadings, and \eqn{\Psi_g} is a diagonal matrix that determines the noise variance for each of the \eqn{p} variables. The family of eight models is formed by introducing different sets of constraints on \eqn{\Lambda} and \eqn{\Psi}.

    The three-letter model names are interpreted as follows. ``C'' indicates that the constraint is imposed, and ``U'' that it is not. The first letter represents constraining \eqn{\Lambda} to be equal across all groups. The second represents constraining \eqn{\Psi} to be equal across all groups. The last letter represents constraining the error variances to be equal \emph{within} groups; that is, isotropic errors.}

    \subsection{Initialization Methods}{
        Because the AECM algorithm cannot guarantee convergence to a global maximum, success of the algorithm in reaching a good fit depends heavily on starting values. Starting values in this case refer to the intitial classification of each observation. This classification may be either ``hard'', in which each observation is assigned to a unique class, or ``soft'', in which each observation is assigned a vector representing the probability of membership in each class.

        The initialization method is specified through the \code{init_method} argument. The default initialization method is \code{kmeans}. \eqn{k}-means is a fast, simple clustering algorithm that returns a hard classification. The second option is \code{pgmm}, which applies the analogous model from the \code{pgmm} package as a starting point. \code{pgmm} is the model that \code{mcgfa} extends, and is similar except that it lacks the capacity to account for and detect outliers. This method begins with a soft classification.

        The user may also manually assign initial values by selecting one of the \code{hard} and \code{soft} options for \code{init_method}. Both the \code{hard} and \code{soft} options require \code{init_class} be a list of length equal to that of \code{rG}. The \code{i}-th element of this list presribes the starting values for models with number of groups corresponding to \code{rG[i]}.

        In the \code{hard} case, each element of the \code{init_class} list must be a vector of length N, where the \code{i}-th observation indicates the initial membership of observation \code{i}.

        In the \code{soft} case, each element of the \code{init_class} list must a matrix with N rows and number of columns corresponding the matching value of \code{rG}. The (i,j)-th element of each matrix indicates the initial probability that observation \code{i} lies in class \code{j}.
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
}

\references{
Meng, X.-L. and Van Dyk, D. (1997). The EM Algorithm - An old folksong sung to a fast new tune. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{59(3)}, 511-567.
Punzo, A. and McNicholas, P. D. (2014). Robust high-dimensional modeling with the contamined Gaussian distribution. \emph{arXiv preprint arXiv:1408.2128v1.}
}

\examples{
\dontrun{
# TOY EXAMPLE: 2 GAUSSIAN COMPONENTS WITH NOISE
# small example for easy visualization
data(gaussnoise)
toy_classes = gaussnoise[,1]
is_noise = gaussnoise[,2]
X = gaussnoise[,-c(1,2)]

toy_fit = mcgfa(X, rG=1:3, rq=1:2)
# plot: outliers => triangle
plot(toy_fit)
# check clustering performance: 2 errors
table(toy_classes[1:175],toy_fit$group[1:175])
# check noise detection performance: 2 false pos, 3 false neg
table(is_noise,toy_fit$isBad)

# REAL DATA EXAMPLE: WINE DATA
data(wine)
# simulate clustering by completely hiding known classes from algorithm
X = wine[,-1]
wine_classes = wine[,1]
wine_fit = mcgfa(X,rG=1:3,rq=1:4)

# check performance:
#  - correctly selected 3 groups
#  - 2 errors
table(wine_classes,wine_fit$group)

# CLASSIFICATION EXAMPLE: OLIVE DATA

# load data
data(olive)
X = olive[,-c(1,2)]
# classes correspond to regions of Italy
regions = olive[,1]

# take fifth of observations to form training data
set.seed(1)
train_ind = sample.int(nrow(X),nrow(X)/5)

known = rep(0,nrow(X))
known[train_ind] = regions[train_ind]

# FULL SUPERVISION (test set completely hidden from classifier)
fit = mcgfa(X[train_ind,],rq=3:4,rG=3)
pred = predict(fit,X[-train_ind,])
# check classification performace: mass confusion
table(regions[-train_ind],pred$hard)

# PARTIAL SUPERVISION (unlabelled data available to classifier)
fit2 = mcgfa(X,rq=3:4,rG=3,known=known)
# check classification performace: only 13 errors (2\% misclassification)
table(regions[-train_ind],fit2$group[-train_ind])
}
}

\keyword{multivariate}
\keyword{cluster}
\keyword{classif}
