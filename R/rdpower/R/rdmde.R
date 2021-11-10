###################################################################
# rdmde: minimum detectable effect calculations for RD designs
# !version 2.1 09-Nov-2021
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' MDE Calculations for RD Designs
#'
#' \code{rdmde()} performs MDE calculations for RD designs.
#'
#'
#' @author
#' Matias Cattaneo, Princeton University. \email{cattaneo@princeton.edu}
#'
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu}
#'
#' Gonzalo Vazquez-Bare, UC Santa Barbara. \email{gvazquez@econ.ucsb.edu}
#'
#' @references
#'
#' Cattaneo, M. D., R. Titiunik and G. Vazquez-Bare. (2019). \href{https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2019_Stata.pdf}{Power Calculations for Regression Discontinuity Designs}. \emph{Stata Journal}, 19(1): 210-245.
#'
#'
#' @param data a matrix (Y,R) containing the outcome variable and the running variable (as column vectors).
#' @param cutoff the RD cutoff (default is 0).
#' @param alpha specifies the significance level for the power function. Default is 0.05.
#' @param beta specifies the desired power. Default is 0.8.
#' @param nsamples sets the total sample size to the left, sample size to the left inside the bandwidth, total sample size to the right and sample size to the right of the cutoff inside the bandwidth to calculate the variance when the running variable is not specified. When not specified, the values are calculated using the running variable.
#' @param sampsi sets the sample size at each side of the cutoff for power calculation. The first number is the sample size to the left of the cutoff and the second number is the sample size to the right. Default values are the sample sizes inside the chosen bandwidth.
#' @param samph sets the bandwidths at each side of the cutoff for power calculation. The first number is the bandwidth to the left of the cutoff and the second number is the bandwidth to the right.  Default values are the bandwidths used by \code{rdrobust}.
#' @param all displays the power using the conventional variance estimator, in addition to the robust bias corrected one.
#' @param bias set bias to the left and right of the cutoff. If not specified, the biases are estimated using \code{rdrobust}.
#' @param variance set variance to the left and right of the cutoff. If not specified, the variances are estimated using \code{rdrobust}.
#' @param init.cond sets the initial condition for the Newton-Raphson algorithm that finds the MDE.  Default is 0.2 times the standard deviation of the outcome below the cutoff.
#' @param covs option for \code{rdrobust()}: specifies additional covariates to be used for estimation and inference.
#' @param covs_drop option for \code{rdrobust()}: if TRUE, it checks for collinear additional covariates and drops them. Default is TRUE.
#' @param deriv option for \code{rdrobust()}: specifies the order of the derivative of the regression functions to be estimated.
#' @param p option for \code{rdrobust()}: specifies the order of the local-polynomial used to construct the point-estimator.
#' @param q option for \code{rdrobust()}: specifies the order of the local-polynomial used to construct the bias-correction.
#' @param h option for \code{rdrobust()}: specifies the values of the main bandwidth to be used on the left and on the right of the cutoff, respectively.
#' @param b option for \code{rdrobust()}: specifies the values of the bias bandwidth $b$ to be used on the left and on the right of the cutoff, respectively.
#' @param rho option for \code{rdrobust()}: specifies the value of \code{rho} so that the bias bandwidth \code{b} equals \code{b=h/rho}.
#' @param kernel option for \code{rdrobust()}: kernel function used to construct the local-polynomial estimators.
#' @param bwselect option for \code{rdrobust()}: specifies the bandwidth selection procedure to be used.
#' @param vce option for \code{rdrobust()}: specifies the procedure used to compute the variance-covariance matrix estimator.
#' @param cluster option for \code{rdrobust()}: indicates the cluster ID variable used for the cluster-robust variance estimation with degrees-of-freedom weights.
#' @param scalepar option for \code{rdrobust()}: specifies scaling factor for RD parameter of interest.
#' @param scaleregul option for \code{rdrobust()}: specifies scaling factor for the regularization terms of bandwidth selectors.
#' @param fuzzy option for \code{rdrobust()}: specifies the treatment status variable used to implement fuzzy RD estimation.
#' @param level option for \code{rdrobust()}: sets the confidence level for confidence intervals.
#' @param weights option for \code{rdrobust()}: is the variable used for optional weighting of the estimation procedure. The unit-specific weights multiply the kernel function.
#' @param masspoints option for \code{rdrobust()}: checks and controls for repeated observations in tue running variable.
#' @param bwcheck option for \code{rdrobust()}: if a positive integer is provided, the preliminary bandwidth used in the calculations is enlarged so that at least \code{bwcheck} unique observations are used.
#' @param bwrestrict option for \code{rdrobust()}: if TRUE, computed bandwidths are restricted to lie withing the range of \code{x}. Default is \code{bwrestrict=TRUE}.
#' @param stdvars option for \code{rdrobust()}: if \code{TRUE}, \code{x} and \code{y} are standardized before computing the bandwidths. Default is \code{stdvars=TRUE}.
#'
#' @return
#' \item{mde}{MDE using robust bias corrected standard error}
#' \item{se.rbc}{robust bias corrected standard error}
#' \item{sampsi.r}{number of observations inside the window to the right of the cutoff}
#' \item{sampsi.l}{number of observations inside the window to the left of the cutoff}
#' \item{samph.r}{bandwidth to the right of the cutoff}
#' \item{samph.l}{bandwidth to the left of the cutoff}
#' \item{alpha}{significance level used in power function}
#' \item{bias.r}{bias to the right of the cutoff}
#' \item{bias.l}{bias to the left of the cutoff}
#' \item{Vr.rb}{Robust bias corrected variance to the right of the cutoff}
#' \item{Vl.rb}{Robust bias corrected variance to the left of the cutoff}
#' \item{N.r}{Total sample size to the right of the cutoff}
#' \item{N.l}{Total sample size to the left of the cutoff}
#' \item{mde.conv}{MDE using conventional inference}
#' \item{se.conv}{conventional standard error}
#'
#'
#' @examples
#' # Toy dataset
#' X <- array(rnorm(2000),dim=c(1000,2))
#' R <- X[,1] + X[,2] + rnorm(1000)
#' Y <- 1 + R -.5*R^2 + .3*R^3 + (R>=0) + rnorm(1000)
#' # MDE calculation
#' tmp <- rdmde(data=cbind(Y,R),init.cond=0.5)
#'
#'
#' @export

rdmde <- function(data = NULL,
                  cutoff = 0,
                  alpha = 0.05,
                  beta = 0.8,
                  nsamples = NULL,
                  sampsi = NULL,
                  samph = NULL,
                  all = FALSE,
                  bias = NULL,
                  variance = NULL,
                  init.cond = NULL,

                  covs = NULL,
                  covs_drop = TRUE,
                  deriv = 0,
                  p = 1,
                  q = NULL,
                  h = NULL,
                  b = NULL,
                  rho = NULL,
                  kernel = 'triangular',
                  bwselect = 'mserd',
                  vce = 'nn',
                  cluster = NULL,
                  scalepar = 1,
                  scaleregul = 1,
                  fuzzy = NULL,
                  level = 95,
                  weights = NULL,
                  masspoints = 'adjust',
                  bwcheck = NULL,
                  bwrestrict = TRUE,
                  stdvars = FALSE){

  #################################################################
  # Options, default values and error checking
  #################################################################

  if (!is.null(data)){
    if (ncol(data)<2){stop('Too few variables specified in data')}
    else if (ncol(data)>2){stop('Too many variables specified in data')}
    else{
      Y <- data[,1]
      R <- data[,2]
    }
  }

  if (!is.null(nsamples)){
    if (!is.null(data)){warning('nsamples ignored when data is specified')}
    else{
      if (!(!is.null(bias) & !is.null(variance) & !is.null(samph) & !is.null(init.cond))){stop('not enough information to calculate power without data')}
      if (length(nsamples)!=4){
        stop('incorrect number of arguments in nsamples')
      } else{
        if (sum(nsamples%%1)!=0){stop('sample sizes in nsamples have to be integers')}
        if (min(nsamples)<=0){stop('sample sizes in nsamples have to be >0')}
        nminus <- nsamples[1]
        n.hnew.l <- nsamples[2]
        nplus <- nsamples[3]
        n.hnew.r <- nsamples[4]
      }
    }
  }

  if (!is.null(sampsi)){
    if (sum(sampsi%%1)!=0){stop('sample sizes in sampsi have to be integers')}
    if (min(sampsi)<=0){stop('sample sizes in sampsi have to be >0')}
    if (length(sampsi)==1){
      ntilde.l <- sampsi
      ntilde.r <- sampsi
    } else if (length(sampsi)==2){
      ntilde.l <- sampsi[1]
      ntilde.r <- sampsi[2]
    } else{
      stop('sampsi incorrectly specified')
    }
  }

  if (!is.null(samph)){
    if (length(samph)==1){
      hnew.l <- samph
      hnew.r <- samph
    } else if (length(samph)==2){
      hnew.l <- samph[1]
      hnew.r <- samph[2]
    } else{
      stop('samph incorrectly specified')
    }
  }

  if (!is.null(bias)){
    if (length(bias)==1){stop('need to specify both Bl and Br')}
    else if (length(bias)==2){
      bias.l <- bias[1]
      bias.r <- bias[2]
    }
    else {stop('bias incorrectly specified')}
  }

  if (!is.null(variance)){
    if (min(variance)<=0){stop('variances have to be >0')}
    if (length(variance)==1){stop('need to specify both Vl and Vr')}
    else if (length(variance)==2){
      Vl.rb <- variance[1]
      Vr.rb <- variance[2]
    }
    else {stop('variance incorrectly specified')}
  }

  if (is.null(q)){ q <- p + 1}


  #################################################################
  # Definition of bias, variance, sample sizes and bandwidths
  #################################################################

  if (!is.null(data)){

    if (is.null(bias) | is.null(variance)){
      aux <- rdrobust::rdrobust(Y,R,c=cutoff,all=TRUE,covs=covs,covs_drop=covs_drop,deriv=deriv,p=p,q=q,h=h,b=b,rho=rho,cluster=cluster,
                     kernel=kernel,bwselect=bwselect,vce=vce,scalepar=scalepar,scaleregul=scaleregul,
                     fuzzy=fuzzy,level=level,weights=weights,masspoints=masspoints,bwcheck=bwcheck,bwrestrict=bwrestrict,stdvars=stdvars)

      h.aux <- aux$bws
      h.l <- h.aux[1,1]
      h.r <- h.aux[1,2]

      if (is.null(bias)){
        bias <- aux$bias
        bias.l <- bias[1]/(h.l^(1+p-deriv))
        bias.r <- bias[2]/(h.r^(1+p-deriv))
      }

      if (is.null(variance)){
        VL.CL <- aux$V_cl_l
        VR.CL <- aux$V_cl_r
        VL.RB <- aux$V_rb_l
        VR.RB <- aux$V_rb_r

        if (!is.null(cluster)){
          N <- length(table(cluster[!is.na(Y) & !is.na(R)]))
        } else{
          N <- sum(!is.na(Y) & !is.na(R))
        }
        pos <- 1 + deriv
        Vl.cl <- N*(h.l^(1+2*deriv))*VL.CL[pos,pos]
        Vr.cl <- N*(h.r^(1+2*deriv))*VR.CL[pos,pos]
        Vl.rb <- N*(h.l^(1+2*deriv))*VL.RB[pos,pos]
        Vr.rb <- N*(h.r^(1+2*deriv))*VR.RB[pos,pos]
      }
    }

    ## Set default new bandwith

    if (is.null(samph)){
      hnew.l <- h.l
      hnew.r <- h.r
    }

    ## Calculate sample sizes

    if (!is.null(cluster)){
      N <- length(table(cluster[!is.na(Y)&!is.na(R)]))
      nplus <- length(table(cluster[R>=cutoff & !is.na(Y) & !is.na(R)]))
      nminus <- length(table(cluster[R<cutoff & !is.na(Y) & !is.na(R)]))
      n.hnew.r <- length(table(cluster[R>=cutoff & R<= cutoff + hnew.r & !is.na(Y) & !is.na(R)]))
      n.hnew.l <- length(table(cluster[R<cutoff & R>= cutoff - hnew.l & !is.na(Y) & !is.na(R)]))
    } else{
      N <- sum(!is.na(Y) & !is.na(R))
      nplus <- sum(R>=cutoff & !is.na(Y) & !is.na(R))
      nminus <- sum(R<cutoff & !is.na(Y) & !is.na(R))
      n.hnew.r <- sum(R>=cutoff & R<= cutoff + hnew.r & !is.na(Y) & !is.na(R))
      n.hnew.l <- sum(R<cutoff & R>= cutoff - hnew.l & !is.na(Y) & !is.na(R))
    }

  }

  if (is.null(sampsi)){
    ntilde.l <- n.hnew.l
    ntilde.r <- n.hnew.r
  }

  ntilde <- nplus*(ntilde.r/n.hnew.r) + nminus*(ntilde.l/n.hnew.l)


  #################################################################
  # Variance and bias adjustment
  #################################################################

  ## Variance adjustment

  V.rbc <- Vl.rb/(hnew.l^(1+2*deriv)) + Vr.rb/(hnew.r^(1+2*deriv))
  se.rbc <- sqrt(V.rbc)

  if (all==TRUE){
    V.conv <- Vl.cl/(hnew.l^(1+2*deriv)) + Vr.cl/(hnew.r^(1+2*deriv))
    se.conv <- sqrt(V.conv)
  } else{
    V.conv <- V.rbc
    se.conv <- se.rbc
  }

  ## Bias adjustment

  bias <- bias.r*hnew.r^(1+p-deriv) + bias.l*hnew.l^(1+p-deriv)

  #################################################################
  # MDE calculation
  #################################################################


  ## Set initial condition for Newton-Raphson

  if (is.null(init.cond)){
    sd0 <- sd(Y[R<cutoff],na.rm=TRUE)
    tau0 <- 0.2*sd0
  } else{
    tau0 <- init.cond
  }

  cat('Calculating MDE...')

  mde.aux <- rdpower.powerNR.mde(ntilde,tau0,se.rbc,qnorm(1-alpha/2),beta)
  mde <- mde.aux$mde

  beta.list <- numeric(4)
  mde.rbc.list <- numeric(4)

  count <- 1
  for (r in c(-0.125,-0.0625,0.0625,0.125)){

    baux <- beta*(1+r)
    beta.list[count] <- baux

    if (baux<1){
      rdpow.aux <- rdpower.powerNR.mde(ntilde,tau0,se.rbc,qnorm(1-alpha/2),baux)
      mde.rbc.list[count] <- rdpow.aux$mde

    }

    count <- count + 1
  }

  if(all==TRUE){
    mde.conv.aux <- rdpower.powerNR.mde(ntilde,tau0+bias,se.conv,qnorm(1-alpha/2),beta)
    mde.conv <- mde.conv.aux$mde

    mde.conv.list <- numeric(4)
    count <- 1
    for (r in c(-0.125,-0.0625,0.0625,0.125)){

      baux <- beta*(1+r)

      if (baux<1){
        rdpow.conv.aux <- rdpower.powerNR.mde(ntilde,tau0+bias,se.conv,qnorm(1-alpha/2),baux)
        mde.conv.list[count] <- rdpow.conv.aux$mde

      }

      count <- count + 1
    }

  }

  cat('MDE obtained')

  #################################################################
  # Descriptive statistics for display
  #################################################################

  if (!is.null(data)){

    # Left panel

    nplus.disp <- sum(R>=cutoff & !is.na(Y) & !is.na(R))
    nminus.disp <- sum(R<cutoff & !is.na(Y) & !is.na(R))

    if (is.null(h)){
      hl <- hnew.l
      hr <- hnew.r
    } else{
      hl <- h.l
      hr <- h.r
    }

    nhr.disp <- sum(R>=cutoff & R<= cutoff + hr & !is.na(Y) & !is.na(R))
    nhl.disp <- sum(R<cutoff & R>= cutoff - hl & !is.na(Y) & !is.na(R))
    
    if (!is.null(cluster)){
      gplus <- length(table(cluster[R>=cutoff & !is.na(Y) & !is.na(R)]))
      gminus <- length(table(cluster[R<cutoff & !is.na(Y) & !is.na(R)]))
      gplus_h_r <- length(table(cluster[R>=cutoff & R<= cutoff + hr & !is.na(Y) & !is.na(R)]))
      gminus_h_l <- length(table(cluster[R<cutoff & R>= cutoff - hl & !is.na(Y) & !is.na(R)]))
    }

    # Right panel

    N.disp <- sum(!is.na(Y) & !is.na(R))

    if (!is.null(bias) & !is.null(variance)){
      bwselect <- NA
      kernel_type <- NA
      vce_type <- NA
    } else{
      bwselect <- aux$bwselect
      kernel_type <- aux$kernel
      vce_type <- aux$vce
    }


  } else{

    # Left panel

    if (!is.null(cluster)){
      gplus <- NA
      gminus <- NA
    }
    nplus.disp <- NA
    nminus.disp <- NA
    nhr.disp <- NA
    nhl.disp <- NA

    hl <- NA
    hr <- NA
    p <- NA

    # Right panel

    N.disp <- NA
    bwselect <- NA
    kernel_type <- NA
    vce_type <- NA
  }

  #################################################################
  # Output
  #################################################################

  cat('\n')
  cat(paste0(format('Number of obs =', width=22), toString(N.disp))); cat("\n")
  cat(paste0(format('BW type       =', width=22), bwselect)); cat("\n")
  cat(paste0(format('Kernel type   =', width=22), kernel_type)); cat("\n")
  cat(paste0(format('VCE method    =', width=22), vce_type)); cat("\n")
  cat(paste0(format('Derivative    =', width=22), toString(deriv))); cat("\n")
  cat('\n\n')

  cat(paste0(format(paste0("Cutoff c = ", toString(round(cutoff, 3))), width=22), format("Left of c", width=16), format("Right of c", width=16))); cat("\n")
  cat(paste0(format("Number of obs",      width=22), format(toString(nminus.disp),     width=16), format(toString(nplus.disp),        width=16))); cat("\n")
  cat(paste0(format("Eff. number of obs", width=22), format(toString(nhl.disp),        width=16), format(toString(nhr.disp),          width=16))); cat("\n")
  cat(paste0(format("BW loc. poly.",      width=22), format(toString(round(hl,3)),     width=16), format(toString(round(hr,3)),       width=16))); cat("\n")
  cat(paste0(format("Order loc. poly.",   width=22), format(toString(p),               width=16), format(toString(p),                 width=16))); cat("\n")
  text_aux <- "New sample"
  if (!is.null(cluster)){
    cat(paste0(format("Number of clusters",    width=22), format(toString(gminus),     width=16), format(toString(gplus),             width=16))); cat("\n")
    cat(paste0(format("Eff. num. of clusters", width=22), format(toString(gminus_h_l), width=16), format(toString(gplus_h_r),         width=16))); cat("\n")
    text_aux<- "New cluster sample"
  }
  
  cat(paste0(format("Sampling BW",    width=22), format(toString(round(hnew.l,3)), width=16), format(toString(round(hnew.r,3)),   width=16))); cat("\n")
  cat(paste0(format(text_aux,         width=22), format(toString(ntilde.l),        width=16), format(toString(ntilde.r),          width=16))); cat("\n")
  cat("\n\n")

  cat(paste0(rep('=',89),collapse='')); cat('\n')
  cat(paste0(format('MDE for power = ', width=25),
             format('beta = ',     width=15),
             format('beta = ',     width=15),
             format('beta = ',     width=15),
             format('beta = ',     width=13),
             format('beta = ',         width=15))); cat('\n')

  cat(paste0(format('', width=25),
             format(toString(round(beta.list[1],3)), width=15),
             format(toString(round(beta.list[2],3)), width=15),
             format(toString(round(beta,3))        , width=15),
             format(toString(round(beta.list[3],3)), width=13),
             format(toString(round(beta.list[4],3)), width=15))); cat('\n')
  cat(paste0(rep('-',89),collapse='')); cat('\n')
  cat(paste0(format('Robust bias-corrected', width=25),
             format(toString(round(mde.rbc.list[1],3)), width=15),
             format(toString(round(mde.rbc.list[2],3)), width=15),
             format(toString(round(mde,3))            , width=15),
             format(toString(round(mde.rbc.list[3],3)), width=13),
             format(toString(round(mde.rbc.list[4],3)), width=15)));

  if (all==TRUE){
    cat('\n')
    cat(paste0(format('Conventional', width=25),
               format(toString(round(mde.conv.list[1],3)), width=15),
               format(toString(round(mde.conv.list[2],3)), width=15),
               format(toString(round(mde.conv,3))        , width=15),
               format(toString(round(mde.conv.list[3],3)), width=13),
               format(toString(round(mde.conv.list[4],3)), width=15))); cat('\n')
    cat(paste0(rep('=',89),collapse='')); cat('\n')

  } else {cat('\n');cat(paste0(rep('=',89),collapse=''));cat('\n\n')}


  #################################################################
  # Return values
  #################################################################

  output <- list(mde = mde,
                se.rbc = se.rbc,
                sampsi.r = ntilde.r,
                sampsi.l = ntilde.l,
                samph.r = hnew.r,
                samph.l = hnew.l,
                N.r = nplus,
                N.l = nminus,
                Nh.l = n.hnew.l,
                Nh.r = n.hnew.r,
                bias.r = bias.r,
                bias.l = bias.l,
                Vr.rb = Vr.rb,
                Vl.rb = Vl.rb,
                alpha = alpha,
                beta = beta)

  if (all==TRUE){
    output <- c(output,
               mde.conv = mde.conv,
               se.conv = se.conv)
  }

  return(output)

}
