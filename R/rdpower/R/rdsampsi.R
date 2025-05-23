###################################################################
# rdsampsi: sample size calculations for RD designs
# !version 2.3 22-May-2025
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' Sample Size Calculations for RD Designs
#'
#' \code{rdsampsi()} performs sample size calculations for RD designs.
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
#' @param tau specifies the treatment effect under the alternative at which the power function is evaluated. The default is half the standard deviation of the outcome for the untreated group.
#' @param alpha specifies the significance level for the power function. Default is 0.05.
#' @param beta specifies the desired power. Default is 0.8.
#' @param nsamples sets the total sample size to the left, sample size to the left inside the bandwidth, total sample size to the right and sample size to the right of the cutoff inside the bandwidth to calculate the variance when the running variable is not specified. When not specified, the values are calculated using the running variable.
#' @param samph sets the bandwidths at each side of the cutoff for power calculation. The first number is the bandwidth to the left of the cutoff and the second number is the bandwidth to the right.  Default values are the bandwidths used by \code{rdrobust}.
#' @param all displays the power using the conventional variance estimator, in addition to the robust bias corrected one.
#' @param bias set bias to the left and right of the cutoff. If not specified, the biases are estimated using \code{rdrobust}.
#' @param variance set variance to the left and right of the cutoff. If not specified, the variances are estimated using \code{rdrobust}.
#' @param nratio specifies the proportion of treated units in the window. Default is the ratio of the standard deviation of the treated to the sum of the standard deviations for treated and controls.
#' @param init.cond sets the initial condition for the Newton-Raphson algorithm that finds the sample size.  Default is the number of observations in the sample with non-missing values of the outcome and running variable.
#' @param plot plots the power function using the conventional and robust bias corrected standard errors from \code{rdrobust}.
#' @param graph.range range of the plot.
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
#' \item{alpha}{significance level}
#' \item{beta}{desired power}
#' \item{tau}{treatment effect under alternative hypothesis}
#' \item{sampsi.h.tot}{total number of observations inside the window}
#' \item{sampsi.h.r}{number of observations inside the window to the right of the cutoff}
#' \item{sampsi.h.l}{number of observations inside the window to the left of the cutoff}
#' \item{N.r}{Total sample size to the right of the cutoff}
#' \item{N.l}{Total sample size to the left of the cutoff}
#' \item{samph.r}{bandwidth to the right of the cutoff}
#' \item{samph.l}{bandwidth to the left of the cutoff}
#' \item{var.r}{Robust bias corrected variance to the right of the cutoff}
#' \item{Var.l}{Robust bias corrected variance to the left of the cutoff}
#' \item{sampsi.h.tot.cl}{implied total number of observations inside the window using conventional s.e.}
#' \item{sampsi.h.r.cl}{number of observations inside the window to the right of the cutoff using conventional s.e.}
#' \item{sampsi.h.l.cl}{number of observations inside the window to the left of the cutoff using conventional s.e.}
#' \item{no.iter}{number of iterations until convergence of the Newton-Raphson algorithm}
#' \item{init.cond}{initial condition of the Newton-Raphson algorithm}
#'
#'
#' @examples
#' # Toy dataset
#' X <- array(rnorm(2000),dim=c(1000,2))
#' R <- X[,1] + X[,2] + rnorm(1000)
#' Y <- 1 + R -.5*R^2 + .3*R^3 + (R>=0) + rnorm(1000)
#' # Sample size to achieve power of 0.8 against tau = 1
#' tmp <- rdsampsi(data=cbind(Y,R),tau=1)
#' # Sample size against tau = 1 including covariates
#' tmp <- rdsampsi(data=cbind(Y,R),tau=1,covs=X)
#'
#'
#' @export

rdsampsi <- function(data = NULL,
                    cutoff = 0,
                    tau = NULL,
                    alpha = 0.05,
                    beta = 0.8,
                    samph = NULL,
                    nsamples = NULL,
                    all = FALSE,
                    bias = NULL,
                    variance = NULL,
                    nratio = NULL,
                    init.cond = NULL,
                    plot = FALSE,
                    graph.range = NULL,

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
    else{
      Y <- data[,1]
      R <- data[,2]
    }
  }

  if (!is.null(nsamples)){
    if (!is.null(data)){warning('nsamples ignored when data is specified')}
    else{
      if (length(nsamples)!=4){
        stop('insufficient arguments in nsamples')
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
      vl <- variance[1]
      vr <- variance[2]
    }
    else {stop('variance incorrectly specified')}
    if (is.null(data)){
      vl.cl <- vl
      vr.cl <- vr
    }
  }

  if (!is.null(variance) & is.null(samph)){stop('need to set samph when variance is specified')}
  if (!is.null(nratio)){
    nratio.cl = nratio
    if (nratio >= 1 | nratio <= 0){
      stop('nratio has to be in (0,1)')
    }
  }

  if (is.null(q)){ q <- p + 1}


  #################################################################
  # Bias and variance
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
        vl <- N*(h.l^(1+2*deriv))*VL.RB[pos,pos]
        vr <- N*(h.r^(1+2*deriv))*VR.RB[pos,pos]
        vl.cl <- N*(h.l^(1+2*deriv))*VL.CL[pos,pos]
        vr.cl <- N*(h.r^(1+2*deriv))*VR.CL[pos,pos]
      }
    }

    if (is.null(samph)){
      hnew.l <- h.l
      hnew.r <- h.r
    }

    if (is.null(vl.cl) | is.null(vr.cl)){
      vl.cl <- vl
      vr.cl <- vr
    }
  }

  # Bias adjustment

  bias <- bias.r*hnew.r^(1+p-deriv) + bias.l*hnew.l^(1+p-deriv)

  # Variance adjustment

  V.rbc <- vl/(hnew.l^(1+2*deriv)) + vr/(hnew.r^(1+2*deriv))
  stilde <- sqrt(V.rbc)
  V.cl <- vl.cl/(hnew.l^(1+2*deriv)) + vr.cl/(hnew.r^(1+2*deriv))
  stilde.cl <- sqrt(V.cl)


  #################################################################
  # Sample size calculation
  #################################################################

  # Critical value

  z <- qnorm(1-alpha/2)

  # Set default value of tau

  if (is.null(tau)){
    sd0 <- sd(Y[R>=cutoff-hnew.l & R<cutoff],na.rm=TRUE)
    tau <- 0.5*sd0
  }

  # Set initial value for Newton-Raphson

  if (is.null(init.cond)){
    init.cond <- sum(!is.na(Y) & !is.na(R))
  }

  # Find m

  cat('Calculating sample size...')

  maux <- rdpower.powerNR(init.cond,tau,stilde,z,beta)
  m <- maux$m

  if (all==TRUE){
    maux1 <- rdpower.powerNR(init.cond,tau+bias,stilde.cl,z,beta)
    m.cl <- maux1$m
  }

  cat('Sample size obtained.')

  # Adjust m to find sample sizes

  if (is.null(nratio)){
    nratio <- sqrt(vr)/(sqrt(vr)+sqrt(vl))
    nratio.cl <- sqrt(vr.cl)/(sqrt(vr.cl)+sqrt(vl.cl))
  }

  if (!is.null(data)){
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

  denom <- nratio*nplus/n.hnew.r + (1-nratio)*nminus/n.hnew.l
  denom.cl <- nratio.cl*nplus/n.hnew.r + (1-nratio.cl)*nminus/n.hnew.l

  M <- m/denom
  Mr <- ceiling(M*nratio)
  Ml <- ceiling(M*(1-nratio))
  M <- Ml + Mr

  if (all==TRUE){
    M.cl <- m.cl/denom.cl
    Mr.cl <- ceiling(M.cl*nratio.cl)
    Ml.cl <- ceiling(M.cl*(1-nratio.cl))
    M.cl <- Ml.cl + Mr.cl
  }

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

    n.hnew.r.disp <- sum(R>=cutoff & R<= cutoff + hr & !is.na(Y) & !is.na(R))
    n.hnew.l.disp <- sum(R<cutoff & R>= cutoff - hl & !is.na(Y) & !is.na(R))

    if (!is.null(cluster)){
      gplus <- length(table(cluster[R>=cutoff & !is.na(Y) & !is.na(R)]))
      gminus <- length(table(cluster[R<cutoff & !is.na(Y) & !is.na(R)]))
      gplus_h_r <- length(table(cluster[R>=cutoff & R<= cutoff + hr & !is.na(Y) & !is.na(R)]))
      gminus_h_l <- length(table(cluster[R<cutoff & R>= cutoff - hl & !is.na(Y) & !is.na(R)]))
    }

    # Right panel

    N.disp <- sum(!is.na(Y) & !is.na(R))

    if (is.null(bias) | is.null(variance)){
      bwselect <- aux$bwselect
      kernel_type <- aux$kernel
      vce_type <- aux$vce
    } else{
      bwselect <- NA
      kernel_type <- NA
      vce_type <- NA
    }


  } else{

    # Left panel

    if (!is.null(cluster)){
      gplus <- NA
      gminus <- NA
    }
    nplus.disp <- NA
    nminus.disp <- NA
    n.hnew.l.disp <- NA
    n.hnew.r.disp <- NA

    hl <- NA
    hr <- NA
    p <- NA

    # Right panel

    N.disp <- NA
    bwselect <- NA
    kernel_type <- NA
    vce_type <- NA
  }

  # Size distortion

  if (all==TRUE){
    se_cl_aux <- stilde.cl / sqrt(m.cl)
    size_dist <- 1 - pnorm(bias/se_cl_aux+z) + pnorm(bias/se_cl_aux-z)

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
  cat(paste0(format('HA:       tau =', width=22), round(tau,3))); cat("\n")
  cat(paste0(format('Power         =', width=22), round(beta,3))); cat("\n")
  if (all==TRUE){
    cat(paste0(format('Size dist.    =', width=22), round(size_dist,3))); cat("\n")
  }
  cat('\n\n')

  cat(paste0(format(paste0("Cutoff c = ", toString(round(cutoff, 3))), width=22), format("Left of c", width=16), format("Right of c", width=16))); cat("\n")
  cat(paste0(format("Number of obs",      width=22), format(toString(nminus.disp),     width=16), format(toString(nplus.disp),        width=16))); cat("\n")
  cat(paste0(format("Eff. number of obs", width=22), format(toString(n.hnew.l.disp),   width=16), format(toString(n.hnew.r.disp),     width=16))); cat("\n")
  cat(paste0(format("BW loc. poly.",      width=22), format(toString(round(hl,3)),     width=16), format(toString(round(hr,3)),       width=16))); cat("\n")
  cat(paste0(format("Order loc. poly.",   width=22), format(toString(p),               width=16), format(toString(p),                 width=16))); cat("\n")
  text_aux <- "New sample"
  if (!is.null(cluster)){
    cat(paste0(format("Number of clusters",    width=22), format(toString(gminus),     width=16), format(toString(gplus),             width=16))); cat("\n")
    cat(paste0(format("Eff. num. of clusters", width=22), format(toString(gminus_h_l), width=16), format(toString(gplus_h_r),         width=16))); cat("\n")
    text_aux<- "New cluster sample"
  }
  cat(paste0(format("Sampling BW",        width=22), format(toString(round(hnew.l,3)), width=16), format(toString(round(hnew.r,3)),   width=16))); cat("\n")
  cat("\n\n")

  cat(paste0(rep('=',89),collapse='')); cat('\n')
  #cat(paste0(format('Chosen sample sizes',   width=33),
  #           format('Sample size in window', width=35),
  #           format('Proportion',            width=15))); cat('\n')

  if (!is.null(cluster)){
    cat(paste0(format('',   width=28),
               format('Number of clusters in window', width=35),
               format('Proportion',            width=15))); cat('\n')
  }
  else{
    cat(paste0(format('',   width=33),
               format('Number of obs in window', width=35),
               format('Proportion',            width=15))); cat('\n')
  }



  cat(paste0(format('',        width=25),
             format('[c-h,c)', width=15),
             format('[c,c+h]', width=15),
             format('Total',   width=15),
             format('[c,c+h]', width=13))); cat('\n')
  cat(paste0(rep('-',89),collapse='')); cat('\n')
  cat(paste0(format('Robust bias-corrected', width=25),
             format(toString(Ml), width=15),
             format(toString(Mr), width=15),
             format(toString(M), width=15),
             format(toString(round(nratio,3)), width=13)))

  if (all==TRUE){
    cat('\n')
    cat(paste0(format('Conventional', width=25),
               format(toString(Ml.cl), width=15),
               format(toString(Mr.cl), width=15),
               format(toString(M.cl), width=15),
               format(toString(round(nratio.cl,3)), width=13))); cat('\n')
    cat(paste0(rep('=',89),collapse='')); cat('\n')
  } else {cat('\n');cat(paste0(rep('=',89),collapse=''));cat('\n\n')}


  #################################################################
  # Power function plot
  #################################################################

  if (plot==TRUE){

    N.plot <- sum(!is.na(Y) & !is.na(R))

    if(is.null(graph.range)){
      left <- 0
      right <- N.plot
    } else{
      left <- graph.range[1]
      right <- graph.range[2]
    }

    plot(function(x) 1 - pnorm(sqrt(x*denom)*tau/stilde+qnorm(1-alpha/2)) + pnorm(sqrt(x*denom)*tau/stilde-qnorm(1-alpha/2)),
          from=left,to=right,xlab='total sample size in window',ylab='power',xlim=c(left,right))
    title('Power function')
    abline(v=M,lty=2,col='gray50')
    abline(h=beta,lty=2,col='gray50')
    if (all==TRUE){
      plot(function(x) 1 - pnorm(sqrt(x*denom.cl)*(tau+bias)/stilde.cl+qnorm(1-alpha/2)) + pnorm(sqrt(x*denom.cl)*(tau+bias)/stilde.cl-qnorm(1-alpha/2)),
            add=TRUE,lty=2)
      legend('bottomleft',legend=c('robust bc','conventional'),lty=c(1,2),bty='n')
    }
  }

  #################################################################
  # Return values
  #################################################################

  output <- list(sampsi.h.tot = M,
                sampsi.h.r = Mr,
                sampsi.h.l = Ml,
                sampsi.tot = m,
                N.r = nplus,
                N.l = nminus,
                Nh.r = n.hnew.r,
                Nh.l = n.hnew.l,
                bias.r = bias.r,
                bias.l = bias.r,
                var.r = vr,
                var.l = vl,
                samph.r = hnew.r,
                samph.l = hnew.l,
                tau = tau,
                beta = beta,
                alpha = alpha,
                init.cond = init.cond,
                no.iter = maux$iter)
  if (all==TRUE){
    output <- c(output,
               sampsi.h.tot.cl = M.cl,
               sampsi.h.r.cl = Mr.cl,
               sampsi.h.l.cl = Ml.cl,
               sampsi.tot.cl = m.cl,
               var.r.cl = vr.cl,
               var.l.cl = vl.cl)
  }

  return(output)

}
