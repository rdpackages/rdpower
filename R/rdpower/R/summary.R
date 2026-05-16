###################################################################
# Summary and print methods for rdpower
# !version 3.0 15-May-2026
# Authors: Matias D. Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' Summarize RDPOWER objects
#'
#' \code{summary()} prints the default robust bias-corrected output. Use
#' \code{all = TRUE} to also display the conventional output.
#'
#' @param object object returned by \code{\link{rdpower}()},
#'   \code{\link{rdsampsi}()}, or \code{\link{rdmde}()}.
#' @param x object returned by \code{\link{rdpower}()},
#'   \code{\link{rdsampsi}()}, or \code{\link{rdmde}()}.
#' @param all displays the conventional output in addition to the robust
#'   bias-corrected output.
#' @param ... additional arguments.
#'
#' @return
#' A summary object containing the original return values, returned invisibly
#' by the print method.
#'
#' @export
summary.rdpower <- function(object, ..., all = FALSE) {
  out <- object
  out$.all <- isTRUE(all)
  class(out) <- c("summary.rdpower", class(object))
  return(out)
}

#' @rdname summary.rdpower
#' @export
print.rdpower <- function(x, ...) {
  print(summary(x, ...))
  invisible(x)
}

#' @rdname summary.rdpower
#' @export
print.summary.rdpower <- function(x, ...) {
  d <- x$.display
  all <- isTRUE(x$.all)

  cat('\n')
  cat(paste0(format('Number of obs =', width=22), toString(d$N.disp))); cat("\n")
  cat(paste0(format('BW type       =', width=22), d$bwselect)); cat("\n")
  cat(paste0(format('Kernel type   =', width=22), d$kernel_type)); cat("\n")
  cat(paste0(format('VCE method    =', width=22), d$vce_type)); cat("\n")
  cat(paste0(format('Derivative    =', width=22), toString(d$deriv))); cat("\n")
  cat(paste0(format('HA:       tau =', width=22), round(x$tau,3))); cat("\n")
  if (all){
    cat(paste0(format('Size dist.    =', width=22), round(d$size.dist,3))); cat("\n")
  }

  cat('\n\n')

  cat(paste0(format(paste0("Cutoff c = ", toString(round(d$cutoff, 3))), width=22), format("Left of c", width=16), format("Right of c", width=16))); cat("\n")
  cat(paste0(format("Number of obs",      width=22), format(toString(d$nminus.disp), width=16), format(toString(d$nplus.disp), width=16))); cat("\n")
  cat(paste0(format("Eff. number of obs", width=22), format(toString(d$nhl.disp),    width=16), format(toString(d$nhr.disp),   width=16))); cat("\n")
  cat(paste0(format("BW loc. poly.",      width=22), format(toString(round(d$hl,3)), width=16), format(toString(round(d$hr,3)), width=16))); cat("\n")
  cat(paste0(format("Order loc. poly.",   width=22), format(toString(d$p),           width=16), format(toString(d$p),           width=16))); cat("\n")

  if (d$clustered){
    cat(paste0(format("Number of clusters",    width=22), format(toString(d$gminus),     width=16), format(toString(d$gplus),     width=16))); cat("\n")
    cat(paste0(format("Eff. num. of clusters", width=22), format(toString(d$gminus_h_l), width=16), format(toString(d$gplus_h_r), width=16))); cat("\n")
  }

  cat(paste0(format("Sampling BW", width=22), format(toString(round(d$hnew.l,3)), width=16), format(toString(round(d$hnew.r,3)), width=16))); cat("\n")
  cat(paste0(format(d$text_aux,    width=22), format(toString(d$ntilde.l),        width=16), format(toString(d$ntilde.r),        width=16))); cat("\n")
  cat("\n\n")

  cat(paste0(rep('=',89),collapse='')); cat('\n')
  cat(paste0(format('Power against:', width=25),
             format('H0: tau = ', width=15),
             format('0.2*tau = ', width=15),
             format('0.5*tau = ', width=15),
             format('0.8*tau = ', width=13),
             format('tau = ', width=15))); cat('\n')

  cat(paste0(format('', width=25),
             format(toString(round(0,3)), width=15),
             format(toString(round(0.2*x$tau,3)), width=15),
             format(toString(round(0.5*x$tau,3)), width=15),
             format(toString(round(0.8*x$tau,3)), width=13),
             format(toString(round(x$tau,3)), width=15))); cat('\n')
  cat(paste0(rep('-',89),collapse='')); cat('\n')
  cat(paste0(format('Robust bias-corrected', width=25),
             format(toString(round(d$power.rbc.list[1],3)), width=15),
             format(toString(round(d$power.rbc.list[2],3)), width=15),
             format(toString(round(d$power.rbc.list[3],3)), width=15),
             format(toString(round(d$power.rbc.list[4],3)), width=13),
             format(toString(round(x$power.rbc,3)), width=15)))

  if (all){
    cat('\n')
    cat(paste0(format('Conventional', width=25),
               format(toString(round(d$power.conv.list[1],3)), width=15),
               format(toString(round(d$power.conv.list[2],3)), width=15),
               format(toString(round(d$power.conv.list[3],3)), width=15),
               format(toString(round(d$power.conv.list[4],3)), width=13),
               format(toString(round(x$power.conv,3)), width=15))); cat('\n')
    cat(paste0(rep('=',89),collapse='')); cat('\n')
  } else {cat('\n');cat(paste0(rep('=',89),collapse=''));cat('\n\n')}

  invisible(x)
}

#' @rdname summary.rdpower
#' @export
summary.rdsampsi <- function(object, ..., all = FALSE) {
  out <- object
  out$.all <- isTRUE(all)
  class(out) <- c("summary.rdsampsi", class(object))
  return(out)
}

#' @rdname summary.rdpower
#' @export
print.rdsampsi <- function(x, ...) {
  print(summary(x, ...))
  invisible(x)
}

#' @rdname summary.rdpower
#' @export
print.summary.rdsampsi <- function(x, ...) {
  d <- x$.display
  all <- isTRUE(x$.all)

  cat('\n')
  cat(paste0(format('Number of obs =', width=22), toString(d$N.disp))); cat("\n")
  cat(paste0(format('BW type       =', width=22), d$bwselect)); cat("\n")
  cat(paste0(format('Kernel type   =', width=22), d$kernel_type)); cat("\n")
  cat(paste0(format('VCE method    =', width=22), d$vce_type)); cat("\n")
  cat(paste0(format('Derivative    =', width=22), toString(d$deriv))); cat("\n")
  cat(paste0(format('HA:       tau =', width=22), round(x$tau,3))); cat("\n")
  cat(paste0(format('Power         =', width=22), round(x$beta,3))); cat("\n")
  if (all){
    cat(paste0(format('Size dist.    =', width=22), round(d$size.dist,3))); cat("\n")
  }
  cat('\n\n')

  cat(paste0(format(paste0("Cutoff c = ", toString(round(d$cutoff, 3))), width=22), format("Left of c", width=16), format("Right of c", width=16))); cat("\n")
  cat(paste0(format("Number of obs",      width=22), format(toString(d$nminus.disp),      width=16), format(toString(d$nplus.disp),        width=16))); cat("\n")
  cat(paste0(format("Eff. number of obs", width=22), format(toString(d$n.hnew.l.disp),    width=16), format(toString(d$n.hnew.r.disp),    width=16))); cat("\n")
  cat(paste0(format("BW loc. poly.",      width=22), format(toString(round(d$hl,3)),      width=16), format(toString(round(d$hr,3)),      width=16))); cat("\n")
  cat(paste0(format("Order loc. poly.",   width=22), format(toString(d$p),                width=16), format(toString(d$p),                width=16))); cat("\n")
  if (d$clustered){
    cat(paste0(format("Number of clusters",    width=22), format(toString(d$gminus),     width=16), format(toString(d$gplus),     width=16))); cat("\n")
    cat(paste0(format("Eff. num. of clusters", width=22), format(toString(d$gminus_h_l), width=16), format(toString(d$gplus_h_r), width=16))); cat("\n")
  }
  cat(paste0(format("Sampling BW", width=22), format(toString(round(d$hnew.l,3)), width=16), format(toString(round(d$hnew.r,3)), width=16))); cat("\n")
  cat("\n\n")

  cat(paste0(rep('=',89),collapse='')); cat('\n')
  if (d$clustered){
    cat(paste0(format('', width=28),
               format('Number of clusters in window', width=35),
               format('Proportion', width=15))); cat('\n')
  } else{
    cat(paste0(format('', width=33),
               format('Number of obs in window', width=35),
               format('Proportion', width=15))); cat('\n')
  }

  cat(paste0(format('', width=25),
             format('[c-h,c)', width=15),
             format('[c,c+h]', width=15),
             format('Total', width=15),
             format('[c,c+h]', width=13))); cat('\n')
  cat(paste0(rep('-',89),collapse='')); cat('\n')
  cat(paste0(format('Robust bias-corrected', width=25),
             format(toString(d$Ml), width=15),
             format(toString(d$Mr), width=15),
             format(toString(d$M), width=15),
             format(toString(round(d$nratio,3)), width=13)))

  if (all){
    cat('\n')
    cat(paste0(format('Conventional', width=25),
               format(toString(d$Ml.cl), width=15),
               format(toString(d$Mr.cl), width=15),
               format(toString(d$M.cl), width=15),
               format(toString(round(d$nratio.cl,3)), width=13))); cat('\n')
    cat(paste0(rep('=',89),collapse='')); cat('\n')
  } else {cat('\n');cat(paste0(rep('=',89),collapse=''));cat('\n\n')}

  invisible(x)
}

#' @rdname summary.rdpower
#' @export
summary.rdmde <- function(object, ..., all = FALSE) {
  out <- object
  out$.all <- isTRUE(all)
  class(out) <- c("summary.rdmde", class(object))
  return(out)
}

#' @rdname summary.rdpower
#' @export
print.rdmde <- function(x, ...) {
  print(summary(x, ...))
  invisible(x)
}

#' @rdname summary.rdpower
#' @export
print.summary.rdmde <- function(x, ...) {
  d <- x$.display
  all <- isTRUE(x$.all)

  cat('\n')
  cat(paste0(format('Number of obs =', width=22), toString(d$N.disp))); cat("\n")
  cat(paste0(format('BW type       =', width=22), d$bwselect)); cat("\n")
  cat(paste0(format('Kernel type   =', width=22), d$kernel_type)); cat("\n")
  cat(paste0(format('VCE method    =', width=22), d$vce_type)); cat("\n")
  cat(paste0(format('Derivative    =', width=22), toString(d$deriv))); cat("\n")
  cat('\n\n')

  cat(paste0(format(paste0("Cutoff c = ", toString(round(d$cutoff, 3))), width=22), format("Left of c", width=16), format("Right of c", width=16))); cat("\n")
  cat(paste0(format("Number of obs",      width=22), format(toString(d$nminus.disp), width=16), format(toString(d$nplus.disp), width=16))); cat("\n")
  cat(paste0(format("Eff. number of obs", width=22), format(toString(d$nhl.disp),    width=16), format(toString(d$nhr.disp),   width=16))); cat("\n")
  cat(paste0(format("BW loc. poly.",      width=22), format(toString(round(d$hl,3)), width=16), format(toString(round(d$hr,3)), width=16))); cat("\n")
  cat(paste0(format("Order loc. poly.",   width=22), format(toString(d$p),           width=16), format(toString(d$p),           width=16))); cat("\n")
  if (d$clustered){
    cat(paste0(format("Number of clusters",    width=22), format(toString(d$gminus),     width=16), format(toString(d$gplus),     width=16))); cat("\n")
    cat(paste0(format("Eff. num. of clusters", width=22), format(toString(d$gminus_h_l), width=16), format(toString(d$gplus_h_r), width=16))); cat("\n")
  }

  cat(paste0(format("Sampling BW", width=22), format(toString(round(d$hnew.l,3)), width=16), format(toString(round(d$hnew.r,3)), width=16))); cat("\n")
  cat(paste0(format(d$text_aux,    width=22), format(toString(d$ntilde.l),        width=16), format(toString(d$ntilde.r),        width=16))); cat("\n")
  cat("\n\n")

  cat(paste0(rep('=',89),collapse='')); cat('\n')
  cat(paste0(format('MDE for power = ', width=25),
             format('beta = ', width=15),
             format('beta = ', width=15),
             format('beta = ', width=15),
             format('beta = ', width=13),
             format('beta = ', width=15))); cat('\n')

  cat(paste0(format('', width=25),
             format(toString(round(d$beta.list[1],3)), width=15),
             format(toString(round(d$beta.list[2],3)), width=15),
             format(toString(round(x$beta,3)), width=15),
             format(toString(round(d$beta.list[3],3)), width=13),
             format(toString(round(d$beta.list[4],3)), width=15))); cat('\n')
  cat(paste0(rep('-',89),collapse='')); cat('\n')
  cat(paste0(format('Robust bias-corrected', width=25),
             format(toString(round(d$mde.rbc.list[1],3)), width=15),
             format(toString(round(d$mde.rbc.list[2],3)), width=15),
             format(toString(round(x$mde,3)),             width=15),
             format(toString(round(d$mde.rbc.list[3],3)), width=13),
             format(toString(round(d$mde.rbc.list[4],3)), width=15)))

  if (all){
    cat('\n')
    cat(paste0(format('Conventional', width=25),
               format(toString(round(d$mde.conv.list[1],3)), width=15),
               format(toString(round(d$mde.conv.list[2],3)), width=15),
               format(toString(round(x$mde.conv,3)),         width=15),
               format(toString(round(d$mde.conv.list[3],3)), width=13),
               format(toString(round(d$mde.conv.list[4],3)), width=15))); cat('\n')
    cat(paste0(rep('=',89),collapse='')); cat('\n')
  } else {cat('\n');cat(paste0(rep('=',89),collapse=''));cat('\n\n')}

  invisible(x)
}
