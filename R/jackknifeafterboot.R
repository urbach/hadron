jab <- function(t, t0, starts, m=1, fn=sd) {
  find.duplicates <- function(xstar, x) {
    duplicated(c(xstar, x))[(length(x) + 1):(2 * length(x))]
  }
  
  jack.boot <- function(indices, xstar, f) {
    if(is.null(dim(xstar)))   apply(xstar[!indices], MARGIN=2L, FUN=f)
    else apply(X=xstar[!indices, ], MARGIN=2L, FUN=f)
  }
  
  ## total number of blocks
  N <- ncol(starts)
  ## number of blocks of blocks
  M <- N - m + 1
  duplicates <- t(apply(X=starts, MARGIN=1L, FUN=find.duplicates, x=c(1:N)))
  if(m > 1) {
    for(i in c(1:M)) {
      duplicates[,i] <- apply(duplicates[,c(i:(i+m-1))], MARGIN=1L, FUN=any)
    }
    duplicates <- duplicates[,c(1:M)]
  }

  jack.boot.values <- apply(X=duplicates, MARGIN=2L, FUN=jack.boot, xstar=t, f=fn)
  phitilde <- (N*t0 - (N-m)*jack.boot.values)/m - t0

  jack.boot.se <- sqrt(m/(N-m)/M *
                       apply(X=phitilde, MARGIN=1L, FUN=function(x) {sum(x^2)})
                       )
  return(jack.boot.se)
}

#' jab.cf
#'
#' apply jackknife after bootstrap to an cf object
#' and compute errors accordingly
#' 
#' @param cf object of type \link{cf}
#' @param m integer. block length
#'
#' @return
#' Returns an object of class `cf`, see \link{cf}, with
#' Jackknife samples added accordingly.
#' 
#' @export
jab.cf <- function(cf, m = 1) {
  stopifnot(inherits(cf, 'cf'))
  stopifnot(inherits(cf, 'cf_boot'))
  stopifnot(cf$cf.tsboot$sim == "fixed")
  stopifnot(cf$resampling_method == 'bootstrap')

  old_seed <- swap_seed(cf$seed)
  ## the resampling block indices
  cf$blockind <- hadron:::boot_ts_array(n=cf$cf.tsboot$n, n.sim=cf$cf.tsboot$n.sim,
                                        R=cf$boot.R, l=cf$boot.l, sim=cf$sim, endcorr=cf$cf.tsboot$endcorr)
  restore_seed(old_seed)

  cf$jack.boot.se <- jab(t=cf$cf.tsboot$t, t0=cf$tsboot.se, starts=cf$blockind$starts, m=m, fn=sd)
  if( has_icf(cf) ){
    # no randomness here so no seed setting required
    cf$ijack.boot.se <- jab(t=cf$icf.tsboot$t, t0=cf$itsboot.se, starts=cf$blockind$starts, m=m, fn=sd)
  }
  
  return(invisible(cf))
}

#' jab.cf.derived
#'
#' apply jackknife after bootstrap to an derived cf object
#' and compute errors accordingly
#'
#' @param cf object of type \link{cf}
#' @param m integer. block length
#'
#' @return
#' Returns a numeric vector with the jackknife estimates of
#' standard error.
#' 
#' @export
jab.cf.derived <- function(cf, m=1) {
  if(cf$cf$cf.tsboot$sim != "fixed") {
    stop("JAB only implemented for 'sim=fixed' at the moment")
  }

  if(is.null(cf$cf$blockind)) {
    old_seed <- swap_seed(cf$seed)
    cf$cf$blockind <- hadron:::boot_ts_array(n=cf$cf$cf.tsboot$n, n.sim=cf$cf$cf.tsboot$n.sim,
                                             R=cf$cf$boot.R, l=cf$cf$boot.l, sim=cf$cf$sim, endcorr=cf$cf$cf.tsboot$endcorr)
    restore_seed(old_seed) 
  }
  jack.boot.se <- jab(t=cf$t[,c(1:length(cf$se))], t0=cf$se, starts=cf$cf$blockind$starts, m=m, fn=sd)

  
  return(jack.boot.se)
}

jab.matrixfit <- function(cf, m=1) {
  if(!any(class(cf) == "matrixfit")) {
    stop("bootstrap.cf requires an object of class 'matrixfit' as input! Aborting!\n")
  }
  cf$jack.boot.se <- jab.cf.derived(cf=cf, m=m)
  return(invisible(cf))
}


jab.effectivemass <- function(cf, m=1) {
  if(!any(class(cf) == "effectivemass")) {
    stop("bootstrap.cf requires an object of class 'matrixfit' as input! Aborting!\n")
  }
  cf$jack.boot.se <- jab.cf.derived(cf=cf, m=m)
  return(invisible(cf))
}

jab.effectivemassfit <- function(cf, m=1) {
  if(!any(class(cf) == "effectivemassfit")) {
    stop("bootstrap.cf requires an object of class 'matrixfit' as input! Aborting!\n")
  }
  cf$jack.boot.se <- jab.cf.derived(cf=cf, m=m)
  return(invisible(cf))
}
