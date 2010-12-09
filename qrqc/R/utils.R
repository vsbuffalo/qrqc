## utils.R - Utility functions for creating parts of graphics, etc

binned2quantilefunc <-
# Assumes names are x values
function(x) {
  y <- as.numeric(x[x!=0])
  x <- as.integer(names(x)[x!=0])
  t <- approxfun(cumsum(y)/sum(y), x, yleft=min(x), yright=max(x))
  return(t)
}

binned2boxplot <-
# Extract useful quantiles for a boxplot from an empirical qunatile
# function.
function(x) {
  f <- binned2quantilefunc(x)
  out <- c(ymin=f(0),
           alt.lower=f(0.1),
           lower=f(0.25),
           middle=f(0.5),
           upper=f(0.75),
           alt.upper=f(0.9),
           ymax=f(1))
  return(out)
}

qualMCLowess <-
# Monte Carlo Lowess, or lowess on binned quality data
function(obj, n=100, f=1/6) {
  set.seed(0)
  vals <- as.numeric(names(obj@qual.freqs[1, -1]))
  binsample <- function(x) {
    sample(vals, n, prob=x/sum(x), replace=TRUE)
  }

  # Use binsample() to sample from the binned quality frequency data, then
  # do some data munging to get it in a usable format for lowess.
  d <- local({
    s <- apply(obj@qual.freqs[, -1], 1, binsample)
    tmp <- t(s)
    nr <- nrow(tmp)
    dim(tmp) <- c(ncol(tmp)*nrow(tmp), 1)
    tmp <- as.data.frame(cbind(1:nr, tmp))    
    colnames(tmp) <- c('position', 'quality')
    return(tmp)
  })
  l <- lowess(d$position, d$quality, f=f)
  return(l)
}
