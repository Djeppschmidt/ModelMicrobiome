### rarefaction tool that allows replacement
### is nearly the same code as rrarefy... just added the argument "replace" to pass to sample

rarefy<-
  function (x, sample, replace) #added replace argument as input
  {
    if (!identical(all.equal(x, round(x)), TRUE))
      stop("function is meaningful only for integers (counts)")
    x <- as.matrix(x)
    if (ncol(x) == 1)
      x <- t(x)
    if (length(sample) > 1 && length(sample) != nrow(x))
      stop(gettextf("length of 'sample' and number of rows of 'x' do not match"))
    sample <- rep(sample, length = nrow(x))
    colnames(x) <- colnames(x, do.NULL = FALSE)
    nm <- colnames(x)
    if (any(rowSums(x) < sample))
      warning("Some row sums < 'sample' and are not rarefied")
    for (i in 1:nrow(x)) {
      if (sum(x[i, ]) <= sample[i])
        next
      row <- sample(rep(nm, times = x[i, ]), sample[i], replace) #use replace argument
      row <- table(row)
      ind <- names(row)
      x[i, ] <- 0
      x[i, ind] <- row
    }
    x
  }
