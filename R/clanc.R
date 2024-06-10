# Everything here comes from https://github.com/naikai/sake or is a modification of it.

#' @param data expression data. (m x n) matrix of class numeric.
#' @param id class id's. n vector of class numeric.
#' @param prior either a vector of length p, "equal", or "class". if "equal",
#'   then equal probabilities will be used for each class. if "class", then the
#'   proportions of each class in the training data will be used as the prior.
#' @param active how many active features to consider? can either be a single
#'   number or a vector containing a range of values to consider.
cvClanc <- function(data, id, priors = "equal", active = 1:10) {
  cvIdx <- balancedFolds(id, 5)
  nrow_data <- nrow(data)
  n <- ncol(data)
  p <- length(unique(id))
  nn <- as.numeric(table(id))
  folds <- length(cvIdx)

  class_priors <- make_class_priors(priors)

  d <- ifelse(is.matrix(active), nrow(active), length(active))

  cv.error <- array(rep(0, d * folds * p), dim = c(d, folds, p))

  cv.err.cnt.cls <- matrix(NA, nrow = d, ncol = p)
  cv.err.prpn.cls <- matrix(NA, nrow = d, ncol = p)
  n.features.cls <- matrix(NA, nrow = d, ncol = p)
  n.features.ttl <- rep(NA, d)

  ID <- model.matrix(~ factor(id) - 1)
  dimnames(ID) <- list(NULL, names(nn))

  ## cross validation
  cat("CV:")
  for(i in seq_len(folds)) {
    cat(i)

    ## form initial statistics
    v <- length(cvIdx[[i]])
    X <- data[, -cvIdx[[i]]]
    Y <- data[, cvIdx[[i]]]
    jd <- id[-cvIdx[[i]]]
    truth <- id[cvIdx[[i]]]
    JD <- ID[-cvIdx[[i]], ]

    mm <- table(jd)
    m.k <- sqrt(1 / mm - 1 / (n - v))

    ## pooled standard deviations
    p.sd <- pooledSDClanc(X, JD)

    ## class- and overall-centroids
    cntrd.k <- scale(X %*% JD, center = FALSE, scale = mm)
    cntrd.o <- drop(X %*% rep(1 / (n - v), n - v))

    ## form statistics and order them
    d.k <- scale((cntrd.k - cntrd.o) / p.sd, center = FALSE, scale = m.k)
    d.k.ord <- orderStatsClanc(d.k = d.k)

    ## select genes, update inactive centroid components
    for(j in seq_len(d)) {
      # FIXME
      aa <- ifelse(is.matrix(active), active[j, ], active[j])
      selected <- selectClanc(d.k = d.k, d.k.ord = d.k.ord, active = aa)
      active.idx <- seq_len(nrow_data)[drop(selected %*% rep(1, p)) != 0]
      cntrds <- cntrd.k[active.idx, ]
      for(k in seq_len(p))
        cntrds[selected[active.idx, k] == 0, k] = cntrd.o[active.idx][selected[active.idx, k] == 0]

      ## classify test sample and assess error status
      for(k in seq_len(v)) {
        dd <- distClanc(data = Y[active.idx, k], cntrds = cntrds, sd = p.sd[active.idx], prior = class_priors)

        if(match(min(dd), dd) != truth[k])
          cv.error[j, i, truth[k]] = cv.error[j, i, truth[k]] + 1
      }
    }
  }

  ## record numbers and proportions of errors
  for(i in seq_len(p)) {
    if(d > 1)
      cv.err.cnt.cls[, i] <- apply(cv.error[, , i], 1, sum)
    else
      cv.err.cnt.cls[, i] <- sum(cv.error[, , i])
    cv.err.prpn.cls[, i] <- cv.err.cnt.cls[, i] / nn[i]
  }
  cv.err.cnt.ttl <- apply(cv.err.cnt.cls, 1, sum)
  cv.err.prpn.ttl <- cv.err.cnt.ttl / n

  cat("\n")
  list(classErrors = cv.err.prpn.cls, overallErrors = cv.err.prpn.ttl, prior = class_priors)
}

make_class_priors <- function(priors, id, data) {
  n_classes <- length(unique(id))
  n_samples <- ncol(data)
  samples_per_class <- as.numeric(table(id))
  class_proportions <- n_samples / samples_per_class

  if (is.numeric(priors)) {
    stopifnot(length(priors) == n_classes, sum(priors) == 1, all(priors >= 0))
    return(priors)
  }

  stopifnot(prior %in% c("equal", "class"))
  if (priors == "equal") return(rep(1 / n_classes, n_classes))
  if (priors == "class") return(class_proportions)
}

buildClanc <- function(data, id, cNames, train, active) {
  m = nrow(data)
  n = ncol(data)
  p = length(unique(id))
  nn = as.numeric(table(id))
  m.k = sqrt(1 / nn - 1 / n)
  geneNames = train$geneNames
  classNames = cNames
  prior = train$prior


  pi.k <- make_class_priors(priors = train$prior, id = id, data = data)

  ## select genes, update inactive centroid components
  selected = selectClanc(d.k = train$tStats, d.k.ord = list("d.k.rnks" = train$tStatsOrderIndices,
    "d.k.o" = train$tStatsOrdered, "no.ties" = train$tStatsNoTies), active = active)
  active.idx = (1:m)[selected %*% rep(1, p) != 0]

  cntrds = train$classCntrds[active.idx, ]
  for(k in 1:p)
    cntrds[selected[active.idx, k] == 0, k] = train$overallCntrd[active.idx][selected[active.idx, k] == 0]

  ## prepare output
  out = list("geneNames" = geneNames[active.idx], "pooledSD" = train$pooledSD[active.idx], "cntrds" = cntrds,
    "prior" = pi.k, "classNames" = classNames, "classID" = seq(p))

  ## assess training error
  error = testClanc(data = data, id = id, geneNames = geneNames, fit = out)
  out$classError = error / nn
  out$overallError = sum(error) / n

  return(out)
}

predictClanc <- function(data, geneNames, fit) {
  n = ncol(data)
  cntrds = fit$cntrds

  active.idx = match(fit$geneNames, geneNames)
  sd = fit$pooledSD
  prior = fit$prior

  if(any(is.na(active.idx)))
    stop("Gene names do not match those in classifier.")

  cls = rep(NA, n)
  for(i in 1:n) {
    dd = distClanc(data = data[active.idx, i], cntrds = cntrds, sd = sd, prior = prior)
    cls[i] = match(min(dd), dd)
  }

  return(cls)
}

selectClanc <- function(d.k, d.k.ord, active) {
  m = nrow(d.k)
  p = ncol(d.k)
  d.k.rnks = d.k.ord$d.k.rnks
  d.k.o = d.k.ord$d.k.o
  no.ties = d.k.ord$no.ties

  if(length(active) != p)
      active = rep(active, p)
  if(any(active >= m))
      selected = matrix(1, nrow = m, ncol = p)
  else {
    if(any((avail <- t(no.ties) %*% rep(1, m)) < active)) {
      active[avail < active] = avail[avail < active]

      warning("Not enough unique genes in at least one class.")
    }

    selected = matrix(0, nrow = m, ncol = p)

    jdx = vector("list", p)
    delta = rep(NA, p)
    for(i in 1:p) {
      jdx[[i]] = d.k.rnks[no.ties[, i], i][1:active[i]]
      delta[i] = d.k.o[no.ties[, i], i][active[i] + 1]
    }
    cls = rep(1:p, active)

    while(any((kdx <- table(unlist(jdx))) > 1)) {
      nms = as.numeric(names(kdx))
      for(j in as.numeric(names(kdx))[kdx > 1]) {
        ldx = cls[unlist(jdx) == j]

        ii = match(max(abs(d.k[j, ldx])), abs(d.k[j, ldx]))

        for(k in ldx[-ii]) {
          active[k] = active[k] + 1
          jdx[[k]] = c(jdx[[k]][-match(j, jdx[[k]])], d.k.rnks[no.ties[, k], k][active[k]])
          delta[k] = d.k.o[no.ties[, k], k][active[k] + 1]
        }
      }
    }

    for(i in 1:p)
      selected[jdx[[i]], i] = 1
  }

  return(selected)
}

testClanc <- function(data, id, geneNames, fit) {
  m = nrow(data)
  n = ncol(data)
  p = ncol(fit$cntrds)
  cntrds = fit$cntrds

  active.idx = match(fit$geneNames, geneNames)
  sd = fit$pooledSD
  prior = fit$prior

  if(any(is.na(active.idx)))
    stop("Gene names do not match those in classifier.")

  error = matrix(0, nrow = n, ncol = p)
  cls = predictClanc(data = data, geneNames = geneNames, fit = fit)

  for(i in 1:n)
    if(cls[i] != id[i])
      error[i, id[i]] = 1

  return(drop(apply(error, 2, sum)))
}

#' @param data expression data. (m x n) matrix of class numeric.
#' @param id class id's. n vector of class numeric.
#' @param geneNames gene names. m vector of class character.
trainClanc <- function(data, id, geneNames, prior = "equal") {
  m = nrow(data)
  n = ncol(data)
  p = length(unique(id))
  nn = as.numeric(table(id))
  m.k = sqrt(1 / nn - 1 / n)

  ID = model.matrix(~ factor(id) - 1)

  ## pooled standard deviations
  p.sd = pooledSDClanc(data, ID)

  ## class and overall centroids
  cntrd.k = scale(data %*% ID, center = F, scale = nn)
  cntrd.o = drop(data %*% rep(1 / n, n))

  ## form statistics and order them
  d.k = scale((cntrd.k - cntrd.o) / p.sd, center = F, scale = m.k)
  d.k.ord = orderStatsClanc(d.k = d.k)

  list(
    pooledSD = p.sd,
    classCntrds = cntrd.k,
    overallCntrd = cntrd.o,
    tStats = d.k,
    tStatsOrdered = d.k.ord$d.k.o,
    tStatsOrderIndices = d.k.ord$d.k.rnks,
    tStatsNoTies = d.k.ord$no.ties,
    geneNames = geneNames,
    prior = prior
  )
}

orderStatsClanc <- function(d.k) {
  m = nrow(d.k)
  p = ncol(d.k)

  d.k.rnks = apply(d.k, 2, function(x) { order(abs(x), decreasing = T) })
  d.k.o = matrix(d.k[as.numeric(d.k.rnks) + m * rep(0:(p - 1), each = m)], nrow = m)
  no.ties = apply(d.k.o, 2, function(x) { !duplicated(x) })

  list(d.k.rnks = d.k.rnks, d.k.o = d.k.o, no.ties = no.ties)
}

pooledSDClanc <- function(X, ID) {
  m = nrow(X)
  n = ncol(X)
  p = ncol(ID)
  nn = drop(t(ID) %*% rep(1, n))

  mn = t(t(X %*% ID) / nn)
  dif2 = (X - mn %*% t(ID)) ^ 2

  sqrt(drop(dif2 %*% rep(1 / (n - p), n)))
}

balancedFolds <- function(y, nfolds) {

  totals = table(y)
  fmax = max(totals)
  nfolds = min(nfolds, fmax)
  folds = as.list(seq(nfolds))
  yids = split(seq(y), y)
  bigmat = matrix(NA, ceiling(fmax / nfolds) * nfolds, length(totals))
  for(i in seq(totals))
    bigmat[seq(totals[i]), i] = sample(yids[[i]])
  smallmat = matrix(bigmat, nrow = nfolds)
  smallmat = permute_rows(t(smallmat))
  res = vector("list", nfolds)
  for(j in 1:nfolds) {
    jj = !is.na(smallmat[, j])
    res[[j]] = smallmat[jj, j]
  }

  res
}

permute_rows <- function(x) {
  t(apply(x, 1, sample))
}

distClanc <- function(data, cntrds, sd, prior) {
  vv = sd ^ 2
  pi.k = prior

  if(length(vv) > 1)
    dd = drop(t(cntrds ^ 2) %*% (1 / vv)) - 2 * drop(t(data * cntrds) %*% (1 / vv)) - 2 * log(pi.k)
  else
    dd = drop(cntrds ^ 2 / vv) - 2 * drop(data * cntrds / vv) - 2 * log(pi.k)

  return(dd)
}

#' Get ClaNC Group info
#'
#' Extract which group each point belongs to
#' @param Build_out
#' @keywords ClaNC
#' @export
#' @examples
#' get_ClaNC_group()
get_ClaNC_group <- function(Build_out){
  idx <- apply(Build_out$cntrds, 1, function(x) which.max(abs(x-Mode(x))))
  summary <- sapply(1:max(idx), function(x) names(idx[idx==x]))
  colnames(summary) <- paste0("Group", 1:ncol(summary))
  summary
}

#' @param data expression data, row(genes), column(samples)
#' @param groups assigned groups for each column(samples) in data
#' @param active.features how many features we want to randomly select to test
#'   cv.errors
#' @param est.num how many times we want to run the test for cv.errors
#' @param select.features how many features(genes) we want to select from each
#'   group for the final results
#' @param skip.est do you want to skip estimate cv errors and go straight to
#'   extract features from each group
#' @examples
#'
#' # Load the example data sets in the "data" directory. These are simulated
#' # microarray data sets, consisting of 10 arrays each from 4 classes; the
#' # first 10 arrays are from class 1, the second 10 from class 2, etc. And each
#' # array has 1500 genes on it. The "train" data set is intended for training
#' # purposes; these are the data you would use in training and building your
#' # classifier. The "test" data are intended for testing a classifier, using
#' # data that was not used in the classifier building process. You may or may
#' # not have test data. If you intend to use all of your samples to build the
#' # classifier and rely on cross-validation to estimate the resulting error
#' # rates, you need not have *any* test data.
#'
#' Y_train = read.delim("data/trainExample.txt", header = T)
#' Y_test = read.delim("data/testExample.txt", header = T)
#' rownames(Y_train) = Y_train[, 1]
#' rownames(Y_test) = Y_test[, 1]
#' Y_train = as.matrix(Y_train[, -1])
#' Y_test = as.matrix(Y_test[, -1])
#'
#' # Create vector for class membership. For these example data, we just need 10
#' # ones followed by 10 twos, then 10 threes, then 10 fours.
#'
#' id = rep(1:4, each = 10)
#' #> Levels: Basal Her2 LumA LumB Normal
#'
#' # The gene names are just the row names from the data matrices, the same for
#' # both the example train and test data sets.
#' gene_names = rownames(Y_train)
#'
#' # We'll just set the class names to be the numbers 1, 2, 3, 4.
#' class_names = 1:4
#'
#' # Specify how many active features(genes) to select from each group Now ready
#' # to use cross-validation to estimate the error rates for classifiers of
#' # different sizes (different numbers of genes used in building the
#' # classifier). THEN View the estimated error rates associated with different
#' # feature-set sizes. By default, cvClanc() will assess the classifiers built
#' # using 1, 2, 3, ..., 10 features. In this example, an estimated
#' # 100% accuracy is attained with as few as 5 features per class.
#'
#' # Assuming we're happy with using 5 features per class as the basis for our
#' # model, we can now train and build our final classifier.
#'
run_ClaNC <- function(data, groups, ColSideColors, active.features=20, est.num=20, select.features=10, file.prefix="ClaNC", skip.est=FALSE){

  Y_train <- as.matrix(data)
  idd <- as.numeric(factor(groups))
  gene_names <- rownames(Y_train)
  class_names <- 1:length(unique(idd))

  if (!skip.est){
    cv_total <- NULL
    cv_out <- NULL
    for (i in 1:est.num){
      cv_out = cvClanc(Y_train, idd, active=1:active.features, prior="class")
      cv_total = rbind(cv_total, cbind(1:active.features, cv_out$overallErrors))
    }
    colnames(cv_total) <- c("NumFeatures", "CV_Error")
    file.prefix <- paste0(file.prefix, ".estnum", est.num, ".act", active.features)
    cv_total <- as.data.frame(cv_total)
    p <- ggplot(data=cv_total, aes(x=factor(NumFeatures), y=CV_Error))
    p <- p + geom_boxplot(aes(fill = factor(NumFeatures)))
    p <- p + xlab("Number of Features")
    pdf(paste0(file.prefix, ".pdf"), height=10, width=16)
    print(p)
    dev.off()
  }

  file.prefix <- paste0(file.prefix, ".sel", select.features)
  train_out = trainClanc(Y_train, idd, gene_names)
  build_out = buildClanc(Y_train, idd, class_names, train_out, active = select.features)
  write.table(build_out$geneNames, paste0(file.prefix, ".genelist.txt"), row.names = F, quote=F)
  clanc.group <- get_ClaNC_group(build_out)
  write.csv(clanc.group, paste0(file.prefix, ".group_features.csv"), quote=F, row.names=F)
  build_out
}
