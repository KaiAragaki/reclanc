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
  selected = select_clanc(d.k = train$tStats, d.k.ord = list("d.k.rnks" = train$tStatsOrderIndices,
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
  p.sd = pooled_sd_clanc(data, ID)

  ## class and overall centroids
  cntrd.k = scale(data %*% ID, center = F, scale = nn)
  cntrd.o = drop(data %*% rep(1 / n, n))

  ## form statistics and order them
  d.k = scale((cntrd.k - cntrd.o) / p.sd, center = F, scale = m.k)
  d.k.ord = order_stats_clanc(d.k = d.k)

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
#' # feature-set sizes. By default, cv_clanc() will assess the classifiers built
#' # using 1, 2, 3, ..., 10 features. In this example, an estimated
#' # 100% accuracy is attained with as few as 5 features per class.
#'
#' # Assuming we're happy with using 5 features per class as the basis for our
#' # model, we can now train and build our final classifier.
#'
run_ClaNC <- function(data,
                      groups,
                      active.features = 20,
                      est.num = 20,
                      select.features = 10,
                      file.prefix = "ClaNC",
                      skip.est = FALSE) {

  Y_train <- as.matrix(data)
  idd <- as.numeric(factor(groups))
  gene_names <- rownames(Y_train)
  class_names <- 1:length(unique(idd))

  if (!skip.est){
    cv_total <- NULL
    cv_out <- NULL
    for (i in 1:est.num){
      cv_out = cv_clanc(Y_train, idd, active=1:active.features, prior="class")
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
