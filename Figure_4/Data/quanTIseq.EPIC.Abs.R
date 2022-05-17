dec1 <- function (mix.mat, arrays = FALSE, signame = "TIL10", tumor = FALSE, 
                  mRNAscale = TRUE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned"){
  deconvolute_quantiseq.default(mix.mat, arrays, signame, tumor, 
                                mRNAscale, method, btotalcells, rmgenes )  
}



deconvolute_quantiseq.default2 <- function (mix.mat, arrays = FALSE, signame = "TIL10", tumor = FALSE, 
          mRNAscale = TRUE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned") 
{
  message("\nRunning quanTIseq deconvolution module\n")
  if (rmgenes == "unassigned" && arrays == TRUE) {
    rmgenes <- "none"
  }
  else if (rmgenes == "unassigned" && arrays == FALSE) {
    rmgenes <- "default"
  }
  listsig <- c("TIL10")
  if (signame %in% listsig) {
    sig.mat.file <- system.file("extdata", "quantiseq", paste0(signame, 
                                                               "_signature.txt"), package = "immunedeconv", mustWork = TRUE)
    mRNA.file <- system.file("extdata", "quantiseq", paste0(signame, 
                                                            "_mRNA_scaling.txt"), package = "immunedeconv", mustWork = TRUE)
    fileab <- system.file("extdata", "quantiseq", paste0(signame, 
                                                         "_TCGA_aberrant_immune_genes.txt"), package = "immunedeconv", 
                          mustWork = TRUE)
    if (rmgenes == "default") {
      filerm <- system.file("extdata", "quantiseq", paste0(signame, 
                                                           "_rmgenes.txt"), package = "immunedeconv", mustWork = TRUE)
    }
    else if (rmgenes == "path") {
      filerm <- system.file("extdata", "quantiseq", paste0(signame, 
                                                           "rmgenes.txt"), package = "immunedeconv", mustWork = TRUE)
    }
  }
  else {
    sig.mat.file <- paste0(signame, "_signature.txt")
    mRNA.file <- paste0(signame, "_mRNA_scaling.txt")
  }
  if (is.numeric(mix.mat[[1, 1]]) != TRUE) {
    stop("Wrong input format for the mixture matrix! Please follow the instructions of the documentation.")
  }
  sig.mat <- read.table(sig.mat.file, header = TRUE, sep = "\t", 
                        row.names = 1)
  if (mRNAscale) {
    mRNA <- read.table(mRNA.file, sep = "\t", header = FALSE, 
                       stringsAsFactors = FALSE)
    colnames(mRNA) <- c("celltype", "scaling")
    mRNA <- as.vector(as.matrix(mRNA$scaling[match(colnames(sig.mat), 
                                                   mRNA$celltype)]))
  }
  else {
    mRNA <- rep(1, ncol(sig.mat))
  }
  message(paste0("Gene expression normalization and re-annotation (arrays: ", 
                 arrays, ")\n"))
  mix.mat <- immunedeconv:::fixMixture(mix.mat, arrays = arrays)
  if (rmgenes != "none") {
    if (signame %in% listsig) {
      lrmgenes <- as.vector(read.table(filerm, header = FALSE, 
                                       sep = "\t")[, 1])
      n1 <- nrow(sig.mat)
      sig.mat <- sig.mat[!rownames(sig.mat) %in% lrmgenes, 
                         , drop = FALSE]
      n2 <- nrow(sig.mat)
      message(paste0("Removing ", n1 - n2, " noisy genes\n"))
    }
  }
  if (tumor) {
    if (signame %in% listsig) {
      abgenes <- as.vector(read.table(fileab, header = FALSE, 
                                      sep = "\t")[, 1])
      n1 <- nrow(sig.mat)
      sig.mat <- sig.mat[!rownames(sig.mat) %in% abgenes, 
                         , drop = FALSE]
      n2 <- nrow(sig.mat)
      message(paste0("Removing ", n1 - n2, " genes with high expression in tumors\n"))
    }
  }
  ns <- nrow(sig.mat)
  us <- length(intersect(rownames(sig.mat), rownames(mix.mat)))
  perc <- round(us * 100/ns, digits = 2)
  message(paste0("Signature genes found in data set: ", us, 
                 "/", ns, " (", perc, "%)\n"))
  message(paste0("Mixture deconvolution (method: ", method, 
                 ")\n"))
  results1 <- immunedeconv:::quanTIseq(sig.mat, mix.mat, scaling = mRNA, method = method)
  if ("Tregs" %in% colnames(sig.mat) && "T.cells.CD4" %in% 
      colnames(sig.mat) && method %in% c("lsei")) {
    minTregs <- 0.02
    i <- which(colnames(sig.mat) == "T.cells.CD4")
    results2 <- immunedeconv:::quanTIseq(sig.mat[, -i], mix.mat, scaling = mRNA[-i], 
                          method = method)
    ind <- which(results1[, "Tregs"] < minTregs)
    if (length(ind) > 0) {
      results1[ind, "Tregs"] <- (results2[ind, "Tregs"] + 
                                   results1[ind, "Tregs"])/2
      results1[ind, "T.cells.CD4"] <- pmax(0, results1[ind, 
                                                       "T.cells.CD4"] - (results2[ind, "Tregs"] + results1[ind, 
                                                                                                           "Tregs"])/2)
    }
  }
  results <- results1
  results <- results/apply(results, 1, sum)
  message("Deconvolution sucessful!")
  DCres <- results
  if (btotalcells == TRUE) {
    celldens <- data.frame(celldensities(DCres))
    celldens <- cbind(rownames(celldens), celldens)
    colnames(celldens)[1] <- "Sample"
    return(celldens)
  }
  else {
    results = data.frame(results)
    results <- cbind(rownames(results), results)
    colnames(results)[1] <- "Sample"
    return(list(Prop=results, Abs=results1))
  }
}

MyEPIC <-function (bulk, reference = NULL, mRNA_cell = NULL, mRNA_cell_sub = NULL, 
          sigGenes = NULL, scaleExprs = TRUE, withOtherCells = TRUE, 
          constrainedSum = TRUE, rangeBasedOptim = FALSE) 
{
  if (!is.matrix(bulk) && !is.data.frame(bulk)) 
    stop("'bulk' needs to be given as a matrix or data.frame")
  with_w <- TRUE
  if (is.null(reference)) {
    reference <- EPIC::TRef
  }
  else if (is.character(reference)) {
    if (reference %in% prebuiltRefNames) {
      reference <- get(reference, pos = "package:EPIC")
    }
    else stop("The reference, '", reference, "' is not part of the allowed ", 
              "references:", paste(prebuiltRefNames, collapse = ", "))
  }
  else if (is.list(reference)) {
    refListNames <- names(reference)
    if ((!all(c("refProfiles", "sigGenes") %in% refListNames)) || 
        (("refProfiles" %in% refListNames) && !is.null(sigGenes))) 
      stop("Reference, when given as a list needs to contain at least the ", 
           "fields 'refProfiles' and 'sigGenes' (sigGenes could also be ", 
           "given as input to EPIC instead)")
    if (!is.matrix(reference$refProfiles) && !is.data.frame(reference$refProfiles)) 
      stop("'reference$refProfiles' needs to be given as a matrix or data.frame")
    if (!("refProfiles.var" %in% refListNames)) {
      warning("'refProfiles.var' not defined; using identical weights ", 
              "for all genes")
      with_w <- FALSE
    }
    else if (!is.matrix(reference$refProfiles.var) && !is.data.frame(reference$refProfiles.var)) {
      stop("'reference$refProfiles.var' needs to be given as a matrix or ", 
           "data.frame when present.")
    }
    else if (!identical(dim(reference$refProfiles.var), dim(reference$refProfiles)) || 
             !identical(dimnames(reference$refProfiles.var), dimnames(reference$refProfiles))) 
      stop("The dimensions and dimnames of 'reference$refProfiles' and ", 
           "'reference$refProfiles.var' need to be the same")
  }
  else {
    stop("Unknown format for 'reference'")
  }
  bulk <- EPIC:::merge_duplicates(bulk, in_type = "bulk samples")
  refProfiles <- EPIC:::merge_duplicates(reference$refProfiles, in_type = "reference profiles")
  if (with_w) {
    refProfiles.var <- EPIC:::merge_duplicates(reference$refProfiles.var, 
                                        warn = F)
  }
  else {
    refProfiles.var <- 0
  }
  nSamples <- NCOL(bulk)
  samplesNames <- colnames(bulk)
  if (is.null(samplesNames)) {
    samplesNames <- 1:nSamples
    colnames(bulk) <- samplesNames
  }
  nRefCells <- NCOL(refProfiles)
  refCellsNames <- colnames(refProfiles)
  bulk_NA <- apply(is.na(bulk), MARGIN = 1, FUN = all)
  if (any(bulk_NA)) {
    warning(sum(bulk_NA), " genes are NA in all bulk samples, removing these.")
    bulk <- bulk[!bulk_NA, ]
  }
  bulkGenes <- rownames(bulk)
  refGenes <- rownames(refProfiles)
  commonGenes <- intersect(bulkGenes, refGenes)
  if (is.null(sigGenes)) 
    sigGenes <- unique(reference$sigGenes)
  sigGenes <- sigGenes[sigGenes %in% commonGenes]
  nSigGenes <- length(sigGenes)
  if (nSigGenes < nRefCells) 
    stop("There are only ", nSigGenes, " signature genes", 
         " matching common genes between bulk and reference profiles,", 
         " but there should be more signature genes than reference cells")
  if (scaleExprs) {
    if (length(commonGenes) < 2000) 
      warning("there are few genes in common between the bulk samples and ", 
              "reference cells:", length(commonGenes), ", so the data scaling ", 
              "might be an issue")
    bulk <- EPIC:::scaleCounts(bulk, sigGenes, commonGenes)$counts
    temp <- EPIC:::scaleCounts(refProfiles, sigGenes, commonGenes)
    refProfiles <- temp$counts
    if (with_w) 
      refProfiles.var <- EPIC:::scaleCounts(refProfiles.var, sigGenes, 
                                     normFact = temp$normFact)$counts
  }
  else {
    bulk <- bulk[sigGenes, ]
    refProfiles <- refProfiles[sigGenes, ]
    if (with_w) 
      refProfiles.var <- refProfiles.var[sigGenes, ]
  }
  if (is.null(mRNA_cell)) 
    mRNA_cell <- EPIC::mRNA_cell_default
  if (!is.null(mRNA_cell_sub)) {
    if (is.null(names(mRNA_cell_sub)) || !is.numeric(mRNA_cell_sub)) 
      stop("When mRNA_cell_sub is given, it needs to be a named numeric vector")
    mRNA_cell[names(mRNA_cell_sub)] <- mRNA_cell_sub
  }
  minFun <- function(x, A, b, w) {
    return(sum((w * (A %*% x - b)^2), na.rm = TRUE))
  }
  minFun.range <- function(x, A, b, A.var) {
    val.max <- (A + A.var) %*% x - b
    val.min <- (A - A.var) %*% x - b
    cErr <- rep(0, length(b))
    outOfRange <- (sign(val.max) * sign(val.min) == 1)
    cErr[outOfRange] <- pmin(abs(val.max[outOfRange]), abs(val.min[outOfRange]))
    return(sum(cErr, na.rm = TRUE))
  }
  if (with_w && !rangeBasedOptim) {
    w <- rowSums(refProfiles/(refProfiles.var + 1e-12), na.rm = TRUE)
    med_w <- stats::median(w[w > 0], na.rm = TRUE)
    w[w > 100 * med_w] <- 100 * med_w
  }
  else w <- 1
  if (withOtherCells) {
    cMin <- 0
  }
  else {
    cMin <- 0.99
  }
  cMax <- 1
  ui <- diag(nRefCells)
  ci <- rep(0, nRefCells)
  if (constrainedSum) {
    ui <- rbind(ui, rep(1, nRefCells), rep(-1, nRefCells))
    ci <- c(ci, cMin, -cMax)
  }
  cInitProp <- (min(1, cMax) - 1e-05)/nRefCells
  tempPropPred <- lapply(1:nSamples, FUN = function(cSample) {
    b <- bulk[, cSample]
    if (!rangeBasedOptim) {
      fit <- stats::constrOptim(theta = rep(cInitProp, 
                                            nRefCells), f = minFun, grad = NULL, ui = ui, 
                                ci = ci, A = refProfiles, b = b, w = w)
    }
    else {
      fit <- stats::constrOptim(theta = rep(cInitProp, 
                                            nRefCells), f = minFun.range, grad = NULL, ui = ui, 
                                ci = ci, A = refProfiles, b = b, A.var = refProfiles.var)
    }
    fit$x <- fit$par
    if (!withOtherCells) 
      fit$x <- fit$x/sum(fit$x, na.rm = TRUE)
    b_estimated <- refProfiles %*% fit$x
    if (nSigGenes > 2) {
      suppressWarnings(corSp.test <- stats::cor.test(b, 
                                                     b_estimated, method = "spearman"))
      corPear.test <- stats::cor.test(b, b_estimated, method = "pearson")
    }
    else {
      corSp.test <- corPear.test <- list()
      corSp.test$estimate <- corSp.test$p.value <- corPear.test$estimate <- corPear.test$p.value <- NA
    }
    regLine <- stats::lm(b_estimated ~ b)
    regLine_through0 <- stats::lm(b_estimated ~ b + 0)
    if (!rangeBasedOptim) {
      rmse_pred <- sqrt(minFun(x = fit$x, A = refProfiles, 
                               b = b, w = w)/nSigGenes)
      rmse_0 <- sqrt(minFun(x = rep(0, nRefCells), A = refProfiles, 
                            b = b, w = w)/nSigGenes)
    }
    else {
      rmse_pred <- sqrt(minFun.range(x = fit$x, A = refProfiles, 
                                     b = b, A.var = refProfiles.var)/nSigGenes)
      rmse_0 <- sqrt(minFun.range(x = rep(0, nRefCells), 
                                  A = refProfiles, b = b, A.var = refProfiles.var)/nSigGenes)
    }
    gof <- data.frame(fit$convergence, ifelse(is.null(fit$message), 
                                              "", fit$message), rmse_pred, rmse_0, corSp.test$estimate, 
                      corSp.test$p.value, corPear.test$estimate, corPear.test$p.value, 
                      regLine$coefficients[2], regLine$coefficients[1], 
                      regLine_through0$coefficients[1], sum(fit$x), stringsAsFactors = FALSE)
    return(list(mRNAProportions = fit$x, fit.gof = gof))
  })
  mRNAProportions <- do.call(rbind, lapply(tempPropPred, function(x) x$mRNAProportions))
  dimnames(mRNAProportions) <- list(samplesNames, refCellsNames)
  fit.gof <- do.call(rbind, lapply(tempPropPred, function(x) x$fit.gof))
  dimnames(fit.gof) <- list(samplesNames, c("convergeCode", 
                                            "convergeMessage", "RMSE_weighted", "Root_mean_squared_geneExpr_weighted", 
                                            "spearmanR", "spearmanP", "pearsonR", "pearsonP", "regline_a_x", 
                                            "regline_b", "regline_a_x_through0", "sum_mRNAProportions"))
  if (any(fit.gof$convergeCode != 0)) 
    warning("The optimization didn't fully converge for some samples:\n", 
            paste(samplesNames[fit.gof$convergeCode != 0], collapse = "; "), 
            "\n - check fit.gof for the convergeCode and convergeMessage")
  if (withOtherCells) 
    mRNAProportions <- cbind(mRNAProportions, otherCells = 1 - 
                               rowSums(mRNAProportions))
  tInds <- match(colnames(mRNAProportions), names(mRNA_cell))
  if (anyNA(tInds)) {
    defaultInd <- match("default", names(mRNA_cell))
    if (is.na(defaultInd)) {
      tStr <- paste(" and no default value is given for this mRNA per cell,", 
                    "so we cannot estimate the cellFractions, only", 
                    "the mRNA proportions")
    }
    else {
      tStr <- paste(" - using the default value of", mRNA_cell[defaultInd], 
                    "for these but this might bias the true cell proportions from", 
                    "all cell types.")
    }
    warning("mRNA_cell value unknown for some cell types: ", 
            paste(colnames(mRNAProportions)[is.na(tInds)], collapse = ", "), 
            tStr)
    tInds[is.na(tInds)] <- defaultInd
  }
  cellFractions <- t(t(mRNAProportions)/mRNA_cell[tInds])
  absCellFractions <- cellFractions
  cellFractions <- cellFractions/rowSums(cellFractions, na.rm = FALSE)
  return(list(mRNAProportions = mRNAProportions, cellFractions = cellFractions, 
              absCellFractions = absCellFractions,
              fit.gof = fit.gof))
}
