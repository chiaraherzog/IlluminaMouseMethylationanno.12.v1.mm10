#' cpg.annotate
#' 
#' Helper function so DMRcate will be compatible with the mouse manifest. Corresponds to DMRcate::cpg.annotate.
#'
#' @param datatype Character string representing the type of data being analysed (array or sequencing).
#' @param object Either:
#' - A matrix of M-values, with unique Illumina probe IDs as rownames and unique sample IDs as column names or,
#' - A GenomicRatioSet, appropriately annotated.
#' @param what Does the data matrix contain Beta or M-values? Not needed if object is a GenomicRatioSet.
#' @param arraytype Is the data matrix sourced from mouse, EPIC or 450K data? 
#' @param analysis.type "differential" for dmrcate() to return DMRs; "variability" to return VMRs; "ANOVA" to return "whole experiment" DMRs, incorporating all possible contrasts from the design matrix using the moderated F-statistics; "diffVar" to return differentially variable methylated regions, using the missMethyl package to generate t-statistics.
#' @param design Study design matrix. Identical context to differential analysis pipeline in limma. Must have an intercept if contrasts=FALSE. Applies only when analysis.type %in% c("differential", "ANOVA", "diffVar").
#' @param contrasts Logical denoting whether a limma-style contrast matrix is specified. Only applicable when datatype="array" and analysis.type %in% c("differential", "diffVar").
#' @param cont.matrix Limma-style contrast matrix for explicit contrasting. For each call to cpg.annotate, only one contrast will be fit. Only applicable when datatype="array" and analysis.type %in% c("differential", "diffVar").
#' @param fdr FDR cutoff (Benjamini-Hochberg) for which CpG sites are individually called as significant. Used to index default thresholding in dmrcate(). Highly recommended as the primary thresholding parameter for calling DMRs. Not used when analysis.type == "variability".
#' @param coef The column index in design corresponding to the phenotype comparison. Corresponds to the comparison of interest in design when contrasts=FALSE, otherwise must be a column name in cont.matrix. Only applicable when analysis.type == "differential".
#' @param varFitcoef The columns of the design matrix containing the comparisons to test for differential variability. If left NULL, will test all columns. Identical context to missMethyl::varFit(). Only applicable when analysis.type %in% "diffVar".
#' @param topVarcoef Column number or column name specifying which coefficient of the linear model fit is of interest. It should be the same coefficient that the differential variability testing was performed on. Default is last column of fit object. Identical context to missMethyl::topVar(). Only applicable when analysis.type %in% "diffVar".
#' @param ... Extra arguments passed to the limma function lmFit() (analysis.type="differential").
#' @returns Annotated CpG object
#' @references
#' Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47.
#' 
#' Feng, H., Conneely, K. N., & Wu, H. (2014). A Bayesian hierarchical model to detect differentially methylated loci from single nucleotide resolution sequencing data. Nucleic Acids Research, 42(8), e69.
#' 
#' Phipson, B., & Oshlack, A. (2014). DiffVar: a new method for detecting differential variability with application to methylation in cancer and aging. Genome Biol, 15(9), 465.
#' 
#' Peters T.J., Buckley M.J., Statham, A., Pidsley R., Samaras K., Lord R.V., Clark S.J. and Molloy P.L. De novo identification of differentially methylated regions in the human genome. Epigenetics & Chromatin 2015, 8:6, doi:10.1186/1756-8935-8-6.
#' @author Tim J. Peters <t.peters@@garvan.org.au>, Chiara Herzog <chiara.herzog@@uibk.ac.at> (modifications) 
#' @export cpg.annotate

cpg.annotate <- function(datatype = c("array", "sequencing"), object, what = c("Beta", 
                                                                                "M"),
                          arraytype = c("EPIC", "450K", "mouse"),
                          analysis.type = c("differential", 
                                            "variability", "ANOVA", "diffVar"),
                          design, contrasts = FALSE, 
                          cont.matrix = NULL, fdr = 0.05, coef, varFitcoef = NULL, 
                          topVarcoef = NULL, ...){
  
  analysis.type <- match.arg(analysis.type)
  what <- match.arg(what)
  arraytype <- match.arg(arraytype)
  if (datatype == "array") {
    stopifnot(class(object)[1] %in% c("matrix", "GenomicRatioSet"))
    if (is(object, "matrix")) {
      if (arraytype == "450K") {
        grset <- minfi::makeGenomicRatioSetFromMatrix(mat = object, 
                                               array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19", 
                                               what = what)
      }
      if (arraytype == "EPIC") {
        grset <- minfi::makeGenomicRatioSetFromMatrix(mat = object, 
                                               array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19", 
                                               what = what)
      }
      
      if(arraytype == "mouse"){
        grset <- minfi::makeGenomicRatioSetFromMatrix(mat = object,
                                               array = "IlluminaMouseMethylation", annotation = "12.v1.mm10",
                                               what = what)
      }
    } else {
      grset <- object
    }
    object <- minfi::getM(grset)
    switch(analysis.type, differential = {
      stopifnot(is.matrix(design))
      if (!contrasts) {
        stopifnot(colnames(design)[1] == "(Intercept)")
      } else {
        stopifnot(!is.null(cont.matrix))
      }
      fit <- limma::lmFit(object, design, ...)
      if (contrasts) {
        stopifnot(coef %in% colnames(cont.matrix))
        fit <- limma::contrasts.fit(fit, cont.matrix)
      }
      fit <- limma::eBayes(fit)
      tt <- limma::topTable(fit, coef = coef, number = nrow(object))
      nsig <- sum(tt$adj.P.Val < fdr)
      if (nsig == 0) {
        message("Your contrast returned no individually significant probes. Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DMRs, but be warned there is an increased risk of Type I errors.")
      }
      if (nsig > 0 & nsig <= 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes; a small but real effect. Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases the risk of Type I errors."))
      }
      if (nsig > 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes. We recommend the default setting of pcutoff in dmrcate()."))
      }
      betafit <- limma::lmFit(minfi::ilogit2(object), design, ...)
      if (contrasts) {
        betafit <- limma::contrasts.fit(betafit, cont.matrix)
      }
      betafit <- limma::eBayes(betafit)
      betatt <- limma::topTable(betafit, coef = coef, number = nrow(object))
      m <- match(rownames(tt), rownames(betatt))
      tt$diff <- betatt$logFC[m]
      m <- match(rownames(object), rownames(tt))
      tt <- tt[m, ]
      anno <- minfi::getAnnotation(grset)
      stat <- tt$t
      annotated <- GenomicRanges::GRanges(as.character(anno$chr), IRanges::IRanges(anno$pos, 
                                                           anno$pos), stat = stat, diff = tt$diff, ind.fdr = tt$adj.P.Val, 
                           is.sig = tt$adj.P.Val < fdr)
      names(annotated) <- rownames(tt)
    }, variability = {
      RSanno <- minfi::getAnnotation(grset)
      wholevar <- var(object)
      weights <- apply(object, 1, var)
      weights <- weights/mean(weights)
      annotated <- GenomicRanges::GRanges(as.character(RSanno$chr), IRanges::IRanges(RSanno$pos, 
                                                             RSanno$pos), stat = weights, diff = rep(0, nrow(object)), 
                           ind.fdr = rep(0, nrow(object)), is.sig = weights > 
                             quantile(weights, 0.95))
      names(annotated) <- rownames(object)
    }, ANOVA = {
      message("You are annotating in ANOVA mode: consider making the value of fdr quite small, e.g. 0.001")
      stopifnot(is.matrix(design))
      fit <- limma::lmFit(object, design, ...)
      fit <- limma::eBayes(fit)
      sqrtFs <- sqrt(fit$F)
      sqrtfdrs <- p.adjust(fit$F.p.value, method = "BH")
      nsig <- sum(sqrtfdrs < fdr)
      if (nsig == 0) {
        message("Your design returned no individually significant probes for ANOVA. Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DMRs, but be warned there is an increased risk of Type I errors.")
      }
      if (nsig > 0 & nsig <= 100) {
        message(paste("Your design returned", nsig, 
                      "individually significant probes for ANOVA; a small but real effect. Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases the risk of Type I errors."))
      }
      if (nsig > 100) {
        message(paste("Your design returned", nsig, 
                      "individually significant probes for ANOVA. We recommend the default setting of pcutoff in dmrcate(). Large numbers (e.g. > 100000) may warrant a smaller value of the argument passed to fdr"))
      }
      anno <- minfi::getAnnotation(grset)
      stat <- sqrtFs
      annotated <- GenomicRanges::GRanges(as.character(anno$chr), IRanges::IRanges(anno$pos, 
                                                           anno$pos), stat = stat, diff = 0, ind.fdr = sqrtfdrs, 
                           is.sig = sqrtfdrs < fdr)
      names(annotated) <- rownames(object)
    }, diffVar = {
      stopifnot(is.matrix(design))
      if (!contrasts) {
        stopifnot(colnames(design)[1] == "(Intercept)")
      } else {
        stopifnot(!is.null(cont.matrix))
      }
      fitvar <- missMethyl::varFit(object, design = design, coef = varFitcoef)
      if (contrasts) {
        fitvar <- missMethyl::contrasts.varFit(fitvar, cont.matrix)
      }
      tt <- missMethyl::topVar(fitvar, coef = topVarcoef, number = nrow(object))
      nsig <- sum(tt$Adj.P.Value < fdr)
      if (nsig == 0) {
        message("Your contrast returned no individually significant probes. Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DVMRs, but be warned there is an increased risk of Type I errors.")
      }
      if (nsig > 0 & nsig <= 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes; a small but real effect. Consider manually setting the value of pcutoff to return more DVMRs, but be warned that doing this increases the risk of Type I errors."))
      }
      if (nsig > 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes. We recommend the default setting of pcutoff in dmrcate()."))
      }
      m <- match(rownames(object), rownames(tt))
      tt <- tt[m, ]
      anno <- minfi::getAnnotation(grset)
      stat <- tt$t
      annotated <- GenomicRanges::GRanges(as.character(anno$chr), IRanges::IRanges(anno$pos, 
                                                           anno$pos), stat = stat, diff = 0, ind.fdr = tt$Adj.P.Value, 
                           is.sig = tt$Adj.P.Value < fdr)
      names(annotated) <- rownames(tt)
    })
    annotated <- sort(annotated)
    return(new("CpGannotated", ranges = annotated))
  }
  if (datatype == "sequencing") {
    stop("Sequencing mode is deprecated for cpg.annotate(). Please use sequencing.annotate().")
  }
  else {
    message("Error: datatype must be one of 'array' or 'sequencing'")
  }
}
