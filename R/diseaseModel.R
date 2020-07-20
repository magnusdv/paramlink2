#' Disease models for linkage analysis
#'
#' Create a disease model in form of a `disModel` object, for use in e.g.
#' [lod()]
#'
#'
#' @param model An existing `disModel` object (to be modified by other
#'   arguments), or one of the integers 1-4:
#'
#'   * 1 = autosomal dominant; fully penetrant, dfreq = 1e-5
#'
#'   * 2 = autosomal recessive; fully penetrant, dfreq = 1e-5
#'
#'   * 3 = X-linked dominant; fully penetrant, dfreq = 1e-5
#'
#'   * 4 = X-linked recessive; fully penetrant, dfreq = 1e-5
#'
#' @param chrom Either "AUTOSOMAL" or "X". Lower case versions are allowed and
#'   will be converted automatically.
#' @param penetrances For autosomal models, a numeric of length 3 corresponding
#'   to `(f0, f1, f2)`, where `fi` is the probability of being affected given
#'   `i` disease alleles. It can also be a matrix with 3 columns, in which case
#'   each row represents a liability class.
#'
#'   For X-linked models, a list of two vectors named `male` and `female`, of
#'   lengths 2 and 3 respectively: `(f0, f1)` for males and `(f0, f1, f2)` for
#'   females. Alternatively, each list entry may be a matrix or data frame (with the same
#'   number of columns) where each row represents a liability class.
#'
#' @param dfreq A number in `[0,1]`: The population frequency of the disease
#'   allele.
#'
#' @return An object of class `disModel`, which is a list with entries `chrom`,
#'   `penetrances` and `dfreq`.
#'
#' @author Magnus Dehli Vigeland
#' @seealso [lod()]
#'
#' @examples
#'
#' # Fully penetrant AD model:
#' m1 = diseaseModel(model = "AD")
#'
#' # The above is equivalent to
#' m2 = diseaseModel(chrom = "Aut", penetrances = c(0,1,1), dfreq = 1e-5)
#' stopifnot(identical(m1, m2))
#'
#' # X-linked recessive model:
#' m3 = diseaseModel(model = "XR", dfreq = 0.01)
#'
#' # Long version of the above:
#' m4 = diseaseModel(chrom = "X", penetrances = list(male = c(0,1), female = c(0,0,1)),
#'               dfreq = 0.01)
#' stopifnot(identical(m3, m4))
#'
#' @export
diseaseModel = function(model = NULL, chrom = NULL, penetrances = NULL, dfreq = NULL) {

  if(is.character(model) && length(model) == 1) {
    model = switch(model,
                   AD = list(chrom = "AUTOSOMAL", penetrances = c(0, 1, 1), dfreq = 1e-05),
                   AR = list(chrom = "AUTOSOMAL", penetrances = c(0, 0, 1), dfreq = 1e-05),
                   XD = list(chrom = "X", penetrances = list(male = c(0, 1), female = c(0, 1, 1)), dfreq = 1e-05),
                   XR = list(chrom = "X", penetrances = list(male = c(0, 1), female = c(0, 0, 1)), dfreq = 1e-05),
                   stop2("Invalid model keyword: ", model, '\n(Valid keywords are AD, AR, XD, XR)'))
  }
  else if(!is.null(model) && !inherits(model, "disModel"))
    stop2("Argument `model` must be either a keyword or an existing `diseaseModel`: ", model)

  # chromosome
  if(!is.null(chrom)) {
    opts = c("AUTOSOMAL", "X")
    idx = pmatch(toupper(chrom), opts)
    if(length(idx) != 1 || is.na(idx))
      stop2("`chrom` must match either AUTOSOMAL or X: ", chrom)
    chrom = opts[idx]
  }
  else
    chrom = model$chrom %||% "AUTOSOMAL"

  # Disease freq
  dfreq = dfreq %||% model$dfreq %||% 1e-05
  if(!isNumber(dfreq, minimum = 0, maximum = 1))
    stop2("Parameter `dfreq` must be a single number between 0 and 1: ", dfreq)

  # Penetrance values
  pen = penetrances %||% model$penetrances
  if (is.null(pen))
    stop2("No penetrance values given")

  if(chrom == "AUTOSOMAL") {
    pen = fixPenetrances(pen)
  }
  else {
    if(!isTRUE(is.list(pen) && setequal(names(pen), c("male", "female"))))
       stop2("For X-linked models, `penetrances` must be a list with elements `male` and `female`")
    pen$male = fixPenetrances(pen$male, maleX = TRUE)
    pen$female = fixPenetrances(pen$female, maleX = FALSE)
  }

  # Return
  structure(list(chrom = chrom, penetrances = pen, dfreq = dfreq), class = "disModel")
}

fixPenetrances = function(p, maleX = FALSE) {
  orig = p
  if(is.vector(p) && !is.list(p)) {
    if(length(p) != 3 - maleX)
      stop2("`penetrances` vector must have length 3 (or 2 for males on X): ", orig)
    dim(p) = c(1, length(p))
  }
  else if(is.data.frame(p))
    p = as.matrix(p)
  else if(!is.matrix(p))
    stop2("`penetrances` must be either a vector or matrix-like")

  if(ncol(p) != 3 - maleX)
    stop2("Illegal number of columns in `penetrances`: ", ncol(p))

  mode(p) = "numeric"
  bad = is.na(p) | p < 0 | p > 1
  if(any(bad))
    stop2("Illegal penetrance value: ", orig[bad])

  # If colnames, sort. Otherwise set colnames
  nms = colnames(p)
  if(!is.null(nms) && all(nms %in% c("f0", "f1", "f2")))
    p = if(maleX) p[, c("f0", "f1")] else p[, c("f0", "f1", "f2")]
  else
    colnames(p) = if(maleX) c("f0", "f1") else c("f0", "f1", "f2")

  rownames(p) = 1:nrow(p)

  p
}


#' @export
print.disModel = function(x, ...) {
  cat(switch(x$chrom, AUTOSOMAL = "Autosomal", X = "X-linked"), "inheritance\n")

  pen = x$penetrances
  if(is.matrix(pen)) {  # autosomal
    if(nrow(pen) == 1)
      cat(sprintf("Penetrance: (f0, f1, f2) = (%s)\n", toString(pen)))
    else {
      cat("Penetrance liability classes:\n")
      print(pen)
    }
  }
  else {
    if(nrow(pen$male) == 1 && nrow(pen$female) == 1)
      cat("Penetrance:",
          sprintf("  Males  : (f0, f1)     = (%s)", toString(x$penetrances$male)),
          sprintf("  Females: (f0, f1, f2) = (%s)", toString(x$penetrances$female)),
          sep = "\n")
    else {
      cat("Penetrance liability classes:\n")
      print(pen)
    }
  }

  cat(sprintf("Disease allele frequency: %.3g\n", x$dfreq))
  invisible(x)
}
