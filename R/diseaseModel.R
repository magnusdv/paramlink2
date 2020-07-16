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
#' @param penetrances If `chrom = "AUTOSOMAL"`: a numeric of length 3 - `(f0,
#'   f1, f2)` - where `fi` is the probability of being affected given `i`
#'   disease alleles.
#'
#'   If `chrom = "X"`: a list of two vectors, of lengths 2 and 3, containing the
#'   penetrances for males and females, respectively: `penetrances = list(male =
#'   c(f0, f1), female = c(f0, f1, f2))`.
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

  if(is.character(model)) {
    model = switch(model,
                   AD = list(chrom = "AUTOSOMAL", penetrances = c(0, 1, 1), dfreq = 1e-05),
                   AR = list(chrom = "AUTOSOMAL", penetrances = c(0, 0, 1), dfreq = 1e-05),
                   XD = list(chrom = "X", penetrances = list(male = c(0, 1), female = c(0, 1, 1)), dfreq = 1e-05),
                   XR = list(chrom = "X", penetrances = list(male = c(0, 1), female = c(0, 0, 1)), dfreq = 1e-05),
                   stop2("Invalid model keyword: ", model, '\n(Valid keywords are AD, AR, XD, XR)'))
  }

  # chromosome
  if(!is.null(chrom))
    chrom = match.arg(toupper(chrom), c("AUTOSOMAL", "X"))
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
  if(is.character(pen))
    pen = switch(match.arg(tolower(pen), c("dominant", "recessive")),
      dominant = if(chrom == "AUTOSOMAL") c(0, 1, 1) else list(male = c(0, 1), female = c(0, 1, 1)),
      recessive = if(chrom == "AUTOSOMAL") c(0, 0, 1) else list(male = c(0, 1), female = c(0, 0, 1)),
      stop2('When `penetrances` is a character, it must be either "dominant" or "recessive": ', pen))

  if(chrom == "AUTOSOMAL") {
    if(!isTRUE(is.numeric(pen) && length(pen) == 3 && !anyNA(pen) && max(pen) <= 1 && min(pen) >= 0))
      stop2("Illegal penetrance vector for autosomal model: ", pen)
  }
  else {
    if(!isTRUE(is.list(pen) && length(pen) == 2))
       stop2("For X-linked models, `penetrances` must be a list of two vectors")
    penM = pen[[1]]
    penF = pen[[2]]
    if(!isTRUE(is.numeric(penM) && length(penM) == 2 && !anyNA(penM) && max(penM) <= 1 && min(penM) >= 0))
      stop2("Illegal male penetrance vector: ", penM)
    if(!isTRUE(is.numeric(penF) && length(penF) == 3 && !anyNA(penF) && max(penF) <= 1 && min(penF) >= 0))
      stop2("Illegal female penetrance vector: ", penF)

    names(pen) <- c("male", "female")
    names(pen$male) <- c("f0_m", "f1_m")
    names(pen$female) <- c("f0_f", "f1_f", "f2_f")
  }

  # Return
  structure(list(chrom = chrom, penetrances = pen, dfreq = dfreq), class = "disModel")
}

#' @export
print.disModel = function(x, ...) {
  txt = switch(x$chrom,
    AUTOSOMAL = c("Autosomal inheritance",
                  sprintf("Penetrance: (f0, f1, f2) = (%s)", toString(x$penetrances)),
                  sprintf("Disease allele frequency: %.3g", x$dfreq)),
    X = c("X-linked inheritance",
          "Penetrance:",
          sprintf("  Males  : (f0, f1) = (%s)", toString(x$penetrances$male)),
          sprintf("  Females: (f0, f1, f2) = (%s)", toString(x$penetrances$female)),
          sprintf("Disease allele frequency: %.3g", x$dfreq))
    )

  cat(paste0(txt, collapse = "\n"))
  invisible(x)
}
