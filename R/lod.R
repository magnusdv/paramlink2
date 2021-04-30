#' Single-point LOD score
#'
#' Calculates the single-point LOD scores of a pedigree for the specified
#' markers. The recombination ratio between the disease and marker loci can be
#' either fixed or optimized.
#'
#' The LOD score of a marker is defined as \deqn{LOD(\rho) = \log[10]
#' \frac{L(\rho)}{L(0.5)}}{% LOD(\rho) = log10 L(\rho)/L(0.5)} where
#' \eqn{L(\rho)} denotes the likelihood of the observed marker genotypes given a
#' recombination ratio \eqn{\rho} between the marker and the disease locus.
#'
#' @param x A `ped` object.
#' @param aff A vector of length `pedsize(x)`, containing the affection status
#'   of each pedigree member (2 = affected; 1 = unaffected; 0 = unknown).
#'   Alternatively, a list of vectors named `affected` and (optionally)
#'   `unaffected`.
#' @param model A `disModel` object, typically created with [diseaseModel()].
#' @param liability NULL (default) or a vector of length `pedsize(x)`
#'   indicating the liability class (a row number of `model$penetrances`) of
#'   each individual.
#' @param markers A vector of marker names or indices referring to markers
#'   attached to `x`. By default all markers are included.
#' @param rho A number between 0 and 0.5 (inclusive); the hypothesised
#'   recombination ratio between the marker and the disease locus.
#' @param maxOnly a logical indicating whether only the maximum LOD score should
#'   be returned. By default this is always done if the number of markers is 1.
#' @param loopBreakers A vector of ID labels indicating loop breakers. (Only
#'   relevant for inbred pedigrees.)
#' @param peelOrder For internal use.
#' @param verbose a logical: verbose output or not.
#'
#' @return If the number of markers is 1, or if `maxOnly = TRUE`, a single
#'   number is returned.
#'
#'   Otherwise a `linkres` object, which is basically a data frame with columns
#'   `CHROM`, `MARKER`, `MB` and `LOD`.
#'
#' @author Magnus Dehli Vigeland
#' @seealso [linkres], [merlinLod()], [diseaseModel()], [lodPeaks()]
#'
#' @examples
#'
#' x = nuclearPed(2)
#' x = setMarkers(x, marker(x))
#' genotype(x, 1, c(1,2:4)) = 1:2
#' genotype(x, 1, 2) = 1
#'
#' aff = c(2,1,2,2)
#' model = diseaseModel(model = "AD")
#'
#' lod(x, aff, model)
#'
#' @importFrom pedprobr likelihood
#' @export
lod = function(x, aff, model, rho = 0, liability = NULL, markers = NULL, maxOnly = NA,
               loopBreakers = NULL, peelOrder = NULL, verbose = FALSE) {

  if(model$chrom == "X")
    stop2("X-linked disease models are not implemented at the moment.")

  if (is.singleton(x))
    stop2("This function is not applicable to singleton objects.")

  if(hasInbredFounders(x))
    stop2("Likelihood of linked markers is not implemented in pedigrees with founder inbreeding.\n",
          "(Note that this is usually not well-defined)")

  if(hasSelfing(x))
    stop2("Likelihood of pedigrees with selfing is not implemented.\n",
          "Contact the maintainer if this is important to you.")

  if(!hasMarkers(x))
    stop2("The pedigree does not have any attached markers")

  if(!isNumber(rho, minimum = 0, maximum = 0.5))
    stop2("`rho` must be a single number in the interval [0, 0.5]")

  if(!is.null(markers))
    x = selectMarkers(x, markers)
  else
    markers = seq_len(nMarkers(x))

  if(is.na(maxOnly))
    maxOnly = nMarkers(x) == 1

  if(hasUnbrokenLoops(x)) {
    x = breakLoops(x, loopBreakers = loopBreakers, verbose = verbose)
    markers = x$MARKERS
  }

  # Convert aff into vector of 0-1-2 status
  aff = fixAff(x, aff, verbose = verbose)

  # Liability classes
  liability = liability %||% rep(1, pedsize(x))

  peel_MD = function(dat, nuc) .peel_MM_AUT(dat, nuc, rho = rho)
  peel0 = function(dat, nuc) .peel_MM_AUT(dat, nuc, rho = 0.5)
  peelOrder = peelingOrder(x)

  peelProcess = pedprobr:::peelingProcess

  lods = vapply(x$MARKERS, function(m) {
      startDat = startdata_MD_AUT(x, m, aff, model = model, liability = liability)

      # Numerator
      numer = peelProcess(x, m, startdata = startDat, peeler = peel_MD, peelOrder = peelOrder)

      # Denominator: Unlinked
      denom = peelProcess(x, m, startdata = startDat, peeler = peel0, peelOrder = peelOrder)

      log10(numer) - log10(denom)
    },
    FUN.VALUE = 0)


  # Quick return if `maxOnly`
  if(maxOnly)
    return(if(all(is.na(lods))) NA else max(lods, na.rm = TRUE))

  # Otherwise, bind with map
  map = getMap(x, na.action = 0, verbose = verbose)
  if(anyNA(map$MARKER))
    map$MARKER = seq_len(nrow(map))

  res = structure(cbind(map, LOD = lods),
                  analysis = "Singlepoint",
                  rho = rho,
                  class = c("linkres", "data.frame"))

  if (verbose)
    summary.linkres(res)

  res
}

