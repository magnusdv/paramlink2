#' Multipoint LOD scores by MERLIN
#'
#' This function is a wrapper for the parametric linkage functionality of the
#' MERLIN software. For this to work, MERLIN must be installed and correctly
#' pointed to in the PATH environment variable.
#'
#' By default the following MERLIN command is run (via a call to [system()])
#' after creating appropriate files in the current working directory:
#'
#' \preformatted{% merlin -p _merlin.ped -d _merlin.dat -m _merlin.map -f
#' _merlin.freq --model _merlin.model --tabulate --markerNames --quiet}
#'
#' The resulting multipoint LOD scores are extracted from the output and
#' returned in R as a [linkres] object.
#'
#' @param x A `ped` object
#' @param aff A vector showing the affection status of each pedigree member (2 =
#'   affected, 1 = unaffected, 0 = unknown).
#' @param model A `disModel` object.
#' @param map A data frame with columns according to MERLIN format for map
#'   files.
#' @param markers A vector of marker names or indices. By default, all markers
#'   are used.
#' @param rho A numeric with values between 0 and 0.5: The recombination
#'   value(s) for which the LOD score is computed. The value of `rho` is
#'   converted to centiMorgans using Haldane's map function and included in the
#'   MERLIN command using the `--position` parameter.
#'
#'   Note: This parameter only works when `markers` has length 1.
#' @param maxOnly a logical indicating whether only the maximum LOD score should
#'   be returned. By default this is always done if the number of markers is 1.
#' @param options A character with additional options to the MERLIN command.
#' @param dir Path to a directory where the MERLIN input files should be
#'   created.
#' @param cleanup A logical indicating if MERLIN files should be removed after
#'   use. Default: TRUE.
#' @param verbose A logical.
#' @param \dots Further arguments passed on to [pedprobr::merlin()].
#'
#' @author Magnus Dehli Vigeland
#' @seealso [lod()]
#'
#' @references \url{http://csg.sph.umich.edu/abecasis/Merlin/}
#'
#' @return If the number of markers is 1, or if `maxOnly = TRUE`, a single
#'   number is returned.
#'
#'   Otherwise a `linkres` object similar to the output of [lod()], but with an
#'   addition column `CM`.
#'
#' @examples
#'
#' \donttest{
#' ### Requires MERLIN to be installed ###
#'
#' x = nuclearPed(2)
#' x = setMarkers(x, marker(x))
#' genotype(x, 1, c(1,2:4)) = 1:2
#' genotype(x, 1, 2) = 1
#'
#' aff = c(2,1,2,2)
#' model = diseaseModel(model = "AD")
#'
#' merlinLod(x, aff, model)
#' }
#'
#' @importFrom pedprobr merlin
#' @importFrom utils read.table
#' @export
merlinLod = function(x, aff, model, map = NULL, markers = NULL,
                     rho = 0, maxOnly = NA, options = "", dir  = ".",
                     cleanup = TRUE, verbose = FALSE, ...) {

  if (!is.null(markers))
    x = selectMarkers(x, markers)

  if(is.na(maxOnly))
    maxOnly = nMarkers(x) == 1

  singlepoint = nMarkers(x) == 1
  if(singlepoint && !isNumber(rho, minimum = 0, maximum = 0.5))
    stop2("`rho` must be a single number in the interval [0, 0.5]")
  if(!singlepoint && isNumber(rho) && rho > 0)
    stop2("Nonzero `rho` is possible only with a single marker")

  # Convert `aff` to status vector
  aff = fixAff(x, aff, verbose = verbose)

  prefix = file.path(dir, "_merlin")

  # Write files here (not by pedprobr::merlin) in order to include aff and model.
  merlinfiles = .writeMerlin(x, aff, model, map = map, prefix = prefix, verbose = verbose)

  options = sprintf("--model %s.model --tabulate --markerNames --quiet", prefix)

  # Singlepoint analysis
  if(singlepoint && rho > 0) {
    pos = as.numeric(map[1, 3]) - 50 * log(1 - 2 * rho)  #Haldane's map: Converting rec.fractions to cM positions.
    options = paste0(options, " --positions:", paste0(pos, collapse = ","))
  }

  # Run MERLIN
  merlin(x, options = options, verbose = verbose, dir = dir, generateFiles = FALSE, ...)

  ## Extract LOD scores (chrom, pos, marker names and LOD)
  res = read.table("merlin-parametric.tbl", sep = "\t", header = TRUE,
                   colClasses = c("numeric", "numeric", "character", "NULL", "numeric", "NULL", "NULL"))
  if (cleanup) {
    unlink("merlin-parametric.tbl")
    unlink(merlinfiles)
    if (verbose)
      cat("MERLIN files removed\n")
  }

  # Quick return if `maxOnly`
  if(maxOnly) {
    lods = res$LOD
    return(if(all(is.na(lods))) NA else max(lods, na.rm = TRUE))
  }

  # Build output data frame
  y = data.frame(CHROM = res$CHR,
                 MARKER = res$LABEL,
                 MB = posMb(x, 1:nMarkers(x)),
                 CM = res$POS * 100,
                 LOD = res$LOD,
                 stringsAsFactors = FALSE)

  res = structure(y,
            analysis = "Multipoint",
            rho = if(singlepoint) rho else NULL,
            class = c("linkres", "data.frame"))

  if (verbose)
    summary.linkres(res)

  res
}
