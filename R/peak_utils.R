#' LOD score peaks
#'
#' Identify and summarise LOD score peaks
#'
#' A peak is defined as a run of at least `width` consecutive markers with LOD
#' score above or equal to `threshold`. If possible, one flanking marker is
#' included on each side of the peak.
#'
#' @param x A [linkres] object, or data frame with columns `CHROM`, `MB`, `LOD`.
#' @param threshold A single number
#' @param width A positive integer
#' @param physmap A matrix or data frame with three columns: Marker
#'   name, chromosome and physical position.
#'
#' @return A list of data frames.
#' @author Magnus Dehli Vigeland
#'
#' @seealso [linkres], [lod()], [merlin()],
#'
#' @examples
#'
#' ## Use built-in dataset `dominantFam`
#' x = dominantFam
#' aff = list(affected = c(2,3,9,10,19,6,16,17))
#' res = lod(x, aff, model = diseaseModel("AD"))
#'
#' # All peaks above LOD = 1.5
#' lodPeaks(res, threshold = 1.5)
#' peakSummary(res, threshold = 1.5)
#'
#' @export
lodPeaks = function(x, threshold, width = 1) {

  peak1chr = function(xchr, threshold, width) { # xchr must be sorted!
    rl = rle(xchr$LOD >= threshold)
    while (1) {
      short = rl$values & (rl$lengths < width)
      if (any(short)) {
        rl$values[short] <- FALSE
        rl = rle(inverse.rle(rl))
      }
      else break
    }
    xchr_nrow = nrow(xchr)
    if (!any(rl$values))
      return(NULL)

    start_ind = c(0, cumsum(rl$lengths))[which(rl$values)]
    stop_ind = start_ind + rl$lengths[rl$values] + 1  # plus 1 to compensate for endpoint[1]

    lapply(1:length(start_ind), function(i) {
      strt = start_ind[i]
      stp = stop_ind[i]
      telomeric = c(if (strt == 0) "start", if (stp > xchr_nrow) "end")
      if (length(telomeric) == 0)
        telomeric = "no"
      strt = max(strt, 1)
      stp = min(stp, xchr_nrow)
      structure(xchr[strt:stp, , drop = FALSE], rownames = NULL, telomeric = telomeric)
    })
  }

  df = x[!is.na(x$LOD), , drop = FALSE]
  chrs = unique.default(df$CHROM)
  res = list()
  for (chr in chrs) {
    dfchr = df[df$CHROM == chr, , drop = FALSE]
    res = c(res, peak1chr(dfchr, threshold = threshold, width = width))
  }
  res

}

#' @rdname lodPeaks
#' @importFrom utils read.table
#' @export
peakSummary = function(x, threshold, width = 1, physmap = NULL) {
  if (inherits(x, "linkres")) {
    if (is.null(threshold))
      stop2("argument `threshold` is missing, with no default")
    x = lodPeaks(x, threshold, width)
  }
  dat = lapply(x, function(df) {
    n = nrow(df)
    from_mb = df$MB[1]
    to_mb = df$MB[n]

    data.frame(
      CHROM = df$CHROM[1],
      FROM_MARKER = df$MARKER[1],
      TO_MARKER = df$MARKER[n],
      FROM_MB = from_mb,
      TO_MB = to_mb,
      LEN = to_mb - from_mb,
      N = n,
      TELOMERIC = attr(df, "telomeric"),
      FROM_LOD = round(df$LOD[1], 3),
      MAX_LOD = round(max(df$LOD),3),
      TO_LOD = round(df$LOD[n],3),
      stringsAsFactors = FALSE
    )
  })

  a = do.call(rbind, dat)
  if (!is.null(physmap)) {
    if (is.matrix(physmap))
      physmap = as.data.frame(physmap)
    if (is.character(physmap))
      physmap = read.table(physmap, header = T, as.is = TRUE)
    if (is.data.frame(physmap)) {
      from_bp = physmap[match(a$FROM_MARKER, physmap[, 1]), 3]
      to_bp = physmap[match(a$TO_MARKER, physmap[, 1]), 3]
      a = cbind(a, FROM_BP = from_bp, TO_BP = to_bp)
      a = a[, c(1:5, 12:13, 6:11)]
    }
  }
  a
}
