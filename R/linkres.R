#' S3 methods for class 'linkres'.
#'
#' Functions for printing, summarizing and plotting the results of a linkage
#' analysis.
#'
#' @param x,object A `linkres` object (normally produced by [lod()] or
#'   [merlinLod()]).
#' @param chrom (Optional) A numeric indicating which chromosomes to be included
#'   in the plot.
#' @param lwd,type,ylim Graphical parameters passed on to `plot()`.
#' @param xlab,ylab Axis labels.
#' @param \dots further arguments.
#'
#' @author Magnus Dehli Vigeland
#' @seealso [lod()], [merlinLod()], [lodPeaks()]
#'
#' @examples
#'
#' # Pedigree with 5 simulated SNP markers
#' x = nuclearPed(3)
#'
#' x = setMarkers(x, alleleMatrix = cbind(
#'   m1 = c("1/1", "1/2", "1/2", "1/2", "1/2"),
#'   m2 = c("1/2", "1/2", "1/2", "1/2", "1/2"),
#'   m3 = c("1/1", "1/2", "1/2", "1/2", "1/1")), sep="/")
#'
#' # Mother and all children affected
#' aff = c(1, 2, 2, 2, 2)
#'
#' # LOD scores under autosomal dominant model
#' lods = lod(x, aff, model = diseaseModel(model = "AD"))
#'
#' summary(lods)
#' as.data.frame(lods)
#' plot(lods)
#'
#' @name linkres
NULL


#' @rdname linkres
#' @export
print.linkres = function(x, ...) {
  x$LOD = round(x$LOD, 4)
  print.data.frame(x)
}

#' @rdname linkres
#' @export
summary.linkres = function(object, ...) {
  lods = object$LOD
  maxlod = max(lods, na.rm = TRUE)
  topMarkers = object$MARKER[maxlod - lods < 1e-04]
  cat("Max LOD score:", maxlod, "\n")
  cat("Achieved at marker(s):", topMarkers, "\n")
}


#' @rdname linkres
#' @importFrom graphics abline axis
#' @export
plot.linkres = function(x, chrom = NULL, type = "l", lwd = NA,
                        ylim = NULL, xlab = NULL, ylab = NULL, ...) {

  x = x[order(x$CHROM, x$MB), , drop = FALSE]

  if (!is.null(chrom))
    x = x[x$CHROM %in% chrom, , drop = FALSE]

  nM = nrow(x)
  if(nM == 0)
    stop2("No values to plot")

  # Default values for x axis ticks
  x_tick = NULL
  x_ticklab = TRUE

  ### X axis
  chr = x$CHROM
  pos = x$MB

  nomap = all(is.na(pos))
  if(nomap) {
    pos = seq_len(nM)
    if(nM < 5) x_tick = x_ticklab = seq_len(nM)
  }

  # Merge multiple chromosomes
  chrBreak = which(chr[-1] != chr[-nM])
  multichr = length(chrBreak) > 0

  if(multichr) {
    # Adjust x-axis positions
    for (b in chrBreak)
      pos[(b + 1):nM] = pos[(b + 1):nM] + x$MB[b]  # NB: by now, map$MB != pos

    x_tick = c(0, pos[chrBreak])
    x_ticklab = chr[c(chrBreak, nM)]
  }

  LOD = x$LOD

  # Fix default parameters
  xlab = xlab %||% if(multichr) "Chromosome" else if(nomap) "Marker" else paste("Chromosome", chr[1])
  ylab = ylab %||% "LOD"
  ylim = ylim %||% c(-1.2, max(c(3, LOD), na.rm = TRUE) + 0.3)

  if(is.na(lwd))
    lwd = switch(attr(x, "analysis"), Singlepoint = 1, Multipoint = 2)

  # Plot (without x axis)
  plot(pos, pmax(LOD, ylim[1]), xaxt = "n", type = type, lwd = lwd,
       ylim = ylim, xlab = xlab, ylab = ylab, ...)

  # Line at y = 0
  abline(h = 0, col = 1, lwd = 1)
  axis(1, at = x_tick, labels = x_ticklab, lwd.ticks = 1 + multichr)
}



#' @rdname linkres
#' @importFrom graphics points par
#' @export
points.linkres = function(x, chrom = NULL, type = "l", lwd = NA, ...) {

  x = x[order(x$CHROM, x$MB), , drop = FALSE]

  if (!is.null(chrom))
    x = x[x$CHROM %in% chrom, , drop = FALSE]

  nM = nrow(x)
  if(nM == 0)
    stop2("No values to plot")

  pos = x$MB
  if(all(is.na(pos)))
    pos = seq_len(nM)

  # Merge multiple chromosomes
  chr = x$CHROM
  chrBreak = which(chr[-1] != chr[-nM])
  multichr = length(chrBreak) > 0

  if(multichr) {
    # Adjust x-axis positions
    for (b in chrBreak)
      pos[(b + 1):nM] = pos[(b + 1):nM] + x$MB[b]  # NB: by now, map$MB != pos
  }

  LOD = x$LOD
  ymin = par("usr")[3] + diff(par("usr")[3:4]) * 0.04 # default expansion
  LODtrunc = pmax(LOD, ymin)

  if(is.na(lwd))
    lwd = switch(attr(x, "analysis"), Singlepoint = 1, Multipoint = 2)

  # Add points to existing plot
  points(pos, LODtrunc, type = type, lwd = lwd, ...)
}
