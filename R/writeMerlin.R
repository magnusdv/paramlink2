########################################################################
# This function is mostly copied from pedprobr:::writePed_merlin(),
# but is modified to accommodate `model`, `aff` and `liability`.
########################################################################

#' @importFrom utils write.table
.writeMerlin = function(x, aff, model, liability, map = NULL, prefix, verbose = TRUE) {

  ### File names with extensions
  ext = c("ped", "map", "dat", "freq", "model")
  fnames = paste(prefix, ext, sep = ".")
  names(fnames) = ext

  ### model file
  dfreq = format(model$dfreq, scientific = F, decimal.mark = ".")

  p = model$penetrances
  if(model$chrom == "AUTOSOMAL")
    pstr = apply(p, 1, function(v)
      paste(format(v, scientific = F, decimal.mark = "."), collapse = ","))
  else
    stop2("X-linked models not implemented yet")

  nliab = length(pstr)
  mod = cbind(DISEASE = c("disease", paste("LIABILITY =", 1:nliab)[-nliab], "OTHERWISE"),
              ALLELE_FREQ = c(dfreq, rep("", nliab)),
              PENETRANCE = c("*", pstr),
              LABEL = c("my_model", rep("", nliab)))

  write.table(mod, file = fnames[["model"]], sep = "\t", quote = FALSE, row.names = FALSE, col.names = F)

  ### ped file
  if(is.pedList(x)) {
    pedmatr = do.call(rbind, lapply(x, as.matrix, include.attrs = FALSE))
    x = x[[1]] # for later stuff
  } else {
    pedmatr = as.matrix(x, include.attrs = FALSE)
  }

  ### add famid, affection status and liability class
  if(length(aff) != nrow(pedmatr))
    stop2("`aff` is incompatible with the pedigree matrix")
  pedmatr = cbind(1, pedmatr[, 1:4], aff, liability, pedmatr[, -(1:4)])

  ### write
  write(t.default(pedmatr), file = fnames[["ped"]], ncolumns = ncol(pedmatr))
  if(verbose) message("File written: ", fnames[["ped"]])

  ### map file
  if(is.null(map))
    map = getMap(x, na.action = 1, verbose = verbose)
  else
    names(map) = toupper(names(map))

  write.table(map, file = fnames[["map"]], col.names = FALSE, row.names = FALSE, quote = FALSE)
  if(verbose) message("File written: ", fnames[["map"]])

  ### dat file
  datmatr = cbind(c("A", "C", rep("M", nrow(map))), c("disease",  "LIABILITY", map$MARKER))
  write.table(datmatr, file = fnames[["dat"]], col.names = FALSE, row.names = FALSE, quote = FALSE)
  if(verbose) message("File written: ", fnames[["dat"]])

  ### freq file
  nalls = nAlleles(x)
  L = sum(nalls) + length(nalls)
  cum = cumsum(c(1, nalls + 1))
  length(cum) = length(nalls)  #remove last

  col1 = rep("A", L)
  col1[cum] = "M"

  col2 = character(L)
  col2[cum] = map$MARKER

  allalleles = unlist(lapply(nalls, seq_len)) # numerical allele names for merlin!
  col2[-cum] = allalleles

  col3 = character(L)
  allfreqs = unlist(lapply(x$MARKERS, afreq))
  col3[-cum] = format(allfreqs, scientifit = FALSE, digits = 6)

  freqmatr = cbind(col1, col2, col3)
  write.table(freqmatr, file = fnames[["freq"]], col.names = FALSE, row.names = FALSE, quote = FALSE)
  if(verbose) message("File written: ", fnames[["freq"]])

  # return file names
  invisible(unname(fnames))
}

