
startdata_MD_AUT = function(x, m, affcode, model, liability, eliminate = 1) {
  nInd = pedsize(x)

  # Build genotype list in internal format
  glist = pedprobr:::.buildGenolist(x, marker = m, eliminate)

  # Return early if impossible
  impossible = attr(glist, "impossible")
  if (impossible)
    return(structure(list(), impossible = TRUE))

  # Genotypes at dis locus (N=1; D=2). Same for all.
  dlist0 = list(pat = c(1,1,2,2), mat = c(1,2,1,2))

  afreq = afreq(m)
  isFounder = logical(nInd)
  isFounder[founders(x, internal = TRUE)] = TRUE

  dfreq = model$dfreq

  ### Penetrances & liability
  penetMat = model$penetrances
  if(!is.matrix(penetMat) && ncol(penetMat) == 3)
    stop2("The `penetrances` model is not a matrix with 3 columns")

  if(!all(liability %in% 1:nrow(penetMat)))
    stop2("Illegal liability class: ", setdiff(liability, 1:nrow(penetMat)))

  dat = lapply(1:nInd, function(i) {

    # If impossible, finish loop quickly (cannot use break in `apply()`)
    if(impossible)
      return(NULL)

    gm = glist[[i]]
    gd = dlist0

    len1 = length(gm$mat)
    len2 = length(gd$mat)
    idx1 = rep(seq_len(len1), each = len2)
    idx2 = rep(seq_len(len2), times = len1)

    pat1 = gm$pat[idx1]
    mat1 = gm$mat[idx1]
    pat2 = gd$pat[idx2]
    mat2 = gd$mat[idx2]

    g = list(pat1 = pat1, mat1 = mat1, pat2 = pat2, mat2 = mat2)

    penet.i = penetMat[liability[i], ]

    # Add probabilities
    prob = startprob_MD_AUT(g, aff = affcode[i], penetrances = penet.i,
                            dfreq = dfreq,  afreq = afreq, founder = isFounder[i])
    g$prob = prob

    keep = prob > 0
    if (!any(keep)) {
      impossible = TRUE
      return(NULL)
    }

    if(!all(keep))
      g[] = lapply(g, function(vec) vec[keep])

    g
  })

  attr(dat, "impossible") = impossible
  dat
}


startprob_MD_AUT = function(g, aff, penetrances, dfreq, afreq, founder) {
  mpat = g$pat1
  mmat = g$mat1
  dpat = g$pat2
  dmat = g$mat2

  # Number of disease alleles: either 0, 1 or 2 for each genotype
  d.no = dpat + dmat - 2

  # Assign prob at disease locus, using the penetrance values
  prob = switch(aff + 1,
                rep.int(1, length(d.no)),      # aff = 0: unknown
                (1 - penetrances)[d.no + 1],   # aff = 1: healthy
                penetrances[d.no + 1])         # aff = 2: affected

  # If founder, multiply with HW probs at both loci
  if (founder) {
    freqM = afreq[mpat] * afreq[mmat] * ((mpat != mmat) + 1)
    freqD = dfreq^d.no * (1 - dfreq)^(2 - d.no)
    prob = prob * freqM * freqD
  }

  as.numeric(prob)
}
