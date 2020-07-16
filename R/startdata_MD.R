
startdata_MD_AUT = function(x, m, affcode, model, eliminate = 1) {
  nInd = pedsize(x)

  glist = pedprobr:::.buildGenolist(x, marker = m, eliminate)
  if (attr(glist, "impossible")) {
    return(structure(list(), impossible = TRUE))
  }

  # Genotypes at dis locus (N=1; D=2). Same for all. TODO: Bedre med N=0, D=1?
  dlist0 = list(pat = c(1,1,2,2), mat = c(1,2,1,2))
  dlist1 = list(pat = c(1,1,2), mat = c(1,2,2))

  afreq = afreq(m)
  isFounder = logical(nInd)
  isFounder[founders(x, internal = TRUE)] = TRUE
  impossible = FALSE

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

    # Add probabilities
    prob = startprob_MD_AUT(g, aff = affcode[i], model = model,
                            afreq = afreq, founder = isFounder[i])
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


startprob_MD_AUT = function(g, aff, model, afreq, founder) {
  mpat = g$pat1
  mmat = g$mat1
  dpat = g$pat2
  dmat = g$mat2

  # Number of disease alleles: either 0, 1 or 2 for each genotype
  d.no = dpat + dmat - 2 # TODO: Burde d-locus ha N=0, D=1?

  # Assign prob at disease locus, using the penetrance values
  prob = switch(aff + 1,
                rep.int(1, length(d.no)),            # aff = 0: unknown
                (1 - model$penetrances)[d.no + 1],   # aff = 1: healthy
                model$penetrances[d.no + 1])         # aff = 2: affected

  # If founder, multiply with HW probs at both loci
  if (founder) {
    freqM = afreq[mpat] * afreq[mmat] * ((mpat != mmat) + 1)
    dfreq = model$dfreq
    freqD = dfreq^d.no * (1 - dfreq)^(2 - d.no)
    prob = prob * freqM * freqD
  }

  as.numeric(prob)
}
