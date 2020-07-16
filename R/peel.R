.peel_MM_AUT = function(dat, sub, rho) {
  fa = sub$father
  mo = sub$mother
  link = sub$link

  faDat = dat[[fa]]
  moDat = dat[[mo]]
  likel = faDat$prob %*% t.default(moDat$prob)
  faLen = nrow(likel)
  moLen = ncol(likel)

  # Loop over the children, except the link if this is a child.
  for (ch in .mysetdiff(sub$children, link)) {
    chDat = dat[[ch]]
    chprob = chDat$prob
    chLen = length(chprob)

    transPat = .transProbMM(faDat, chDat[c('pat1', 'pat2')], rho = rho)
    transMat = .transProbMM(moDat, chDat[c('mat1', 'mat2')], rho = rho)
    dim(transMat) = c(chLen, moLen)

    transMat_rep = transMat[rep(seq_len(chLen), faLen), ] # as.numeric(do.call(rbind, rep(list(transMat), faLen)))
    mm = .colSums((transPat * chprob) * transMat_rep, chLen, faLen * moLen)
    likel = likel * mm
  }

  if (link == 0)
    return(sum(likel))

  # Goal is to update the probabilities of the link individual ("pivot"):
  pivDat = dat[[link]]

  if (link == fa)
    res = .rowSums(likel, faLen, moLen)
  else if (link == mo)
    res = .colSums(likel, faLen, moLen)
  else { # link is a child
    pivp = pivDat$prob
    pivLen = length(pivp)

    TRarray = array(0, dim = c(faLen, moLen, pivLen))
    transPat = .transProbMM(faDat, pivDat[c('pat1', 'pat2')], rho = rho)
    transMat = .transProbMM(moDat, pivDat[c('mat1', 'mat2')], rho = rho)
    dim(transPat) = c(pivLen, faLen)
    dim(transMat) = c(pivLen, moLen)

    for (i in seq_len(faLen)) {
      transpat = transPat[, i]
      for (j in seq_len(moLen))
        TRarray[i, j, ] = transpat * transMat[, j]
    }

    arr = as.vector(TRarray) * as.vector(likel)
    dim(arr) = dim(TRarray)
    res = .colSums(arr, faLen * moLen, pivLen)  #sum for each entry of haps[[link]]
    res = res * pivp
  }

  # Update the probabilities
  pivDat$prob = res

  # Reduce if possible
  if(any(res == 0)) {
    pivDat$pat1 = pivDat$pat1[res > 0]
    pivDat$mat1 = pivDat$mat1[res > 0]
    pivDat$pat2 = pivDat$pat2[res > 0]
    pivDat$mat2 = pivDat$mat2[res > 0]
    pivDat$prob = pivDat$prob[res > 0]
  }

  dat[[link]] = pivDat

  if (sum(res) == 0)
    attr(dat, "impossible") = TRUE

  dat
}


.transProbMM = function(par, gamete, rho) {
  # parent = list(pat1, mat1, pat2, mat2); vectors of same length
  # gamete = list(pat1, pat2) (or mat); vecs of same length
  lenPar = length(par[[1]])
  lenGam = length(gamete[[1]])

  p1 = rep(par$pat1, each = lenGam)
  m1 = rep(par$mat1, each = lenGam)
  p2 = rep(par$pat2, each = lenGam)
  m2 = rep(par$mat2, each = lenGam)

  gam1 = rep(gamete[[1]], lenPar)
  gam2 = rep(gamete[[2]], lenPar)

  res =
    (p1 == gam1 & p2 == gam2) * (1 - rho) +
    (m1 == gam1 & m2 == gam2) * (1 - rho) +
    (p1 == gam1 & m2 == gam2) * rho +
    (m1 == gam1 & p2 == gam2) * rho

  res/2
}
