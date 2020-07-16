stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

# Test that input is a single integer.
isCount = function(x, minimum = 1, maximum = NA) {
  isTRUE(length(x) == 1 &&
           (is.integer(x) || (is.numeric(x) && x == as.integer(x))) &&
           x >= minimum && (is.na(maximum) || x <= maximum))
}

# Test that input is a single number, with optional range constraints
isNumber = function(x, minimum = NA, maximum = NA) {
  isTRUE(length(x) == 1 &&
           is.numeric(x) &&
           (is.na(minimum) || x >= minimum) &&
           (is.na(maximum) || x <= maximum))
}

`%||%` = function(x, y) {
  if(is.null(x)) y else x
}

.mysetdiff = function(x, y) {
  unique.default(x[match(x, y, 0L) == 0L])
}

# Fast intersection. NB: assumes no duplicates!
.myintersect = function (x, y) {
  y[match(x, y, 0L)]
}

printDat = function(dat) {
  print(data.frame(
    M = paste(dat$pat1, dat$mat1, sep="/"),
    D = paste(c("N","D")[dat$pat2], c("N","D")[dat$mat2], sep="/"),
    prob = dat$prob))
}

# Convert various forms of `aff` input to affection status vector
fixAff = function(x, aff, verbose = FALSE) {
  origaff = aff

  n = pedsize(x)
 labs = labels(x)

  # If vector of ID labels, interpret as affecteds
  if(is.vector(aff) && !is.list(aff) && !anyDuplicated(aff) && length(aff) < n)
    aff = list(affected = aff, unaffected = setdiff(labs, aff))

  # If list of `affected` and `unaffected`, convert to status vector
  if(is.list(aff)) {
    validnames = c("affected", "unaffected", "unknown")
    if(!all(names(aff) %in% validnames))
      stop2("Invalid name of `aff`: ", setdiff(names(aff), validnames))

    a = aff$affected   %||% setdiff(labs, c(aff$unaffected, aff$unknown))
    u = aff$unaffected %||% setdiff(labs, c(a, aff$unknown))
    unk = aff$unknown  %||% setdiff(labs, c(a, u))

    if(dup <- anyDuplicated.default(c(a,u,unk)))
      stop2("Individual with multiple affection statuses: ", c(a,u,unk)[dup])

    affvec = rep(NA_integer_, n)
    affvec[internalID(x, a)]  = 2
    affvec[internalID(x, u)]  = 1
    affvec[internalID(x, unk)]= 0
    aff = affvec

    # Check for doubles (do after affvec, letting internalID catches unknown IDs)
    if(anyNA(aff))
      stop2("Individual with missing affection status: ", labs[is.na(aff)])
  }

  # Check
  ok = is.vector(aff) && length(aff) == n && is.numeric(aff) && all(aff %in% 0:2)
  if(!ok)
    stop2("Invalid `aff`: ", origaff)

  if(verbose) {
    cat("Affected:", toString(labs[aff == 2]), "\n")
    cat("Healthy :", toString(labs[aff == 1]), "\n")
    cat("Unknown :", toString(labs[aff == 0]), "\n")
  }

  aff
}
