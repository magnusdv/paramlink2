
test_that("lod() catches bad `aff` input", {
  x = nuclearPed(1)
  x = setMarkers(x, marker(x))

  mod = diseaseModel("AD")

  expect_error(lod(x, aff = 1:4, mod = mod), "Invalid `aff`")
  expect_error(lod(x, aff = 0:1, mod = mod), "Unknown ID label: 0")
  expect_error(lod(x, aff = list(aff = 1), mod = mod), "Invalid name of `aff`")
  expect_error(lod(x, aff = list(affected = 1, unaffected = 1), mod = mod),
               "Individual with multiple affection statuses")
})


test_that("lod() with liability classes", {
  x = nuclearPed(3)
  x = setMarkers(x, marker(x, geno = c("1/2", "1/1", "1/2","1/2","1/2")))

  mod = diseaseModel("AD", pen = data.frame(f0 = 0, f1 = 1:0, f2 = 1:0))

  # Switch off last sib
  res = lod(x, aff = c(2,1,2,2,1), mod = mod, liab = c(1,1,1,1,2))

  expect_equal(round(res - log10(2), 3), 0)

  # Same: Set last sib to unknown
  expect_equal(res, lod(x, aff = c(2,1,2,2,0), mod = diseaseModel("AD")))
})
