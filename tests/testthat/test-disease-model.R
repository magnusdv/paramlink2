context("disease models")

test_that("diseaseModel() catches errors", {

  expect_error(diseaseModel("A"), "Invalid model keyword: A")
  expect_error(diseaseModel(), "No penetrance values given")
  expect_error(diseaseModel(1:2), "Argument `model` must be either a keyword or an existing `diseaseModel`")
  expect_error(diseaseModel(T), "Argument `model` must be either a keyword or an existing `diseaseModel`")

  expect_error(diseaseModel(chrom = 1:10), "`chrom` must match either AUTOSOMAL or X")
  expect_error(diseaseModel(chrom = "chromX"), "`chrom` must match either AUTOSOMAL or X")

  expect_error(diseaseModel(dfreq = -1), "Parameter `dfreq` must be a single number between 0 and 1")
  expect_error(diseaseModel(dfreq = 1:2), "Parameter `dfreq` must be a single number between 0 and 1")

  expect_error(diseaseModel(pen = 1), "`penetrances` vector must have length 3 (or 2 for males on X)", fixed = T)
  expect_error(diseaseModel(pen = 1:4), "`penetrances` vector must have length 3 (or 2 for males on X)", fixed = T)
  expect_error(diseaseModel(pen = c(0,0,2)), "Illegal penetrance value")

  expect_error(diseaseModel("XR", pen = 1),
               "For X-linked models, `penetrances` must be a list with elements `male` and `female`")
  expect_error(diseaseModel("XR", pen = list(0:1, c(0,1,1))),
               "For X-linked models, `penetrances` must be a list with elements `male` and `female`")
  expect_error(diseaseModel("XR", pen = list(male=1,female=2)),
               "`penetrances` vector must have length 3 (or 2 for males on X)", fixed=T)
  expect_error(diseaseModel("XR", pen = list(male=0:1,female=c(1,1,NA))),
               "Illegal penetrance value: NA")
})
