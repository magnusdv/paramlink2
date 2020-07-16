context("lod function")

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
