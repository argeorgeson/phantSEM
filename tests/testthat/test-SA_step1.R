library(lavaan)

test_that("SA_step1 returns a list of length 8", {
  mod_obs <-
    " Y2 ~ X + M2
    M2 ~ X"

  mod_phant <-
    "Y2 ~ X + M1 + Y1 + M2
   M2 ~ X + M1 + Y1"

  cov <- matrix(c(1, .07, -.19, .07, 1, -.22, -.19, -.22, 1), nrow = 3, byrow = 3)
  rownames(cov) <- c("X", "M2", "Y2")
  colnames(cov) <- rownames(cov)

  lavoutput <- sem(model = mod_obs, sample.cov = cov, sample.nobs = 50)
  expect_equal(length(SA_step1(lavoutput, mod_obs, mod_phant)), 8)
})



test_that("SA_step1 has same results for sem() and lavaan(auto.var=TRUE)", {
  mod_obs <-
    " Y2 ~ X + M2
    M2 ~ X"

  mod_phant <-
    "Y2 ~ X + M1 + Y1 + M2
   M2 ~ X + M1 + Y1"

  cov <- matrix(c(1, .07, -.19, .07, 1, -.22, -.19, -.22, 1), nrow = 3, byrow = 3)
  rownames(cov) <- c("X", "M2", "Y2")
  colnames(cov) <- rownames(cov)

  lavoutputsem <- sem(model = mod_obs, sample.cov = cov, sample.nobs = 50)
  lavoutputlavaan <- lavaan(model = mod_obs, sample.cov = cov, sample.nobs = 50, auto.var = TRUE)
  expect_equal(SA_step1(lavoutputsem, mod_obs, mod_phant), SA_step1(lavoutputlavaan, mod_obs, mod_phant))
})


test_that("check that matrix template output is correct dimensions", {
  mod_obs <-
    " Y2 ~ X + M2
    M2 ~ X"

  mod_phant <-
    "Y2 ~ X + M1 + Y1 + M2
   M2 ~ X + M1 + Y1"

  cov <- matrix(c(1, .07, -.19, .07, 1, -.22, -.19, -.22, 1), nrow = 3, byrow = 3)
  rownames(cov) <- c("X", "M2", "Y2")
  colnames(cov) <- rownames(cov)

  lavoutput <- sem(model = mod_obs, sample.cov = cov, sample.nobs = 50)

  expect_equal(nrow(SA_step1(lavoutput, mod_obs, mod_phant)[[1]]), 5)
  expect_equal(nrow(SA_step1(lavoutput, mod_obs, mod_phant)[[3]]), 5)
  expect_equal(nrow(SA_step1(lavoutput, mod_obs, mod_phant)[[4]]), 5)
  expect_equal(nrow(SA_step1(lavoutput, mod_obs, mod_phant)[[7]]), 25)
})

test_that("check that symmetry is correct in output by comparing output to transpose", {
  mod_obs <-
    " Y2 ~ X + M2
    M2 ~ X"

  mod_phant <-
    "Y2 ~ X + M1 + Y1 + M2
   M2 ~ X + M1 + Y1"

  cov <- matrix(c(1, .07, -.19, .07, 1, -.22, -.19, -.22, 1), nrow = 3, byrow = 3)
  rownames(cov) <- c("X", "M2", "Y2")
  colnames(cov) <- rownames(cov)

  lavoutput <- sem(model = mod_obs, sample.cov = cov, sample.nobs = 50)
  mat <- SA_step1(lavoutput, mod_obs, mod_phant)[[4]]
  expect_equal(mat, t(mat))

  mat2 <- SA_step1(lavoutput, mod_obs, mod_phant)[[3]]
  expect_false(isTRUE(all.equal(mat2, t(mat2))))
})


test_that("check that non-nested models error out", {
  mod_obs <-
    " Y2 ~ X + M2
    M2 ~ X"

  mod_phant <-
    "Y2 ~  M1 + Y1 + M2
   M2 ~  M1 + Y1"

  cov <- matrix(c(1, .07, -.19, .07, 1, -.22, -.19, -.22, 1), nrow = 3, byrow = 3)
  rownames(cov) <- c("X", "M2", "Y2")
  colnames(cov) <- rownames(cov)

  lavoutput <- sem(model = mod_obs, sample.cov = cov, sample.nobs = 50)

  expect_error(SA_step1(lavoutput, mod_obs, mod_phant))
})
