library(lavaan)

test_that("Check number of output objects", {
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
  step1 <- SA_step1(lavoutput, mod_obs, mod_phant)

  phantom_assignment <- list(
    "CovM1M2" = seq(0, .4, .1),
    "CovM1X" = 0,
    "CovM1Y2" = seq(-.3, .3, .1),
    "CovY1M1" = "CovY2M2",
    "CovY1M2" = "CovM1Y2",
    "CovY1X" = 0,
    "CovY1Y2" = "CovM1M2",
    "VarM1" = 1,
    "VarY1" = 1
  )

  expect_equal(length(SA_step2(phantom_assignment, step1)), 5)
})


test_that("Check if function works when given one set of test values to be used for two phantom params", {
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
  step1 <- SA_step1(lavoutput, mod_obs, mod_phant)

  phantom_assignment <- list(
    "CovM1M2" = seq(0, .4, .1),
    "CovM1X" = 0,
    "CovM1Y2" = "CovM2Y2",
    "CovY1M1" = "CovY2M2",
    "CovY1M2" = "CovM2Y2",
    "CovY1X" = 0,
    "CovY1Y2" = "CovM1M2",
    "VarM1" = 1,
    "VarY1" = 1
  )

  expect_equal(length(SA_step2(phantom_assignment, step1)), 5)
})


test_that("Check that phantom variables equal to other phantom variables will error out if the phantom variable name is not written correctly", {
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
  step1 <- SA_step1(lavoutput, mod_obs, mod_phant)



  phantom_assignment2 <- list(
    "CovM1M2" = seq(0, .4, .1),
    "CovM1X" = 0,
    "CovM1Y2" = seq(0, .3, .1),
    "CovY1M1" = "CovM2Y2",
    "CovY1M2" = "CovY2M1",
    "CovY1X" = 0,
    "CovY1Y2" = "CovM2M1",
    "VarM1" = 1,
    "VarY1" = 1
  )


  expect_error(SA_step2(phantom_assignment2, step1))
})


test_that("Check if covariance names can be written with the variables in either order (e.g., CovY1Y2 and CovY2Y1)", {
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
  step1 <- SA_step1(lavoutput, mod_obs, mod_phant)

  phantom_assignment1 <- list(
    "CovM1M2" = seq(0, .4, .1),
    "CovM1X" = "CovXM2",
    "CovM1Y2" = "CovM2Y2",
    "CovY1M1" = "CovY2M2",
    "CovY1M2" = "CovM2Y2",
    "CovY1X" = "CovY2X",
    "CovY1Y2" = "CovM1M2",
    "VarM1" = 1,
    "VarY1" = 1
  )

  phantom_assignment2 <- list(
    "CovM1M2" = seq(0, .4, .1),
    "CovM1X" = "CovM2X",
    "CovM1Y2" = "CovY2M2",
    "CovY1M1" = "CovM2Y2",
    "CovY1M2" = "CovY2M2",
    "CovY1X" = "CovXY2",
    "CovY1Y2" = "CovM1M2",
    "VarM1" = 1,
    "VarY1" = 1
  )


  expect_equal(SA_step2(phantom_assignment1, step1)[[4]], SA_step2(phantom_assignment2, step1)[[4]])
  expect_equal(SA_step2(phantom_assignment1, step1), SA_step2(phantom_assignment2, step1))
})


test_that("Check that function works if numeric values are given for all phantom parameters", {
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
  step1 <- SA_step1(lavoutput, mod_obs, mod_phant)

  phantom_assignment <- list(
    "CovM1M2" = .4,
    "CovM1X" = 0,
    "CovM1Y2" = .5,
    "CovY1M1" = .3,
    "CovY1M2" = .1,
    "CovY1X" = 0,
    "CovY1Y2" = .4,
    "VarM1" = 1,
    "VarY1" = 1
  )

  expect_equal(length(SA_step2(phantom_assignment, step1)), 5)
})


test_that("Check that when phantom parameter is not included that the default test values are used", {
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
  step1 <- SA_step1(lavoutput, mod_obs, mod_phant)

  phantom_assignment <- list(
    "CovM1M2" = .4,
    "CovM1X" = 0,
    "CovM1Y2" = .5,
    "CovY1M1" = .3,
    "CovY1M2" = .1,
    "CovY1X" = 0,
    # "CovY1Y2"= .4,
    "VarM1" = 1,
    "VarY1" = 1
  )

  expect_equal(SA_step2(phantom_assignment, step1)[[3]][1], data.frame(CovY1Y2 = c(seq(-.3, .3, .1))))
})


# phantom_assignment <-list("CovM1M2"=seq(-.1,.1,.05),
#                          "CovM1X"= 0,
#                          "CovM1Y2"= .5,
#                          "CovY1M1"= .3,
#                          "CovY1M2"= .1,
#                          "CovY1X"= 0,
#                          #"CovY1Y2"= .4,
#                          "VarM1"= 1,
#                          "VarY1" = 1)




test_that("Check that function works when provided with observed variable names for all phantom variables", {
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
  step1 <- SA_step1(lavoutput, mod_obs, mod_phant)

  phantom_assignment <- list(
    "CovM1M2" = "CovM2Y2",
    "CovM1X" = "CovM2X",
    "CovM1Y2" = "CovM2Y2",
    "CovY1M1" = "CovM2Y2",
    "CovY1M2" = "CovM2Y2",
    "CovY1X" = "CovY2X",
    "CovY1Y2" = "CovM2Y2",
    "VarM1" = "VarM2",
    "VarY1" = "VarY2"
  )

  expect_equal(length(SA_step2(phantom_assignment, step1)), 5)
})
