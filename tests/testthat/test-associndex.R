# Test 1.- La función da error si no existe open
test_that("associndex falla si no existe Open", {
  int_no_open <- Amoladeras_int[Amoladeras_int$Canopy != "Open", ]
  expect_error(
    associndex(int_no_open, Amoladeras_cover),
    "does not contain a node named Open"
  )
})
# Test 2.- Devuelve null en la rama no allsp
test_that("associndex devuelve NULL en combinación unused", {
  expect_warning(
    res <- associndex(Amoladeras_int, Amoladeras_cover,
                      expand = "no",
                      rm_sp_no_cover = "allsp")
  )
  expect_null(res)
})

#Test 3.- rama yes sp devuelve dataframe
test_that("associndex modo rec devuelve dataframe", {

  res <- associndex(Amoladeras_int, Amoladeras_cover,
                    expand = "yes",
                    rm_sp_no_cover = "allsp")

  expect_s3_class(res, "data.frame")
  expect_gt(nrow(res), 0)

})

# Test 4.- Rama no onlycanopy devuelve dataframe
test_that("associndex modo fac devuelve dataframe", {

  res <- associndex(Amoladeras_int, Amoladeras_cover,
                    expand = "no",
                    rm_sp_no_cover = "onlycanopy")

  expect_s3_class(res, "data.frame")
  expect_gt(nrow(res), 0)

})

# Test 5.- Rama yes onlycanopy devuelve dataframe
test_that("associndex modo comp devuelve dataframe", {

  res <- associndex(Amoladeras_int, Amoladeras_cover,
                    expand = "yes",
                    rm_sp_no_cover = "onlycanopy")

  expect_s3_class(res, "data.frame")
  expect_gt(nrow(res), 0)

})

# Test 6.- la función calcula el thresold automáticamente
test_that("associndex calcula threshold automáticamente", {

  expect_message(
    associndex(Amoladeras_int, Amoladeras_cover,
               expand = "yes",
               rm_sp_no_cover = "allsp",
               threshold_density = "Weibull"),
    regexp = "threshold density has been set"
  )

})




