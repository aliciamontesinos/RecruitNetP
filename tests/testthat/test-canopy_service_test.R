
res <- canopy_service_test(Amoladeras_int, Amoladeras_cover)

test_that("canopy_service_test  data.frame with correct columns", {

  expect_s3_class(res, "data.frame")

  expect_equal(
    colnames(res),
    c("Canopy", "Fc", "Ac", "Fro", "Ao",
      "testability", "Significance", "Test_type", "Canopy_effect")
  )

  expect_equal(nrow(res), 23)
})


#------------------------------------

test_that("Aggregated values correct for Artemisia_barrelieri", {

  fila <- res[res$Canopy == "Artemisia_barrelieri", ]

  expect_equal(fila$Fc, 18)
  expect_equal(fila$Fc, sum(Amoladeras_int[Amoladeras_int$Canopy=="Artemisia_barrelieri" ,"Frequency"]))
  expect_equal(fila$Fro, 5111)
  expect_equal(fila$Fro, sum(Amoladeras_int[Amoladeras_int$Canopy=="Open" ,"Frequency"]))
  expect_equal(fila$Ac, 1,tolerance = 1)
  expect_equal(fila$Ao, 6927,tolerance = 1)


})


#--------------------------------------

test_that("Clasification Canopy_effect correct for known species", {

  expect_equal(
    res$Canopy_effect[res$Canopy == "Artemisia_barrelieri"],
    "Facilitative"
  )

  expect_equal(
    res$Canopy_effect[res$Canopy == "Artemisia_campestris"],
    "Neutral"
  )
})



