
test_that(".AsStockList() wraps a bare stock object using its Name slot", {
  s <- methods::new("stock")
  s@Name <- "MyStock"
  out <- .AsStockList(s)
  expect_equal(names(out), "MyStock")
  expect_identical(out$MyStock, s)
})

test_that(".AsStockList() passes an already-named list through unchanged", {
  s1 <- methods::new("stock")
  lst <- list(A = s1)
  expect_identical(.AsStockList(lst), lst)
})

test_that(".AsStockList() rejects an unnamed list", {
  s1 <- methods::new("stock")
  expect_error(.AsStockList(list(s1)))
})

test_that(".OneData() takes one data object or a list of one", {
  d <- MSEtool::Data()
  expect_identical(.OneData(d), d)
  expect_identical(.OneData(list(A = d)), d)
  expect_error(.OneData(list(d, d)), "one data set")
  expect_error(.OneData("x"), "data")
})

test_that(".RenameSPFrom() renames stock names in values and weight names", {
  renames <- c(Old1 = "A", Old2 = "B")
  expect_equal(.RenameSPFrom("Old1", renames), "A")
  expect_equal(.RenameSPFrom(c(Old1 = 0.5, Old2 = 0.5), renames), c(A = 0.5, B = 0.5))
  expect_equal(.RenameSPFrom("Other", renames), "Other")
  expect_equal(.RenameSPFrom(2, renames), 2)
  expect_null(.RenameSPFrom(NULL, renames))
})
