
test_that(".PerStockArg() broadcasts a bare value to every stock", {
  out <- .PerStockArg("catch", c("A", "B", "C"), "catch")
  expect_equal(out, list(A = "catch", B = "catch", C = "catch"))
})

test_that(".PerStockArg() broadcasts a bare list value (not per-stock) to every stock", {
  LWT <- list(Chist = 1, Index = 2)
  out <- .PerStockArg(LWT, c("A", "B"), list())
  expect_equal(out$A, LWT)
  expect_equal(out$B, LWT)
})

test_that(".PerStockArg() takes per-stock values from a named list keyed by stock name", {
  out <- .PerStockArg(list(A = "effort"), c("A", "B"), "catch")
  expect_equal(out$A, "effort")
  expect_equal(out$B, "catch") # gap filled with default
})

test_that(".PerStockArg() per-stock list can hold list-valued entries", {
  arg <- list(A = list(Chist = 5), B = list(Chist = 10))
  out <- .PerStockArg(arg, c("A", "B"), list())
  expect_equal(out$A, list(Chist = 5))
  expect_equal(out$B, list(Chist = 10))
})

test_that(".PerStockArg() ignores stocks argument order and subsets to stock_names", {
  out <- .PerStockArg(list(B = "effort", Z = "ignored"), c("A", "B"), "catch")
  expect_equal(names(out), c("A", "B"))
  expect_equal(out$A, "catch")
  expect_equal(out$B, "effort")
})

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

test_that(".AsStockKeyedList() fills NULL input with one NULL per stock", {
  out <- .AsStockKeyedList(NULL, c("A", "B"), "fleet")
  expect_equal(names(out), c("A", "B"))
  expect_true(is.null(out$A) && is.null(out$B))
})

test_that(".AsStockKeyedList() keys a bare object to the first stock name", {
  f <- methods::new("fleet")
  out <- .AsStockKeyedList(f, c("A", "B"), "fleet")
  expect_equal(names(out), "A")
  expect_identical(out$A, f)
})
