library(testthat)

window = IRanges(c(1,11), c(10, 20))

context("Test average_in_window")

# test will be added later
test_that("Test average_in_window", {
    ir = IRanges(1, 10)
	mtch = as.matrix(findOverlaps(window, ir))
	x = average_in_window(window, ir, mtch, 1, mean_mode = "w0")
	x = as.vector(x)
    expect_that(x, equals(1))

    ir = IRanges(1, 20)
    mtch = as.matrix(findOverlaps(window, ir))
	x = average_in_window(window, ir, mtch, 1, mean_mode = "w0")
	x = as.vector(x)
    expect_that(x, equals(c(1, 1)))

    ir = IRanges(6, 15)
	mtch = as.matrix(findOverlaps(window, ir))
	x = average_in_window(window, ir, mtch, 1, mean_mode = "w0")
	x = as.vector(x)
    expect_that(x, equals(c(0.5, 0.5)))

    ir = IRanges(c(1, 11), c(5, 15))
	mtch = as.matrix(findOverlaps(window, ir))
	x = average_in_window(window, ir, mtch, c(1, 1), mean_mode = "w0")
	x = as.vector(x)
    expect_that(x, equals(c(1/2, 1/2)))

    ir = IRanges(c(1, 11), c(5, 15))
	mtch = as.matrix(findOverlaps(window, ir))
	x = average_in_window(window, ir, mtch, c(2, 3), mean_mode = "w0")
	x = as.vector(x)
    expect_that(x, equals(c(1, 1.5)))

    x = average_in_window(window, ir, mtch, list(c(1, 1), c(2, 3)), mean_mode = "w0")
	x = as.vector(x)
    expect_that(x, equals(c(0.5, 0.5, 1, 1.5)))

    ir = IRanges(c(1, 2), c(3, 5))
    mtch = as.matrix(findOverlaps(window, ir))
	x = average_in_window(window, ir, mtch, c(1, 1), mean_mode = "w0")
	x = as.vector(x)
    expect_that(x, equals(0.5))

    ir = IRanges(21, 22)
    mtch = as.matrix(findOverlaps(window, ir))
	x = average_in_window(window, ir, mtch, c(1, 1), mean_mode = "w0")
	x = as.vector(x)
    expect_that(length(x), equals(0))

})
