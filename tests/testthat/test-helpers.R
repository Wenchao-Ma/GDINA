test_that("bdiagMatrix combines rectangular blocks and fill values", {

	out <- bdiagMatrix(list(matrix(1:4, 2), matrix(5:6, 1)), fill = -1)
	expected <- matrix(c(1, 2, -1,
											 3, 4, -1,
											 -1, -1, 5,
											 -1, -1, 6), nrow = 3)

	expect_equal(out, expected)

})


test_that("cjoint pads shorter inputs and strips generated names", {

	out <- cjoint(a = 1:2, b = matrix(3:8, nrow = 3), fill = 99)
	expected <- data.frame(c(1, 2, 99), c(3, 4, 5), c(6, 7, 8))
	names(expected) <- NULL

	expect_equal(out, expected)
	expect_null(cjoint())

})


test_that("rowMatch returns matching count and row numbers", {

	df <- data.frame(a = c(1, 2, 1), b = c("x", "y", "x"))

	expect_equal(rowMatch(df, c(1, "x")), list(count = 2L, row.no = c(1L, 3L)))
	expect_equal(rowMatch(df, c(9, "z")), list(count = 0L, row.no = integer(0)))

})


test_that("unique_only returns sorted values occurring exactly once", {

	expect_equal(unique_only(c(4, 2, 3, 5, 4, 4, 4)), c(2, 3, 5))
	expect_equal(unique_only(c("b", "a", "b", "c")), c("a", "c"))

})


test_that("unrestrQ expands category rows within each sequential item", {

	Qc <- matrix(c(1, 1, 1, 0, 0,
								 1, 2, 0, 1, 0,
								 2, 1, 0, 0, 1,
								 2, 2, 1, 0, 1), ncol = 5, byrow = TRUE)
	expected <- matrix(c(1, 1, 1, 1, 0,
											 1, 2, 1, 1, 0,
											 2, 1, 1, 0, 1,
											 2, 2, 1, 0, 1), ncol = 5, byrow = TRUE)

	expect_equal(unrestrQ(Qc), expected)

})


test_that("ClassRate computes attribute and vector classification rates", {

	att1 <- matrix(c(1, 0,
									 1, 1,
									 0, 0), ncol = 2, byrow = TRUE)
	att2 <- matrix(c(1, 1,
									 1, 1,
									 0, 0), ncol = 2, byrow = TRUE)

	out <- ClassRate(att1, att2)

	expect_equal(out$PCA, 5 / 6)
	expect_equal(out$PCV, c(1, 2 / 3))
	expect_error(ClassRate(att1, att2[-1, ]), "same dimensions")

	perfect <- ClassRate(att1, att1)
	expect_equal(perfect$PCA, 1)
	expect_equal(perfect$PCV, c(1, 1))

})
