test_that("raises if cases arg is missing", {
	expect_error(epimutations(), "Cases argument is missing.")
})

test_that("raises if ExpressionSet cases is empty", {
	expect_error(
		epimutations(cases = Biobase::ExpressionSet()),
		"If controls is missing, cases must contain at least 3 samples."
	)
})

test_that("raises if GenomicRatioSet and less than (1 case, 2 control) samples", {
	data("genomicratioset")
	expect_error(
		epimutations(cases = new("GenomicRatioSet")),
		"If controls is missing, cases must contain at least 3 samples."
	)
	expect_error(
		epimutations(cases = genomicratioset[, 1]),
		"If controls is missing, cases must contain at least 3 samples."
	)
	expect_error(
		epimutations(cases = genomicratioset[, 2:3]),
		"If controls is missing, cases must contain at least 3 samples."
	)
	expect_error(
		epimutations(cases = new("GenomicRatioSet"), controls = genomicratioset),
		"Cases must contain at least 1 sample."
	)
	expect_error(
		epimutations(cases = genomicratioset, controls = new("GenomicRatioSet")),
		"Controls must contain at least 2 samples."
	)
	expect_error(
		epimutations(cases = genomicratioset, controls = genomicratioset[, 3]),
		"Controls must contain at least 2 samples."
	)
})


test_that("returns zero rows if bumphunter finds no bumps", {
	data("genomicratioset")
	actual <- nrow(epimutations(cases = genomicratioset[1,]))
	expect_equal(actual, 0)
})

test_that("returns with toy genomicratioset", {
	data("genomicratioset")
	out <- epimutations(cases = genomicratioset, num.cpgs = 3, method = "iso.forest")
	expect_equal(nrow(out), 0)
})



test_that("raises if all post-bumphunter methods are specified", {
	data("genomicratioset")
	expect_error(epimutations(cases = genomicratioset, method=c("manova", "mlm", "iso.forest", "Mahdist.MCD"), "Method must be 'manova', 'mlm','iso.forest','Mahdist.MCD'"))
})



test_that("raises if subject ids are not avaialble", {
	data("genomicratioset")
	nullids <- genomicratioset
	colnames(nullids) <- NULL
	expect_error(epimutations(cases = nullids, num.cpgs = 3, method = "iso.forest"))
})



test_that("raises if cases are in controls", {
	data("genomicratioset")
	expect_error(epimutations(cases = genomicratioset, controls = genomicratioset, num.cpgs = 3, method = "iso.forest"))
})




