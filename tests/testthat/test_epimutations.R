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

test_that("returns >= 1 epimutation for all methods with toy genomicratioset", {
	data("genomicratioset")
	expect_gte(nrow(epimutations(cases = genomicratioset, method = "manova")), 1)
	expect_gte(nrow(epimutations(cases = genomicratioset, method = "mlm")), 1)
	expect_gte(nrow(epimutations(cases = genomicratioset, method = "iso.forest")), 1)
	expect_gte(nrow(epimutations(cases = genomicratioset, method = "Mahdist.MCD")), 1)
})

test_that("returns >= 1 epimutation for aref-eshghi with toy genomicratioset", {
	data("genomicratioset")
	expect_gte(nrow(epimutations(cases = genomicratioset, method = "manova",
								 epis_kept="aref-eshghi")), 1)
})


test_that("raises if all post-bumphunter methods are specified", {
	data("genomicratioset")
	expect_error(epimutations(cases = genomicratioset, method=c("manova", "mlm", "iso.forest", "Mahdist.MCD"), "Method must be 'manova', 'mlm','iso.forest','Mahdist.MCD'"))
})



test_that("raises if subject ids are not available", {
	data("genomicratioset")
	nullids <- genomicratioset
	colnames(nullids) <- NULL
	expect_error(epimutations(cases = nullids, num.cpgs = 3))
})



test_that("raises if cases are in controls", {
	data("genomicratioset")
	expect_error(epimutations(cases = genomicratioset, controls = genomicratioset, num.cpgs = 3))
})



test_that("raises num cpgs is 0", {
	data("genomicratioset")
	expect_error(epimutations(cases = genomicratioset, num.cpgs = 0, method = "manova"), "minimum number of cpgs must be 3")
	expect_error(epimutations(cases = genomicratioset, num.cpgs = 0, method = "mlm"), "minimum number of cpgs must be 3")
	expect_error(epimutations(cases = genomicratioset, num.cpgs = 0, method = "iso.forest"), "minimum number of cpgs must be 3")
	expect_error(epimutations(cases = genomicratioset, num.cpgs = 0, method = "Mahdist.MCD"), "minimum number of cpgs must be 3")
})




test_that("returns output when NA's are present", {
	data("genomicratioset")
	subjectwithna <- genomicratioset
	assay(subjectwithna)[,1] <- NA
	expect_warning(class(epimutations(cases = subjectwithna)))
})


