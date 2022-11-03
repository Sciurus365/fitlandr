#' Normalizer for the input data
#'
#' This is an adapted version of [SparseVFC::norm_vecs()] with `na.rm = TRUE`.
#' @param x The matrix to be normalized. Each row of x represent a vector.
#' @noRd
normalize_vecs <- function(x) {
	x_noNA <- stats::na.omit(x)
	n_noNA <- nrow(x_noNA)
	xmean <- colMeans(x_noNA)
	x_noNA_centered <- x_noNA - matrix(rep(xmean, n_noNA), nrow = n_noNA, byrow = TRUE)
	xscale <- sqrt(sum(x_noNA_centered^2)/n_noNA)
	x_noNA_normalized <- x_noNA_centered / xscale

	x <- (x - matrix(rep(xmean, nrow(x)), nrow = nrow(x), byrow = TRUE))/xscale

	return(
		structure(
			x,
			x_noNA = structure(
				x_noNA_normalized,
				mean = xmean,
				scale = xscale
			),
			mean = xmean,
			scale = xscale
		)
	)
}

# several helpers
normalize_v <- function(v, ref) {
	(v-attr(ref, "mean"))/attr(ref, "scale")
}

normalize_x <- function(x, ref) {
	(x-attr(ref, "mean")[1])/attr(ref, "scale")
}

normalize_y <- function(y, ref) {
	(y-attr(ref, "mean")[2])/attr(ref, "scale")
}

denormalize_v <- function(v, ref) {
	v*attr(ref, "scale")+attr(ref, "mean")
}

scale_up <- function(v, ref) {
	v*attr(ref, "scale")
}

scale_up2 <- function(v, ref) {
	v*attr(ref, "scale")^2
}

denormalize_x <- function(x, ref) {
	x*attr(ref, "scale")+attr(ref, "mean")[1]
}

denormalize_y <- function(y, ref) {
	y*attr(ref, "scale")+attr(ref, "mean")[2]
}
