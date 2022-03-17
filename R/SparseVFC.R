#' Sparse Vector Field Consensus
#'
#' Translated from codes in https://github.com/jiayi-ma/VFC
#' Jiayi Ma, Ji Zhao, Jinwen Tian, Xiang Bai, and Zhuowen Tu. "Regularized Vector Field Learning with Sparse Approximation for Mismatch Removal", Pattern Recognition, 46(12), pp. 3519-3532, 2013.
#' @export
SparseVFC <- function(X, Y, MaxIter = 500, gamma = 0.9, beta = 0.1,
											lambda = 3, theta = 0.75, a = 10, ecr = 1e-5,
											minP = 1e-5, silent = TRUE){
	if(!silent) message("Start mismatch removal...\n")
	N <- nrow(Y)
	D <- ncol(Y)

	# Construct kernel matrix K
	M <- 16
	tmp_X <- unique(X)
	idx <- sample.int(nrow(tmp_X))
	idx <- idx[1:min(M, nrow(tmp_X))]
	ctrl_pts <- tmp_X[idx,]
	K <- con_K(ctrl_pts, ctrl_pts, beta)
	U <- con_K(X, ctrl_pts, beta)
	M <- nrow(ctrl_pts)

	# Initialization
	V <- matrix(0, nrow = N, ncol = D)
	iter <- 1
	tecr <- 1
	C <- matrix(0, nrow = M, ncol = D)
	E <- 1
	sigma2 <- sum((Y-V)^2)/(N*D)

	while(iter < MaxIter & tecr > ecr){
		# E-step.
		E_old <- E
		temp_PE <- get_P(Y, V, sigma2, gamma, a)
		P <- temp_PE$P
		E <- temp_PE$E
		tecr <- abs((E-E_old)/E)
		if(!silent) message(
			sprintf('iterate: %dth, gamma: %f, the energy change rate: %f, sigma2=%f\n', iter, gamma, tecr, sigma2)
		)

		# M-step. Solve linear system for C.
		P <- pmax(P, minP)
		C <- solve((t(U) * matrix(rep(P, M), nrow = M, byrow = TRUE)) %*% U + lambda*sigma2*K) %*%
			((t(U) * matrix(rep(P, M), nrow = M, byrow = TRUE)) %*% Y)

		# Update V and sigma^2
		V <- U %*% C
		Sp <- sum(P)
		sigma2 <- sum(P*rowSums((Y-V)^2)) /(Sp*D)

		# Update gamma
		numcorr <- sum(P > theta)
		gamma <- numcorr/nrow(X)
		if(gamma > 0.95) gamma <- 0.95
		if(gamma < 0.05) gamma <- 0.05

		iter <- iter+1
	}

	if(!silent) message("Removing outliers succesfully completed.")

	result <- list(
		X = ctrl_pts,
		Y = Y,
		beta = beta,
		V = V,
		C = C,
		P = P,
		VFCIndex = which(P > theta),
		sigma2 = sigma2
	)

	class(result) <- c("VFC")
	return(
		result
	)
}

#' Construct the kernel K
#'
#' K[i, j] = k(x[i,], y[j,]) = exp(-beta*||x[i,]-y[j,]||^2)
con_K <- function(x, y, beta){
	n <- nrow(x)
	m <- nrow(y)
	if(ncol(x) != ncol(y)) stop("ncol(x) != ncol(y)")
	d <- ncol(x)

	if(n == 1){
		K <- matrix(sum((x - y)^2)^(1/2))
	}else{
		K <- as.matrix(purrr::quietly(pdist::pdist)(x, y)$result)
	}
	K <- exp(-beta*(K^2))

	return(K)
}

#' Estimate the posterior probability and part of the energy
#'
get_P <- function(Y, V, sigma2, gamma, a){
	D <- ncol(Y)
	temp1 <- exp(-rowSums((Y-V)^2)/(2*sigma2))
	temp2 <- (2*pi*sigma2)^(D/2)*(1-gamma)/(gamma*a)
	P <- temp1/(temp1+temp2)
	E <- sum(P * rowSums((Y-V)^2)/(2*sigma2)) + sum(P)*log(sigma2)*D/2
	return(list(P = P,E = E))
}


norm2 <- function(x, y){
	n <- nrow(x)
	m <- nrow(y)

	xm <- colMeans(x)
	ym <- colMeans(y)

	x <- x - matrix(rep(xm, n), nrow = n, byrow = TRUE)
	y <- y - matrix(rep(ym, m), nrow = m, byrow = TRUE)

	xscale <- sqrt(sum(x^2)/n)
	yscale <- sqrt(sum(y^2)/m)

	X <- x/xscale
	Y <- y/yscale

	return(
		list(
			X = X, Y = Y, xscale = xscale, yscale = yscale
		)
	)
}
