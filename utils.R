library(Rcpp)
library(RcppArmadillo)

sourceCpp("compute_T.cpp")

compute_T_w_R_Gaussian <- function( n, p, m, dat ){
	X = cbind( rep(1, n), dat$X )
	mdls = list()
	sigma_hat = rep(0, p)
	for( i in 1:p ){
		mdls[[i]] = lm( dat$Y[, i] ~ dat$X )
		sigma_hat[i] = sqrt( sum( mdls[[i]]$residuals^2 ) / n )
	}

	fitted = sapply( 1:p, function(i) mdls[[i]]$fitted.values )

	partial_log_f = array(0, c(p, m+2, n))	# p*n
	for( i in 1:p ){
		partial_log_f[i, ,] = rbind( t( X * mdls[[i]]$residuals ) / sigma_hat[i]^2, mdls[[i]]$residuals^2 / sigma_hat[i]^3 - 1 / sigma_hat[i] )
	}

	J_inv = array( 0, c(p, m+2, m+2) )
	XtX_inv = solve( t( X ) %*% X / n )
	for( i in 1:p ){
		J_inv[i, m+2, m+2] = sigma_hat[i]^2 / 2
		J_inv[i, 1:(m+1), 1:(m+1)] = XtX_inv * sigma_hat[i]^2
	}
	A_inv_partial_log_f = array(0, c(n, m+2, p))
	for(i in 1:p){
		A_inv_partial_log_f[, , i] = t( J_inv[i, ,] %*% partial_log_f[i, ,] )
	}

	R = array(0, c(p, n, n)) # R^i_(k, k')(beta_hat, sigma_hat)
	for( i in 1:p ){
		R[i, , ] = exp( - outer(dat$Y[, i], dat$Y[, i], FUN="-") * outer( fitted[, i], fitted[, i], FUN = "-" ) / sigma_hat[i]^2 )
	}

	res = compute_T_cpp_w_R( n, p, m, X, as.matrix(dat$Y), fitted, sigma_hat, A_inv_partial_log_f, R )	
	T.sigma = apply( res$T_sd, 2, sd ); T = sqrt(n) * res$T / T.sigma
	return( list( T = T, T_unnorm=res$T, T.sigma=T.sigma ) )
}