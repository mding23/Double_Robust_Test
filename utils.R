library(Rcpp)
library(RcppArmadillo)

source( "ut_utils.R" )
# Laplace error case
sourceCpp( "compute_T_Laplace.cpp" )

compute_T_w_truncation_w_R_laplace <- function( n, p, m, dat, trunc_type = 0, pct=NA, return_T_sd=FALSE ){
	X = cbind( rep(1, n), dat$X )
	# diff = array( NA, c(p, n, n) ); 
	diff_abs_sum = array( NA, c(p, n, n) )
	mdls = list()
	phi_hat = rep(0, p)
	for( i in 1:p ){
		mdls[[i]] = lad( dat$Y[, i] ~ dat$X, method="EM" )
		phi_hat[i] = sqrt(2) * mean( abs( mdls[[i]]$residuals ) )
	}

	fitted = sapply( 1:p, function(i) mdls[[i]]$fitted.values )

	XtX_inv = solve( t( X ) %*% X / n )
	A_inv_partial_log_f = array(0, c(n, m+2, p))
	for(i in 1:p){
		A_inv_partial_log_f[, , i] = cbind( sqrt(2) * phi_hat[i] * ( X * psi( mdls[[i]]$residuals ) ) %*% XtX_inv, sqrt(2) * abs( mdls[[i]]$residuals ) - phi_hat[i] )
	}

	diff_sign_sum= array( NA, c(p, n, n) ); 

	for( i in 1:p ){
		tmp = outer( fitted[, i], dat$Y[, i], FUN="-" ); atmp = abs(tmp); stmp = sign( tmp )
		# diff[i, ,] = tmp
		diff_abs_sum[i, ,] = ( atmp - diag( atmp ) + t( atmp - diag( atmp ) ) ) * ( sqrt(2) / phi_hat[i]^2 )
		diff_sign_sum[i, ,] = ( stmp - diag( stmp ) ) * ( - sqrt(2) / phi_hat[i] )
	}

	R = array(0, c(p, n, n)) # R^i_(k, k')(beta_hat, sigma_hat)
	for( i in 1:p ){
		R[i, , ] = exp( - phi_hat[i] * diff_abs_sum[i, ,] )
		R[i, , ] = trunc( R=R[i, ,], n=n, trunc_type=trunc_type, pct=pct )
	}

	res = compute_T_cpp_w_truncation_w_R_laplace( n, p, m, X, as.matrix(dat$Y), A_inv_partial_log_f, R, diff_abs_sum, diff_sign_sum )
	T.sigma = apply( res$T_sd, 2, sd ); T = sqrt(n) * res$T / T.sigma;
	if( ! return_T_sd ){
		res$T_sd = NA
	}
	return( list( T_unnorm=res$T, T = T, T.sigma=T.sigma, T_sd=res$T_sd ) )
}

sourceCpp("compute_T.cpp")

compute_T_w_truncation_w_R_Gaussian <- function( n, p, m, dat, trunc_type = 0, pct=NA ){
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
		R[i, , ] = trunc( R=R[i, ,], n=n, trunc_type=trunc_type, pct=pct )
	}

	res = compute_T_cpp_w_truncation_w_R( n, p, m, X, as.matrix(dat$Y), fitted, sigma_hat, A_inv_partial_log_f, R )	
	T.sigma = apply( res$T_sd, 2, sd ); T = sqrt(n) * res$T / T.sigma
	return( list( T = T, T_unnorm=res$T, T.sigma=T.sigma ) )
}

sourceCpp("compute_T_t_dist.cpp")

compute_T_w_truncation_w_R_t_dist <- function( n, p, m, nus, dat, trunc_type = 0, pct=NA ){
	X = cbind( rep(1, n), dat$X )
	mdls = list()
	sigma_hat_sq = rep(0, p)
	for( i in 1:p ){
		mdls[[i]] = lm( dat$Y[, i] ~ dat$X )
		sigma_hat_sq[i] = ( ( nus[i] - 2 ) / nus[i] ) * sum( mdls[[i]]$residuals^2 ) / n
	}
	fitted = sapply( 1:p, function(i) mdls[[i]]$fitted.values )

	partial_log_f = array(0, c(p, m+2, n))	# p*n
	for( i in 1:p ){
		partial_log_f[i, ,] = rbind( t( X * mdls[[i]]$residuals ), ( ( nus[i] - 2 ) / nus[i] ) * mdls[[i]]$residuals^2 - sigma_hat_sq[i] )
	}
	XtX_inv = solve( t( X ) %*% X / n )
	A_inv_partial_log_f = array(0, c(n, m+2, p))
	for(i in 1:p){
		A_inv_partial_log_f[, 1:(m + 1), i] = t( XtX_inv %*% partial_log_f[i, 1:(m + 1), ] )
		A_inv_partial_log_f[, m + 2, i] = partial_log_f[i, m + 2, ]
	}

	R = array(0, c(p, n, n)) # R^i_(k, k')(beta_hat, sigma_hat)
	for( i in 1:p ){
		diff_mat = outer( fitted[, i], dat$Y[, i], FUN="-" )
		R[i, , ] = outer( nus[i] * sigma_hat_sq[i] + diag( diff_mat ) ^ 2, nus[i] * sigma_hat_sq[i] + diag( diff_mat ) ^ 2, FUN="*" ) / ( ( nus[i] * sigma_hat_sq[i] + diff_mat^2 ) * ( nus[i] * sigma_hat_sq[i] + t( diff_mat )^2 ) )
		R[i, , ] = R[i, , ] ^ ( ( nus[i] + 1 ) / 2 )
		R[i, , ] = trunc( R=R[i, ,], n=n, trunc_type=trunc_type, pct=pct )
	}

	res = compute_T_cpp_w_truncation_w_R_t_dist( n, p, m, X, as.matrix(dat$Y), fitted, nus, sigma_hat_sq, A_inv_partial_log_f, R )	
	T.sigma = apply( res$T_sd, 2, sd ); T = sqrt(n) * res$T / T.sigma
	return( list( T = T, T_unnorm=res$T, T.sigma=T.sigma ) )
}

sourceCpp( "compute_T_naive.cpp" )

compute_T_w_truncation_w_R_Gaussian_naive <- function( n, p, m, dat, trunc_type = 0, pct=NA ){
	X = cbind( rep(1, n), dat$X )
	mdls = list()
	sigma_hat = rep(0, p)
	for( i in 1:p ){
		mdls[[i]] = lm( dat$Y[, i] ~ dat$X )
		sigma_hat[i] = sqrt( sum( mdls[[i]]$residuals^2 ) / n )
	}

	fitted = sapply( 1:p, function(i) mdls[[i]]$fitted.values )

	R = array(0, c(p, n, n)) # R^i_(k, k')(beta_hat, sigma_hat)
	for( i in 1:p ){
		R[i, , ] = exp( - outer(dat$Y[, i], dat$Y[, i], FUN="-") * outer( fitted[, i], fitted[, i], FUN = "-" ) / sigma_hat[i]^2 )
		R[i, , ] = trunc( R=R[i, ,], n=n, trunc_type=trunc_type, pct=pct )
	}

	res = compute_T_cpp_w_truncation_w_R_naive( n, p, m, X, as.matrix(dat$Y), R )	
	T.sigma = apply( res$T_sd, 2, sd ); T = sqrt(n) * res$T / T.sigma
	return( list( T = T, T_unnorm=res$T, T.sigma=T.sigma ) )
}

sourceCpp("compute_T_screen.cpp")

compute_T_w_truncation_w_R_Gaussian_screen <- function( n, p, m, dat, trunc_type = 0, pct=NA ){
	X = cbind( rep(1, n), dat$X ) # n * (m+1)
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
		R[i, , ] = trunc( R=R[i, ,], n=n, trunc_type=trunc_type, pct=pct )
	}

	M_vec = rep(0, p); kappa_vec = rep(0, p); screen_vec = rep(0, p)
	for( i in 1:p ){
		tmp_comp = outer( dat$Y[, i], dat$Y[, i], FUN=">" )
		tmp_mat = tmp_comp - t( tmp_comp ) * R[i, ,]
		M_vec[i] = sum( tmp_mat ) / ( n * (n - 1) )
		tmp_prod = t( tmp_comp ) * R[i, ,] * outer( dat$Y[, i], dat$Y[, i], FUN="-" )
		nabla_M = array(0, c(m+2, 1));
		nabla_M[m+2] = - 2 * sum( tmp_prod * outer( fitted[, i], fitted[, i], FUN = "-" ) ) / ( n * (n - 1) * ( sigma_hat[i] ^ 3 ) )
		nabla_M[ 1:(m+1) ] = ( apply( tmp_prod, 1, sum ) - apply( tmp_prod, 2, sum ) ) %*% X / ( n * (n - 1) * ( sigma_hat[i] ^ 2 ) )
		kappa_vec[i] = sd( apply( tmp_mat, 1, sum ) / (n - 1) + A_inv_partial_log_f[, , i] %*% nabla_M )
		screen_vec[i] = as.numeric( abs( sqrt(n) * M_vec[i] / kappa_vec[i] ) >= sqrt( 2 * log( p ) ) ) # 1: misspecified, 0:correctly specified
	}

	res = compute_T_cpp_w_truncation_w_R_screen( n, p, m, X, as.matrix(dat$Y), fitted, sigma_hat, A_inv_partial_log_f, R, screen_vec )	
	T.sigma = apply( res$T_sd, 2, sd ); T = sqrt(n) * res$T / T.sigma
	return( list( T = T, T_unnorm=res$T, T.sigma=T.sigma, M_vec=M_vec, kappa_vec=kappa_vec, screen_vec=screen_vec ) )
}

sourceCpp("compute_T_t_dist_screen.cpp")

compute_T_w_truncation_w_R_t_dist_screen <- function( n, p, m, nus, dat, trunc_type = 0, pct=NA ){
	X = cbind( rep(1, n), dat$X )
	mdls = list()
	sigma_hat_sq = rep(0, p)
	for( i in 1:p ){
		mdls[[i]] = lm( dat$Y[, i] ~ dat$X )
		sigma_hat_sq[i] = ( ( nus[i] - 2 ) / nus[i] ) * sum( mdls[[i]]$residuals^2 ) / n
	}
	fitted = sapply( 1:p, function(i) mdls[[i]]$fitted.values )

	partial_log_f = array(0, c(p, m+2, n))	# p*n
	for( i in 1:p ){
		partial_log_f[i, ,] = rbind( t( X * mdls[[i]]$residuals ), ( ( nus[i] - 2 ) / nus[i] ) * mdls[[i]]$residuals^2 - sigma_hat_sq[i] )
	}
	XtX_inv = solve( t( X ) %*% X / n )
	A_inv_partial_log_f = array(0, c(n, m+2, p))
	for(i in 1:p){
		A_inv_partial_log_f[, 1:(m + 1), i] = t( XtX_inv %*% partial_log_f[i, 1:(m + 1), ] )
		A_inv_partial_log_f[, m + 2, i] = partial_log_f[i, m + 2, ]
	}

	R = array(0, c(p, n, n)) # R^i_(k, k')(beta_hat, sigma_hat)
	for( i in 1:p ){
		diff_mat = outer( fitted[, i], dat$Y[, i], FUN="-" )
		R[i, , ] = outer( nus[i] * sigma_hat_sq[i] + diag( diff_mat ) ^ 2, nus[i] * sigma_hat_sq[i] + diag( diff_mat ) ^ 2, FUN="*" ) / ( ( nus[i] * sigma_hat_sq[i] + diff_mat^2 ) * ( nus[i] * sigma_hat_sq[i] + t( diff_mat )^2 ) )
		R[i, , ] = R[i, , ] ^ ( ( nus[i] + 1 ) / 2 )
		R[i, , ] = trunc( R=R[i, ,], n=n, trunc_type=trunc_type, pct=pct )
	}

	M_vec = rep(0, p); kappa_vec = rep(0, p); screen_vec = rep(0, p)
	for( i in 1:p ){
		tmp_comp = outer( dat$Y[, i], dat$Y[, i], FUN=">" )
		tmp_mat = tmp_comp - t( tmp_comp ) * R[i, ,]
		M_vec[i] = sum( tmp_mat ) / ( n * (n - 1) )

		tmp_prod = - t( tmp_comp ) * R[i, ,]
		nabla_M = array(0, c(m+2, 1));
		st_mat = outer( fitted[, i], dat$Y[, i], FUN = "-" )
		st_inv = 1 / ( nus[i] * sigma_hat_sq[i] + ( st_mat ^ 2 ) )
		nabla_M[m+2] = ( (nus[i] + 1) * nus[i] / 2 ) * sum( tmp_prod * ( outer( diag( st_inv ), diag(st_inv), FUN="+" ) - st_inv - t( st_inv ) ) ) / ( n * (n - 1) )
		tmp_st_st_inv = st_mat * st_inv; 
		tmp_st_st_inv_diff = diag( tmp_st_st_inv ) - tmp_st_st_inv
		nabla_M[ 1:(m+1) ] = ( nus[i] + 1 ) * ( apply( tmp_prod * tmp_st_st_inv_diff, 1, sum ) + apply( tmp_prod * t(tmp_st_st_inv_diff), 2, sum ) ) %*% X / ( n * (n - 1) )
		kappa_vec[i] = sd( apply( tmp_mat, 1, sum ) / (n - 1) + A_inv_partial_log_f[, , i] %*% nabla_M )
		screen_vec[i] = as.numeric( abs( sqrt(n) * M_vec[i] / kappa_vec[i] ) >= sqrt( 2 * log( p ) ) ) # 1: misspecified, 0:correctly specified
	}

	res = compute_T_cpp_w_truncation_w_R_t_dist_screen( n, p, m, X, as.matrix(dat$Y), fitted, nus, sigma_hat_sq, A_inv_partial_log_f, R, screen_vec )	
	T.sigma = apply( res$T_sd, 2, sd ); T = sqrt(n) * res$T / T.sigma
	return( list( T = T, T_unnorm=res$T, T.sigma=T.sigma, M_vec=M_vec, kappa_vec=kappa_vec, screen_vec=screen_vec ) )
}