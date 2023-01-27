psi <- function( x, tau=0.5 ){
	return( tau - ( x < 0 ) )
}

trunc <- function( R, n, trunc_type, pct=NA ){
	if( trunc_type == 0 ){
		return(R)
	}
	if( trunc_type == 1 ){
		pct = 1 - 0.01 * ( 200 / n )
		q = quantile( R, pct )
		R = replace( R, R > q, q )		
	}
	if( trunc_type == 2 ){
		q = quantile( R, pct )
		R = replace( R, R > q, q )			
	}
	return(R)
}

G <- function(t, df=0){
	return( 2 * ( 1 - pnorm(t) ) )
	# return( 2 * ( 1 - pt(q=t, df=df) ) )
}

fdr <- function( T, n, p, alpha ){
	c_p = sqrt( 4 * log( p ) - 2 * log( log( p ) ) ); n_H = (p - 1) * p / 2
	b_p = 2 * sqrt( log(p) )
	abs_T_sort = sort( abs(T) )
	tmp = G( abs_T_sort, df = n - 1 ) * n_H / ( n_H:1 ) <= alpha
	idx = which( tmp )[1]
	if( is.na(idx) ){
		t_hat = b_p
		reach_threshold = 1
	}else{
		if( c_p < abs_T_sort[idx] ){
			t_hat = b_p
			reach_threshold = 1
		}else{
			t_hat = abs_T_sort[idx]
			reach_threshold = 0
		}
	}
	return( list( t_hat=t_hat, reach_threshold=reach_threshold ) )
}

eval_fdr <- function( T, t_hat, true_corr){
	est_corr = abs(T) >= t_hat
	FDP = sum( est_corr * ( 1 - true_corr ) ) / max( 1, sum(est_corr) )
	FN = sum( ( 1 - est_corr ) * true_corr )
	power = 1 - FN / sum( true_corr )
	return( list( FDP=FDP, FN=FN, power=power ) )
}

ij2idx <- function(p, i, j){
	if( ( i == j ) | ( i < 1 ) | ( i > p ) | ( j < 1 ) | ( j > p) ){
		return(0)
	}
	n_H = p * (p - 1) / 2
	mat = array( 0, c(p, p) )
	mat[ lower.tri( mat, diag=FALSE ) ] = 1:n_H
	return( mat[j, i] )
}

idx2ij <- function( p, idx ){
	n_H = p * (p - 1) / 2
	if( ( idx < 1 ) | ( idx > n_H ) ){
		return( list( i=0, j=0 ) )
	}
	mat = array( 0, c(p, p) )
	mat[ lower.tri( mat, diag=FALSE ) ] = 1:n_H
	tmp = which( mat == idx, arr.ind = TRUE )
	return( list( i=tmp[2], j=tmp[1] ) )
}