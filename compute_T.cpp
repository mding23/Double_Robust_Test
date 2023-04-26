#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List compute_T_cpp_w_R( int const n, int const p, int const m, const arma::mat& X, const arma::mat& Y, const arma::mat& fitted, const arma::vec& sigmaHat, const arma::cube& A_inv_partial_log_f, const arma::cube& R ){
	// for Y continuous case
	int n_H = p * (p-1) / 2;
	arma::vec T(n_H), Y_diff(p), Y_gt(p), fitted_diff(p);
	T.zeros(); Y_diff.zeros(); Y_gt.zeros(); fitted_diff.zeros();
	arma::mat T_sd( n, n_H ); T_sd.zeros(); 
	double rho_ij_st;
	arma::colvec v(m+2); v.zeros();
	arma::mat E_ij_i(m+2, n_H), E_ij_j(m+2, n_H);
	E_ij_i.zeros(); E_ij_j.zeros();
	int k;
	for( int s=0; s<n; s++ ){
		for( int t=0; t<n; t++ ){
			if( s == t ){
				continue;
			}
			for( int i=0; i<p; i++ ){
				Y_diff[i] = Y(s, i) - Y(t, i);
				// Y_gt[i] = Y_diff[i] > 0;
				fitted_diff[i] = ( fitted(s, i) - fitted(t, i) ) / sigmaHat[i];
			}
			v.subvec(0, m) = trans( X.row(t) - X.row(s) );
			k = -1;
			for( int i=0; i < (p-1); i++ ){
				for( int j=i+1; j < p; j++ ){
					k += 1;
					rho_ij_st = ( ( Y_diff[i] > 0 ) ? -1 : R(i, s, t) ) * ( ( Y_diff[j] > 0 ) ? -1 : R(j, s, t) );
					T[k] += rho_ij_st; T_sd(s, k) += rho_ij_st; T_sd(t, k) += rho_ij_st;
					v[m+1] = 2 * fitted_diff[i];
					E_ij_i.col(k) += ( ( Y_diff[i] < 0 ) ? ( R(i, s, t) * Y_diff[i] * ( ( Y_diff[j] > 0 ) ? -1 : R(j, s, t) ) / ( pow( sigmaHat[i], 2 ) ) ) : 0 ) * v;
					v[m+1] = 2 * fitted_diff[j];
					E_ij_j.col(k) += ( ( Y_diff[j] < 0 ) ? ( R(j, s, t) * Y_diff[j] * ( ( Y_diff[i] > 0 ) ? -1 : R(i, s, t) ) / ( pow( sigmaHat[j], 2 ) ) ) : 0 ) * v;	
				}
			}
		}
	}
	k = -1;
	for( int i=0; i < (p-1); i++ ){
		for( int j=i+1; j < p; j++ ){
			k += 1;
			T_sd.col(k) += ( A_inv_partial_log_f.slice(i) * E_ij_i.col(k) + A_inv_partial_log_f.slice(j) * E_ij_j.col(k) ) / n;
		}
	}
	T = T / ( n * (n-1) );
	T_sd = T_sd / ( n - 1 );
	return List::create( Named("T") = T, Named("T_sd") = T_sd );
}

/*** R
*/

