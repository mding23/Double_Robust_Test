#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List compute_T_cpp_w_truncation_w_R_t_dist_screen( int const n, int const p, int const m, const arma::mat& X, const arma::mat& Y, const arma::mat& fitted, const arma::vec& nus, const arma::vec& sigmaHat_sq, const arma::cube& A_inv_partial_log_f, const arma::cube& R, const arma::vec& screen_vec ){
	// for Y continuous case
	int n_H = p * (p-1) / 2;
	arma::vec T(n_H), Y_diff(p), res_ss(p), res_tt(p), res_st(p), res_ts(p), nu_sig_sq(p), ss_inv(p), tt_inv(p), st_inv(p), ts_inv(p);
	T.zeros(); Y_diff.zeros();
	arma::mat T_sd( n, n_H ); T_sd.zeros(); 
	double rho_ij_st;
	arma::colvec v(m+2);
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
				res_ss[i] = fitted(s, i) - Y(s, i);
				res_tt[i] = fitted(t, i) - Y(t, i);
				res_st[i] = fitted(s, i) - Y(t, i);
				res_ts[i] = fitted(t, i) - Y(s, i);
				nu_sig_sq[i] = nus[i] * sigmaHat_sq[i];
				ss_inv[i] = 1 / ( nu_sig_sq[i] + pow( res_ss[i], 2 ) );
				tt_inv[i] = 1 / ( nu_sig_sq[i] + pow( res_tt[i], 2 ) );
				st_inv[i] = 1 / ( nu_sig_sq[i] + pow( res_st[i], 2 ) );
				ts_inv[i] = 1 / ( nu_sig_sq[i] + pow( res_ts[i], 2 ) );
			}
			k = -1;
			for( int i=0; i < (p-1); i++ ){
				for( int j=i+1; j < p; j++ ){
					k += 1;
					rho_ij_st = ( ( Y_diff[i] > 0 ) ? -1 : R(i, s, t) ) * ( ( Y_diff[j] > 0 ) ? -1 : R(j, s, t) );
					T[k] += rho_ij_st; T_sd(s, k) += rho_ij_st; T_sd(t, k) += rho_ij_st;
					v.zeros();
					if( ( Y_diff[i] < 0 ) & ( screen_vec[j] > 0.5 ) ){
						v.subvec(0, m) = ( nus[i] + 1 ) * trans( ( res_ss[i] * ss_inv[i] - res_st[i] * st_inv[i] ) * X.row(s) + ( res_tt[i] * tt_inv[i] - res_ts[i] * ts_inv[i] ) * X.row(t) );
						v[m+1] = ( ( nus[i] + 1 ) * nus[i] / 2 ) * ( ss_inv[i] + tt_inv[i] - st_inv[i] - ts_inv[i] );
						E_ij_i.col(k) -= R(i, s, t) * ( ( Y_diff[j] > 0 ) ? 1 : - R(j, s, t) ) * v;
					}
					v.zeros();
					if( ( Y_diff[j] < 0 ) & ( screen_vec[i] > 0.5 ) ){
						v.subvec(0, m) = ( nus[j] + 1 ) * trans( ( res_ss[j] * ss_inv[j] - res_st[j] * st_inv[j] ) * X.row(s) + ( res_tt[j] * tt_inv[j] - res_ts[j] * ts_inv[j] ) * X.row(t) );
						v[m+1] = ( ( nus[j] + 1 ) * nus[j] / 2 ) * ( ss_inv[j] + tt_inv[j] - st_inv[j] - ts_inv[j] );
						E_ij_j.col(k) -= R(j, s, t) * ( ( Y_diff[i] > 0 ) ? 1 : - R(i, s, t) ) * v;
					}
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

