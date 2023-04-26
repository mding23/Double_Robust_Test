source("utils.R")

load("indep_data.RData")
load("correlatad_data.RData")

m = 2; p = 2; n = 200

indep_res = compute_T_w_R_Gaussian( n=n, p=p, m=m, dat=indep_data )
# > indep_res
# $T
#        [,1]
# [1,] 0.5265

# $T_unnorm
#            [,1]
# [1,] 0.05050615

# $T.sigma
# [1] 1.356628


correlated_res = compute_T_w_R_Gaussian( n=n, p=p, m=m, dat=correlatad_data )
# > correlated_res
# $T
#          [,1]
# [1,] 5.663558

# $T_unnorm
#           [,1]
# [1,] 0.4375647

# $T.sigma
# [1] 1.092617



