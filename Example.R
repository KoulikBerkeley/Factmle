# Packages for main function
library(rARPACK)
library(psych)
library(Matrix)
library(psych)

# packages for simulation
library(SoDA)

## Tridiagonal example

# problem_parameters

dim_tridiag = 15
Lambda_rank = 3
sample_size = 10000

# simulation details

Psi_inverse = triDiag(5,1,1, dim_tridiag) 
Psi = solve(Psi_inverse)
Lambda = mvrnorm(n = 15, mu = rep(0, Lambda_rank), Sigma = diag(1, Lambda_rank))

data_tridiag = mvrnorm(n = sample_size, mu = rep(0, dim_tridiag), Sigma = Psi + Lambda%*%t(Lambda))

# fitting the data

rank_fit = 3
out_tridiag = FA_tridiag(data_tridiag, curr_rank = 3)

Psi_out = out_tridiag$Psi
Psi_inverse_out = out_tridiag$Psi_inverse

## Block-diagonal example

# Problem parameters
dim_blk = 10
Lambda_rank = 3
sample_size = 10000

# simulating the data

mat1 = matrix(mvrnorm(25, mu = 0, Sigma = 1), nrow = 5)
mat1 = mat1 %*% t(mat1)

mat2 = matrix(mvrnorm(25, mu = 0, Sigma = 1), nrow = 5)
mat2 = mat2 %*% t(mat2)

Psi = bdiag(mat1, mat2)
Lambda = mvrnorm(n = 10, mu = rep(0, Lambda_rank), Sigma = diag(1, Lambda_rank))

data_blkdiag = mvrnorm(n = sample_size, mu = rep(0, dim_blk), Sigma = Psi + Lambda%*%t(Lambda))

# fitting the data

out_blk = FA_blk(data_blkdiag, curr_rank = 3, blk_number = 2)


