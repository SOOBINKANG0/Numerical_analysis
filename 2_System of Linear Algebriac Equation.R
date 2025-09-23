## textbook: Numerical Methods in Engineering with python 3
## I tried to write by myself but I was supported by Chatgpt to get un hint

## package
install.packages("pracma")
library(pracma)

## 1. Gauss elimination

## Example 2.1
A = matrix(c(4, -2, 1, -2, 4, -2, 1, -2, 4), byrow = T, nrow = 3)
B = matrix(c(11, -16 , 17), byrow = T, ncol = 1)

Sol = solve(A) %*% B # 2, -1, 3

## 안쪽루프1 = i(행)
if(A[2,1] != 0){
  lambda = -A[1,1]/A[2,1]
  A[2, ] <- lambda*A[2,] + A[1,]
  B[2, ] <- lambda*B[2,] + B[1,]
}

if(A[3,1] != 0){
  lambda = -A[1,1]/A[3,1]
  A[3, ] <- lambda*A[3,] + A[1,]
  B[3, ] <- lambda*B[3,] + B[1,]
}

## 바깥쪽 루프 = j(열)
lambda = -A[2,2]/A[3,2]
if(A[2,2] != 0){
  A[3, ] <- lambda*A[3,] + A[2,]
  B[3, ] <- lambda*B[3,] + B[2,]
}

## back substitution

B[3,] = B[3,]/A[3,3]
B[2,] = (B[2,]-A[2,3]*B[3,])/A[2,2]
B[1,] = (B[1,]-A[1,2]*B[2,]-A[1,3]*B[3,])/A[1,1]
print(B)

##
A = matrix(c(4, -2, 1, -2, 4, -2, 1, -2, 4), byrow = T, nrow = 3)
B = matrix(c(11, -16 , 17), byrow = T, ncol = 1)

## 일반화
Gauss_eli = function(A, B){
  
  n = nrow(A)
  x = numeric(n)
  count = 0
  
  ## Elimination part(행: i, 열: k, A파트) 
  for(k in 1:(n-1)){ # 바깥루프, 열 # 마지막 행은 할필요가 없음
    
    for(i in (k+1):n){ # 안쪽루프, 행
      
      if(A[i,k] != 0){
        lambda = A[i,k] / A[k,k]
        
        A[i, ] = A[i, ] - lambda*A[k,] # k가 기준이 된다!(k_11부터...)
        
        B[i]  = B[i] - lambda*B[k]    
        
      }
      
    }
    
  }
  ## back substitution(B 파트)
  x[n] = B[n] / A[n,n]
  
  for(k in (n-1):1){
    sum_of_product = sum(A[k, (k + 1):n] * x[(k + 1):n])
    
    x[k] = (B[k] - sum_of_product)/A[k,k] 
  } 
  
  return(x)
}

Gauss_eli(A,B)

##Example 2.5
A2 = matrix(c(1,4,1,1,6,-1,2,-1,2), byrow = T, ncol = 3)
B2 = matrix(c(7,13,5), byrow = T, ncol = 1)

BBB = matrix(c(4, 3, 1, 4,
               2, 1, 0 ,1,
               1, 0, 1, 2,
               4, 1, 2, -1), byrow = T, nrow = 4)

CCC = matrix(c(1,3,2,4), nrow = 4)

solve(BBB) %*% CCC
Gauss_eli(BBB, CCC)

####################################################################

## 2. LU decomposition
AA = matrix(c(1, 4, 1, 1, 6, -1, 2, -1, 2), byrow = T, ncol = 3)
BB = matrix(c(4,-2, 2,-2, 2, -4, 2, -4,11), byrow = T, ncol = 3)

## 
L = lu(AA)$L; U = lu(AA)$U
L2 = lu(BB)$L; U2 = lu(BB)$U

## elimination(decomposition)
## 행:i(안쪽 루프), 열:k(바깥 루프)

LU_decom = function(A){
  
  n = ncol(A)
  res = list(L = matrix(0,ncol = ncol(A), nrow = nrow(A)), 
             U = matrix(0,ncol = ncol(A), nrow = nrow(A)))
  
  for(k in 1:(n-1)){
    
    for(i in (k+1):n){
      
      if(A[i,k] != 0){
        lambda = A[i,k]/A[k,k]
        A[i,] = A[i,] - lambda*A[k,]
        
        res$L[i,k] <- lambda
      } else {
        next
      }
      
    }
  }
  
  diag(res$L) <- 1
  res$U <- A
  return(res)
}

result = LU_decom(A= BB)
rm(A, A2, AA, B, B2, BB, L, L2, U, U2, lambda, BBB,CCC, Sol, x, result)
####################################################################
sym = function(x){
  
  # 1.squred matrix 점검
  if(ncol(x) == nrow(x)){
    
    # 2. 실제로 대칭인지 점검
    # i = 행, k = 열
    for(i in 1:(nrow(x)-1)){
      
      for(k in (i+1):(ncol(x))){
        
        if(i == k){
          next
        }
        
        if(x[i,k] != x[k,i]){
          print("it is not Symmetric")
          return()
        }
        
      }
      
    }
    print("it is symmetric matrix")
    
  } else {
    print("this matrix cannot do Cholesky decomposition because it is not a square matrix")
  }
}

## 3. Choleski decomposition(Pivot 연산이 필요없어서 수치적 안정)
## example 2.6

A = matrix(c(4, -2, 2, -2, 2, -4, 2, -4, 11), byrow = T, nrow = 3)
B = matrix(c(5, 2, 1, 4,
             2, 1, 0 ,1,
             1, 0, 1, 2,
             4, 1, 2, 5), byrow = T, nrow = 4)
t(chol(A)) 

chol(B)

## 1) line by line
x = matrix(0, nrow = 3, ncol = 3)

x[1,1] <- sqrt(A[1,1]); x[2,1] <- A[2,1]/x[1,1]; x[3,1] <- A[3,1]/x[1,1]; 

x[2,2] <- sqrt(A[2,2] - x[2,1]^2); x[3,2] <- (A[3,2]-x[2,1]*x[3,1])/x[2,2];

x[3,3] <- sqrt(A[3,3]- x[3,2]^2 - x[3,1]^2)

## 2) generalization
Chol_decom = function(x){
  
  ## 저장공간
  n = ncol(x)
  res = list(L = matrix(0,ncol = ncol(x), nrow = nrow(x)))
  
  ## decomposition 행:i 열:j
  for(j in 1:n){ # 열
    
    for(i in j:n){ # 행 #j = 2일때는 2,1 은 해선 안됌
      
      if(i == j){ 
        squared_sum = pracma::dot()
        res$L[i,j] <- sqrt(x[i,j] - squared_sum)
        
      } else {
        res$L[i,j] <- 
      }
      
    }
  }
  
  ##
  return(res)
  
}
## 4. swap and partial pivoting
## Swap
## 1) line by line
exam_piv = matrix(c(-1, 3, 2,
                    4, 2, 4,
                    2, 0.5, 1), byrow = T, ncol = 3)

exam_piv[c(2,3,1),]

## 2)generalize
swap = function(X){
  
  idx = order(X[,1], decreasing = T)
  X = X[idx,]
  
  return(X)
}

## partial pivoting
## 1)


## 5. Gauss-seidel method(indirect or iterative method)

## 6. Final result
## -> use 'solve' from LAPACK