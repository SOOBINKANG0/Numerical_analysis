## incremental Search method
## example 4.1
rm(list=ls())
increment = function(x, y, dx){
  
  start = x
  end = y
  
  grid = seq(from = x, to = y, by = dx)
  
  for(i in 1:length(grid)){
    
    f = function(x){
      return(x^3-10*x^2+5)
    }

    res1 = f(grid[i])
    res2 = f(grid[i+1])
    
    if(sign(res1)*sign(res2) < 0){
        break
    } else {
      
    }
    
  }
  return(c(grid[i], res2))
  rm(res1, res2, grid, f, increment)
}

increment(x = 0, y = 1, dx = 0.0001) # [1]  0.734600000 -0.001261529

## bisection
## example 4.2

bisection = function(a, b, tol){
  
  f = function(x){
    return(x^3-10*x^2+5)
  }
  
  
  ##
  while(1){
    
    mid = (a+b)/2
    
    if(sign(f(a))*sign(f(mid)) > 0){
      a <- mid
    } else {
      b <- mid
    }
    
    if(abs(b-a) < tol){
      break
    }
  }
  return((a+b)/2)  
}

bisection(a = 0, b =1 ,tol = 0.0001) # [1] 0.7345886

## newton-rahpson method
## example 4.4
library(pracma)

f = function(x){
  return((x^3-10*x^2+5))
}

f_prime = function(x){
  return(3*x^2-20*x)
}

##

newton = function(x, f, iter, tol){
  
  bin = numeric(iter)
  bin[1] <- x
  iter_count = 0
  
  for(i in 1:iter){
    
    bin[i+1] = bin[i] - f(bin[i])/f_prime(bin[i])
    
    if(abs(bin[i+1]-bin[i]) < tol){
      break
    }
  
    iter_count = iter_count +1
      
  }
  return(print(paste0("x is", bin[i+1], " and iteration is", iter_count)))
}
given_f = expression(x^3-10*x^2+5)

newton(x = 0.6, f = given_f, iter = 100, tol = 0.0001)#
 
## secant method
plot(x = seq(from = -3, to = 3, by = 0.01),
     y = f(seq(from = -3, to = 3, by = 0.01)), type = "l", 
     main = "대략적인 해값", xlab = "-3~+3", ylab = "f(x)"); abline(h = 0)

## example
## secant(0.6, 0.8, maxiter = 100, fun = f)

my_secant = function(a, b, f = f, max_iter, tol){
  
  x = numeric(max_iter+2)
  x[1] <- a
  x[2] <- b
  iter_count = 0
  
  ## main loop
  for(i in 1:max_iter){
    x[i+2] = x[i+1] - f(x[i+1])/
      (
        (f(x[i+1])-f(x[i])) / (x[i+1] - x[i])
      )
    
    if(abs((x[i+1]- x[i])/x[i]) < tol){
       break
    }
    
    iter_count = iter_count +1
  }
  
  return(cat(paste0("X is ",x[i+2],"iteration is", iter_count)))
}

my_secant(a= 0.6, b= 0.8, f= f, max_iter = 100, tol = 0.0001)

## Ridder's method
## example 4.5
## 어렵다... GPT의 도움 받음음

my_func <- function(x){
  return(x^3 - 10*x^2 + 5)
}

my_ridders <- function(f, a, b, max_iter, tol){
  
  x1 <- a; x2 <- b
  f1 <- f(x1); f2 <- f(x2)
  
  # 이전 근사값을 저장할 변수
  x_root <- NULL 
  
  for(i in 1:max_iter){
    
    x3 <- (x1 + x2) / 2
    f3 <- f(x3)
    
    # 분모가 0이 되는 예외 처리
    if (f3^2 - f1*f2 == 0) return(x3)
    s <- sqrt(f3^2 - f1*f2)
    
    x4 <- x3 + (x3 - x1) * sign(f1 - f2) * f3 / s
    
    if (!is.null(x_root)) {
      if (abs(x4 - x_root) < tol) {
        return(x4)
      }
    }
    
    x_root <- x4
    f4 <- f(x_root)
    
    if (f3 * f4 < 0) {
      x1 <- x3
      f1 <- f3
      x2 <- x4
      f2 <- f4
    } else if (f1 * f4 < 0) {
      x2 <- x4
      f2 <- f4
    } else {
      x1 <- x4
      f1 <- f4
    }
  }

  return(x_root)
}


my_ridders(a = 0.6, b = 0.8, f = my_func, max_iter = 100, tol = 0.001)


######################
## newton 수렴속도  ##
######################
alpha <- 0.734603507789303

newton <- function(x, f, f_prime, iter, tol, alpha) {

  bin <- numeric(iter + 1)
  bin[1] <- x
  
  iter_count <- 0
  for (i in 1:iter) {
    iter_count <- i
    bin[i+1] <- bin[i] - f(bin[i]) / f_prime(bin[i])
    
    if (abs(bin[i+1] - bin[i]) < tol) {
      break
    }
  }
  
  results <- data.frame(iteration = 0:iter_count)
  results$x_n <- bin[1:(iter_count + 1)]
  
  results$error <- abs(results$x_n - alpha)
  
  p <- numeric(iter_count + 1)
  p[] <- NA 
  
  for (i in 3:(iter_count + 1)) {
    if (results$error[i-1] < 0.0001 || results$error[i-2] < 0.0001) next
    p[i] <- log(results$error[i] / results$error[i-1]) / log(results$error[i-1] / results$error[i-2])
  }
  
  results$p <- p
  return(results)
}
res = newton(x = 0.6, f = f, f_prime = f_prime, alpha = alpha, iter = 100, tol = 0.0001)
res = as.data.frame(res)
write.csv(res, "res.csv")

#####################
## secant 수렴속도 ##
#####################
alpha = 0.734603507789303
my_secant <- function(a, b, f, max_iter, tol, alpha) {
  
  x <- numeric(max_iter + 2)
  x[1] <- a
  x[2] <- b
  
  iter_count <- 0
  
  for (i in 1:max_iter) {
    iter_count <- i
    
    fx_i1 <- f(x[i+1])
    fx_i <- f(x[i])
    
    x[i+2] <- x[i+1] - fx_i1 * (x[i+1] - x[i]) / (fx_i1 - fx_i)
    
    if (abs((x[i+1] - x[i]) / (x[i])) < tol) {
      break
    }
    
  }
  
  results <- data.frame(iteration = 0:(iter_count + 1))
  results$x_n <- x[1:(iter_count + 2)]
  
  results$error <- abs(results$x_n - alpha)
  p <- numeric(iter_count + 2)
  p[] <- NA
  
  for (i in 3:(iter_count + 2)) {
    if (results$error[i-1] < 0.0001 || results$error[i-2] < 0.0001) next
    p[i] <- log(results$error[i] / results$error[i-1]) / log(results$error[i-1] / results$error[i-2])
  }
  results$p <- p
  
  return(results)
}


res2 = my_secant(a = 0.6, b = 0.8, f = f, alpha = alpha, max_iter = 100, tol = 0.0001)
res2 = as.data.frame(res2)
write.csv(res2, "res2.csv")

##########################
#### System equations ####
##########################
install.packages("nleqslv")
library(nleqslv)

## 1) packages
sys_eq = function(vars){
  
  x = vars[1]
  y = vars[2]
  
  res = numeric(length(vars))
  
  res[1] = x^2+y^2-3
  res[2] = x*y-1
  
  return(res)
}

sol = nleqslv(x = c(3,0.5), fn = sys_eq)
print(sol$x) #[1] 1.618034 0.618034

## 2) 구현
origin_f <- function(vars) {
  x <- vars[1]
  y <- vars[2]
  
  results <- numeric(2)
  results[1] <- x^2+y^2-3 # f1
  results[2] <- x*y-1            # f2
  
  return(results)
}


Jacobi_f = function(vars){
  
  x = vars[1]
  y = vars[2]
  
  J = matrix(NA, ncol = 2, nrow = 2)
  
  J[1,1] <-2*x   ; J[1,2] <- 2*y
  J[2,1] <- y ; J[2,2] <- x
    
  return(J)
}

my_newton = function(vars, max_iter, tol){
  
  x0 <- vars
  
  for(i in 1:max_iter){
    
    origin_val = origin_f(x0)
    Jacobi_val = Jacobi_f(x0)
    
    delta_x <- solve(Jacobi_val, -origin_val)
    x_next <- x0 + delta_x
    
    if(sqrt(sum(delta_x^2)) < tol){
      break
    }
 
    x0 <- x_next
  }
  return(x_next)
}

my_newton(vars = c(3,0.5), max_iter = 100, tol = 0.0001)

#########################
## zero of polynomials ##
#########################

# Evaluation of polynomial

## hard coding
eval1 = function(x){
  return(x^3-2*x^2-8*x +27) # n(n-1)/2
}

eval1(x=3) # 12

## generalized
eval2 = function(a, x){ # a is coef, x is evaluation point
  
  p = 0
  
  for(i in 1:4){
    p = p + a[i]*x^(i-1)
  }
  
  return(p)
}

eval2(a = c(27,-8,-2, 1),x = 3)

## 