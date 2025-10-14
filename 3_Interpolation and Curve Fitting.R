## 3. Interpolation and curve fitting

install.packages("polynom")
library(polynom)
library(tidyverse)
library(pracma)

## 1. 라그랑지 방법
lagrange_method = function(xdata, ydata, x){ 
  
  n = length(xdata) #데이터 포인트 개수
  L_values = numeric(n) #다항식 결과과
      
    for(i in 1:n){
      
      temp = 1
      current_tem = numeric(n)
      
      for(j in 1:n){
        if(i != j){
          current_tem[j] = (x - xdata[j])/ (xdata[i] - xdata[j])
        } else {
          current_tem[j] <- 1
        }
        
        temp = temp * current_tem[j]
      }

      L_values[i] = temp
    }
  
  result = as.numeric(ydata %*% L_values)
  return(result)
}

## example 3.1 
x = c(0,2,3)
y = c(7,11,28)
lagrangeInterp(x= x, y = y, 1) # 4
lagrange_method(xdata = x, ydata = y, x = 1) #4

## 2. Newton's Method
## example 3.2
x = c(-2, 1, 4 , -1, 3,  -4)
y = c(-1, 2, 59,  4, 24, -53)

## 1)coef mat
my_coef = function(x,y){
  
  n  = length(y)
  bin = matrix(NA, nrow = n, ncol = n)
  bin[,1] <- y
    
  for(j in 2:n){ #열
      
    for(i in j:n){ # 행
    
      bin[i,j] <- (bin[i,j-1]-bin[j-1,j-1])/(x[i]-x[j-1])
      
      # inductino
      # bin[2,2] <- (bin[2,1]-bin[1,1])/(x[2]-x[1]) # 1
      # bin[3,2] <- (bin[3,1]-bin[1,1])/(x[3]-x[1]) # 10
      # bin[4,2] <- (bin[4,1]-bin[1,1])/(x[4]-x[1]) # 5
      # bin[5,2] <- (bin[5,1]-bin[1,1])/(x[5]-x[1]) # 5
      # bin[6,2] <- (bin[6,1]-bin[1,1])/(x[6]-x[1]) # 26
      
      # bin[3,3] <- (bin[3,2]-bin[2,2])/(x[3]-x[2]) #3
      # bin[4,3] <- (bin[4,2]-bin[2,2])/(x[4]-x[2]) #-2
      # bin[5,3] <- (bin[5,2]-bin[2,2])/(x[5]-x[2]) #2
      # bin[6,3] <- (bin[6,2]-bin[2,2])/(x[6]-x[2]) #-5

    }
  }
  return(bin)
}

res_coef = my_coef(x=x,y=y)

## 2) evaluation
my_eval = function(x, coef, x_new){
  
  a = as.numeric(diag(coef))
  n = length(x)-1 # 다항식 차수
  p = a[n+1]
  
  for(i in n:1){
  p <- a[i] + (x_new - x[i])*p
  }
  
  return(p)
}

my_eval(x = x, coef = res_coef, x_new = 1)

plot(x=x,y=y)

x_grid = seq(from = min(x), to = max(x), by = 0.01)
y_grid = my_eval(x = x, coef = res_coef, x_new = x_grid)
lines(x_grid, y_grid, lwd = 3)

## 3. Neville's Method
## 자주 쓰지 않는 알고리즘이므로 GPT로 만듦

## exercise 14.
data = data.frame(x = c(-2, -0.1, -1.5, 0.5, -0.6, 2.2, 1.0, 1.8),
                  y = c(2.2796, 1.0025, 1.6467, 1.0635, 1.0920, 2.6291, 1.2661, 1.9896))

## GPT
neville <- function(xData, yData, x) {

  m <- length(xData)
  y <- yData
  
  for (k in 1:(m - 1)) {
    
    # Neville의 재귀 공식을 벡터 연산으로 계산
    # P_i, P_{i+1}, ..., P_{i+k} 를 계산하는 과정
    indices <- 1:(m - k)
    
    y[indices] <- ((x - xData[indices + k]) * y[indices] + (xData[indices] - x) * y[indices + 1]) /
                  (xData[indices] - xData[indices + k])
  }
  
  # 최종 결과는 항상 첫 번째 요소에 저장됨
  return(y[1])
}

neville(xData = data$x, yData = data$y, 1.1)
neville(xData = data$x, yData = data$y, 1.2)
neville(xData = data$x, yData = data$y, 1.3)

## 4. polynomial fit
## 1)example 3.11, find a*exp(bx) = ln(a) + bx

x = c(1.2, 2.8, 4.3, 5.4, 6.8, 7.9)
y = c(7.5,16.1,38.9,67.0,146.6,266.2)
ln_y = log(y, base= exp(1))

res = round(lm(ln_y~x)$coefficients,4) 
# 1.3321+0.5366x
# ln(a) = 1.3321, exp(1.3321) = 3.788992
par(mfrow = c(1,2))
plot(3.788992*exp(0.5366*x)~x, type = "l"); plot(1.3321+0.5336*x~x, type = 'l')

## 2)polynomial fitting
## 생략

## 5. cubic spline
## 이걸 잘해야할듯