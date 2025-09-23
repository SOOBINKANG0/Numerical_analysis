## 3. Interpolation and curve fitting

install.packages("polynom")
library(polynom)
library(tidyverse)
library(pracma)

## Polynomial function
P = polynomial(c(6,-1,-7,1,1))
f = as.function(P)
approx #linearly interp

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
##


## Example 3.4
x = c(0.15, 2.3, 3.15, 4.85, 6.25, 7.95)
y = c(4.79867, 4.49013, 4.2243, 3.47313, 2.66674, 1.51909)



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