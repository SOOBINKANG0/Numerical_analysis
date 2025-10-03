## incremental Search method
## example 4.1

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
## example 4.7
library(pracma)

f = function(x){
  return((x^4 -6.4*x^3 + 6.45*x^2+20.538*x-31.752))
}

f_prime = function(x){
  return(4*x^3 - 6.4*3*x^2 + 6.45*2*x + 20.538)
}

##

newton = function(x, f, iter, tol){
  
  bin = numeric(iter)
  bin[1] <- x
  
  for(i in 1:iter){
    
    bin[i+1] = bin[i] - f(bin[i])/f_prime(bin[i])
    
    if(abs(bin[i+1]-bin[i]) < tol){
      break
    }
    
  }
  return(bin[i+1])
}

given_f = expression(x^4 -6.4*x^3 + 6.45*x^2+20.538*x-31.752)
newton(x = 0, f = given_f, iter = 1000, tol = 0.0001)# [1] 2.099926
 
## secant method
## example