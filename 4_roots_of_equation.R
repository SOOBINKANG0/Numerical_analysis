## 4.2 incremental Search method
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

## 4.3 bisection
## example 4.2



