inner = function(f,g,grid){
  h = f*g
  return(sum((0.5*h[1:(length(grid)-1)] + 0.5*h[2:(length(grid))])*(grid[2] - grid[1])))
}
