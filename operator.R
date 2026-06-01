operator=function(basis_temp,eigen,f,grid){
  temp = 0 
  for (k in 1:length(eigen)){
    temp = temp + eigen[k]*inner(basis_temp[,k],f,grid)*basis_temp[,k]
  }
  return(temp)
}