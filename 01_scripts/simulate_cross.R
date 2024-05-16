# simulate_cross function sourced from 
# Konstantin Divilov, 2023 https://github.com/kdivilov/Aquaculture_2023/blob/main/GWAS.R

simulate_cross = function(x){
  if(x[1]==0 && x[2]==0){
    return(0)
  }
  if(x[1]==2 && x[2]==2){
    return(2)
  }
  if(x[1]==0 && x[2]==2){
    return(1)
  }
  if(x[1]==2 && x[2]==0){
    return(1)
  }
  if(x[1]==1 && x[2]==1){
    return(1)
  }
  if(x[1]==0 && x[2]==1){
    return(0.5)
  }
  if(x[1]==1 && x[2]==0){
    return(0.5)
  }
  if(x[1]==2 && x[2]==1){
    return(1.5)
  }
  if(x[1]==1 && x[2]==2){
    return(1.5)
  }
}

