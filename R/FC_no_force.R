#SAD generators

NoForce = function(matlist, Jmax){

  SV <- c()
  for (i in 1:length(matlist)) {
    SV <- append(SV, length(matlist[[i]]))
  }

  S <- SV

  J <- c()
  for ( i in 1:length(Jmax)){
    J <- append(J, Jmax[i] * (S[i] / (100 + S[i])))
  }


oldJ <- c()
for ( i in 1:length(matlist)){
  old <- sum(matlist[[i]])
  oldJ <- append(oldJ, old)
}

newamt <- matlist


for (i in 1:length(matlist)) {
  if (oldJ[[i]]  > J[[i]]) {
   # print("above J nulcline")
   # newamt[[i]] <- newamt[[i]] + ((1-(oldJ[[i]]/J[[i]]))*oldJ[[i]])*pers
   # newamt[[i]] <-newamt[[i]]* 1 - (J[[i]]-oldJ[[i]])/sum(unlist(matlist))
  } else {
   # print("BELOW J nulcline")
   # newamt[[i]] <-newamt[[i]] * (J[[i]]-oldJ[[i]])/sum(unlist(matlist))
   # newamt[[i]] <-newamt[[i]] + newamt[[i]]*.05

  }

}


return(newamt)

}





