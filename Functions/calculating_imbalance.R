imbalance_metric <- function (tree) {
  x <-balance(tree)
  for (i in (1:nrow(x))) {
    if (x[i,1] < x[i,2]) {
      x[i,1:2] <- x[i,2:1]
    }
      }
 return(mean(x[,1]/(x[,1]+ x[,2])))}

