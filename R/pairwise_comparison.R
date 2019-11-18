pairwise_comparison <- function( x, out='delta' ){
  if ( length(x) >2 ){
    cat(red$bold$underline('Warning! There is more than 2 values, will only use first 2!\n'))
  }
  if (out=='delta'){
    return( abs( x[1] - x[2]) )
  } else if (out=='strcat') {
    return( paste(x[1],x[2],sep='_') )
  } else if (out=='isequal'){
    return( equals(x[1],x[2]))
  } else {
    cat(red(underline(out), ' is not a valid option. One of: {"delta","strcat","isequal"}'))
  }
}