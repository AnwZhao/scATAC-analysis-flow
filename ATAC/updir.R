updir -> function(path){
  a <- unlist(strsplit(path, split="/"))
  b <- paste(a[-(length(a))], collapse="/")
  return(b)
}