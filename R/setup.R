get_c_trans <- function(location) {
  return(read_lines(paste(location,"/f_protein_code_from_gencodev43.txt",sep="")))
}

get_gtf <- function(location) {
  return(read.csv(paste(location,"/gtf_withInfo.csv",sep=""), header = T)[,-c(1)])
}
