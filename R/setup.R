c_trans <- function() {
  return(read_lines("../f_protein_code_from_gencodev43.txt"))
}

gtf <- function() {
  return(read.csv("../gtf_withInfo.csv", header = T)[,-c(1)])
}
