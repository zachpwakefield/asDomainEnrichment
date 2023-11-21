c_trans <- function() {
  return(read_lines(paste(system.file(package="domainEnrichment"),"/f_protein_code_from_gencodev43.txt",sep="")))
}

gtf <- function() {
  return(read.csv(paste(system.file(package="domainEnrichment"),"/gtf_withInfo.csv",sep=""), header = T)[,-c(1)])
}
