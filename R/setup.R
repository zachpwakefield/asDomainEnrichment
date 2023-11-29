# get_c_trans
# input: directory containing f_protein_code_from_gencodev43.txt, default should be package directory or: system.file(package="domainEnrichment")
# output: f_protein_code_from_gencodev43.txt lines

get_c_trans <- function(location) {
  return(read_lines(paste(location,"/f_protein_code_from_gencodev43.txt",sep="")))
}

# get_gtf
# input: directory containing gtf_withInfo, default should be package directory or: system.file(package="domainEnrichment")
# output: gtf_withInfo.csv csv

get_gtf <- function(location) {
  return(read.csv(paste(location,"/gtf_withInfo.csv",sep=""), header = T)[,-c(1)])
}
