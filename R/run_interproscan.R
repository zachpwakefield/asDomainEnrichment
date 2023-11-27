run_interproscan <- function(bash_location, interproscan_location, cluster = F, project, memory, name_job, email, output_location) {
  if (cluster) {
    temp_in <- readr::read_lines(paste0(bash_location, "/run_cluster_interproscan.sh"))
    temp_in[3] <- paste(temp_in[3], project)
    temp_in[4] <- paste(temp_in[4], memory)
    temp_in[6] <- paste(temp_in[6], name_job)
    temp_in[8] <- paste(temp_in[8], email)
    temp_in[10] <- paste0(temp_in[10], output_location)
    temp_in[11] <- paste0(temp_in[11], interproscan_location)
    write_lines(temp_in, paste0(bash_location, "/run_cluster_interproscan_withVars.sh"))
    system(paste("qsub ", bash_location, "/run_cluster_interproscan_withVars.sh", sep=""))
  } else {
    system(paste("bash ", bash_location, "/run_interproscan.sh",
          " ", output_location,
          " ", interproscan_location,
          sep=""))
  }

}

