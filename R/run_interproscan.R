run_interproscan <- function(bash_location, interproscan_location, cluster = F, project, memory, name_job, email, output_location) {
  if (cluster) {
    system(paste("qsub ", bash_location, "/run_cluster_interproscan.sh",
          " ", project,
          " ", memory,
          " ", name_job,
          " ", email,
          " ", output_location,
          " ", interproscan_location, sep=""))
  } else {
    system(paste("bash ", bash_location, "/run_interproscan.sh",
          " ", output_location,
          " ", interproscan_location,
          sep=""))
  }

}

