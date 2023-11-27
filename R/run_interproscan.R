run_interproscan <- function(bash_location, interproscan_location, cluster = T, project, memory, name_job, email, output_location) {
  if (bu_scc) {
    paste("qsub ", bash_location, "run_cluster_interproscan.sh",
          " ", project,
          " ", memory,
          " ", name_job,
          " ", email,
          " ", output_location,
          " ", interproscan_location, sep="")
  } else {
    paste("bash ", bash_location, "run_interproscan.sh",
          " ", output_location,
          " ", interproscan_location,
          sep="")
  }

}

