# make directories
output_dirs <- c(
  "results",
  "results/objects",
  "results/dimplots",
  "results/cluster_markers")

for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# run analysis in order
source("scripts/SMC-Has3-clustering.R")
source("scripts/SMC-Has3-annotation.R")
source("scripts/SMC-Has3-dimplot.R")
source("scripts/SMC-Has3-cluster-markers.R")
