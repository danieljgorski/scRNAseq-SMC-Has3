# make directories
output_dirs <- c(
  "results",
  "results/objects",
  "results/dimplots",
  "results/cluster-markers"
)

for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# run analysis in order
source("scripts/1-SMC-Has3-clustering.R")
source("scripts/2-SMC-Has3-annotation.R")
source("scripts/3-SMC-Has3-dimplots.R")
source("scripts/4-SMC-Has3-cluster-markers.R")
