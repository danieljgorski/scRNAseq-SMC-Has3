# make directories
output_dirs <- c(
  "results",
  "results/objects",
  "results/dimplots",
  "results/cluster-markers",
  "results/genes-of-interest",
  "results/composition",
  "results/differential-gene-expression",
  "results/volcano-plots"
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
source("scripts/5-SMC-Has3-genes-of-interest.R")
source("scripts/6-SMC-Has3-composition.R")
source("scripts/7-SMC-Has3-differential-gene-expression.R")
source("scripts/8-SMC-Has3-volcano-plots.R")
