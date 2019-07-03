# Overview of files and analyses in the **`vector-occurrence-sdm`** folder

Scripts in this folder:

- `cross-val-gam-tgb.R`: Runs cross-validation as described
in the manuscript to obtain GAMs that perform best w.r.t. to
AUC on hold out test data (fold 5). The cross-validation is
a spatial block-wise CV (see manuscript). The data used
for each species is obtained by the target-group background (TGB)
approach as described in the manuscript.

- `eval-cv-gam-tgb.R`: Following `cross-val-gam-tgb.R` this script
evaluates the results and produces the predicted maps (based on
latest available covariate layers) based on the
best performing model for each species respectively. These are
saved as PDF as well as `.tif` raster files in **`prediction-maps`**

- `analyses-publication-sdm.R`: Produces tables and figures that are
used in the manuscript/publication (manuscript is developed at Overleaf)
