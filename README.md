# Species distribution models of chagas vectors

This folder contains code used to fit species ditribution models for chagas vectors.

**Disclaimer**: The analysis is not reproducible as it relies on
a specific location of covariate layers and some confidential data that
could not be made public.


## Folder structure
Overview of project files and folders:

- Preliminary notes:
    + The initial steps of this project requiere access to the **`mastergrids`** folder at Maleria Atlas Project and will thus not be fully reproducible

- **`mastergrids`**:
An **`R`** package that facilitates the import of environmental
raster data from the **`mastergrids`** folder at BDI MAP and some
utility functions to transform rasters to data frames and vice versa
(the data import won't work outside of the BDI/without connection to mastergrids). Also contains two functions `grid_to_df` and `df_to_grid` which
convert RasterLayer/RasterBrick object to a data frame and vice versa.

- **`tcruziutils`**: An **`R`** package that facilitates all steps of the
analysis. Most functions are specific to this project and should not be used
for general purpose projects. It can be loaded at the beginning of the script
using `devtools::load_all("<path to project>/tcruziutils")` or installed via
`devtools::install("<path to project>/tcruziutils")` and then loaded as usual
with `library(tcruziutils)`

- **`infection data`** (not included):
Contains data on infection prevalence and presence in vectors and humans
    + **`External vector database`**: Additional data set (compiled by another
    research group that contains *presence only* data on different
    vector species)

- **`endemic zone`** (not included):
Contains a shape-file that defines the endemic zone of the
disease ("mask"). Can be used to crop environmental grid data and other spatial
objects. The spatial extent defined by this mask is stored in the
**`tcruziutils`** package as `raster::extent` and `sp::bbox` objects for
convenience (see `tcruziutils::extent_tcruzi` and `tcruziutils::bbox_tcruzi`).

- **`polygon boundaries`** (not included):
Shapefiles containing polygon boundaries on administrative district levels 1 and 2. Used to extract location/area and polygon information based on the GAUL code. Cropped according to endemic zone and stored as `shp_admin_1.Rds` and `shp_admin_2.Rds` in the **`preprocessing`** folder

- **`preprocessing`**:
This folder contains the main pre-processing scripts and stores the pre-processed data sets that will be used for modeling

    + `import.R`: Initial data import of the presence/prevalence. For
    observations recorded on a polygon level, also adds the centroid
    and area information of the respective polygon to the data (see function
    `add_spat_dat` in **`tcruziutils`** package)

    + `prep-for-modeling.Rmd`: Builds on the `import.R` script and pre-processes the initially imported data for modeling purposes. This includes:
      - application of inclusion/exclusion criteria
      - imputation of (some) missing data
      - addition of covariate layers to the observed presence/absence data
      (based on coordinates; see `raster::extract` and
      `tcruziutils::add_grid_data`)
      - Splits the data set in *train*/*test* data (and possibly *evaluation*
      data). Additionally creates block-wise cross-validation scheme,
      stored as `fold` variable in the original data.
      - also produces some additional visualizations (stored in
      `tcruzi/figures` as `pdf` and `png`)

- **`descripts`**: Some preliminary pre-processing (after initial import) that
informs steps in `prep-for-modeling.Rmd`
    + `descriptive-analysis.Rmd`: Initial descriptive analysis after import


- **`modeling`**:

     - **`vector_occurrence-sdm`**:
     Folder containing Species Distribution Models (SDM) based on presense only data:
     - See README therein
