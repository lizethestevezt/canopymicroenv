# =============================================================================
# csv_processing.R
# =============================================================================
# Pipeline for processing field epiphyte observations into a standardised
# format suitable for canopy-zone analysis.
#
# WORKFLOW OVERVIEW:
#   1. Raw GeoJSON files (from field mapping apps) are converted to CSV via
#      a Python script (geojson-to-csv.py).
#   2. Those CSVs are reformatted into EpiphytesDatabase format via a second
#      Python script (convert_observations.py).
#   3. All processed CSVs are combined into a single dataframe (combined.csv).
#   4. process_csv() loops over each observation row, resolves elevation
#      (from field or AWS), fetches canopy height (from field or GEE), and
#      produces one record per individual observation.
#   5. get_sites() extracts unique site metadata for spatial queries.
#
# COMBINED CSV COLUMN SCHEMA:
#   Source, Area_or_Site, lat, lon, Elevation_m,
#   FieldID, Abundance, Height_m, CanopyHeight_m,
#   pictures, note, AI_ID, FinalID, Genus, species
#
# DEPENDENCIES:
#   R packages : readr, tidyr, dplyr, elevatr, sf, parzer, rgee, terra, RWmisc
#   Python scripts:
#     - scripts/geojson-csv-sql-conversion-tools/python/geojson-to-csv.py
#     - scripts/config_processing/convert_observations.py
#   Internal R scripts:
#     - scripts/get_data/canopy_height.R   (provides get_canopy_height())
#     - scripts/config_processing/paths.R  (project path constants)
#
# USAGE:
#   Source this file, edit the `sites` list in Section 4, then run the
#   combined and downstream blocks.
# =============================================================================

reticulate::use_condaenv("base", required = TRUE)

library(readr)
library(tidyr)
library(dplyr)
library(elevatr)
library(sf)
library(parzer)
library(rgee)
library(terra)

source("scripts/get_data/canopy_height.R")
source("scripts/config_processing/paths.R")

ee$Authenticate(auth_mode = "notebook")
ee$Initialize(project = "ee-lizethestevezt")





# =============================================================================
# SECTION 2: GeoJSON -> EpiphytesDatabase conversion (calls Python scripts)
# =============================================================================



# =============================================================================
# SECTION 3: Core processing functions
# =============================================================================









# =============================================================================
# SECTION 4: Run pipeline
