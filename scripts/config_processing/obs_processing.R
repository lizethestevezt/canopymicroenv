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
#   3. All processed CSVs are combined into a single dataframe.
#   4. process_csv() loops over each observation row, resolves coordinates and
#      elevation, fetches canopy height, and expands Johansson zone abundance
#      into individual records.
#   5. get_sites() extracts unique site metadata for spatial queries.
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
#   Source this file and edit the `sites` list at the bottom, then run the
#   `combined` and downstream blocks.
# =============================================================================

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


# =============================================================================
# SECTION 1: Helper functions
# =============================================================================

# -----------------------------------------------------------------------------
# get_crown_depth()
# Returns the vertical crown depth (m) for a given location.
# Currently a placeholder returning a fixed value; replace with real logic
# (e.g. a lookup table or remote-sensing model) when available.
#
# Args:
#   lat  : Latitude in decimal degrees
#   lon  : Longitude in decimal degrees
#
# Returns: numeric — crown depth in metres
# -----------------------------------------------------------------------------
get_crown_depth <- function(lat, lon) {
  return(2.3)  # placeholder – replace with real logic
}


# -----------------------------------------------------------------------------
# get_hMed()
# Calculates the median height (m) for a given Johansson zone, based on
# canopy structure measurements.
#
# Johansson zones (JZ1–JZ5) divide the vertical profile of a tree from
# ground to canopy top:
#   JZ1 : ground/root zone           (~0.5 m)
#   JZ2 : lower trunk                (ground to crown base)
#   JZ3 : lower crown third
#   JZ4 : middle crown third
#   JZ5 : upper crown third
#
# Args:
#   zone         : Character, one of "JZ1"–"JZ5"
#   height_at_pt : Canopy height at the observation point (m)
#   crowndepth   : Vertical depth of the crown (m)
#   crown_base   : Height of the crown base above ground (m)
#
# Returns: numeric — median height of the zone in metres
# -----------------------------------------------------------------------------
get_hMed <- function(zone, height_at_pt, crowndepth, crown_base) {
  
  if (zone == "JZ1") return(0.5)
  
  if (zone == "JZ2") {
    return(mean(c(1, crown_base)))
  }
  
  if (zone == "JZ3") {
    return(mean(c(crown_base,
                  crown_base + crowndepth / 3)))
  }
  
  if (zone == "JZ4") {
    return(mean(c(crown_base + crowndepth / 3,
                  crown_base + 2 * crowndepth / 3)))
  }
  
  if (zone == "JZ5") {
    return(mean(c(crown_base + 2 * crowndepth / 3,
                  height_at_pt)))
  }
}


# -----------------------------------------------------------------------------
# normalize_dms()
# Standardises degree-minute-second (DMS) coordinate strings by replacing
# variant Unicode degree/minute/second symbols with ASCII equivalents.
# Needed because field data often contains mixed encodings.
#
# Args:
#   x : Character — raw DMS string (e.g. "05° 54′ 17.9″ S")
#
# Returns: character — cleaned DMS string
# -----------------------------------------------------------------------------
normalize_dms <- function(x) {
  x <- trimws(x)
  x <- gsub("°", "°", x)
  x <- gsub("′|'", "'", x)
  x <- gsub("″|"|"", "\"", x)
  x <- gsub("\\s+([NSEW])$", " \\1", x)
  x
}


# -----------------------------------------------------------------------------
# parse_one()
# Converts a flexible date/year string into a POSIXlt timestamp.
# Handles three input formats:
#   "YYYY"       → start: Jan 1 00:00:00 / end: Dec 31 23:59:59
#   "YYYY-MM"    → start: 1st of month   / end: last day of month 23:59:59
#   "YYYY-MM-DD" → start: 00:00:00       / end: 23:59:59
#
# Args:
#   x        : Character or numeric — date value from the CSV
#   is_start : Logical — TRUE for expedition start, FALSE for end
#   tz       : Timezone string (default "UTC")
#
# Returns: POSIXlt timestamp
# -----------------------------------------------------------------------------
parse_one <- function(x, is_start = TRUE, tz = "UTC") {
  
  x <- as.character(x)
  
  if (grepl("^\\d{4}$", x)) {
    year <- as.integer(x)
    if (is_start) {
      as.POSIXlt(sprintf("%04d-01-01 00:00:00", year), tz = tz)
    } else {
      as.POSIXlt(sprintf("%04d-12-31 23:59:59", year), tz = tz)
    }
    
  } else if (grepl("^\\d{4}-\\d{2}$", x)) {
    year  <- as.integer(substr(x, 1, 4))
    month <- as.integer(substr(x, 6, 7))
    if (is_start) {
      as.POSIXlt(sprintf("%04d-%02d-01 00:00:00", year, month), tz = tz)
    } else {
      last_day <- last_day_of_month(year, month)
      as.POSIXlt(
        sprintf("%04d-%02d-%02d 23:59:59", year, month, last_day),
        tz = tz
      )
    }
    
  } else {
    as.POSIXlt(
      paste0(x, if (is_start) " 00:00:00" else " 23:59:59"),
      tz = tz
    )
  }
}


# =============================================================================
# SECTION 2: GeoJSON → EpiphytesDatabase conversion (calls Python scripts)
# =============================================================================

# -----------------------------------------------------------------------------
# run_conversion()
# Runs the two-step Python conversion pipeline for a single field site:
#   Step 1 — geojson-to-csv.py    : GeoJSON → raw CSV
#   Step 2 — convert_observations.py : raw CSV → EpiphytesDatabase format
#
# Args:
#   site_name            : Character — site label used in filenames and the
#                          --site flag passed to convert_observations.py
#   geojson_path         : Character — path to the input .geojson file
#   geojson_to_csv_script: Character — path to geojson-to-csv.py
#   convert_obs_script   : Character — path to convert_observations.py
#   csv_dir              : Character — directory for intermediate raw CSVs
#   output_dir           : Character — directory for final processed CSVs
#
# Returns: character — path to the processed output CSV
# -----------------------------------------------------------------------------
run_conversion <- function(site_name,
                           geojson_path,
                           geojson_to_csv_script = "scripts/geojson-csv-sql-conversion-tools/python/geojson-to-csv.py",
                           convert_obs_script    = "scripts/config_processing/convert_observations.py",
                           csv_dir               = "geojson_to_csv/csv",
                           output_dir            = "data/csv") {
  
  raw_csv <- file.path(csv_dir,    paste0(site_name, ".csv"))
  out_csv <- file.path(output_dir, paste0("Processed", site_name, ".csv"))
  
  # Step 1: GeoJSON → raw CSV
  message("  [1/2] GeoJSON → CSV : ", basename(geojson_path))
  system2("python", args = c(
    geojson_to_csv_script,
    "--input",  geojson_path,
    "--output", raw_csv
  ))
  
  # Step 2: raw CSV → EpiphytesDatabase format
  message("  [2/2] Formatting    : ", basename(raw_csv))
  system2("python", args = c(
    convert_obs_script,
    "--input",  raw_csv,
    "--output", out_csv,
    "--site",   shQuote(site_name)
  ))
  
  message("  Done: ", out_csv)
  return(out_csv)
}


# =============================================================================
# SECTION 3: Core processing functions
# =============================================================================

# -----------------------------------------------------------------------------
# process_csv()
# Processes a single EpiphytesDatabase-format CSV into analysis-ready records.
# For each observation row, resolves coordinates and elevation, retrieves
# canopy height, then expands Johansson zone abundances into one record
# per zone per observation.
#
# Rows are skipped (next) if:
#   - Canopy height is NA
#   - Crown depth >= canopy height (structurally invalid)
#   - Abundance in a given zone is 0
#
# Processing stops (break) if Source is NA (treats blank rows as end of data).
#
# Args:
#   csv_path : Character — path to a processed EpiphytesDatabase CSV
#
# Returns: data.frame with columns:
#   Spp, lat, lon, elevMed, zone, hMed, abund, hCanopy, tme_start, tme_end
# -----------------------------------------------------------------------------
process_csv <- function(csv_path) {
  
  df <- read_csv(csv_path, na = c("", "NA", "N/A"))
  
  results <- list()
  k <- 1
  
  for (i in seq_len(nrow(df))) {
    
    # Stop at first empty Source (treats trailing blank rows as end of data)
    if (is.na(df$Source[i])) break
    
    # ---- Species ----
    source <- df$Source[i]
    spp    <- paste(df$Genus[i], df$species[i])
    
    # ---- Coordinates ----
    lat_dms <- normalize_dms(df$lat[i])
    lon_dms <- normalize_dms(df$lon[i])
    
    coords <- RWmisc::dms2dd(lon = lon_dms, lat = lat_dms)
    
    lat  <- coords[[2]]
    lon  <- coords[[1]]
    lat0 <- lat + 0.00001   # small offset used for raster point queries
    lon0 <- lon + 0.00001
    
    # ---- Elevation ----
    if (is.na(df$Elevation_m[i])) {
      # Fetch elevation from AWS terrain tiles when not recorded in the CSV
      pt_df  <- data.frame(x = c(lon, lon0), y = c(lat, lat0))
      pt_sf  <- st_as_sf(x = pt_df, coords = c("x", "y"), crs = 4326)
      pt_rast <- rast(pt_sf, nrow = 2, ncol = 2)
      elev   <- get_elev_point(locations = pt_rast, prj = 4326, src = "aws")$elevation[1]
    } else {
      # CSV may store elevation as a range string e.g. "850-890"
      elev <- as.numeric(strsplit(df$Elevation_m[i], "-")[[1]])
    }
    
    elevMed <- mean(elev)
    
    # ---- Canopy structure ----
    canopy_height <- get_canopy_height(lon, lat)
    if (is.na(canopy_height)) next
    
    crowndepth <- get_crown_depth(lat, lon)
    crown_base <- canopy_height - crowndepth
    
    if (crowndepth >= canopy_height) next
    
    # ---- Expedition time window ----
    tme_start <- parse_one(df$Exp_Start[i], TRUE)
    tme_end   <- parse_one(df$Exp_End[i],   FALSE)
    
    # ---- Johansson zones ----
    for (zone in paste0("JZ", 1:5)) {
      
      abund <- df[[zone]][i]
      if (abund == 0) next
      
      hMed <- get_hMed(zone, canopy_height, crowndepth, crown_base)
      
      results[[k]] <- data.frame(
        Spp       = spp,
        lat       = lat,
        lon       = lon,
        elevMed   = elevMed,
        zone      = zone,
        hMed      = hMed,
        abund     = abund,
        hCanopy   = canopy_height,
        tme_start = tme_start,
        tme_end   = tme_end
      )
      
      k <- k + 1
    }
  }
  
  observations <- do.call(rbind, results)
  return(observations)
}


# -----------------------------------------------------------------------------
# get_sites()
# Extracts unique site metadata from a processed EpiphytesDatabase CSV.
# Used to define spatial bounding boxes for downstream remote-sensing queries.
#
# Args:
#   csv_path : Character — path to a processed EpiphytesDatabase CSV
#
# Returns: data.frame with columns:
#   Site, lat_min, lon_min, lat_max, lon_max, tme_start, tme_end
# -----------------------------------------------------------------------------
get_sites <- function(csv_path) {
  
  df <- read_csv(csv_path, na = c("", "NA", "N/A")) |>
    dplyr::filter(!is.na(Source), !is.na(Area_or_Site))
  
  sites_df <- df |>
    dplyr::distinct(Area_or_Site, lat, lon, Exp_Start, Exp_End, .keep_all = TRUE)
  
  results <- vector("list", nrow(sites_df))
  
  for (i in seq_len(nrow(sites_df))) {
    
    lat_dms <- normalize_dms(sites_df$lat[i])
    lon_dms <- normalize_dms(sites_df$lon[i])
    
    coords <- RWmisc::dms2dd(lon = lon_dms, lat = lat_dms)
    
    lat_min <- coords[[2]]
    lon_min <- coords[[1]]
    
    results[[i]] <- data.frame(
      Site      = sites_df$Area_or_Site[i],
      lat_min   = lat_min,
      lon_min   = lon_min,
      lat_max   = lat_min + 0.00001,
      lon_max   = lon_min + 0.00001,
      tme_start = parse_one(sites_df$Exp_Start[i], TRUE),
      tme_end   = parse_one(sites_df$Exp_End[i],   FALSE)
    )
  }
  
  dplyr::bind_rows(results)
}


# =============================================================================
# SECTION 4: Run pipeline
# =============================================================================
# Edit `sites` to add or remove field sites. Each entry needs:
#   name   : used as the site label and in output filenames
#   geojson: path to the raw GeoJSON file exported from the mapping app
#
# Output CSVs are written to data/csv/ as "Processed{name}.csv"
# =============================================================================

sites <- list(
  list(name = "Maquipucuna",   geojson = "geojson_to_csv/raw/Maquipucuna.geojson"),
  list(name = "MiradorMindo",  geojson = "geojson_to_csv/raw/MiradorMindo.geojson"),
  list(name = "Mashpi",        geojson = "geojson_to_csv/raw/Mashpi.geojson"),
  list(name = "MindoTarabita", geojson = "geojson_to_csv/raw/MindoTarabita.geojson")
)

# --- Step 1: Convert all GeoJSON files and combine into one dataframe ---
message("=== Converting GeoJSON files ===")

combined <- lapply(sites, function(s) {
  message("Processing site: ", s$name)
  path <- run_conversion(s$name, s$geojson)
  read_csv(path, na = c("", "NA", "N/A"))
}) |>
  dplyr::bind_rows()

write_csv(combined, "data/csv/combined.csv")
message("Combined CSV saved to data/csv/combined.csv (", nrow(combined), " rows)")

# --- Step 2: Process combined CSV through the main pipeline ---
message("\n=== Processing observations ===")
observations <- process_csv("data/csv/combined.csv")

# --- Step 3: Extract site metadata ---
message("\n=== Extracting site metadata ===")
sites_metadata <- get_sites("data/csv/combined.csv")