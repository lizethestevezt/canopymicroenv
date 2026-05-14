# patches
original_wf_request <- ecmwfr::wf_request
assignInNamespace(
  "wf_request",
  function(request, user = "ecmwfr", transfer = TRUE, path = tempdir(),
           time_out = 3600, retry = 120, job_name, verbose = TRUE) {
    original_wf_request(
      request  = request,
      user     = user,
      transfer = transfer,
      path     = path,
      time_out = time_out,
      retry    = retry,
      verbose  = verbose
    )
  },
  ns = "ecmwfr"
)

assignInNamespace(
  "reflectance_calc",
  function(alb, lai, x, plotprogress = TRUE, maxiter = 50, tol = 0.001, bwgt = 0.5) {
    e1  <- terra::intersect(terra::ext(lai), terra::ext(alb))
    e   <- terra::intersect(e1, terra::ext(x))
    lai <- terra::crop(lai, e)
    alb <- terra::crop(alb, e)
    x   <- terra::crop(x, e)
    all_same <- terra::compareGeom(lai, alb, x)
    if (all_same) {
      tst  <- exp(-mean(as.vector(lai), na.rm = TRUE))
      lref <- (x * 0 + 0.5) * (1 - bwgt) + bwgt * alb  # fix: wgt -> bwgt
      gref <- x * 0 + 0.15
      mxdif <- tol * 10
      paim  <- as.matrix(lai, wide = TRUE)
      xm    <- as.matrix(x,   wide = TRUE)
      albm  <- as.matrix(alb, wide = TRUE)
      itr   <- 1
      while (mxdif > tol) {
        if (tst < 0.5) {
          lref2 <- microclimdata:::.rast(microclimdata:::find_lref(paim, as.matrix(gref, wide = TRUE), xm, albm), x)
          lref2 <- microclimdata:::.fillna(lref2, x, zerotoNA = FALSE)
          gref2 <- microclimdata:::.rast(microclimdata:::find_gref(as.matrix(lref2, wide = TRUE), paim, xm, albm), x)
          gref2 <- microclimdata:::.fillna(gref2, x, zerotoNA = FALSE)
        } else {
          gref2 <- microclimdata:::.rast(microclimdata:::find_gref(as.matrix(lref, wide = TRUE), paim, xm, albm), x)
          gref2 <- microclimdata:::.fillna(gref2, x, zerotoNA = FALSE)
          lref2 <- microclimdata:::.rast(microclimdata:::find_lref(paim, as.matrix(gref, wide = TRUE), xm, albm), x)
          lref2 <- microclimdata:::.fillna(lref2, x, zerotoNA = FALSE)
        }
        gref  <- bwgt * gref + (1 - bwgt) * gref2
        lref  <- bwgt * lref + (1 - bwgt) * lref2
        mxdif1 <- mean(abs(as.vector(gref) - as.vector(gref2)), na.rm = TRUE)
        mxdif2 <- mean(abs(as.vector(lref) - as.vector(lref2)), na.rm = TRUE)
        mxdif  <- max(mxdif1, mxdif2)
        itr <- itr + 1
        if (itr > maxiter) mxdif <- 0
      }
    } else {
      stop("Geometries of input rasters do not match")
    }
    return(list(gref = gref, lref = lref))
  },
  ns = "microclimdata"
)

library(rgee)
library(readr)
library(mcera5)
library(microclimf)
library(microclimdata)
library(terra)
# ── specify python environment to use ── 
reticulate::use_python("/Users/lizethestevezt/miniforge3/bin/python", required = TRUE)

# ── logging ──────────────────────────────────────────────────────────────────
log_file <- file.path("/Users/lizethestevezt/canopymicroenv/logs",
                      sprintf("run_micropoint_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
log_msg <- function(msg) {
  stamped <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  message(stamped)
  cat(stamped, "\n", file = log_file, append = TRUE)
}
log_msg("Script started")

# ── source functions, import credentials and site information ── 
source("scripts/functions.R")
mycredentials <- readRDS("/Users/lizethestevezt/canopymicroenv/credentials.rds")
site <- make_site("data/csv/combinedv3.csv")

# 0 make the SpatialRaster, TemplateRaster, TimePeriodObjects and import directories
source("scripts/paths.R") # here I need a small function that checks if all directories
                          # exist, so creates them if necessary


raster <- terra::rast(nrows = 2, ncols = 2,
                             xmin  = site$lon_min, xmax = site$lon_max,
                             ymin  = site$lat_min, ymax = site$lat_max,
                             crs   = "EPSG:4326")
raster_utm <- terra::project(raster, "EPSG:32717")
r_dtm <- terra::rast(nrows = 50, ncols = 50, xmin  = site$lon_min, xmax = site$lon_max,
                     ymin  = site$lat_min, ymax = site$lat_max, crs = "EPSG:4326")
r_buffered <- terra::rast(
  nrows = 2, ncols = 2,
  xmin  = site$lon_min - 0.5,
  xmax  = site$lon_max + 0.5,
  ymin  = site$lat_min - 0.5,
  ymax  = site$lat_max + 0.5,
  crs   = "EPSG:4326"
)

TemplateRaster <- terra::rast(
  extent     = terra::ext(site$lon_min, site$lon_max,
                          site$lat_min, site$lat_max),
  resolution = 0.0001,
  crs        = "EPSG:4326")
tme_start <- as.POSIXlt(site$tme_start, tz = "UTC")
tme_end <- as.POSIXlt(site$tme_end, tz = "UTC")

tme <- as.POSIXlt(seq(from = tme_start, to = tme_end, by = "day"))
tmeMonth <- as.POSIXlt(seq(
  from = as.POSIXlt(format(tme_start, "%Y-%m-01 00:00:00"), tz = "UTC"),
  to   = as.POSIXlt(format(seq(as.POSIXlt(format(tme_start, "%Y-%m-01"), tz="UTC"), 
                               by = "month", length.out = 2)[2] - 3600, 
                           "%Y-%m-%d %H:00:00"), tz = "UTC"), by   = "day"))
tmeHourly <- as.POSIXlt(seq(
  from = as.POSIXlt(format(tme_start, "%Y-%m-01 00:00:00"), tz = "UTC"),
  to   = as.POSIXlt(format(seq(as.POSIXlt(format(tme_start, "%Y-%m-01"), tz="UTC"),
                               by = "month", length.out = 2)[2] - 3600,
                           "%Y-%m-%d %H:00:00"), tz = "UTC"), by   = "hour"))
  
# 1) get input information
  # weather data
    weatherdata  <- get_weather(site = site, credentials = mycredentials, 
                            r = r_buffered, tme = tmeHourly, dir = ERA5_DIR) 
  # dtm data, not masking cause we're not in coastal areas
    dtmdata <- get_dtm(r = r_dtm, dir = DTM_DIR, mask = FALSE)
    dtmdata <- terra::project(dtmdata, "EPSG:4326")
  # landcover data, needed for both vegetation and soil data, earth engine needed 
    ee_Authenticate()
    ee_Initialize(project = "ee-lizethestevezt", user = "lizethestevezt@gmail.com", drive = TRUE)
    landcoverdata <- get_landcover(site = site, r = raster,  out_dir = LCOVER_DIR)
  # vegetation data used in the model 
    vegetationdata <- get_vegetation(r = raster, lcover = landcoverdata, tme = tmeHourly, 
                                 lat = mean(c(site$lat_min, site$lat_max)),
                                 lon = mean(c(site$lon_min, site$lon_max)))
    for (n in names(vegetationdata)) {
      if (class(vegetationdata[[n]])[1] == "PackedSpatRaster")
        vegetationdata[[n]] <- terra::unwrap(vegetationdata[[n]])
    }
  # soil because it's important (:
    soildata <- get_soil(r = raster, template = TemplateRaster, 
                     tme = tmeMonth, credentials = mycredentials, 
                     landcover = landcoverdata, dir = SOIL_DIR, albedodir = ALB_DIR)
    for (n in names(soildata)) {
      if (class(soildata[[n]])[1] == "PackedSpatRaster")
        soildata[[n]] <- terra::unwrap(soildata[[n]])
    }

# 2) Use inputs in point model within a loop 
    tme_model <- as.POSIXlt(tmeHourly, tz = "UTC")
    heights   <- seq(from = site$hObs_min, to = site$hObs_max, by = 0.2)
    models    <- list()
    
    log_msg(sprintf("Running point model: %.1f – %.1f m (%d steps)",
                    min(heights), max(heights), length(heights)))
    
    for (h in heights) {
      message("")  # blank line before each step
      log_msg(sprintf("  runpointmodela @ %.1f m", h))
      models[[sprintf("h%.1f", h)]] <- microclimf::runpointmodela(
        climarrayr = weatherdata,
        tme        = tme_model,
        reqhgt     = h,
        dtm        = dtmdata,
        vegp       = vegetationdata,
        soilc      = soildata
      )
      message("")  # blank line after progress bar clears
    }
    
    log_msg(sprintf("Saving %d models to pointmodel1.rds", length(models)))
    saveRDS(models, file.path(BASE_DIR, "data/processed/pointmodelv1.rds"))
    log_msg("Done")
