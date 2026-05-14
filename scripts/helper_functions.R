# =============================================================================
# Helper functions
# =============================================================================
# Column spec shared by read_csv() calls throughout — keeps types consistent
# when binding rows across sites and prevents pictures from coercing to logical.
COMBINED_COL_TYPES <- cols(
  Source         = col_character(),
  Area_or_Site   = col_character(),
  lat            = col_double(),
  lon            = col_double(),
  Elevation_m    = col_character(),   # may be a range string e.g. "850-890"
  FieldID        = col_character(),
  Abundance      = col_double(),
  Height_m       = col_double(),
  CanopyHeight_m = col_double(),
  pictures       = col_character(),
  note           = col_character(),
  AI_ID          = col_character(),
  FinalID        = col_character(),
  Genus          = col_character(),
  species        = col_character()
)
# -----------------------------------------------------------------------------
# normalize_dms()
# Standardises degree-minute-second (DMS) coordinate strings by replacing
# variant Unicode degree/minute/second symbols with ASCII equivalents.
# Needed because legacy CSVs may use DMS coordinates.
# Field observation CSVs already have decimal degrees and bypass this.
#
# Args:
#   x : Character — raw DMS string (e.g. "05° 54′ 17.9″ S")
#
# Returns: character — cleaned DMS string
# -----------------------------------------------------------------------------
normalize_dms <- function(x) {
  x <- trimws(x)
  x <- gsub("\u00b0", "\u00b0", x)           # degree sign
  x <- gsub("\u2032|\u0027", "\u0027", x)    # prime / apostrophe
  x <- gsub("\u2033|\u201c|\u201d", "\"", x) # double prime / curly quotes
  x <- gsub("\\s+([NSEW])$", " \\1", x)
  x
}

# --------------------------
# log_msg()
# --------------------------
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...))
  message(msg)
  cat(msg, "\n", file = log_file, append = TRUE)
}

# -----------------------------------------------------------------------------
# parse_one()
# Converts a flexible date/year string into a POSIXlt timestamp.
# Handles three input formats:
#   "YYYY"       -> start: Jan 1 00:00:00  / end: Dec 31 23:59:59
#   "YYYY-MM"    -> start: 1st of month    / end: last day of month 23:59:59
#   "YYYY-MM-DD" -> start: 00:00:00        / end: 23:59:59
#
# Args:
#   x        : Character or numeric - date value from the CSV
#   is_start : Logical - TRUE for expedition start, FALSE for end
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

# --------------------------
# bbox()
# --------------------------

bbox <- function (r) # to get the bounding box of a SpatRaster
{
  e <- ext(r)
  xy <- data.frame(x = c(e$xmin, e$xmin, e$xmax, e$xmax), 
                   y = c(e$ymin, e$ymax, e$ymin, e$ymax))
  xy <- sf::st_as_sf(xy, coords = c("x", "y"), crs = crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(sf::st_coordinates(ll))
  out <- c(min(ll$X), min(ll$Y), max(ll$X), max(ll$Y))
  return(out)
}

# --------------------------
# .assert_lonlat()
# --------------------------

.assert_lonlat <- function(lon, lat) {
  if (!is.numeric(lon) || !is.numeric(lat)) {
    stop("lon and lat must be numeric")
  }
  if (abs(lat) > 90 && abs(lon) <= 90) {
    stop("Coordinates look swapped: lat > 90 but lon <= 90")
  }
  if (abs(lon) > 180 || abs(lat) > 90) {
    stop("Invalid lon/lat values")
  }
}

# --------------------------
# get_canopy_height()
# extract height (using Earth Engine)
# --------------------------

get_canopy_height <- function(lon, lat) {
  
  .assert_lonlat(lon, lat)
  
  canopy_img <- ee$Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1")$select("b1")
  
  pt <- ee$Geometry$Point(list(lon, lat))  # explicit list, not c()
  
  val <- canopy_img$sample(region = pt, scale = 10, geometries = FALSE)$first()
  
  if (is.null(val)) return(NA_real_)
  
  h <- val$get("b1")$getInfo()
  
  if (length(h) == 0 || is.null(h) || is.na(h)) {
    return(NA_real_)
  }
  
  as.numeric(h)
}




########### FROM MEMO

# functions
library(fields)
library(viridis)

pink_viridis <- colorRampPalette(c("#3777FF", "#FFB5C2","#FFE9CE", "#FFE156", "#FFBE86"))

ricker <- function(N, r, K){
  y<-N*exp(r*(1-(N/K)))
  return(y)
}

stochRicker <- function(N,r,K){
  y <- rpois(1, lambda = ricker(N,r,K))
  return(y)
}



popDemo <- function(timesteps = 30, xDim = 20 , yDim = 50, repRate = 1.8, carCap = 200, Ninit = 1, stochastic = FALSE, Visualize = FALSE, sleeptime = 0.2){
  #initialisation:
  # landscape:
  landscape <- array(data = 0, dim = c(yDim, xDim, timesteps))
  for (i in 1:timesteps){
    for (r in 1:yDim){
      landscape[r, ,i] <- runif(n=xDim, min = 0, max = 0.01*r)
    }
  }
  #abundance:
  abundance <- array (data = Ninit, dim= c(yDim, xDim, timesteps))
  totalabundance <- vector(mode = "numeric", length = timesteps)
  for (r in 1:yDim){
    for (c in 1:xDim){
      if (landscape[r, c, 1]>0) abundance[r, c, 1] <- runif(n=1, min=1, max=1+carCap*landscape[r,c,1])
    }
  }
  # eco loop:
  for (index in 1:(timesteps-1)){
    for (r in 1:yDim){
      for (c in 1:xDim){
        if (stochastic==FALSE) abundance[r, c, index+1] <- ricker(N=abundance[r, c, index], 
                                                                  r = repRate, K = carCap*landscape[r, c, index]) # pop dynamics
        if (stochastic==TRUE) abundance[r, c, index+1] <- stochRicker(N=abundance[r, c, index], 
                                                                      r = repRate, K = carCap*landscape[r,c, index]) # pop dynamics
      }
    }
    # Output
    totalabundance[index+1] <- sum(abundance[, , index+1])
    Output <- list(landscape, abundance, totalabundance)
    if (Visualize==TRUE){
      par(mfrow=c(1,3))
      image.plot(landscape[, , index], col = pink_viridis(25))
      mtext(paste("Landscape (t = ", index+1, ")", sep = ""))
      image.plot(abundance[, , index], col = pink_viridis(25), axis.args = list(c(0: carCap)))
      mtext(paste("Abundance (t = ", index+1, ")", sep = ""))
      plot(totalabundance[1:index], type = "b", col= "pink", ylab = "Time step", main = "Total abundance")
      abline (h = carCap*xDim*yDim, col = "red")
      Sys.sleep(sleeptime)
    }
  }
  return(Output)
}

indDisp <- function(maxDispersal, r, c, indN, Disp){
  for (ind in 1:indN){
    yDir <- sample(-1:1, size = 1, replace = TRUE)
    xDir <- sample(-1:1, size = 1, replace = TRUE)
    rDist <- sample(0:maxDispersal, size = 1, replace = TRUE)
    Disp[r+rDist*yDir+maxDispersal, c+rDist*xDir+maxDispersal] <- 
      Disp[r+rDist*yDir+maxDispersal, c+rDist*xDir+maxDispersal] + 1
  }
  return(Disp)
}

metapop <- function(timesteps = 30, xDim = 10 , yDim = 100, repRate = 1.05, carCap = 200, Ninit = 1, stochastic = FALSE, Visualize = TRUE, sleeptime = 0.2, maxDisp = 4){
  #initialisation:
  # landscape:
  landscape <- array(data = 0, dim = c(yDim, xDim, timesteps))
  for (i in 1:timesteps){
    for (r in 1:yDim){
      landscape[r, ,i] <- runif(n=xDim, min = 0, max = 0.01*r)
    }
  }
  #abundance:
  abundance <- array (data = Ninit, dim= c(yDim, xDim, timesteps))
  totalabundance <- vector(mode = "numeric", length = timesteps)
  #dispersal 
  dispersalmatrix <- array(data = 0, dim = c((yDim+2*maxDisp), (xDim+2*maxDisp), timesteps))
  
  # modeling the landscape: 
  for (r in 1:yDim){
    for (c in 1:xDim){
      if (landscape[r, c, 1]>0) abundance[r, c, 1] <- runif(n=1, min=1, max=1+carCap*landscape[r,c,1])
    }
  }
  # eco loop:
  for (index in 1:(timesteps-1)){
    for (r in 1:yDim){
      for (c in 1:xDim){
        if (stochastic==FALSE){
          offspring <- ricker(N=abundance[r, c, index], 
                              r = repRate, K = carCap*landscape[r, c, index]) # pop dynamics
        }
        if (stochastic==TRUE){
          offspring <- stochRicker(N=abundance[r, c, index], 
                                   r = repRate, K = carCap*landscape[r,c, index]) # pop dynamics
        }
        if (offspring>0){
          dispersalmatrix[, , index] <- indDisp(maxDispersal = maxDisp,r = r, c = c, 
                                                indN = offspring, Disp = dispersalmatrix[, , index])
        }
      }
    }
    for (r in 1:yDim){
      for (c in 1:xDim){
        abundance[r, c, index+1] <- dispersalmatrix[r + maxDisp, c + maxDisp, index]  
      }
    }
    # Output
    totalabundance[index+1] <- sum(abundance[, , index+1])
    Output <- list(landscape, abundance, totalabundance)
    
    #Visualization
    if (Visualize==TRUE){
      par(mfrow=c(1,3))
      image.plot(landscape[, , index], col = pink_viridis(25))
      mtext(paste("Landscape (t = ", index+1, ")", sep = ""))
      image.plot(abundance[, , index], col = pink_viridis(25), axis.args = list(c(0: carCap)))
      mtext(paste("Abundance (t = ", index+1, ")", sep = ""))
      plot(totalabundance[1:index], type = "b", col= "pink", ylab = "Time step", main = "Total abundance", las = 1)
      abline (h = carCap*xDim*yDim, col = "red")
      Sys.sleep(sleeptime)
    }
  }
  return(Output)
}


my.OU <- function(n , theta = 0.5, alpha = 1, sigma2 = 0.1){
  sd = sqrt(sigma2/2*alpha)
  rnorm(n, mean = theta, sd)
}

simulate_network <- function(resourceMin = 2, resourceMax = 20, # range of resource species
                             consumerMin = 2, consumerMax = 20, # range of consumer species
                             networksN = 1, # number of networks
                             resTraitMin = 0, resTraitMax = 1, # range of mean trait values for resources
                             consTraitMin = 0, consTraitMax = 1, # range of mean trait values for consumers
                             resourceTolerance = 0.1, consumerTolerance = 0.1, #range of tolerance
                             alpha = 0.01){ # pull for OU
  # define a list for result storage
  simulations <- vector("list", networksN)
  for (k in 1:networksN){
    # define spp richness in each simulation -> network size
    resourceN <- sample(resourceMin:resourceMax, 1)
    consumerN <- sample(consumerMin:consumerMax, 1)
    
    # plugin an OU trait distribution
    resourceTraits <- my.OU(n = resourceN, theta = 10, alpha, sigma2 = 0.1)
    consumerTraits <- my.OU(n = consumerN, theta = 10, alpha, sigma2 = 0.1)
    
    # tolerance in trait spread
    delta_resource <- runif(resourceN, 0.05, resourceTolerance)
    delta_consumer <- runif(consumerN, 0.05, consumerTolerance)
    
    # store the results
    web <- matrix(0, nrow = resourceN, ncol = consumerN)
    
    # trait matching rule
    for (i in 1:resourceN){ # loop for resource spp
      for (j in 1:consumerN){ # loop for consumer spp 
        trait_distance <- abs(resourceTraits[i] - consumerTraits[j]) # first part of the equation
        tolerance <- 0.5 * (delta_resource[i] + delta_consumer[j]) # tolerance threshold
        if (trait_distance < tolerance) web[i, j] <- 1
      }
    }
    #store values
    linksN <- sum(web)
    simulations[[k]] <- list(
      resourceN = resourceN,
      consumerN = consumerN,
      resourceTraits = resourceTraits,
      consumerTraits = consumerTraits,
      delta_resource = delta_resource,
      delta_consumer = delta_consumer,
      web = web,
      linksN = linksN 
    )
  }
  return(simulations)
} 

merge_era5_steptype_files <- function(pathin, pathout) {
  library(ncdf4)
  
  nc_files <- list.files(pathin, pattern = "stepType.*\\.nc$", full.names = TRUE)
  if (length(nc_files) == 0) stop("No stepType .nc files found in ", pathin)
  
  # Variable rename map (new CDS names -> what microclimdata expects)
  rename_map <- c(
    avg_snlwrf = "msnlwrf",
    avg_sdlwrf = "msdwlwrf"
  )
  
  # Open all source files
  datasets <- lapply(nc_files, nc_open)
  
  # Use first file as dimension template
  src <- datasets[[1]]
  
  # Collect all dimensions from source
  out_dims <- lapply(src$dim, function(d) {
    ncdim_def(d$name, d$units, d$vals, unlim = d$unlim)
  })
  names(out_dims) <- names(src$dim)
  
  # Collect all variables from all files, skipping duplicates and dim vars
  seen_vars  <- names(src$dim)  # skip dimension variables
  out_vars   <- list()
  var_data   <- list()
  
  for (ds in datasets) {
    for (vname in names(ds$var)) {
      if (vname %in% seen_vars) next
      seen_vars <- c(seen_vars, vname)
      
      var      <- ds$var[[vname]]
      out_name <- ifelse(vname %in% names(rename_map), rename_map[vname], vname)
      
      # Get dimension objects for this variable
      var_dims <- lapply(var$dim, function(d) out_dims[[d$name]])
      
      out_vars[[out_name]] <- ncvar_def(
        name  = out_name,
        units = var$units,
        dim   = var_dims,
        missval = var$missval
      )
      var_data[[out_name]] <- ncvar_get(ds, vname)
      
      message("  ", vname, " -> ", out_name)
    }
  }
  
  # Create output file and write
  nc_out <- nc_create(pathout, vars = out_vars)
  for (vname in names(var_data)) {
    ncvar_put(nc_out, vname, var_data[[vname]])
  }
  nc_close(nc_out)
  lapply(datasets, nc_close)
  
  message("Merged file written to ", pathout)
  return(pathout)
}


create_vegpoint_fixed <- function(landcover, vhgt, lai, refldata, 
                                  lctype = "ESA", lat = NA, long = NA) {
  if (class(lat) == "logical") {
    e  <- terra::ext(landcover)
    xy <- data.frame(x = (e$xmin + e$xmax) / 2, y = (e$ymin + e$ymax) / 2)
  } else {
    xy <- data.frame(x = long, y = lat)
    xy <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 4326)
    xy <- sf::st_transform(xy, terra::crs(landcover))
    xy <- data.frame(x = sf::st_coordinates(xy)[1], y = sf::st_coordinates(xy)[2])
  }
  
  reso <- terra::res(landcover)[1]
  xmin <- floor((xy$x - reso) / reso) * reso
  xmax <- xmin + 2 * reso
  ymin <- floor((xy$y - reso) / reso) * reso
  ymax <- ymin + 2 * reso
  e    <- terra::ext(xmin, xmax, ymin, ymax)
  
  lcc     <- terra::crop(landcover, e)
  vhgtc   <- terra::crop(vhgt,      e)
  laic    <- terra::crop(lai,        e)
  refldata$lref <- terra::crop(refldata$lref, e)
  
  # Fixed: create_veggrid instead of create_vegpgrid
  vegp  <- microclimdata::create_veggrid(lcc, vhgtc, laic, refldata, lctype)
  
  x     <- terra::rast(vegp$x)
  gsmax <- terra::rast(vegp$gsmax)
  clump <- terra::rast(vegp$clump)
  leafd <- terra::rast(vegp$leafd)
  leaft <- terra::rast(vegp$leaft)
  
  vegpp <- list(
    h     = terra::extract(vhgtc,        xy)[, 2],
    pai   = terra::extract(laic,         xy)[, 2],
    x     = terra::extract(x,            xy)[, 2],
    clump = terra::extract(clump,        xy)[, 2],
    lref  = terra::extract(refldata$lref, xy)[, 2],
    ltra  = terra::extract(leaft,        xy)[, 2],
    leafd = terra::extract(leafd,        xy)[, 2],
    gsmax = terra::extract(gsmax,        xy)[, 2],
    q50   = 100
  )
  
  class(vegpp) <- "vegparams"
  return(vegpp)
}
