# canopymicroenv

> **Microclimate niche modelling for epiphytic orchids along an elevational gradient in NW Ecuador**
> Master's thesis project — University of Bonn, 2026

---

## What this is

This repository contains the full analysis pipeline for my master's thesis on the **vertical stratification of epiphytic Maxillariinae orchids** and their microenvironmental correlates across five sites in the Chocó Andino cloud forest of northwestern Ecuador (Maquipucuna, Mindo Mirador, Mashpi, MindoTarabita, Yanayacu).

The core idea: instead of using coarse climate data to describe orchid habitat, this pipeline models the **exact microclimate at the height and location where each individual was observed** — down to 0.2 m resolution in the canopy vertical gradient. These microclimate profiles are then used to characterise the realised microclimatic niche of each species and explore how vertical space partitioning emerges in a diverse epiphyte community.

> ⚠️ **This pipeline is a work in progress.** It's running, but not all steps are fully implemented yet. See [Status](#status) below.

---

## Study system

- **Focal group:** Maxillariinae orchids (tribe Maxillarieae, subtribe Maxillariinae)
- **Sites:** 5 cloud forest sites, ~800–1800 m elevation, NW Ecuador
- **Observations:** 141 individuals across all sites, with observed height in canopy (`hObs`), canopy height, elevation, genus, and photo references
- **Approach:** Field observations → microclimate modelling → niche characterisation → vertical stratification analysis

---

## Pipeline overview

```
GeoJSON field exports
        ↓
  convert_observations.py        # parse field notes → structured CSV
        ↓
  csv_processing.R               # clean, validate, extract hObs
        ↓
  make_site()                    # bounding box + time window from observations
        ↓
  main.R  ── ERA5 climate data   # hourly gridded weather
          ── DTM                 # digital elevation model
          ── ESA landcover       # via Google Earth Engine
          ── Vegetation params   # from habitat classification
          ── Soil params         # SoilGrids + MODIS LAI + albedo
        ↓
  runpointmodela()               # microclimf point model, per height step
        ↓
  [pending] runmicro()           # grid model
        ↓
  [pending] per-observation extraction   # join microclimate → hObs
        ↓
  [pending] statistical analysis        # vertical stratification ~ microclimate
```

---

## Repository structure

```
canopymicroenv/
├── scripts/
│   ├── functions.R                     # all pipeline functions
│   ├── helper_functions.R              # small utilities (logging, parsing, etc.)
│   ├── main.R                          # executable: data acquisition + model run
│   ├── paths.R                         # directory constants
│   ├── convert_observations.py         # GeoJSON → CSV (Python)
│   ├── build_photo_lookup.py           # photo reference table builder
│   └── config_processing/
│       └── csv_processing.R            # observation cleaning
├── geojson_to_csv/
│   ├── raw/                            # raw GeoJSON exports from field app
│   └── csv/                            # per-site processed CSVs
├── data/
│   ├── csv/                            # combined observation dataset
│   ├── raw/                            # ERA5, DTM, LAI, etc. (not tracked)
│   └── processed/                      # model outputs (not tracked)
├── logs/                               # run logs (not tracked)
└── output/                             # figures and results (not tracked)
```

---

## Status

| Step | Status |
|---|---|
| Field data collection (5 sites) | ✅ Complete |
| GeoJSON → CSV conversion | ✅ Complete |
| Observation cleaning + `hObs` extraction | ✅ Complete |
| ERA5 climate data download | ✅ Complete |
| DTM, landcover, vegetation, soil parameters | ✅ Complete |
| Point model height loop (`runpointmodela`) | ✅ Running |
| Species identification | 🔄 In progress |
| Grid model (`runmicro`) | ⏳ Not yet implemented |
| Per-observation microclimate extraction | ⏳ Not yet implemented |
| Statistical analysis | ⏳ Not yet started |

---

## Key dependencies

- [`microclimf`](https://github.com/ilyamaclean/microclimf) — mechanistic microclimate modelling (Maclean 2026)
- [`microclimdata`](https://github.com/ilyamaclean/microclimdata) — automated input data acquisition
- [`mcera5`](https://github.com/dklinges9/mcera5) — ERA5 climate data download
- [`rgee`](https://github.com/r-spatial/rgee) — Google Earth Engine interface from R
- `terra`, `dplyr`, `readr`

---

## Notes

- Johansson canopy zones (JZ1–JZ5) were considered and dropped in favour of directly observed individual heights (`hObs`). All analysis uses actual measurement data.
- Two monkey-patches are applied at runtime in `main.R` to fix known bugs in `ecmwfr` and `microclimdata` without modifying package source.
- `credentials.rds` (CDS API, NASA Earthdata, Google credentials) is excluded from version control. You'll need your own.

---

*Lizeth Estevez Tobar — Universidad de Bonn, 2026*
*Supervisor: Juliano Sarmento Cabral*