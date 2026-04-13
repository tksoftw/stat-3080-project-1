# analysis.R  —  Capital Bikeshare & Crime Proximity Analysis
# R equivalent of analyze.ipynb

# ── Constants ──────────────────────────────────────────────────────────────────
EARTH_RADIUS_KM    <- 6371
KM_PER_DEG         <- 2 * pi * EARTH_RADIUS_KM / 360
DC_BLOCK_RADIUS_KM <- 62 / 1000   # ~sqrt(12 000 m² / π)
BLOCKS             <- 5
CIRCLE_RADIUS_KM   <- DC_BLOCK_RADIUS_KM * BLOCKS
CELL_DEG           <- CIRCLE_RADIUS_KM / KM_PER_DEG

# ── Load data ──────────────────────────────────────────────────────────────────
df_bike  <- read.csv("Capital_Bike_Share_Locations.csv")[ , c("LATITUDE","LONGITUDE")]
df_crime <- read.csv("Crime_Incidents_in_2026.csv"     )[ , c("LATITUDE","LONGITUDE","OFFENSE")]

cat("Bike stations:", nrow(df_bike),  "\n")
cat("Crime records:", nrow(df_crime), "\n\n")

# ── Haversine distance (km) ────────────────────────────────────────────────────
distance_earth <- function(lat1, lon1, lat2, lon2) {
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  a    <- sin(dlat/2)^2 + cos(lat1*pi/180) * cos(lat2*pi/180) * sin(dlon/2)^2
  2 * EARTH_RADIUS_KM * asin(sqrt(a))
}

# ── Build spatial grid (environment as hash map) ───────────────────────────────
grid <- new.env(hash = TRUE, parent = emptyenv())

for (k in seq_len(nrow(df_bike))) {
  lat <- df_bike$LATITUDE[k]
  lon <- df_bike$LONGITUDE[k]
  key <- paste0(floor(lat / CELL_DEG), ",", floor(lon / CELL_DEG))
  if (exists(key, envir = grid)) {
    grid[[key]] <- rbind(grid[[key]], c(lat, lon))
  } else {
    grid[[key]] <- matrix(c(lat, lon), nrow = 1)
  }
}

# ── Proximity check ────────────────────────────────────────────────────────────
point_in_any_circle <- function(lat, lon) {
  i <- floor(lat / CELL_DEG)
  j <- floor(lon / CELL_DEG)
  for (di in -1:1) {
    for (dj in -1:1) {
      key <- paste0(i + di, ",", j + dj)
      if (exists(key, envir = grid)) {
        stations <- grid[[key]]
        dists    <- distance_earth(lat, lon, stations[, 1], stations[, 2])
        if (any(dists <= CIRCLE_RADIUS_KM)) return(TRUE)
      }
    }
  }
  FALSE
}

cat("Computing proximity for", nrow(df_crime), "crimes...\n")
df_crime$bike_nearby <- mapply(point_in_any_circle,
                               df_crime$LATITUDE, df_crime$LONGITUDE)

# ── Crime classification ───────────────────────────────────────────────────────
VIOLENT  <- c("ROBBERY", "HOMICIDE", "ASSAULT W/DANGEROUS WEAPON", "SEX ABUSE")
PROPERTY <- c("MOTOR VEHICLE THEFT", "THEFT/OTHER", "THEFT F/AUTO", "BURGLARY")

df_crime$is_violent  <- df_crime$OFFENSE %in% VIOLENT
df_crime$is_property <- df_crime$OFFENSE %in% PROPERTY

# ── Part 1: Summary statistics ─────────────────────────────────────────────────
cat("\n=== OFFENSE BREAKDOWN ===\n")
offense_tbl <- sort(table(df_crime$OFFENSE), decreasing = TRUE)
print(offense_tbl)

near <- df_crime[ df_crime$bike_nearby, ]
far  <- df_crime[!df_crime$bike_nearby, ]

cat("\n=== PROXIMITY BREAKDOWN ===\n")
cat("Near bikeshare station : ", nrow(near), "\n")
cat("Far from any station   : ", nrow(far),  "\n")

prop_near <- mean(near$is_property) * 100
prop_far  <- mean(far$is_property)  * 100
viol_near <- mean(near$is_violent)  * 100
viol_far  <- mean(far$is_violent)   * 100

cat("\nProperty crime rate — Near:", round(prop_near, 1), "%   Far:", round(prop_far, 1), "%\n")
cat("Violent  crime rate — Near:", round(viol_near, 1), "%   Far:", round(viol_far, 1), "%\n")

# Contingency tables
cat("\n=== CONTINGENCY TABLE: Property Crime ===\n")
ct_prop <- table(Proximity = ifelse(df_crime$bike_nearby, "Near", "Far"),
                 Property  = ifelse(df_crime$is_property,  "Yes",  "No"))
print(ct_prop)

cat("\n=== CONTINGENCY TABLE: Violent Crime ===\n")
ct_viol <- table(Proximity = ifelse(df_crime$bike_nearby, "Near", "Far"),
                 Violent   = ifelse(df_crime$is_violent,  "Yes",  "No"))
print(ct_viol)

# ── Part 2: Chi-squared tests ──────────────────────────────────────────────────
alpha <- 0.05

cat("\n=== CHI-SQUARED TEST: Property Crime ===\n")
test_prop <- chisq.test(ct_prop)
print(test_prop)

cat("\n=== CHI-SQUARED TEST: Violent Crime ===\n")
test_viol <- chisq.test(ct_viol)
print(test_viol)

# ── Conclusion ─────────────────────────────────────────────────────────────────
cat("\n=== CONCLUSION ===\n")

if (test_prop$p.value < alpha && prop_near > prop_far) {
  cat("SET 1: Property crime IS significantly HIGHER near stations. Supports Ha.\n")
} else if (test_prop$p.value < alpha) {
  cat("SET 1: Property crime IS significantly LOWER near stations. Opposite of Ha.\n")
} else {
  cat("SET 1: No significant difference in property crime rates.\n")
}

if (test_viol$p.value < alpha && viol_near < viol_far) {
  cat("SET 2: Violent crime IS significantly LOWER near stations. Supports Ha.\n")
} else if (test_viol$p.value < alpha) {
  cat("SET 2: Violent crime IS significantly HIGHER near stations. Opposite of Ha.\n")
} else {
  cat("SET 2: No significant difference in violent crime rates.\n")
}
