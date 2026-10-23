# =====================================================================
# 00_prepare_data.R  -  run ONCE before launching the app
# Reads the Excel workbook, cleans it, downloads province boundaries,
# and saves small .rds files that app.R loads quickly.
# =====================================================================
# install.packages(c("readxl","dplyr","tidyr","sf","geodata"))
library(readxl); library(dplyr); library(tidyr); library(sf)

cat("Working directory:", getwd(), "\n")
xl <- "Data (1).xlsx"          # <- put the workbook in the project folder; edit name if needed
if (!file.exists(xl)) {
  cat("Excel files found here:", list.files(pattern = "\\.xlsx?$"), "\n")
  stop("Cannot find '", xl, "' in ", getwd(), ". Edit the xl <- line to match your file name exactly.")
}
dir.create("data", showWarnings = FALSE)
cat("STEP 1: reading time series...\n")

# ---- helper: make province names comparable across sources ----------
key <- function(x) gsub("[^a-z]", "", tolower(x))   # "North-West" -> "northwest"

# ---- 1. Time series (cumulative counts) -----------------------------
read_series <- function(sheet, value) {
  read_excel(xl, sheet = sheet) %>%
    select(-Total) %>%
    filter(!is.na(Gauteng)) %>%            # drops the all-NA reporting gaps (27 Mar, 7 Apr in Cases)
    mutate(Date = as.Date(Date)) %>%
    pivot_longer(-Date, names_to = "province", values_to = value)
}
cases <- read_series("Confirmed Cases", "cases")
deaths <- read_series("Deaths", "deaths")
recov <- read_series("Recoveries", "recoveries")

# Start where all three series exist (recoveries begin 17 Apr 2020)
dates <- seq(as.Date("2020-04-17"), max(cases$Date), by = "day")

ts <- cases %>%
  full_join(deaths, by = c("Date", "province")) %>%
  full_join(recov,  by = c("Date", "province")) %>%
  complete(province, Date = dates) %>%         # make one row per province per day
  arrange(province, Date) %>%
  group_by(province) %>%
  fill(cases, deaths, recoveries, .direction = "down") %>%   # cumulative => carry last value forward
  ungroup() %>%
  filter(Date %in% dates)

stopifnot(!anyNA(ts))
saveRDS(ts, "data/ts.rds")

# ---- 2. Population (ASSUMPTION - verify!) ---------------------------
# Stats SA mid-year population estimates 2020 (P0302). These are NOT in your
# workbook. Check each number against the Stats SA release before citing.
pop <- tribble(
  ~province,        ~pop,
  "Eastern Cape",    6734001,
  "Free State",      2928903,
  "Gauteng",        15488137,
  "KwaZulu-Natal",  11531628,
  "Limpopo",         5852553,
  "Mpumalanga",      4679786,
  "Northern Cape",   1292786,
  "North West",      4108816,
  "Western Cape",    7005741
)
write.csv(pop, "data/population_2020.csv", row.names = FALSE)
saveRDS(pop, "data/pop.rds")

cat("STEP 3: downloading province boundaries (needs internet)...\n")
# ---- 3. Province boundaries (GADM via geodata) ----------------------
# Needs internet once. Alternative: read your own shapefile with st_read().
v <- geodata::gadm(country = "ZAF", level = 1, path = "data")
prov <- st_as_sf(v) %>% select(NAME_1) %>% st_make_valid()
prov$province <- unique(ts$province)[match(key(prov$NAME_1), key(unique(ts$province)))]
stopifnot(nrow(prov) == 9, !anyNA(prov$province))   # all 9 names must match
prov <- prov %>% select(province) %>% st_transform(4326)
saveRDS(prov, "data/prov.rds")

cat("STEP 4: reading health centres...\n")
# ---- 4. Community Health Centres (points) ---------------------------
chc <- read_excel(xl, sheet = "Community Health Centers")
cat("CHCs with missing coordinates (dropped):\n"); print(chc[!complete.cases(chc), ])
chc <- chc %>% filter(!is.na(Longitude), !is.na(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
saveRDS(chc, "data/chc.rds")

# ---- 5. QA: does each CHC's recorded province match its location? ---
chk <- st_join(chc, prov, suffix = c("_listed", "_spatial"))
cat("\nCHCs whose listed province differs from the polygon they fall in:",
    sum(chk$Province != chk$province, na.rm = TRUE), "\n")
cat("CHCs falling outside all polygons:", sum(is.na(chk$province)), "\n")
cat("\nDone. Files written to:", file.path(getwd(), "data"), "\n")
print(list.files("data"))
