# =====================================================================
# Spatial Statistics Dashboard - COVID-19 & Community Health Centres, South Africa
# Methods: (1) Global Moran's I  (2) Local Moran's I (LISA)  (3) CSR tests on CHC point pattern
# NOTE: written without being run in R by the author - test it and fix small errors.
# =====================================================================
library(dplyr); library(ggplot2)
library(sf); library(leaflet); library(spdep); library(spatstat)
# Load these LAST: spatstat defines its own card(), which would mask bslib::card()
library(shiny); library(bslib)

# ---------------------------------------------------------------- DATA
needed <- file.path("data", c("ts.rds", "pop.rds", "prov.rds", "chc.rds"))
if (!all(file.exists(needed)))
  stop("Missing: ", paste(needed[!file.exists(needed)], collapse = ", "),
       ". Working directory is ", getwd(),
       ". Run 00_prepare_data.R from the app folder first.")
ts   <- readRDS("data/ts.rds")
pop  <- readRDS("data/pop.rds")
prov <- readRDS("data/prov.rds")          # sf, EPSG:4326, column 'province'
chc  <- readRDS("data/chc.rds")           # sf points, EPSG:4326

# Project to metres so neighbours, centroids and distances are valid
crs_sa <- "+proj=aea +lat_1=-24 +lat_2=-33 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m +no_defs"
prov_p <- st_transform(prov, crs_sa)
chc_p  <- st_transform(chc, crs_sa)
cent   <- st_coordinates(st_centroid(st_geometry(prov_p)))

chc_n <- st_drop_geometry(chc) %>% count(province = Province, name = "n_chc")

# Panel of indicators for every province x date
panel <- ts %>%
  left_join(pop, by = "province") %>%
  left_join(chc_n, by = "province") %>%
  mutate(cases_raw   = cases,
         cases_100k  = cases / pop * 1e5,
         deaths_100k = deaths / pop * 1e5,
         cfr         = ifelse(cases > 0, deaths / cases * 100, NA_real_),
         rec_ratio   = ifelse(cases > 0, recoveries / cases * 100, NA_real_),
         chc_100k    = n_chc / pop * 1e5)

indicators <- c("Cumulative cases per 100k" = "cases_100k",
                "Cumulative deaths per 100k" = "deaths_100k",
                "Case fatality ratio (%)" = "cfr",
                "Recovery ratio (%)" = "rec_ratio",
                "CHCs per 100k people (static)" = "chc_100k",
                "Cumulative cases, raw count (not population-adjusted)" = "cases_raw")

# ---- point pattern objects (units = km) -----------------------------
g_km <- st_geometry(prov_p) * 0.001
st_crs(g_km) <- NA
prov_km <- st_sf(province = prov_p$province, geometry = g_km)

make_owin <- function(sfc_km) {
  u <- st_union(sfc_km) %>% st_simplify(dTolerance = 3, preserveTopology = TRUE) %>% st_cast("POLYGON")
  a <- as.numeric(st_area(u)); u <- u[a >= 0.01 * max(a)]       # drop tiny islands
  orient <- function(m, ccw) {                                    # spatstat: exterior ccw, holes cw
    m <- m[-nrow(m), , drop = FALSE]
    s <- 0.5 * sum(m[, 1] * c(m[-1, 2], m[1, 2]) - c(m[-1, 1], m[1, 1]) * m[, 2])
    if ((s > 0) != ccw) m <- m[nrow(m):1, , drop = FALSE]
    list(x = m[, 1], y = m[, 2])
  }
  rings <- list()
  for (i in seq_along(u)) for (r in seq_along(u[[i]]))
    rings[[length(rings) + 1]] <- orient(u[[i]][[r]], ccw = (r == 1))
  tryCatch(owin(poly = rings), error = function(e) NULL)
}
win <- make_owin(g_km)
win_ok <- !is.null(win)
if (!win_ok) win <- owin(range(cent[, 1] / 1000) + c(-500, 500), range(cent[, 2] / 1000) + c(-500, 500))
xy <- st_coordinates(chc_p) / 1000
X  <- suppressWarnings(ppp(xy[, 1], xy[, 2], window = win))
n_dropped <- nrow(xy) - npoints(X)

# ------------------------------------------------------------ HELPERS
nb_edges <- function(nb, cent) {
  do.call(rbind, lapply(seq_along(nb), function(i) {
    j <- nb[[i]]; j <- j[j > 0]
    if (length(j)) data.frame(x = cent[i, 1], y = cent[i, 2], xend = cent[j, 1], yend = cent[j, 2])
  }))
}
moran_text <- function(I, EI, p, alpha) {
  if (p < alpha) {
    if (I > EI) paste0("p < ", alpha, ": reject H0. I is above E[I], so similar values tend to be neighbours (positive spatial autocorrelation / clustering).")
    else paste0("p < ", alpha, ": reject H0. I is below E[I], so neighbours tend to have dissimilar values (negative spatial autocorrelation).")
  } else paste0("p >= ", alpha, ": fail to reject H0. No statistical evidence against a spatially random pattern. With only 9 provinces the test has low power, so this is not proof of randomness.")
}
cluster_cols <- c("High-High" = "#d7191c", "Low-Low" = "#2c7bb6", "High-Low" = "#fdae61",
                  "Low-High" = "#abd9e9", "Not significant" = "#d9d9d9")

latest <- panel[panel$Date == max(panel$Date), ]
css <- "
.hero{background:linear-gradient(120deg,#1f4e79,#2e75b6);color:#fff;padding:2.2rem 2rem;border-radius:12px;margin-bottom:1.2rem}
.hero h1{font-weight:700;margin-bottom:.4rem} .hero p{font-size:1.1rem;max-width:760px;opacity:.95}
.intro{background:#f2f5f9;border-left:5px solid #e07b39;padding:.7rem 1rem;border-radius:6px;margin-bottom:1rem}
.method-card .card-header{font-weight:600;color:#1f4e79}
.navbar{box-shadow:0 1px 4px rgba(0,0,0,.15)}
"
intro <- function(...) div(class = "intro", ...)

# ----------------------------------------------------------------- UI
ui <- page_navbar(
  title = "Spatial Statistics: COVID-19 & Health Centres in South Africa",
  id = "nav", selected = "home", header = tags$style(HTML(css)),
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#1f4e79", base_font = "Segoe UI, Helvetica, Arial, sans-serif"),
  sidebar = sidebar(id = "side", width = 320, open = FALSE,
    conditionalPanel("input.nav != 'pp'",
      h6("Lattice data (provinces)"),
      selectInput("ind", "Indicator", indicators),
      sliderInput("date", "Date (cumulative to)", min = min(ts$Date), max = max(ts$Date),
                  value = max(ts$Date), timeFormat = "%d %b %Y", animate = TRUE),
      hr(), h6("Spatial weights W"),
      selectInput("nbtype", "Neighbourhood", c("Queen contiguity" = "queen", "Rook contiguity" = "rook",
                                               "k nearest neighbours" = "knn")),
      conditionalPanel("input.nbtype == 'knn'", sliderInput("k", "k", 1, 5, 2, step = 1)),
      selectInput("wstyle", "Weights style", c("Row-standardised (W)" = "W", "Binary (B)" = "B")),
      selectInput("alt", "Alternative hypothesis", c("Two-sided" = "two.sided",
                  "Positive autocorrelation (greater)" = "greater", "Negative (less)" = "less")),
      sliderInput("alpha", "Significance level", 0.01, 0.10, 0.05, step = 0.01)
    ),
    conditionalPanel("input.nav == 'pp'",
      h6("Point pattern (CHC locations)"),
      radioButtons("ppview", "View", c("Intensity (kernel density)" = "kde",
                                       "Quadrat test" = "quad",
                                       "Distance functions (Monte Carlo)" = "fun")),
      conditionalPanel("input.ppview == 'kde'", sliderInput("bw", "Bandwidth sigma (km)", 20, 250, 80, step = 10),
                       checkboxInput("showpts", "Show points", TRUE)),
      conditionalPanel("input.ppview == 'quad'", sliderInput("nx", "Quadrats: columns", 2, 8, 4),
                       sliderInput("ny", "Quadrats: rows", 2, 8, 3)),
      conditionalPanel("input.ppview == 'fun'",
        selectInput("fn", "Function", c("K  (Ripley)" = "K", "L  (variance-stabilised K)" = "L",
                                        "G  (nearest-neighbour)" = "G", "F  (empty space)" = "F")),
        sliderInput("nsim", "CSR simulations", 19, 199, 39, step = 20),
        actionButton("run", "Run simulation", class = "btn-primary"))
    )
  ),
  nav_panel("Home", value = "home",
    div(class = "hero", h1("COVID-19 and Health Centres in South Africa"),
        p("An interactive spatial statistics dashboard. Explore where COVID-19 burden clusters across provinces, and whether Community Health Centres are spread randomly across the country."),
        actionButton("go_ov", "Explore the map", class = "btn-light btn-lg")),
    layout_columns(
      value_box("Confirmed cases", format(sum(latest$cases), big.mark = " "), p(format(max(latest$Date), "%d %b %Y")), theme = "primary"),
      value_box("Deaths", format(sum(latest$deaths), big.mark = " "), theme = "danger"),
      value_box("Recoveries", format(sum(latest$recoveries), big.mark = " "), theme = "success"),
      value_box("Health centres", nrow(chc), p("with coordinates"), theme = "secondary")),
    h4("Choose a spatial method", style = "margin-top:1rem"),
    layout_columns(col_widths = c(4, 4, 4),
      card(class = "method-card", card_header("1. Global Moran's I"),
           p("Are provinces with similar COVID-19 values also neighbours? One number for the whole country."),
           p(em("Data: provinces (lattice)")), actionButton("go_gm", "Open", class = "btn-primary")),
      card(class = "method-card", card_header("2. Local Moran's I (LISA)"),
           p("Which provinces form hot spots, cold spots or spatial outliers?"),
           p(em("Data: provinces (lattice)")), actionButton("go_lm", "Open", class = "btn-primary")),
      card(class = "method-card", card_header("3. Point pattern (CSR)"),
           p("Are health centre locations random, clustered or regular? Density maps, quadrat test, K/L/G/F functions."),
           p(em("Data: CHC locations (points)")), actionButton("go_pp", "Open", class = "btn-primary")))),
  nav_panel("Map explorer", value = "ov", uiOutput("vb"),
    layout_columns(col_widths = c(8, 4),
      card(card_header("Map"), checkboxInput("showchc", "Show Community Health Centres", FALSE),
           leafletOutput("map", height = 430)),
      card(card_header("Through time"), plotOutput("tsplot", height = 430)))),
  nav_panel("1. Global Moran's I", value = "gm",
    intro(strong("Question: "), "is the pattern of the chosen indicator across provinces spatially random? ", strong("H0: "), "I = E[I] (no spatial autocorrelation)."),
    layout_columns(col_widths = c(12), uiOutput("gm_vb")),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Moran scatterplot"), plotOutput("gm_scatter", height = 360)),
      card(card_header("Permutation test (Monte Carlo)"), sliderInput("mcsim", "Permutations", 99, 4999, 999, step = 100),
           plotOutput("gm_mc", height = 300))),
    navset_card_tab(
      nav_panel("Neighbour structure", plotOutput("nbplot", height = 340)),
      nav_panel("Weights matrix W", tableOutput("wtab")),
      nav_panel("Moran's I through time", plotOutput("gm_time", height = 300)))),
  nav_panel("2. Local Moran's I (LISA)", value = "lm",
    intro(strong("Question: "), "which individual provinces are significant clusters or outliers? ", strong("H0 (per province): "), "no local spatial association."),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Cluster map"), leafletOutput("lmap", height = 430)),
      card(card_header("Local statistics"), selectInput("padj", "p-value adjustment",
             c("none", "bonferroni", "fdr" = "BH")), tableOutput("ltab"))),
    card(card_header("How to read this"), uiOutput("ltxt"))),
  nav_panel("3. Point pattern (CSR)", value = "pp",
    intro(strong("Question: "), "are health centre locations consistent with Complete Spatial Randomness? ", strong("H0: "), "homogeneous Poisson process."), uiOutput("pp_head"), uiOutput("pp_ui")),
  nav_panel("About", value = "ab", card(uiOutput("about")))
)

# ------------------------------------------------------------- SERVER
server <- function(input, output, session) {

  observeEvent(input$go_ov, nav_select("nav", "ov"))
  observeEvent(input$go_gm, nav_select("nav", "gm"))
  observeEvent(input$go_lm, nav_select("nav", "lm"))
  observeEvent(input$go_pp, nav_select("nav", "pp"))
  observeEvent(input$nav, sidebar_toggle("side", open = !(input$nav %in% c("home", "ab"))))

  dat <- reactive({
    d <- panel[panel$Date == input$date, ]
    out <- prov_p %>% left_join(d, by = "province")
    out$value <- out[[input$ind]]
    out
  })
  v <- reactive({ x <- dat()$value; validate(need(!anyNA(x) && sd(x) > 0, "Indicator is missing or constant on this date.")); x })

  nb <- reactive({
    switch(input$nbtype,
      queen = poly2nb(prov_p, queen = TRUE, snap = 1000),
      rook  = poly2nb(prov_p, queen = FALSE, snap = 1000),
      knn   = knn2nb(knearneigh(cent, k = input$k)))
  })
  lw <- reactive(nb2listw(nb(), style = input$wstyle, zero.policy = TRUE))

  # ---------- Overview
  output$vb <- renderUI({
    d <- panel[panel$Date == input$date, ]
    layout_columns(
      value_box("Confirmed cases", format(sum(d$cases), big.mark = " "), theme = "primary"),
      value_box("Deaths", format(sum(d$deaths), big.mark = " "), theme = "danger"),
      value_box("Recoveries", format(sum(d$recoveries), big.mark = " "), theme = "success"),
      value_box("Health centres mapped", nrow(chc), theme = "secondary"))
  })
  output$map <- renderLeaflet({
    d <- st_transform(dat(), 4326); pal <- colorNumeric("YlOrRd", d$value)
    m <- leaflet(d) %>% addProviderTiles("CartoDB.Positron") %>%
      addPolygons(fillColor = ~pal(value), fillOpacity = 0.8, color = "white", weight = 1,
                  label = ~sprintf("%s: %s", province, round(value, 2))) %>%
      addLegend(pal = pal, values = ~value, title = names(indicators)[indicators == input$ind])
    if (isTRUE(input$showchc))
      m <- m %>% addCircleMarkers(data = chc, radius = 3, stroke = FALSE, fillOpacity = 0.7,
                                 color = "#1f4e79", label = ~Name)
    m
  })
  output$tsplot <- renderPlot({
    ggplot(panel, aes(Date, .data[[input$ind]], colour = province)) + geom_line() +
      geom_vline(xintercept = input$date, linetype = 2) +
      labs(y = names(indicators)[indicators == input$ind], colour = NULL) + theme_minimal(base_size = 13)
  })

  # ---------- 1. Global Moran
  gm <- reactive(moran.test(v(), lw(), alternative = input$alt, zero.policy = TRUE))
  output$gm_vb <- renderUI({
    m <- gm(); I <- unname(m$estimate[1]); EI <- unname(m$estimate[2])
    tagList(
      layout_columns(
        value_box("Moran's I", sprintf("%.3f", I), theme = "primary"),
        value_box("E[I] = -1/(n-1)", sprintf("%.3f", EI), theme = "secondary"),
        value_box("z-score", sprintf("%.2f", unname(m$statistic)), theme = "secondary"),
        value_box("p-value", sprintf("%.4f", m$p.value), theme = if (m$p.value < input$alpha) "success" else "warning")),
      p(strong("H0:"), " the spatial pattern is random (I = E[I]).  ", strong("Conclusion: "),
        moran_text(I, EI, m$p.value, input$alpha)))
  })
  output$gm_scatter <- renderPlot({
    z <- as.numeric(scale(v())); lz <- lag.listw(lw(), z, zero.policy = TRUE)
    df <- data.frame(z, lz, province = dat()$province)
    ggplot(df, aes(z, lz)) + geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
      geom_smooth(method = "lm", se = FALSE, colour = "#d7191c") + geom_point(size = 3, colour = "#1f4e79") +
      geom_text(aes(label = province), vjust = -0.9, size = 3.3) +
      labs(x = "Standardised value (z)", y = "Spatial lag of z (neighbour average)",
           subtitle = "Slope of the line = Moran's I when W is row-standardised") + theme_minimal(base_size = 13)
  })
  output$gm_mc <- renderPlot({
    set.seed(1)
    mc <- moran.mc(v(), lw(), nsim = input$mcsim, alternative = input$alt, zero.policy = TRUE)
    sims <- mc$res[seq_len(input$mcsim)]
    ggplot(data.frame(sims), aes(sims)) + geom_histogram(bins = 30, fill = "grey70", colour = "white") +
      geom_vline(xintercept = unname(mc$statistic), colour = "#d7191c", linewidth = 1.2) +
      labs(x = "Moran's I under random re-allocation of values", y = "Count",
           subtitle = sprintf("Observed I (red) = %.3f; Monte Carlo p = %.4f", mc$statistic, mc$p.value)) +
      theme_minimal(base_size = 13)
  })
  output$nbplot <- renderPlot({
    e <- nb_edges(nb(), cent)
    ggplot() + geom_sf(data = prov_p, fill = "grey95", colour = "grey50") +
      geom_segment(data = e, aes(x = x, y = y, xend = xend, yend = yend), colour = "#d7191c") +
      geom_point(data = data.frame(x = cent[, 1], y = cent[, 2]), aes(x, y)) + theme_void()
  })
  output$wtab <- renderTable({
    Wm <- round(listw2mat(lw()), 3); dimnames(Wm) <- list(prov_p$province, substr(prov_p$province, 1, 4))
    cbind(Province = rownames(Wm), as.data.frame(Wm), `Row sum` = round(rowSums(Wm), 3))
  })
  output$gm_time <- renderPlot({
    l <- lw(); dts <- sort(unique(panel$Date))
    res <- lapply(dts, function(dt) {
      d <- panel[panel$Date == dt, ]; x <- d[[input$ind]][match(prov_p$province, d$province)]
      if (anyNA(x) || sd(x) == 0) return(NULL)
      m <- moran.test(x, l, alternative = input$alt, zero.policy = TRUE)
      data.frame(Date = dt, I = unname(m$estimate[1]), EI = unname(m$estimate[2]), p = m$p.value)
    })
    r <- bind_rows(res); validate(need(nrow(r) > 0, "No valid dates."))
    r$sig <- ifelse(r$p < input$alpha, "Significant", "Not significant")
    ggplot(r, aes(Date, I)) + geom_hline(aes(yintercept = EI), linetype = 2, colour = "grey40") +
      geom_line(colour = "grey60") + geom_point(aes(colour = sig), size = 1.6) +
      scale_colour_manual(values = c("Significant" = "#d7191c", "Not significant" = "grey50"), name = NULL) +
      labs(y = "Moran's I", subtitle = "Dashed line = E[I]") + theme_minimal(base_size = 13)
  })

  # ---------- 2. Local Moran
  lisa <- reactive({
    x <- v(); l <- lw(); z <- as.numeric(scale(x)); lz <- lag.listw(l, z, zero.policy = TRUE)
    lm <- localmoran(x, l, alternative = "two.sided", zero.policy = TRUE)
    p <- p.adjust(lm[, 5], method = input$padj)
    quad <- ifelse(z > 0 & lz > 0, "High-High", ifelse(z < 0 & lz < 0, "Low-Low",
            ifelse(z > 0 & lz < 0, "High-Low", "Low-High")))
    cl <- ifelse(p < input$alpha, quad, "Not significant")
    data.frame(Province = prov_p$province, Value = x, Ii = lm[, 1], z = lm[, 4], p = p, Quadrant = quad, Cluster = cl)
  })
  output$lmap <- renderLeaflet({
    r <- lisa(); d <- st_transform(prov_p, 4326); d$Cluster <- factor(r$Cluster, levels = names(cluster_cols))
    pal <- colorFactor(cluster_cols, levels = names(cluster_cols))
    leaflet(d) %>% addProviderTiles("CartoDB.Positron") %>%
      addPolygons(fillColor = ~pal(Cluster), fillOpacity = 0.85, color = "white", weight = 1,
                  label = ~sprintf("%s: %s", province, Cluster)) %>%
      addLegend(colors = unname(cluster_cols), labels = names(cluster_cols), title = "LISA cluster")
  })
  output$ltab <- renderTable({ r <- lisa(); r$Value <- round(r$Value, 2); r$Ii <- round(r$Ii, 3)
    r$z <- round(r$z, 2); r$p <- round(r$p, 4); r[, c("Province", "Value", "Ii", "z", "p", "Cluster")] })
  output$ltxt <- renderUI(tagList(
    p("Local Moran's I decomposes the global statistic into one value per province. A province is a ",
      strong("High-High"), " cluster if it and its neighbours are above average, ", strong("Low-Low"),
      " if both are below, and a spatial outlier (", strong("High-Low / Low-High"),
      ") if it differs from its neighbours. Only provinces with p below the chosen level are coloured."),
    p(strong("Caution: "), "with n = 9 each province has only 1-5 neighbours, so local p-values are very unstable, and 9 simultaneous tests inflate false positives. Compare p-value adjustments.")))

  # ---------- 3. Point pattern
  output$pp_head <- renderUI(tagList(
    p(sprintf("%d health centres with valid coordinates are analysed (%d fell outside the simplified window).", npoints(X), n_dropped)),
    if (!win_ok) p(strong("Warning: "), "province window could not be built; a bounding box is used instead.")))
  output$pp_ui <- renderUI({
    switch(input$ppview,
      kde = card(card_header("Kernel intensity estimate"), plotOutput("kde", height = 520),
                 p("Intensity lambda(s) = expected CHCs per unit area. Bandwidth sigma controls smoothing: small = bumpy, large = smooth.")),
      quad = layout_columns(col_widths = c(7, 5), card(card_header("Quadrat counts"), plotOutput("quadplot", height = 460)),
                            card(card_header("Chi-square test of CSR"), uiOutput("quadtxt"))),
      fun = layout_columns(col_widths = c(8, 4), card(card_header("Observed vs CSR envelope"), plotOutput("funplot", height = 460)),
                           card(card_header("How to read"), uiOutput("funtxt"))))
  })
  output$kde <- renderPlot({
    d <- as.data.frame(density(X, sigma = input$bw))
    p <- ggplot() + geom_raster(data = d, aes(x, y, fill = value * 1e4)) +
      geom_sf(data = prov_km, fill = NA, colour = "white", linewidth = 0.3) +
      scale_fill_viridis_c(name = "CHCs per\n10 000 km^2") + theme_void(base_size = 13)
    if (isTRUE(input$showpts)) p <- p + geom_point(data = as.data.frame(X), aes(x, y), size = 0.6, colour = "white", alpha = 0.7)
    p
  })
  qt <- reactive(quadrat.test(X, nx = input$nx, ny = input$ny, method = "Chisq"))
  output$quadplot <- renderPlot(plot(qt(), main = ""))
  output$quadtxt <- renderUI({
    t <- qt(); cnt <- as.vector(quadratcount(X, nx = input$nx, ny = input$ny))
    tagList(
      p(strong("H0:"), " CHCs follow a homogeneous Poisson process (CSR)."),
      p(sprintf("X-squared = %.1f, df = %d, p = %.4g", unname(t$statistic), unname(t$parameter), t$p.value)),
      p(sprintf("Approx. variance-to-mean ratio of counts = %.2f (CSR is about 1; >1 clustered).", var(cnt) / mean(cnt))),
      p(em("Quadrats cut by the border have smaller area; the chi-square test accounts for this, the VMR does not. Expected counts below 5 weaken the test; try fewer quadrats."), p(em("Results depend on quadrat size (scale).")))
    )
  })
  fsim <- eventReactive(input$run, {
    withProgress(message = "Simulating CSR patterns...", {
      set.seed(42)
      switch(input$fn,
        K = envelope(X, Kest, nsim = input$nsim, correction = "border", verbose = FALSE),
        L = envelope(X, Lest, nsim = input$nsim, correction = "border", verbose = FALSE),
        G = envelope(X, Gest, nsim = input$nsim, correction = "km", verbose = FALSE),
        F = envelope(X, Fest, nsim = input$nsim, correction = "km", verbose = FALSE))
    })
  })
  output$funplot <- renderPlot({
    e <- fsim(); fn <- isolate(input$fn)
    if (fn == "L") plot(e, . - r ~ r, main = "L(r) - r", legend = FALSE) else plot(e, main = paste0(fn, "(r)"), legend = FALSE)
  })
  output$funtxt <- renderUI(tagList(
    p("Black = observed. Grey band = range of ", isolate(input$nsim), " simulated CSR patterns."),
    tags$ul(
      tags$li(strong("K / L: "), "above the band = more clustering than CSR at that distance; below = more regular."),
      tags$li(strong("G: "), "above the band = shorter nearest-neighbour distances (clustering)."),
      tags$li(strong("F: "), "reversed: below the band = large empty gaps (clustering).")),
    p(strong("Caveat: "), "CHCs are sited where people live, so the process is almost certainly not homogeneous. Rejection of CSR here mostly reflects varying intensity (first-order), not attraction between centres (second-order). Distances are in km; edge effects handled by border/Kaplan-Meier corrections."),
    p(em("Click 'Run simulation' in the sidebar after changing settings."))))

  output$about <- renderUI(tagList(
    h4("About this dashboard"),
    p("Data: provincial cumulative COVID-19 cases, deaths and recoveries (5 Mar - 4 Oct 2020) and 287 Community Health Centre locations (one without coordinates)."),
    p("Population denominators are Stats SA 2020 mid-year estimates entered manually (data/population_2020.csv) and should be verified."),
    p("Boundaries: GADM level 1. Methods: Global and Local Moran's I (spdep), CSR tests, kernel intensity and K/L/G/F functions (spatstat)."),
    p(strong("Limitations: "), "9 provinces is a very small lattice (low power, MAUP, ecological fallacy); the CHC list may be incomplete; reporting and testing differ by province.")))
}

shinyApp(ui, server)
