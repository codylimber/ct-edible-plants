#' Phenology Site Engine
#'
#' Generalized engine for building phenology calendar websites from iNaturalist
#' or GBIF observation data. Reads params.yml for all site-specific configuration.
#'
#' Usage: source("R/phenology_engine.R") from a QMD script after loading params.

# ============================================================================
# Configuration
# ============================================================================

#' Load and validate site parameters from params.yml
#'
#' @param path Path to params.yml (default: "params.yml" in project root)
#' @return Named list of parameters with defaults applied
load_params <- function(path = "params.yml") {
  params <- yaml::read_yaml(path)

  # Apply defaults
  defaults <- list(
    data_source = "inat",
    years = c(2015, 2025),
    min_obs = 20,
    quality = "research",
    mode = "auto",
    species_list = NULL,
    state_abbrev = NULL,
    county_year = 2020,
    color_scale = "inferno",
    plot_height = 1200,
    rarity = list(rare_below = 50, uncommon_below = 200),
    active_months = NULL
  )

  for (nm in names(defaults)) {
    if (is.null(params[[nm]])) params[[nm]] <- defaults[[nm]]
  }

  # Auto-detect terminology if not provided
  if (is.null(params$terminology)) {
    params$terminology <- detect_terminology(params$taxon_name)
  }

  # Expand years to sequence
  params$year_seq <- seq(params$years[1], params$years[2])

  params
}

#' Auto-detect terminology based on taxon name
#'
#' @param taxon_name Common taxon group name (e.g., "Odonata", "Lepidoptera")
#' @return Named list with period_name, multi_period, group labels
detect_terminology <- function(taxon_name) {
  tn <- tolower(taxon_name)

  if (grepl("odonata|dragonfl|damselfl", tn)) {
    list(period_name = "flight", multi_period = "flights",
         group_singular = "species", group_plural = "species",
         whats_active_verb = "flying", show_multi_period = TRUE)
  } else if (grepl("lepidoptera|butterfly|butterfl|papilion|moth", tn)) {
    list(period_name = "flight", multi_period = "broods",
         group_singular = "species", group_plural = "species",
         whats_active_verb = "flying", show_multi_period = TRUE)
  } else if (grepl("plant|flora|wildflower|flower", tn)) {
    list(period_name = "bloom", multi_period = "flushes",
         group_singular = "species", group_plural = "species",
         whats_active_verb = "blooming", show_multi_period = TRUE)
  } else if (grepl("fungi|mushroom", tn)) {
    list(period_name = "fruiting", multi_period = "flushes",
         group_singular = "species", group_plural = "species",
         whats_active_verb = "fruiting", show_multi_period = TRUE)
  } else if (grepl("bird|aves", tn)) {
    list(period_name = "active period", multi_period = "waves",
         group_singular = "species", group_plural = "species",
         whats_active_verb = "active", show_multi_period = FALSE)
  } else {
    list(period_name = "active period", multi_period = "periods",
         group_singular = "species", group_plural = "species",
         whats_active_verb = "active", show_multi_period = FALSE)
  }
}


# ============================================================================
# Data fetching
# ============================================================================

#' Fetch all observations for a taxon in a place, paginated by year + month
#'
#' Used in auto-discovery mode (fetch everything in a higher taxon).
#'
#' @param params Site parameters list
#' @param cache_dir Optional cache directory for RDS files
#' @return Data frame of raw observations
fetch_observations_auto <- function(params, cache_dir = NULL) {
  cache_file <- NULL
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    cache_file <- file.path(cache_dir, paste0("taxon_", params$taxon_id, "_all.rds"))
    if (file.exists(cache_file)) {
      cat("Loading cached observations...\n")
      return(readRDS(cache_file))
    }
  }

  years <- params$year_seq
  months <- 1:12

  grid <- tidyr::expand_grid(yr = years, mo = months)
  obs_list <- purrr::pmap(grid, \(yr, mo) {
      cat("Fetching", yr, month.abb[mo], "... ")
      Sys.sleep(1)
      tryCatch({
        result <- rinat::get_inat_obs(
          taxon_id   = params$taxon_id,
          place_id   = params$place_id,
          quality    = params$quality,
          geo        = TRUE,
          maxresults = 10000,
          meta       = FALSE,
          year       = yr,
          month      = mo
        )
        cat(nrow(result), "obs\n")
        result
      },
      error = function(e) {
        cat("(no results:", conditionMessage(e), ")\n")
        NULL
      })
    }) |>
    purrr::compact()

  obs_raw <- dplyr::bind_rows(obs_list)

  # Check for months that hit the API cap
  month_counts <- obs_raw |>
    dplyr::mutate(yr = lubridate::year(as.Date(datetime)),
                  mo = lubridate::month(as.Date(datetime))) |>
    dplyr::count(yr, mo) |>
    dplyr::filter(n >= 9999)
  if (nrow(month_counts) > 0) {
    warning("Some year-months hit the 10k cap -- results may be incomplete:\n",
            paste(month_counts$yr, month_counts$mo, sep = "-", collapse = ", "))
  }

  if (!is.null(cache_file)) saveRDS(obs_raw, cache_file)
  cat("Total observations fetched:", nrow(obs_raw), "\n")
  obs_raw
}

#' Fetch observations for a curated species list
#'
#' Used when params$mode == "curated" and a species_list CSV is provided.
#'
#' @param params Site parameters list
#' @param cache_dir Optional cache directory
#' @return Data frame of raw observations with .taxon_id_query column
fetch_observations_curated <- function(params, cache_dir = NULL) {
  species_df <- readr::read_csv(params$species_list, show_col_types = FALSE)

  if (!is.null(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  obs_list <- purrr::pmap(
    list(species_df$inat_taxon_id, species_df$common_name),
    \(tid, name) {
      if (!is.null(cache_dir)) {
        cache_file <- file.path(cache_dir, paste0("taxon_", tid, ".rds"))
        if (file.exists(cache_file)) {
          cat("Cached:", name, "(taxon", tid, ")\n")
          return(readRDS(cache_file))
        }
      }

      cat("Fetching", name, "(taxon", tid, ") ... ")
      result <- fetch_species_obs(tid, params$place_id, params$year_seq,
                                  params$quality)
      n <- if (is.null(result)) 0L else nrow(result)
      cat(n, "obs\n")

      if (!is.null(result) && n > 0) {
        result <- result |> dplyr::mutate(.taxon_id_query = tid)
        if (!is.null(cache_dir)) {
          saveRDS(result, file.path(cache_dir, paste0("taxon_", tid, ".rds")))
        }
      }
      result
    }
  ) |> purrr::compact()

  dplyr::bind_rows(obs_list)
}

#' Fetch observations for a single taxon, with smart pagination
#'
#' Tries a single API call first; paginates by year if results exceed 9k.
#'
#' @param taxon_id iNaturalist taxon ID
#' @param place_id iNaturalist place ID
#' @param years Integer vector of years
#' @param quality Quality grade filter
#' @return Data frame of observations, or NULL
fetch_species_obs <- function(taxon_id, place_id, years, quality = "research") {
  all_at_once <- tryCatch({
    Sys.sleep(1)
    rinat::get_inat_obs(
      taxon_id   = taxon_id,
      place_id   = place_id,
      quality    = quality,
      geo        = TRUE,
      maxresults = 10000,
      meta       = FALSE
    )
  }, error = function(e) NULL)

  if (!is.null(all_at_once) && nrow(all_at_once) < 9000) return(all_at_once)

  obs_list <- purrr::map(years, \(yr) {
    Sys.sleep(1)
    tryCatch({
      rinat::get_inat_obs(
        taxon_id = taxon_id, place_id = place_id,
        quality = quality, geo = TRUE, maxresults = 10000,
        meta = FALSE, year = yr
      )
    }, error = function(e) NULL)
  }) |> purrr::compact()

  if (length(obs_list) == 0) return(all_at_once)
  dplyr::bind_rows(obs_list)
}

#' Main fetch dispatcher -- routes by data_source and mode
#'
#' @param params Site parameters list
#' @param cache_dir Optional cache directory
#' @return Data frame of raw observations
fetch_observations <- function(params, cache_dir = NULL) {
  if (params$data_source == "gbif") {
    fetch_observations_gbif(params, cache_dir)
  } else if (params$mode == "curated") {
    fetch_observations_curated(params, cache_dir)
  } else {
    fetch_observations_auto(params, cache_dir)
  }
}


# ============================================================================
# GBIF data fetching
# ============================================================================

#' Standardize GBIF occurrence data to match iNat column format
#'
#' Maps GBIF column names to the names expected by clean_observations().
#'
#' @param df Data frame from rgbif::occ_search()
#' @return Data frame with columns: datetime, scientific_name, common_name,
#'   latitude, longitude, taxon_id
standardize_gbif <- function(df) {
  df |>
    dplyr::transmute(
      datetime       = eventDate,
      scientific_name = species,
      common_name    = dplyr::if_else(
        is.na(vernacularName), "", vernacularName
      ),
      latitude       = decimalLatitude,
      longitude      = decimalLongitude,
      taxon_id       = speciesKey
    )
}

#' Fetch GBIF observations, paginated by year + month
#'
#' occ_search() returns max 100k records per call. We paginate by year+month
#' to stay within limits. Results are standardized to match the iNat column
#' format so downstream cleaning/phenology code works unchanged.
#'
#' @param params Site parameters list (needs gbif_taxon_key, country_code)
#' @param cache_dir Optional cache directory for RDS files
#' @return Data frame of raw observations in iNat-compatible format
fetch_observations_gbif <- function(params, cache_dir = NULL) {
  cache_file <- NULL
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    cache_file <- file.path(cache_dir,
                            paste0("gbif_", params$gbif_taxon_key, "_all.rds"))
    if (file.exists(cache_file)) {
      cat("Loading cached GBIF observations...\n")
      return(readRDS(cache_file))
    }
  }

  years <- params$year_seq
  months <- 1:12
  gbif_limit <- 10000  # cap per call; peak months may exceed but phenology is robust

  grid <- tidyr::expand_grid(yr = years, mo = months)
  obs_list <- purrr::pmap(grid, \(yr, mo) {
      cat("GBIF:", yr, month.abb[mo], "... ")
      Sys.sleep(0.5)  # be polite to GBIF API
      tryCatch({
        result <- rgbif::occ_search(
          taxonKey       = params$gbif_taxon_key,
          country        = params$country_code,
          hasCoordinate  = TRUE,
          year           = yr,
          month          = mo,
          occurrenceStatus = "PRESENT",
          limit          = gbif_limit
        )
        if (is.null(result$data) || nrow(result$data) == 0) {
          cat("0 obs\n")
          return(NULL)
        }
        n <- nrow(result$data)
        cat(n, "obs")
        if (n >= gbif_limit) cat(" [HIT LIMIT]")
        cat("\n")
        standardize_gbif(result$data)
      },
      error = function(e) {
        cat("(error:", conditionMessage(e), ")\n")
        NULL
      })
    }) |>
    purrr::compact()

  obs_raw <- dplyr::bind_rows(obs_list)

  if (!is.null(cache_file)) saveRDS(obs_raw, cache_file)
  cat("Total GBIF observations fetched:", nrow(obs_raw), "\n")
  obs_raw
}


# ============================================================================
# Data cleaning
# ============================================================================

#' Clean raw observations: parse dates, extract species, resolve common names
#'
#' @param obs_raw Raw observation data frame
#' @param min_obs Minimum observations to include a species
#' @return List with obs (filtered), species_counts, common_names
clean_observations <- function(obs_raw, min_obs = 20) {
  obs <- obs_raw |>
    dplyr::mutate(
      date = as.Date(datetime),
      week = lubridate::isoweek(date),
      species = stringr::word(stringr::str_trim(scientific_name), 1, 2),
      common_name = stringr::str_trim(common_name)
    ) |>
    dplyr::filter(!is.na(date), !is.na(species), species != "")

  # Most common name per species
  common_names <- obs |>
    dplyr::filter(common_name != "") |>
    dplyr::count(species, common_name, sort = TRUE) |>
    dplyr::group_by(species) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(species, common_name)

  # Species with enough observations
  species_counts <- obs |>
    dplyr::count(species, name = "total_obs") |>
    dplyr::filter(total_obs >= min_obs) |>
    dplyr::left_join(common_names, by = "species") |>
    dplyr::mutate(
      display_name = dplyr::if_else(
        is.na(common_name) | common_name == "",
        species,
        paste0(common_name, " (", species, ")")
      )
    ) |>
    dplyr::arrange(dplyr::desc(total_obs))

  obs_filtered <- obs |>
    dplyr::semi_join(species_counts, by = "species")

  list(obs = obs_filtered, species_counts = species_counts,
       common_names = common_names)
}


# ============================================================================
# County assignment (US states only)
# ============================================================================

#' Assign observations to counties via spatial join
#'
#' @param obs Cleaned observation data frame (must have latitude, longitude)
#' @param state_abbrev Two-letter US state abbreviation (e.g., "CT")
#' @param county_year Census year for county boundaries (default 2020)
#' @return List with obs_county (observations with county column),
#'   location_hidden_count, county_boundaries
assign_counties <- function(obs, state_abbrev, county_year = 2020) {
  counties_sf <- tigris::counties(state = state_abbrev, cb = TRUE,
                                  year = county_year) |>
    sf::st_transform(4326) |>
    dplyr::select(county = NAME, geometry)

  obs_with_coords <- obs |>
    dplyr::filter(!is.na(latitude), !is.na(longitude))

  obs_sf <- sf::st_as_sf(obs_with_coords,
                          coords = c("longitude", "latitude"), crs = 4326)

  obs_county <- sf::st_join(obs_sf, counties_sf, join = sf::st_within) |>
    sf::st_drop_geometry()

  location_hidden_count <- nrow(obs) - nrow(obs_with_coords)

  list(obs_county = obs_county,
       location_hidden_count = location_hidden_count,
       counties_sf = counties_sf)
}


# ============================================================================
# Rarity and conservation status
# ============================================================================

#' Assign rarity tiers based on observation counts
#'
#' @param species_counts Data frame with total_obs column
#' @param rare_below Threshold for "Rare" (default 50)
#' @param uncommon_below Threshold for "Uncommon" (default 200)
#' @return species_counts with rarity column added
assign_rarity <- function(species_counts, rare_below = 50,
                          uncommon_below = 200) {
  species_counts |>
    dplyr::mutate(
      rarity = dplyr::case_when(
        total_obs < rare_below  ~ "Rare",
        total_obs <= uncommon_below ~ "Uncommon",
        TRUE ~ "Common"
      ),
      rarity = factor(rarity, levels = c("Common", "Uncommon", "Rare"))
    )
}

#' Fetch conservation status from iNat API for a single taxon
#'
#' @param tid iNaturalist taxon ID
#' @param place_pattern Regex pattern for place matching in NatureServe
#'   (e.g., "Connecticut|CT")
#' @return Tibble with taxon_id, conservation, concern_level
fetch_conservation_single <- function(tid, place_pattern = NULL) {
  url <- paste0("https://api.inaturalist.org/v1/taxa/", tid)
  resp <- tryCatch(jsonlite::fromJSON(url, flatten = TRUE),
                   error = function(e) NULL)
  if (is.null(resp)) {
    return(tibble::tibble(taxon_id = tid, conservation = NA_character_))
  }

  statuses <- tryCatch(resp$results$conservation_statuses[[1]],
                       error = function(e) NULL)
  if (is.null(statuses) || !is.data.frame(statuses) || nrow(statuses) == 0) {
    return(tibble::tibble(taxon_id = tid, conservation = "Not assessed"))
  }
  if (!"authority" %in% names(statuses)) {
    return(tibble::tibble(taxon_id = tid, conservation = "Not assessed"))
  }

  # Prefer IUCN global
  iucn <- statuses |> dplyr::filter(authority == "IUCN Red List") |> head(1)
  if (nrow(iucn) > 0) {
    return(tibble::tibble(taxon_id = tid,
                          conservation = paste("IUCN:", iucn$status)))
  }

  # NatureServe -- try place-specific first, then global
  has_place <- "place.display_name" %in% names(statuses)
  if (has_place && !is.null(place_pattern)) {
    ns_local <- statuses |>
      dplyr::filter(
        stringr::str_detect(authority, "NatureServe"),
        stringr::str_detect(
          tidyr::replace_na(place.display_name, ""), place_pattern
        )
      ) |> head(1)
    if (nrow(ns_local) > 0) {
      return(tibble::tibble(taxon_id = tid,
                            conservation = paste("NS-local:", ns_local$status)))
    }
  }

  if (has_place) {
    ns_global <- statuses |>
      dplyr::filter(
        stringr::str_detect(authority, "NatureServe"),
        is.na(place.display_name) | place.display_name == ""
      ) |> head(1)
    if (nrow(ns_global) > 0) {
      return(tibble::tibble(taxon_id = tid,
                            conservation = paste("NS:", ns_global$status)))
    }
  }

  tibble::tibble(taxon_id = tid, conservation = "Not assessed")
}

#' Fetch conservation status for all species and add concern levels
#'
#' For iNat data, queries the iNat API directly using taxon IDs.
#' For GBIF data, looks up iNat taxon IDs by species name first, then
#' queries the iNat API for conservation status.
#'
#' @param obs_filtered Filtered observation data (to extract taxon IDs)
#' @param place_pattern Regex for NatureServe place matching
#' @param data_source "inat" or "gbif" — determines ID lookup strategy
#' @return Data frame with taxon_id, conservation, concern_level
fetch_conservation_all <- function(obs_filtered, place_pattern = NULL,
                                   data_source = "inat") {
  if (data_source == "gbif") {
    taxon_ids <- resolve_gbif_to_inat_ids(obs_filtered)
  } else {
    taxon_ids <- obs_filtered |>
      dplyr::count(species, taxon_id, sort = TRUE) |>
      dplyr::group_by(species) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::select(species, taxon_id)
  }

  # Drop species where we couldn't resolve an iNat taxon ID
  taxon_ids_valid <- taxon_ids |> dplyr::filter(!is.na(taxon_id))

  cat("Fetching conservation status for", nrow(taxon_ids_valid), "species...\n")
  cons_df <- purrr::map(taxon_ids_valid$taxon_id, \(tid) {
    Sys.sleep(0.5)
    fetch_conservation_single(tid, place_pattern)
  }, .progress = TRUE) |>
    dplyr::bind_rows()

  # Add back unresolved species as "Unknown"
  unresolved <- taxon_ids |>
    dplyr::filter(is.na(taxon_id)) |>
    dplyr::mutate(conservation = NA_character_)
  if (nrow(unresolved) > 0) {
    cons_df <- dplyr::bind_rows(
      cons_df,
      unresolved |> dplyr::select(taxon_id, conservation)
    )
  }

  cons_df <- cons_df |>
    dplyr::mutate(
      concern_level = dplyr::case_when(
        stringr::str_detect(conservation, "EN|S1|Endangered|CR") ~ "Endangered",
        stringr::str_detect(conservation, "VU|S2|Threatened|NT") ~ "Vulnerable",
        stringr::str_detect(conservation, "LC|S[3-5]|Secure|Not assessed") ~
          "Least Concern",
        is.na(conservation) ~ "Unknown",
        TRUE ~ "Least Concern"
      ),
      concern_level = factor(
        concern_level,
        levels = c("Least Concern", "Vulnerable", "Endangered", "Unknown")
      )
    )

  list(taxon_ids = taxon_ids, conservation = cons_df)
}

#' Resolve GBIF species names to iNaturalist taxon IDs
#'
#' For GBIF-sourced data, taxon_id contains GBIF speciesKey values. This
#' function queries the iNat taxa autocomplete API by species name to find
#' matching iNat taxon IDs for conservation status lookups.
#'
#' @param obs_filtered Cleaned observation data with species column
#' @return Tibble with species and taxon_id (iNat ID, NA if unresolved)
resolve_gbif_to_inat_ids <- function(obs_filtered) {
  species_list <- obs_filtered |>
    dplyr::distinct(species) |>
    dplyr::pull(species)

  cat("Resolving", length(species_list), "species to iNat taxon IDs...\n")
  purrr::map_dfr(species_list, \(sp) {
    Sys.sleep(0.5)
    url <- paste0("https://api.inaturalist.org/v1/taxa/autocomplete?q=",
                  utils::URLencode(sp), "&rank=species&per_page=1")
    resp <- tryCatch(jsonlite::fromJSON(url, flatten = TRUE),
                     error = function(e) NULL)
    if (!is.null(resp) && resp$total_results > 0) {
      tibble::tibble(species = sp, taxon_id = resp$results$id[1])
    } else {
      tibble::tibble(species = sp, taxon_id = NA_integer_)
    }
  }, .progress = TRUE)
}


# ============================================================================
# Phenology computation
# ============================================================================

#' Compute weekly phenology: relative abundance, peak weeks, active periods
#'
#' @param obs_filtered Cleaned, filtered observation data frame
#' @param species_counts Species counts data frame
#' @param params Site parameters (for terminology)
#' @return List with weekly, species_summary
compute_phenology <- function(obs_filtered, species_counts, params) {
  # Weekly observation counts

  weekly <- obs_filtered |>
    dplyr::count(species, week) |>
    tidyr::complete(species, week = 1:52, fill = list(n = 0)) |>
    dplyr::group_by(species) |>
    dplyr::mutate(rel_abundance = n / max(n)) |>
    dplyr::ungroup()

  # Peak week per species
  peak_weeks <- weekly |>
    dplyr::group_by(species) |>
    dplyr::slice_max(rel_abundance, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(species, peak_week = week)

  # Active period (weeks >= 25% of peak)
  active_periods <- weekly |>
    dplyr::filter(rel_abundance >= 0.25) |>
    dplyr::group_by(species) |>
    dplyr::summarize(
      first_week = min(week),
      last_week  = max(week),
      .groups = "drop"
    )

  # Multi-period detection
  period_info <- detect_multi_periods(weekly)

  # Assemble species summary
  species_counts_dedup <- species_counts |>
    dplyr::distinct(species, .keep_all = TRUE)

  species_summary <- species_counts_dedup |>
    dplyr::left_join(peak_weeks, by = "species") |>
    dplyr::left_join(active_periods, by = "species") |>
    dplyr::left_join(period_info, by = "species") |>
    dplyr::arrange(peak_week, first_week)

  list(weekly = weekly, species_summary = species_summary)
}

#' Detect species with multiple active periods (flights, broods, flushes)
#'
#' A species has multiple periods if there's a gap of >= 3 weeks with < 10%
#' peak abundance between active periods.
#'
#' @param weekly Weekly phenology data frame
#' @return Data frame with species and period_count
detect_multi_periods <- function(weekly) {
  weekly |>
    dplyr::group_by(species) |>
    dplyr::mutate(active = rel_abundance >= 0.10) |>
    dplyr::summarize(
      period_count = {
        r <- rle(active)
        active_runs <- which(r$values)
        if (length(active_runs) <= 1) {
          1L
        } else {
          gap_lengths <- r$lengths[active_runs[-length(active_runs)] + 1]
          1L + sum(gap_lengths >= 3)
        }
      },
      .groups = "drop"
    )
}


# ============================================================================
# Plotting
# ============================================================================

#' Build the interactive heatmap calendar
#'
#' @param weekly Weekly phenology data
#' @param species_summary Species summary data
#' @param params Site parameters
#' @return A plotly object
build_calendar_plot <- function(weekly, species_summary, params) {
  term <- params$terminology
  week_dates <- week_date_lookup()

  # Build display labels
  species_summary <- species_summary |>
    dplyr::mutate(
      rarity_icon = dplyr::case_when(
        rarity == "Rare"     ~ "\U0001f534",
        rarity == "Uncommon" ~ "\U0001f7e1",
        TRUE ~ ""
      ),
      period_label = dplyr::if_else(
        period_count > 1,
        paste0(" [", period_count, " ", term$multi_period, "]"),
        ""
      ),
      plot_label = paste0(
        dplyr::if_else(rarity_icon != "", paste0(rarity_icon, " "), ""),
        display_name, period_label
      )
    )

  display_order <- species_summary |>
    dplyr::arrange(dplyr::desc(peak_week), dplyr::desc(first_week)) |>
    dplyr::pull(plot_label)

  plot_data <- weekly |>
    dplyr::left_join(week_dates, by = "week") |>
    dplyr::left_join(
      species_summary |>
        dplyr::select(species, plot_label, common_name, rarity, period_count),
      by = "species"
    ) |>
    dplyr::mutate(
      plot_label = factor(plot_label, levels = display_order),
      hover_text = paste0(
        common_name, " (", species, ")\n",
        "Week of ", format(approx_date, "%b %d"), "\n",
        "Observations: ", n, "\n",
        "Relative abundance: ", round(rel_abundance * 100), "%\n",
        "Rarity: ", rarity,
        dplyr::if_else(
          period_count > 1,
          paste0("\n", stringr::str_to_title(term$multi_period), ": ",
                 period_count),
          ""
        )
      )
    )

  month_starts <- as.Date(paste0("2024-", 1:12, "-01"))
  week_starts  <- as.Date("2024-01-01") + (0:51) * 7

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = approx_date, y = plot_label,
                 fill = rel_abundance, text = hover_text)
  ) +
    ggplot2::geom_vline(xintercept = week_starts, color = "grey80",
                        linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = month_starts, color = "grey50",
                        linewidth = 0.5) +
    ggplot2::geom_tile(height = 0.8) +
    ggplot2::scale_fill_viridis_c(
      option = params$color_scale,
      name = "Relative\nabundance",
      limits = c(0, 1)
    ) +
    ggplot2::scale_x_date(
      date_labels = "%b",
      date_breaks = "1 month",
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      title = paste(params$place_name, params$taxon_common_name,
                    stringr::str_to_title(term$period_name), "Periods"),
      subtitle = paste0(
        if (params$data_source == "gbif") "GBIF" else "iNaturalist",
        if (params$data_source != "gbif") paste0(" ", params$quality, "-grade") else "",
        " observations, ",
        params$years[1], "\u2013", params$years[2],
        " (species with \u2265 ", params$min_obs, " obs)\n",
        "\U0001f534 = rare  \U0001f7e1 = uncommon"
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  # Return both ggplot and plotly
  list(
    gg = p,
    plotly = plotly::ggplotly(p, tooltip = "text",
                              height = params$plot_height) |>
      plotly::layout(yaxis = list(tickfont = list(size = 9)))
  )
}


# ============================================================================
# Summary tables
# ============================================================================

#' Build the species summary datatable
#'
#' @param species_summary Species summary data frame
#' @param params Site parameters
#' @return DT datatable object
build_summary_table <- function(species_summary, params) {
  term <- params$terminology
  week_dates <- week_date_lookup()

  summary_table <- species_summary |>
    dplyr::left_join(week_dates |> dplyr::rename(peak_date = approx_date),
                     by = c("peak_week" = "week")) |>
    dplyr::left_join(week_dates |> dplyr::rename(start_date = approx_date),
                     by = c("first_week" = "week")) |>
    dplyr::left_join(week_dates |> dplyr::rename(end_date = approx_date),
                     by = c("last_week" = "week")) |>
    dplyr::mutate(
      peak     = format(peak_date, "%b %d"),
      active   = paste(format(start_date, "%b %d"), "\u2013",
                        format(end_date, "%b %d")),
      duration = last_week - first_week + 1
    ) |>
    dplyr::select(
      species, common_name, total_obs, rarity,
      dplyr::any_of(c("concern_level", "conservation")),
      period_count, peak, active, duration_weeks = duration
    ) |>
    dplyr::arrange(species)

  list(
    table = summary_table,
    dt = DT::datatable(
      summary_table,
      filter = "top",
      options = list(pageLength = 20, autoWidth = TRUE),
      caption = paste(params$taxon_common_name, term$period_name,
                      "periods in", params$place_name)
    )
  )
}


# ============================================================================
# "What's active?" query
# ============================================================================

#' Query which species are active for a given date
#'
#' @param target_date Date to query (default: today)
#' @param weekly Weekly phenology data
#' @param species_summary Species summary data
#' @param params Site parameters
#' @return Data frame of active species sorted by relative abundance
whats_active <- function(target_date = Sys.Date(), weekly, species_summary,
                         params) {
  target_week <- lubridate::isoweek(target_date)

  weekly |>
    dplyr::filter(week == target_week, rel_abundance > 0) |>
    dplyr::left_join(
      species_summary |>
        dplyr::select(species, common_name, rarity,
                      dplyr::any_of("concern_level"),
                      period_count, peak_week),
      by = "species"
    ) |>
    dplyr::mutate(
      distance_from_peak = abs(target_week - peak_week),
      status = dplyr::case_when(
        rel_abundance >= 0.75 ~ "Peak",
        rel_abundance >= 0.25 ~ "Active",
        target_week <= peak_week ~ "Early",
        TRUE ~ "Late"
      )
    ) |>
    dplyr::arrange(dplyr::desc(rel_abundance)) |>
    dplyr::select(common_name, species, status, rel_abundance, rarity,
                  dplyr::any_of("concern_level"), period_count)
}

#' Build monthly overview tables
#'
#' @param weekly Weekly phenology data
#' @param species_summary Species summary data
#' @param params Site parameters
#' @param months Months to show (default: Apr-Oct)
#' @return Data frame with month column
build_monthly_overview <- function(weekly, species_summary, params,
                                   months = 4:10) {
  month_samples <- as.Date(paste0("2024-", months, "-01"))

  purrr::map_dfr(month_samples, \(d) {
    result <- whats_active(d, weekly, species_summary, params)
    if (nrow(result) > 0) {
      result |>
        dplyr::mutate(month = format(d, "%B")) |>
        dplyr::select(month, common_name, species, status, rarity)
    }
  })
}


# ============================================================================
# JSON export
# ============================================================================

#' Export phenology data as JSON
#'
#' @param weekly Weekly phenology data
#' @param species_summary Species summary data
#' @param params Site parameters
#' @param obs_raw Raw observation data (for total count)
#' @param county_data Optional county data list from assign_counties()
#' @param out_path Output file path
export_json <- function(weekly, species_summary, params, obs_raw,
                        county_data = NULL, out_path) {
  weekly_sparse <- weekly |>
    dplyr::filter(n > 0) |>
    dplyr::select(species, week, n, rel_abundance)

  species_export <- species_summary |>
    dplyr::select(species, common_name,
                  dplyr::any_of(c("taxon_id", "total_obs", "rarity",
                                  "period_count", "peak_week",
                                  "first_week", "last_week")))

  export_data <- list(
    metadata = list(
      last_updated = as.character(Sys.Date()),
      data_source = params$data_source,
      taxon = params$taxon_name,
      place = params$place_name,
      date_range = paste0(params$years[1], "\u2013", params$years[2]),
      total_observations = nrow(obs_raw),
      species_count = nrow(species_summary),
      min_obs_threshold = params$min_obs
    ),
    species = species_export,
    weekly = weekly_sparse
  )

  # Add county data if available

  if (!is.null(county_data)) {
    county_weekly <- county_data$obs_county |>
      dplyr::filter(!is.na(county)) |>
      dplyr::count(species, county, week) |>
      tidyr::complete(species, county, week = 1:52, fill = list(n = 0)) |>
      dplyr::group_by(species, county) |>
      dplyr::mutate(rel_abundance = n / max(max(n), 1)) |>
      dplyr::ungroup() |>
      dplyr::filter(n > 0) |>
      dplyr::select(species, county, week, n, rel_abundance)

    export_data$metadata$location_hidden_count <-
      county_data$location_hidden_count
    export_data$metadata$counties <-
      sort(unique(county_data$obs_county$county[
        !is.na(county_data$obs_county$county)
      ]))
    export_data$county_weekly <- county_weekly
  }

  jsonlite::write_json(export_data, out_path, auto_unbox = TRUE, pretty = FALSE)
  cat("Exported JSON:", file.size(out_path), "bytes\n")
}


# ============================================================================
# Reference helpers
# ============================================================================

#' Week-to-date lookup table (using 2024 as reference year)
week_date_lookup <- function() {
  tibble::tibble(
    week = 1:52,
    approx_date = as.Date("2024-01-01") + (week - 1) * 7
  )
}

#' Format a week number as approximate date string
format_week <- function(week, fmt = "%b %d") {
  dates <- as.Date("2024-01-01") + (week - 1) * 7
  format(dates, fmt)
}


# ============================================================================
# Build info
# ============================================================================

#' Write BUILD_INFO.txt
#'
#' @param out_dir Output directory
#' @param params Site parameters
#' @param obs_raw Raw observations (for count)
#' @param species_count Number of species plotted
write_build_info <- function(out_dir, params, obs_raw, species_count) {
  writeLines(
    c(
      paste("site:", params$site_title),
      paste("taxon:", params$taxon_name, "(id:", params$taxon_id, ")"),
      paste("place:", params$place_name, "(id:", params$place_id, ")"),
      paste("git_hash:", tryCatch(
        system("git rev-parse --short HEAD", intern = TRUE),
        error = function(e) "unknown"
      )),
      paste("date:", Sys.time()),
      paste("observations:", nrow(obs_raw)),
      paste("species_plotted:", species_count)
    ),
    file.path(out_dir, "BUILD_INFO.txt")
  )
}