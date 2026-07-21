# Configuration
config <- list(
  survey_file   = Sys.getenv("SURVEY_FILE"),
  metadata_file = Sys.getenv("METADATA_FILE"),
  trial_folder  = Sys.getenv("TRIAL_FOLDER"),
  output_file   = Sys.getenv("OUTPUT_FILE", unset = file.path("outputs", "extracted_accelerometer_features.csv")),
  survey_names = list(
    baseline = c(
      "Mood Survey - Baseline",
      "Wellbeing Survey - Baseline"
    ),
    mid = c(
      "Mood Survey - Mid - Intervention A",
      "Mood Survey - Mid - Intervention B",
      "Mood Survey - Mid - Intervention C",
      "Mood Survey - Mid - Control"
    ),
    post = c(
      "Mood Survey - Post - Intervention A",
      "Mood Survey - Post - Intervention B",
      "Mood Survey - Post - Intervention C",
      "Mood Survey - Post - Control",
      "Wellbeing Survey - Post - Control",
      "Wellbeing Survey - Post - Intervention A",
      "Wellbeing Survey - Post - Intervention B",
      "Wellbeing Survey - Post - Intervention C",
      "Feedback Survey - Post - Control",
      "Feedback Survey - Post - Intervention A",
      "Feedback Survey - Post - Intervention B",
      "Feedback Survey - Post - Intervention C"
    )
  )
)

  # Load required packages
  load_packages <- function() {
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load(
      mhealthtools, arrow, data.table, dplyr, ggplot2, tidyr, tidyverse,
      haven, future, furrr, signal, zoo, pracma, e1071, seewave, statcomp,
      lubridate
    )
  }
  
  # Fast base-R time-domain summary
  time_domain_summary_base <- function(x) {
    x <- x[!is.na(x)]
    nm <- c("mean","median","mode","mx","mn","sd",
            "skewness","kurtosis","Q25","Q75","range","IQR",
            "energy","rmsmag","rugo")
    if (length(x) == 0) return(setNames(rep(NA_real_, length(nm)), nm))
    mu   <- mean(x); med <- median(x)
    ux   <- unique(x); mode_val <- ux[which.max(tabulate(match(x, ux)))]
    mx   <- max(x); mn <- min(x)
    sdev <- if (length(x) > 1) sd(x) else NA_real_
    if (requireNamespace("e1071", quietly=TRUE) && length(x) >= 4) {
      skew <- e1071::skewness(x, type=2)
      kurt <- e1071::kurtosis(x, type=2)
    } else {
      skew <- NA_real_; kurt <- NA_real_
    }
    q <- quantile(x, c(0.25, 0.75), names=FALSE)
    q25 <- q[1]; q75 <- q[2]; rng <- mx - mn; iqr <- q75 - q25
    energy <- sum(x^2); rmsmag <- sqrt(mean(x^2))
    rug <- seewave::rugo(x)
    setNames(c(mu, med, mode_val, mx, mn, sdev, skew, kurt,
               q25, q75, rng, iqr, energy, rmsmag, rug), nm)
  }
  
  # Prepare surveys function with dynamic timepoints
  prepare_surveys <- function(t0, t1) {
    survey_data <- read.csv(file.path(config$survey_file))
    
    # last survey at t0 (fixed column names)
    survey_t0 <- survey_data %>%
      dplyr::filter(survey_name %in% config$survey_names[[t0]]) %>%
      group_by(participant_id) %>%
      slice_max(order_by = recorded_at, n = 1) %>%
      rename(survey_name_t0 = survey_name) %>%
      mutate(
        last_recorded_at_t0 = as.POSIXct(
          recorded_at, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"
        ) %>% with_tz(.data$time_zone)
      ) %>%
      select(participant_id, survey_name_t0, last_recorded_at_t0)
    
    # last survey at t1 (fixed column names)
    survey_t1 <- survey_data %>%
      dplyr::filter(survey_name %in% config$survey_names[[t1]]) %>%
      group_by(participant_id) %>%
      slice_max(order_by = recorded_at, n = 1) %>%
      rename(survey_name_t1 = survey_name) %>%
      mutate(
        last_recorded_at_t1 = as.POSIXct(
          recorded_at, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"
        ) %>% with_tz(.data$time_zone)
      ) %>%
      select(participant_id, survey_name_t1, last_recorded_at_t1)
    
    list(
      survey_t0 = survey_t0,
      survey_t1 = survey_t1
    )
  }
  
  # Read metadata
  read_metadata <- function() {
    arrow::read_parquet(
      file.path(config$metadata_file),
      col_select = c("submission_id","participant_id","platform","localtime","time_zone")
    )
  }
  
  # Process one participant folder
  process_file <- function(participant_folder, metadata, surveys) {
    start.time <- Sys.time()
    parquet_files <- list.files(participant_folder, pattern = "\\.parquet$", full.names = TRUE)
    if (length(parquet_files) == 0) return(NULL)
    
    participant_data <- parquet_files %>%
      lapply(arrow::read_parquet) %>%
      bind_rows() %>%
      left_join(metadata, by = "submission_id") %>%
      distinct() %>%
      mutate(
        platform = toupper(platform),
        localtime = lubridate::as_datetime(localtime, tz = "UTC")
      ) %>%
      group_by(submission_id) %>%
      arrange(localtime, .by_group = TRUE) %>%
      mutate(
        # Create best available timestamp early
        localtime_and_timestamp = dplyr::if_else(
          platform == "ANDROID",
          localtime[1] + lubridate::dseconds(row_number() - 1L),
          localtime + lubridate::dseconds(as.numeric(time_offset))
        ),
        # Convert Android acceleration from m/s^2 to g.
        across(
          c(x_axis, y_axis, z_axis),
          ~ dplyr::if_else(platform == "ANDROID", .x / 9.81, .x)
        )
      ) %>%
      ungroup()
    
    message("Loaded parquet files for folder: ", participant_folder)
    
    participant_data <- participant_data %>%
      dplyr::filter(
        between(x_axis, -20, 20),
        between(y_axis, -20, 20),
        between(z_axis, -20, 20)
      ) %>%
      left_join(surveys$survey_t0, by = "participant_id") %>%
      left_join(surveys$survey_t1, by = "participant_id") %>%
      dplyr::filter(
        is.na(last_recorded_at_t1) |
          dplyr::between(
            localtime_and_timestamp,
            last_recorded_at_t1 - lubridate::days(7),
            last_recorded_at_t1
          )
      ) %>%
      as.data.table()
    
    # --- Detect sampling rate & (if >1 Hz) downsample to 1 Hz ----------------------
    setorder(participant_data, localtime_and_timestamp)
    
    # inter-sample gaps (s)
    participant_data[, dt_diff := as.numeric(difftime(
      localtime_and_timestamp,
      data.table::shift(localtime_and_timestamp),
      units = "secs"
    ))]
    
    # guard against non-positive/NA gaps
    dt_pos <- participant_data$dt_diff[is.finite(participant_data$dt_diff) & participant_data$dt_diff > 0]
    orig_sr_hz <- if (length(dt_pos) > 0) 1 / median(dt_pos, na.rm = TRUE) else NA_real_
    participant_data[, original_sampling_rate_hz := orig_sr_hz]
    
    # default: did not downsample (<= 1 Hz or unknown)
    down_status <- "retained_original_<=1Hz"
    
    if (is.finite(orig_sr_hz) && orig_sr_hz > 1) {
      # collapse to one sample per whole second (keep first observation within each second)
      tz_ <- attr(participant_data$localtime_and_timestamp, "tzone")
      participant_data[, second_ts := as.POSIXct(
        floor(as.numeric(localtime_and_timestamp)),
        origin = "1970-01-01",
        tz = ifelse(is.null(tz_), "UTC", tz_)
      )]
      
      # make sure rows are ordered within each second, then take the first
      data.table::setorder(participant_data, second_ts, localtime_and_timestamp)
      
      participant_data <- participant_data[
        ,
        .(
          localtime_and_timestamp = first(second_ts),
          x_axis       = first(x_axis),
          y_axis       = first(y_axis),
          z_axis       = first(z_axis),
          participant_id = first(participant_id),
          submission_id  = first(submission_id),
          platform       = first(platform),
          time_zone      = first(time_zone)
        ),
        by = second_ts
      ]
      
      participant_data[, second_ts := NULL]
      down_status <- "downsampled_1Hz"  # label back to old behaviour
    }
    
    # recompute post-downsampling coverage over the inclusive span of seconds
    participant_data[, ts_sec := as.integer(floor(as.numeric(localtime_and_timestamp)))]
    n_secs   <- data.table::uniqueN(participant_data$ts_sec)
    sec_span <- (max(participant_data$ts_sec) - min(participant_data$ts_sec) + 1L)
    coverage <- if (is.finite(sec_span) && sec_span > 0L) n_secs / sec_span else NA_real_
    
    participant_data[, `:=`(
      downsampled_sampling_rate_hz = if (is.finite(orig_sr_hz) && orig_sr_hz > 1) 1.0 else orig_sr_hz,
      coverage_1hz                 = coverage,
      downsample_status            = down_status
    )]
    participant_data[, ts_sec := NULL]
    
    # Check again after filtering and downsampling
    message("Filtered and processed participant data includes ", nrow(participant_data), " rows.")
    if (nrow(participant_data) == 0) return(NULL)
    
    # per-second index (used for windowing/completion)
    participant_data[, time_from_first_s := round(as.numeric(
      difftime(localtime_and_timestamp, first(localtime_and_timestamp), units = "secs")
    ))]
    
    # mark seconds that actually come from sensor data
    participant_data[, has_raw := 1L]
    
    # define origin_time AFTER (down)sampling
    origin_time <- min(participant_data$localtime_and_timestamp, na.rm = TRUE)
    
    participant_data_completed <- participant_data %>%
      tidyr::complete(time_from_first_s = seq(
        min(time_from_first_s), max(time_from_first_s), by = 1
      )) %>%
      mutate(
        localtime_and_timestamp = origin_time + lubridate::seconds(time_from_first_s),
        has_raw = ifelse(is.na(has_raw), 0L, has_raw)   # 0 = created by complete()
      )
    
    message("applying filters...")
    # Build dt and mark validity
    dt <- as.data.table(participant_data_completed)
    dt[, valid := !is.na(x_axis) & !is.na(y_axis) & !is.na(z_axis)]
    
    # Butterworth low-pass for gravity estimate
    fs <- 1
    cutoff <- 0.3
    wn <- cutoff / (fs / 2)
    bf <- signal::butter(n = 4, W = wn, type = "low")
    
    # protect against short valid runs
    min_run_len <- 10  # seconds; adjust as needed
    
    dt[, run := data.table::rleid(valid)]
    dt[, run_len := .N, by = run]
    
    # Only filter sufficiently long valid runs
    dt[
      valid == TRUE & run_len >= min_run_len,
      c("x_filt", "y_filt", "z_filt") := .(
        signal::filtfilt(bf, x_axis),
        signal::filtfilt(bf, y_axis),
        signal::filtfilt(bf, z_axis)
      ),
      by = run
    ]
    
    # Short runs remain NA in x_filt/y_filt/z_filt
    
    setorder(dt, localtime_and_timestamp)
    
    dt[, `:=`(
      magnitude_raw = sqrt(x_axis^2 + y_axis^2 + z_axis^2),
      jerk_raw = sqrt((x_axis - data.table::shift(x_axis))^2 +
                        (y_axis - data.table::shift(y_axis))^2 +
                        (z_axis - data.table::shift(z_axis))^2),
      magnitude_filt = sqrt(x_filt^2 + y_filt^2 + z_filt^2),
      jerk = sqrt((x_filt - data.table::shift(x_filt))^2 +
                    (y_filt - data.table::shift(y_filt))^2 +
                    (z_filt - data.table::shift(z_filt))^2),
      hour = lubridate::hour(localtime_and_timestamp)
    )]
    
    # Assign gx/gy/gz 
    dt[, `:=`(
      gx = x_filt,
      gy = y_filt,
      gz = z_filt
    )]
    
    # Assign dx/dy/dz (these rely on gx/gy/gz)
    dt[, `:=`(
      dx = x_axis - gx,
      dy = y_axis - gy,
      dz = z_axis - gz
    )]
    
    # Assign columns that depend on dx/dy/dz
    dt[, `:=`(
      dyn_mag = sqrt(dx^2 + dy^2 + dz^2),
      sma     = abs(dx) + abs(dy) + abs(dz)
    )]
    
    # ENMO from raw data, only when valid
    dt[, enmo := ifelse(
      valid,
      pmax(0, sqrt(x_axis^2 + y_axis^2 + z_axis^2) - 1),
      NA_real_
    )]
    
    # Gravity unit vector 
    dt[, g_norm := sqrt(gx^2 + gy^2 + gz^2)]
    dt[, `:=`(
      ghx = data.table::fifelse(g_norm > 0, gx/g_norm, NA_real_),
      ghy = data.table::fifelse(g_norm > 0, gy/g_norm, NA_real_),
      ghz = data.table::fifelse(g_norm > 0, gz/g_norm, NA_real_)
    )]
    
    # Vertical/horizontal split + dynamic jerk
    dt[, `:=`(
      dyn_vert = dx*ghx + dy*ghy + dz*ghz,
      jerk_dyn = sqrt((dx - data.table::shift(dx))^2 +
                        (dy - data.table::shift(dy))^2 +
                        (dz - data.table::shift(dz))^2)
    )]
    
    dt[, dyn_horz := sqrt(pmax(0, dyn_mag^2 - dyn_vert^2))]
    
    # 5-min windows based on time_from_first_s
    windowed_data <- dt %>%
      mutate(window_id = floor(time_from_first_s / 300)) %>%
      as.data.table()
    
    windowed_data[, `:=`(
      is_daytime = hour >= 7 & hour <= 22,
      is_moving  = dyn_mag > 0.05
    )]
    
    # ---- 7-day window accounting (valid / invalid / no-data) -------------
    # 7-day anchor: last_recorded_at_t1 defines week_end
    week_end <- unique(participant_data$last_recorded_at_t1)
    if (length(week_end) != 1 || is.na(week_end)) {
      week_end <- max(participant_data$localtime_and_timestamp, na.rm = TRUE)
    }
    week_start <- week_end - lubridate::days(7)
    
    dt[, t_from_week_start := as.numeric(
      difftime(localtime_and_timestamp, week_start, units = "secs")
    )]
    
    dt_7d <- dt[t_from_week_start >= 0 & t_from_week_start < 7 * 24 * 60 * 60]
    
    dt_7d[, window_id_7d := floor(t_from_week_start / 300)]
    
    win_all_7d <- dt_7d[, .(
      n_total = .N,
      n_valid = sum(valid),
      n_raw   = sum(has_raw)
    ), by = .(participant_id, window_id_7d)]
    
    min_valid_prop <- 0.5
    min_valid_secs <- 150
    
    win_all_7d[, window_type_7d := data.table::fifelse(
      n_raw == 0,
      "no_data",
      data.table::fifelse(
        n_valid >= min_valid_secs & (n_valid / n_total) >= min_valid_prop,
        "valid",
        "invalid"
      )
    )]
    
    n_week_windows <- ceiling(7 * 24 * 60 / 5)  # 2016
    
    n_windows_with_rows <- win_all_7d[, data.table::uniqueN(window_id_7d)]
    n_windows_no_row    <- n_week_windows - n_windows_with_rows
    
    per_participant_windows_7d <- win_all_7d[, .(
      n_windows_valid_7d   = sum(window_type_7d == "valid"),
      n_windows_invalid_7d = sum(window_type_7d == "invalid"),
      n_windows_no_data_7d = sum(window_type_7d == "no_data") + n_windows_no_row
    ), by = participant_id]
    
    per_participant_windows_7d[, n_windows_total_7d :=
                                 n_windows_valid_7d + n_windows_invalid_7d + n_windows_no_data_7d
    ]
    # ----------------------------------------------------------------------
    
    # Ensure window order and apply same validity rule (≥50% valid in a window)
    setorder(windowed_data, window_id, localtime_and_timestamp)
    
    win_stats <- windowed_data[, .(
      n_total = .N,
      n_valid = sum(valid),
      dyn_mag_mean = mean(dyn_mag, na.rm = TRUE)
    ), by = window_id][ n_valid / n_total >= 0.5 ]
    
    # Flag low-movement windows; no NAs
    win_stats[, low_win := !is.na(dyn_mag_mean) & dyn_mag_mean < 0.02]
    
    # Participant-level summaries of low windows
    if (nrow(win_stats) == 0) {
      prop_low_windows   <- NA_real_
      longest_low_streak <- 0
    } else {
      prop_low_windows <- mean(win_stats$low_win)
      r <- rle(win_stats$low_win)
      longest_low_streak <- if (length(r$lengths) && any(r$values)) {
        max(r$lengths[r$values])
      } else {
        0
      }
    }
    
    # Features to summarise
    cols_to_summarise <- c(
      "magnitude_raw", "magnitude_filt",
      "dyn_mag", "sma", "enmo",
      "dyn_vert", "dyn_horz",
      "x_axis", "y_axis", "z_axis",
      "jerk_raw", "jerk", "jerk_dyn"
    )
    
    # function to calculate features per subset
    calc_features <- function(data, condition_prefix = "") {
      feats <- data[
        , {
          stats <- unlist(lapply(cols_to_summarise, function(col) {
            v <- get(col)
            s <- time_domain_summary_base(v)
            setNames(s, paste0(condition_prefix, col, "_", names(s)))
          }), recursive = FALSE)
          c(list(n_total = .N, n_valid = sum(valid)), stats)
        }, by = window_id][n_valid / n_total >= 0.5]
      
      feats %>%
        summarise(
          n_windows     = n(),
          total_n_valid = sum(n_valid),
          across(
            where(is.numeric) & !any_of(c("window_id","n_total","n_valid")),
            ~ mean(.x, na.rm = TRUE)
          )
        ) %>%
        mutate(participant_id = unique(participant_data$participant_id))
    }
    
    # Calculate all conditions clearly
    summary_all     <- calc_features(windowed_data, "")
    summary_moving  <- calc_features(windowed_data[is_moving == TRUE], "moving_")
    summary_daytime <- calc_features(windowed_data[is_daytime == TRUE], "daytime_")
    
    summary_all     <- summary_all     %>% dplyr::rename_with(~ paste0("all_", .x),     -participant_id)
    summary_moving  <- summary_moving  %>% dplyr::rename_with(~ paste0("moving_", .x),  -participant_id)
    summary_daytime <- summary_daytime %>% dplyr::rename_with(~ paste0("daytime_", .x), -participant_id)
    
    # Combine summaries by participant_id
    participant_summary <- purrr::reduce(
      list(summary_all, summary_moving, summary_daytime),
      dplyr::left_join,
      by = "participant_id"
    )
    
    participant_summary <- participant_summary %>%
      dplyr::left_join(per_participant_windows_7d, by = "participant_id") %>%
      mutate(
        prop_low_windows            = prop_low_windows,
        longest_low_streak_windows  = as.numeric(longest_low_streak),
        original_sampling_rate_hz   = unique(participant_data$original_sampling_rate_hz),
        downsampled_sampling_rate_hz= unique(participant_data$downsampled_sampling_rate_hz),
        coverage_1hz_with_missing   = unique(participant_data$coverage_1hz),
        downsample_status           = unique(participant_data$downsample_status),
        platform                    = unique(participant_data$platform)
      )
    
    finish.time <- Sys.time()
    message(finish.time - start.time)
    participant_summary
  }
  
  # Main execution function
  run_pipeline <- function() {
    load_packages()
    surveys <- prepare_surveys("baseline", "mid")
    metadata <- read_metadata()
    participant_folders <- list.dirs(config$trial_folder, recursive = FALSE)
    
    all_results <- purrr::map_dfr(
      participant_folders,
      process_file,
      metadata = metadata,
      surveys = surveys
    )
    
    all_results
  }
  
  # Execute and time the pipeline
  start.time <- Sys.time()
  all_accel_feats <- run_pipeline()
  finish.time <- Sys.time()
  message("Total time: ", finish.time - start.time)
  
  # Write extracted features
  dir.create(dirname(config$output_file), recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    all_accel_feats,
    file = config$output_file,
    row.names = FALSE
  )
