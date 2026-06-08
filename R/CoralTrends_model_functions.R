CoralTrends_prepare_tow_level_data <- function(data, manta_tow) {
  manta_tow <- manta_tow |>
    ## filter(Region == domain) |>
    droplevels() |>
    mutate(oLIVE_CORAL = factor(LIVE_CORAL,
      levels = c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
      ordered = TRUE),
      nREEF_NAME = as.numeric(as.factor(REEF_NAME)),
      nLIVE_CORAL = as.numeric(oLIVE_CORAL),
      REEF_YEAR = interaction(REEF_NAME, Year),
      fYEAR = REPORT_YEAR     #fYEAR needs to be there for downstream dashboard table reasons
    )
  data <- data |>
    mutate(data = map(.x = lab,
      .f = ~ {
        .x <- manta_tow
        fname <- paste0(data_path, "modelled/", lab,
          "_HC_beta_ _ _ _Cover_ _raw_data.rds")
        saveRDS(manta_tow,
          file = fname)
        return(fname)
      })) |>
  mutate(data_group = map(.x = data,
    .f = ~ {
      .x
    }))
  return(data)
}


CoralTrends_model_formula <- function(data) {
  form <- bf(Cover ~ Year + (1|REEF_NAME/REEF_YEAR),
    phi~0+Year,
    family = Beta(link = "logit")
  )
  data <- data |>
    mutate(form = list(form))
  return(data)
}

CoralTrends_model_priors <- function(data) {
  priors <- prior(normal(0, 3), class = "b") +
    prior(normal(0, 3), class = "Intercept") +
    prior(gamma(2, 1), class = "sd")
  data <- data |>
    mutate(priors = list(priors))
  return(data)
}

CoralTrends_fit_model <- function(data) {
  data <- data |>
    mutate(mod = pmap(.l = c(data_group, list(form), list(priors), label),
      .f = ~ {
        manta_tow <- readRDS(..1)
        form <- ..2
        ## environment(form) <- environment()  ## not needed for brms
        priors <- ..3
        label <- ..4
        fname <- paste0(data_path, "modelled/", label, ".rds")
        if (refit) {
          ## cat("form: ", form,  "\n\n\n", append = TRUE)
          ## cat("manta_tow", capture.output(print(manta_tow)), "\n", sep = "n", append=TRUE)
          ## cat("manta_tow", capture.output(nrow(manta_tow)), "\n", sep = "n")
          ## cat("manta_tow:", class(manta_tow), "\n\n\n", append = TRUE)
          mod <- brm(form,
            data = manta_tow,
            iter = 1e4,
            warmup = 5e3,
            thin = 5,
            chains = 4, cores = 4,
            control = list(adapt_delta = 0.95)
          )
          saveRDS(mod, file = fname)
        } else {
          Sys.setFileTime(fname, Sys.time())
          cat(paste0("\t - Skipped refitting (user specified)\n"),
            append = TRUE)
        }
        fname
      }
    ))
  return(data)
}

CoralTrends_prepare_for_contrasts <- function(data) {
  data <- data |>
    mutate(posteriors = map(.x = label,
      .f = ~ {
        list(year_group_posteriors = NA,
          year_group_sum = NA,
          year_posteriors = NA,
          year_sum = NA,
          yearcomp_posteriors = NA,
          yearcomp_sum = NA,
          all_yearcomp_posteriors = NA,
          all_yearcomp_sum = NA
          )
      }))
  data
}

CoralTrends_cellmeans_posteriors <- function(data) {
  data <- data |>
    mutate(posteriors = pmap(.l = list(mod, domain, label, posteriors),
      .f = ~ {
        mod <- readRDS(..1)
        domain <- ..2
        label <- ..3
        posteriors <- ..4
        mod_em_draws <- mod |>
          emmeans(~Year, type = "link") |>
          tidybayes::gather_emmeans_draws() |>
          mutate(.value = plogis(.value)) |>
          dplyr::rename(REPORT_YEAR = Year) |>
          mutate(Region = domain,
            GROUP = "HARD CORAL",
            FAMILY = "beta") |>
          dplyr::select(-.chain, -.iteration)
        fname <- paste0(data_path, "modelled/", label,
          "_year_posteriors.rds")
        saveRDS(mod_em_draws, file = fname)
        posteriors[["year_posteriors"]] <- fname
        posteriors
      }))
  return(data)
}

CoralTrends_cellmeans_sum <- function(data) {
  data <- data |>
    mutate(posteriors = pmap(.l = list(label, posteriors),
      .f = ~ {
        label <- ..1
        posteriors <- ..2
        mod_em_draws <- readRDS(posteriors$year_posteriors)
        mod_em_sum <- mod_em_draws |>
          group_by(REPORT_YEAR, Region, GROUP, FAMILY) |>
          posterior::summarise_draws(median, mean, HDInterval::hdi) |>
          dplyr::select(-variable) |>
          mutate(fYEAR = REPORT_YEAR, DATE = NA)   ## to fit in with the dashboard
        fname <- paste0(data_path, "modelled/", label,
          "_year_sum.rds")
        saveRDS(mod_em_sum, file = fname)
        posteriors[["year_sum"]] <- fname
        posteriors
      }
    ))
  return(data)
}

CoralTrends_contrast_posteriors <- function(data) {
  data <- data |>
    mutate(posteriors = pmap(.l = list(label, posteriors, domain),
      .f = ~ {
        label <- ..1
        posteriors <- ..2
        domain <- ..3
        mod_em_draws <- readRDS(posteriors$year_posteriors)
        report_year_levels <- mod_em_draws |> ungroup() |>
          filter(.draw == 1) |> pull(REPORT_YEAR) |> sort() |>
          rev() |> as.character()
        yearcomp_posteriors <- mod_em_draws |>
          ungroup() |>
          group_by(.draw) |>
          ## filter(.draw == 1) |>
          arrange(desc(REPORT_YEAR)) |>
          mutate(REPORT_YEAR = factor(REPORT_YEAR, levels = report_year_levels)) |>
          ## summarise(
          reframe(
            frac = exp(as.vector(as.vector(log(.value)) %*%
                                   t(cbind(1, -1 * model.matrix(~REPORT_YEAR)[-1, -1])))),
            value = as.vector(as.vector(.value) %*%
                                t(cbind(1, -1 * model.matrix(~REPORT_YEAR)[-1, -1]))),
            YearComp = paste0(first(REPORT_YEAR), "-", REPORT_YEAR[-1]),
            .groups = "keep"
          ) |>
          dplyr::select(YearComp, .draw, value, frac) |>
          arrange(desc(YearComp)) |>
          mutate(Region = domain,
            GROUP = "HARD CORAL",
            FAMILY = "beta")
        fname <- paste0(data_path, "modelled/", label,
          "_yearcomp_posteriors.rds")
        saveRDS(yearcomp_posteriors, file = fname)
        posteriors[["yearcomp_posteriors"]] <- fname
        posteriors
      }
    ))
  return(data)
}

CoralTrends_contrast_sum <- function(data) {
  data <- data |>
    mutate(posteriors = pmap(.l = list(label, posteriors),
      .f = ~ {
        label <- ..1
        posteriors <- ..2
        yearcomp_posteriors <- readRDS(posteriors$yearcomp_posteriors)
        yearcomp_sum <- yearcomp_posteriors |>
          group_by(YearComp, Region, GROUP, FAMILY) |>
          posterior::summarise_draws(median, mean, HDInterval::hdi,
            Pl = ~ mean(.x < 0), Pg = ~ mean(.x > 0)) |>
          ungroup() |>
          group_by(YearComp, Region, GROUP, FAMILY) |>
          mutate(across(c(Pl, Pg), ~ first(.x))) |>
          ungroup() |>
          arrange(desc(YearComp))
        fname <- paste0(data_path, "modelled/", label,
          "_yearcomp_sum.rds")
        saveRDS(yearcomp_sum, file = fname)
        posteriors[["yearcomp_sum"]] <- fname
        posteriors
      }
    ))
  return(data)
}

CoralTrends_all_contrast_posteriors <- function(data) {
  data <- data |>
    mutate(posteriors = pmap(.l = list(label, posteriors, domain),
      .f = ~ {
        label <- ..1
        posteriors <- ..2
        domain <- ..3
        mod_em_draws <- readRDS(posteriors$year_posteriors)

        years <- levels(mod_em_draws$REPORT_YEAR)
        xmat <- emmeans:::tukey.emmc(rev(years))
        yearcomp_posteriors <- mod_em_draws |>
          ungroup() |>
          group_by(.draw) |>
          arrange(desc(REPORT_YEAR)) |>     #must be used with tukey.emmc(rev(years))
          ## summarise(
          reframe(
            frac = exp(as.vector(as.vector(log(.value)) %*% as.matrix(xmat))),
            value = as.vector(as.vector(.value) %*% as.matrix(xmat)),
            YearComp = names(xmat),
            .groups = "keep"
          ) |>
          dplyr::select(YearComp, .draw, value, frac) |>
          mutate(Region = domain,
            GROUP = "HARD CORAL",
            FAMILY = "beta")
        fname <- paste0(data_path, "modelled/", label,
          "_all_yearcomp_posteriors.rds")
        saveRDS(yearcomp_posteriors, file = fname)
        posteriors[["all_yearcomp_posteriors"]] <- fname
        posteriors
      }
    ))
  return(data)
}

CoralTrends_all_contrast_sum <- function(data) {
  data <- data |>
    mutate(posteriors = pmap(.l = list(label, posteriors),
      .f = ~ {
        label <- ..1
        posteriors <- ..2
        yearcomp_posteriors <- readRDS(posteriors$all_yearcomp_posteriors)

        yearcomp_sum <- yearcomp_posteriors |>
          group_by(YearComp, Region, GROUP, FAMILY) |>
          posterior::summarise_draws(median, mean, HDInterval::hdi,
            Pl = ~ mean(.x < 0), Pg = ~ mean(.x > 0)) |>
          ungroup() |>
          group_by(YearComp, Region, GROUP, FAMILY) |>
          mutate(across(c(Pl, Pg), ~ first(.x))) |>
          ungroup() |>
          arrange(desc(YearComp))
        fname <- paste0(data_path, "modelled/", label,
          "_all_yearcomp_sum.rds")
        saveRDS(yearcomp_sum, file = fname)
        posteriors[["all_yearcomp_posteriors"]] <- fname
        posteriors
      }
    ))
  return(data)
}


CoralTrends_raw_summary <- function(data) {
  data <- data |>
    mutate(raw_sum = pmap(.l = list(data, label),
      .f = ~ {
        data <- readRDS(..1)
        label <- ..2
        data_sum <- data |>
          ungroup() |>
          group_by(AIMS_REEF_NAME, REEF_NAME, REPORT_YEAR,
            Zone) |>
          summarise(
            mean = mean(Cover, na.rm = TRUE),
            median = median(Cover, na.rm = TRUE),
            .groups = "drop"
          ) |>
          group_by(REPORT_YEAR, Zone) |>
          summarise(
            mean = mean(mean, na.rm = TRUE),
            median = median(median, na.rm = TRUE),
            .groups = "drop"
          ) |> mutate(fYEAR = REPORT_YEAR)  ## needed for dashboard compatibility
        fname <- paste0(data_path, "modelled/", label,
          "_raw_sums.rds")
        saveRDS(data_sum, file = fname)
        fname
      }
    ))
  return(data)
}

CoralTrends_raw_summary_plots <- function(data) {
  data <- data |>
    mutate(gg = pmap(.l = list(data_group, raw_sum, fig_label),
      .f = ~ {
        data_group <- readRDS(..1)
        raw_sum <- readRDS(..2)
        fig_label <- ..3
        dg <- data_group |>
          ungroup() |>
          group_by(REPORT_YEAR, AIMS_REEF_NAME) |>
          summarise(value = mean(Cover), .groups = "keep") |>
          ungroup()
        p <-
          dg |>
          ggplot(aes(x = REPORT_YEAR, y = value)) +
          geom_line(show.legend = FALSE, colour = "red") +
          geom_line(data = raw_sum,
            aes(y = mean,
              x = as.numeric(as.character(REPORT_YEAR)))) +
          geom_line(data = raw_sum,
            aes(y = median,
              x = as.numeric(as.character(REPORT_YEAR))),
            linetype = "dashed") +
          scale_y_continuous("Percent cover") +
          scale_x_continuous("") +
          facet_wrap(~AIMS_REEF_NAME) +
          theme_bw(base_size = 5)
        fname <- paste0(fig_path, "gg_raw_sum_", fig_label,
          ".png")
        ggsave(fname, plot = p, width = 12, height = 7)
        fname
      }
    ))
 return(data)
}

CoralTrends_partial_plots <- function(data) {
  data <- data |>
    mutate(gg = pmap(.l = list(posteriors, raw_sum, fig_label),
      .f = ~ {
        year_sum <- readRDS(..1$year_sum)
        raw_sum <- readRDS(..2)
        fig_label <- ..3
        family_type <- unique(year_sum$FAMILY)

        p <-
          year_sum |>
          ggplot(aes(x = as.numeric(as.character(REPORT_YEAR)), y = median)) +
          geom_ribbon(aes(ymin = lower, ymax = upper, x = as.numeric(as.character(REPORT_YEAR))),
            fill = "orange", alpha = 0.3) +
          geom_line(aes(group = family_type, colour = family_type)) +
          geom_point()  +
          geom_line(data = raw_sum,
            aes(y = mean, x = as.numeric(as.character(REPORT_YEAR)),
              colour = "Raw mean")) +
          geom_line(data = raw_sum,
            aes(y = median, x = as.numeric(as.character(REPORT_YEAR)),
              colour = "Raw median")) +
          scale_y_continuous("Percent cover") +
          scale_x_continuous("") +
          facet_wrap(~FAMILY, scales = "free_y") +
          theme_bw() +
          theme(axis.title.x = element_blank(),
            strip.background = element_rect(fill = "lightblue")
          )
        fname <- paste0(fig_path, "gg_", fig_label,
          ".png")
        ggsave(fname, plot = p, width = 12, height = 7)
        fname
      }
    ))
 return(data)
}

CoralTrends_longterm_posteriors <- function(data) {
  data <- data |>
    mutate(posteriors = pmap(.l = list(label, posteriors, domain),
      .f = ~ {
        label <- ..1
        posteriors <- ..2
        domain <- ..3
        mod_em_draws <- readRDS(posteriors$year_posteriors)
        report_year_levels <- mod_em_draws |> ungroup() |>
          filter(.draw == 1) |> pull(REPORT_YEAR) |> sort() |>
          rev() |> as.character()
        longterm.mat <- matrix(rep(1/length(report_year_levels),
          length(report_year_levels)),
          nrow = 1)
        years <- levels(mod_em_draws$REPORT_YEAR)
        xmat <- emmeans:::eff.emmc(rev(years))
        yearcomp_posteriors <- mod_em_draws |>
          ungroup() |>
          group_by(.draw) |>
          ## filter(.draw == 1) |>
          arrange(desc(REPORT_YEAR)) |>
          mutate(REPORT_YEAR = factor(REPORT_YEAR, levels = report_year_levels)) |>
          ## summarise(
          reframe(
            longterm = as.vector(.value %*% t(longterm.mat)),
            frac = exp(as.vector(as.vector(log(.value)) %*% as.matrix(xmat))),
            value = as.vector(as.vector(.value) %*% as.matrix(xmat)),
            YearComp = names(xmat),

            ## YearComp = paste0("Longterm", "-", REPORT_YEAR[-1]),
            .groups = "keep"
          ) |>
          dplyr::select(YearComp, .draw, value, frac, longterm) |>
          arrange(desc(YearComp)) |>
          mutate(Region = domain,
            GROUP = "HARD CORAL",
            FAMILY = "beta")
        lngtrm <- yearcomp_posteriors |>
          distinct(.draw, .keep_all = TRUE) |>
          mutate(YearComp = "Long-term", value = longterm) |>
          dplyr::select(YearComp, Region, GROUP, FAMILY, .draw, value)
        longterm_posteriors <- yearcomp_posteriors |>
          dplyr::select(YearComp, Region, GROUP, FAMILY, .draw, value, frac) |>
          bind_rows(lngtrm) |>
          arrange(.draw, desc(YearComp == "Long-term"))
        fname <- paste0(data_path, "modelled/", label,
          "_long-term_posteriors.rds")
        saveRDS(longterm_posteriors, file = fname)
        posteriors[["longterm_posteriors"]] <- fname
        posteriors
      }
    ))
  return(data)
}

CoralTrends_longterm_contrast_sum <- function(data) {
  data <- data |>
    mutate(posteriors = pmap(.l = list(label, posteriors),
      .f = ~ {
        label <- ..1
        posteriors <- ..2
        longterm_posteriors <- readRDS(posteriors$longterm_posteriors)

        longterm_sum <- longterm_posteriors |>
          group_by(YearComp, Region, GROUP, FAMILY) |>
          posterior::summarise_draws(median, mean, HDInterval::hdi,
            Pl = ~ mean(.x < 0), Pg = ~ mean(.x > 0)) |>
          ungroup() |>
          group_by(YearComp, Region, GROUP, FAMILY) |>
          mutate(across(c(Pl, Pg), ~ first(.x))) |>
          mutate(across(c(Pl, Pg), ~ ifelse(YearComp == "Long-term", NA, .x))) |>
          ungroup() |>
          arrange(desc(YearComp == "Long-term"), desc(YearComp)) |>
          filter(!(YearComp == "Long-term" & variable == "frac"))
        fname <- paste0(data_path, "modelled/", label,
          "_longterm_sum.rds")
        saveRDS(longterm_sum, file = fname)
        posteriors[["longterm_posteriors"]] <- fname
        posteriors
      }
    ))
  return(data)
}
