################################################################
## The following function converts a list with x and y into a ##
## data.frame                                                 ##
##   parameters:                                              ##
##      xy:    a list with elements x and y                   ##
##   returns:  a data.frame with fields x and y               ##
################################################################
xy2df <- function(xy) {
  data.frame(x = xy$x, y = xy$y)
}

CoralTrends_minimap_base_ <- function(qld_sf) {
  base <- ggplot() +
    geom_blank(aes(x = 270, y = -20)) +
    geom_sf(data = qld_sf, fill = NA, colour = "black") +
    coord_sf() +
    theme_classic() +
    theme(panel.background = element_rect(fill = NA),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.background = element_blank(),
      panel.spacing = unit(0, "pt"),
      plot.margin = unit(c(0, 0, 0, 0), "pt")
    )
  return(base)
}

CoralTrends_minimap_ <- function(region = "GBR", gbr_3Zone, qld_sf, base) {
  gbr_3Zone_all <- gbr_3Zone |>
    st_union()
  ## basemap
  i <- case_when(
    region == "GBR" ~ 1:3,
    region == "Northern GBR" ~ 1,
    region == "Central GBR" ~ 2,
    region == "Southern GBR" ~ 3
  )
  p <- base +
    geom_sf(data = gbr_3Zone[i,], fill = hues[4], color = NA) +
    geom_sf(data = gbr_3Zone_all, fill = NA, color = "black", size = 0.2) +
    geom_sf(data = qld_sf, fill = "white", color = hues[4])
  return(p)
}

CoralTrends_generate_single_banner <- function(region = "GBR") {
  gbr_3Zone <- readRDS(file = paste0(spatial_path, "gbr_3Zone.rds"))
  aus <- oz:::ozRegion()
  qld <- rbind(
    xy2df(aus$lines[[3]]),
    xy2df(aus$lines[[13]]),
    xy2df(aus$lines[[12]])[nrow(xy2df(aus$lines[[12]])):1,],
    xy2df(aus$lines[[11]]))

  qld_sf <- sf::st_polygon(list(as.matrix(qld))) |>
    st_sfc(crs = st_crs(gbr_3Zone))

  base <- CoralTrends_minimap_base_(qld_sf)
  minimap <- CoralTrends_minimap_(region = region, gbr_3Zone, qld_sf, base)
  ## data <- data |>
  ##   mutate(minimap = list(minimap))
  ## return(data)
  return(minimap)
}

###################################################################################
## The plots will have tick marks along the temporal (x) axis.                   ##
## For the forseable future, these should represent every five years.            ##
## The starting year will always be 1985, but the final year will keep evolving. ##
## Lets take the variable final_year (defined in parameters/CoralTrends.conf     ##
## and round it to the nearest 5 and use that.                                   ##
###################################################################################
mceiling <- function(x,base){
        base*ceiling(x/base)
}

corporate_theme <-
  theme(
    strip.background = element_rect(
      fill = "white",
      color = "white",
      size = 0.5
    ),
    ## panel.background = element_rect(color = "black"),
    plot.margin = margin(t = 2, r = 7, b = 0, l = 0),
    axis.title.y = element_text(size = rel(1.5),
      margin = margin(r = 0.75, unit = "lines")),
    axis.text.x = element_text(size = rel(1.5)),
    axis.text.y = element_text(size = rel(1.5)),
    panel.grid.minor = element_line(size = 0.1, color = "gray70"),
    panel.grid.major = element_line(size = 0.1, color = "gray70"),
    panel.grid.minor.x = element_line(size = 0, color = "white", linetype = NULL),
    panel.grid.major.x = element_line(size = 0, color = "white", linetype = NULL),
    strip.text = element_text(margin = margin(t = 4, b = 1, r = 5, unit = "lines"),
      size = 0, lineheight = 0, colour = "white", face = "bold",
      family = "Arial",
      hjust = 0.50, vjust = -1))

hues <- RColorBrewer::brewer.pal(4, "Blues")

minimap_theme <-
  theme(
    strip.background = element_rect(
      fill = hues[2],
      color = "black",
      size = 0.5
    ),
    panel.background = element_rect(color = "black"),
    plot.margin = margin(t = 2, r = 7, b = 0, l = 0),
    axis.title.y = element_text(size = rel(1.5),
      margin = margin(r = 0.75, unit = "lines")),
    axis.text.x = element_text(size = rel(1.5)),
    axis.text.y = element_text(size = rel(1.5)),
    panel.grid.minor = element_line(size = 0.1, color = "gray70"),
    panel.grid.major = element_line(size = 0.1, color = "gray70"),
    panel.grid.minor.x = element_line(size = 0, color = "white", linetype = NULL),
    panel.grid.major.x = element_line(size = 0, color = "white", linetype = NULL),
    strip.text = element_text(margin = margin(t = 3, b = 4, r = 4, unit = "lines"),
      size = 20, lineheight = 0.5,
      colour = "black", face = "bold",
      family = "Arial",
      hjust = 0.3, vjust = -1))

CoralTrends_generate_single_trend <- function(data, region = "Northern GBR", final_year = NA) {
  dat <- readRDS(data$posteriors[[1]]$year_sum)
  if (is.na(final_year)) {
    final_year <- dat |>
      mutate(Year = as.numeric(as.character(REPORT_YEAR))) |>
      pull(Year) |>
      max()
  }
  final_year_seq <- mceiling(final_year,5)
  headings_lookup <- CoralTrends_get_headings_lookup__()

  uncertainty_lookup <- tribble(
     ~ribbon, ~ribbon_str, ~description,
     TRUE,    "ribbon",    "Ribbon uncertainty",
     FALSE,   "bars",      "Bars uncertainty"
  )
  n_lookup <- tribble(
    ~with_n, ~n_str,      ~description,
    TRUE,    "with n",    "with n",
    FALSE,   "without n", "without n"
  )
  axis_lookup <- tribble(
    ~fixed_axis, ~axis_str,      ~description,
    TRUE,        "fixed y",      "same y axis range",
    FALSE,       "free y",       "individual y axis range"
  )
  theme_lookup <- tribble(
    ~theme_type, ~theme_str,        ~description,
    "corporate", "corporate theme", "corporate theme",
    "minimap",   "minimap theme",   "minimap theme"
  )

  trend_plots <- uncertainty_lookup |>
    crossing(n_lookup, .name_repair = "universal") |>
    crossing(axis_lookup, .name_repair = "universal") |>
    crossing(theme_lookup, .name_repair = "universal") |>
    unite("description", contains("description"), sep = ", ") |>
    unite("file_str", ribbon_str, n_str, axis_str, theme_str, sep = "_") |>
    suppressMessages()

  trend_plots <- trend_plots |>
    mutate(plot = pmap(.l = list(file_str, ribbon, with_n, fixed_axis, theme_type),
      .f = ~ {
        file_str <- ..1
        ribbon <- ..2
        with_n <- ..3
        fixed_axis <- ..4
        theme_type <- ..5
        p <-
          dat |>
          left_join(headings_lookup, by = "Region") |>
          ggplot(aes(median, x = as.numeric(as.character(REPORT_YEAR))))
        if (ribbon) {
          p <- p +
            geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#9ccbed") +
            geom_point(colour = "#004785")
        } else {
          p <- p +
            geom_pointrange(aes(ymin = lower, ymax = upper), colour = "#004785")
        }
        p <- p +
          geom_line(color = "#004785") +
          scale_x_continuous('',
            breaks = seq(1985, final_year_seq, by = 5),
            limit = c(1985, final_year)) +
          theme_classic(base_family = "Arial", base_size = 12)
           ## theme_classic()
        if (with_n) {
          n_year <- readRDS(data$data_group[[1]]) |>
            dplyr::select(REPORT_YEAR, AIMS_REEF_NAME) |>
            dplyr::distinct() |>
            dplyr::mutate(REPORT_YEAR = factor(REPORT_YEAR)) |>
            group_by(REPORT_YEAR) |>
            count()
          dat <- dat |>
            left_join(n_year, by = "REPORT_YEAR")
          p <- p +
            ## geom_text(data = dat, aes(y = Inf, label = n), vjust = 1, size = rel(0.5))
            geom_text(data = dat,
              aes(y = upper, label = n), vjust = -0.3, size = rel(0.5))
        }
        if (fixed_axis) {
          p <- p +
            scale_y_continuous(expression(Coral~cover~('%')),
              labels = function(x) x*100,
              limits = c(0,0.49),
              expand = c(0,0))
        } else {
          p <- p +
            scale_y_continuous(expression(Coral~cover~("%")),
              labels = function(x) x*100)
        }
        if (theme_type == "corporate") {
          p <- p +
            facet_wrap(~Heading, nrow=1, scales='fixed',
              labeller=label_bquote(rows = "")) +
            corporate_theme
        } else {
          p <- p +
            ## facet_wrap(~Heading, nrow=1, scales='fixed',
            ##   labeller = labeller(Heading = setNames(paste0("\n", unique(dat$Heading),"\n"), unique(dat$Heading)))) +
            facet_wrap(~Heading, nrow=1, scales='fixed') +
            minimap_theme
        }
        ## p <- p + switch(theme_type,
        ##   "corporate" = corporate_theme,
        ##   "minimap" = minimap_theme
        ## )
        p
      }
    ))

  data <- data |> mutate(trend_plot = list(trend_plots))
  return(data)
}

## The following function defines a lookup that maps Region names to
## alternative names that are used in various graphs
CoralTrends_get_headings_lookup__ <- function() {
  headings_lookup <- tribble(
    ~Region,        ~Heading,                      ~Subheading,
    "Northern GBR", "Northern Great Barrier Reef", "Cape York to Cooktown",
    "Central GBR",  "Central Great Barrier Reef",  "Cooktown to Proserpine",
    "Southern GBR", "Southern Great Barrier Reef", "Proserpine to Gladstone",
    "GBR", "Great Barrier Reef", "Cape York to Gladstone",
    )
  return(headings_lookup)
}
## The following function defines a lookup that maps Region names to
## alternative names that are used in various graphs
CoralTrends_get_fonts_lookup__ <- function() {
  fonts_lookup <- tribble(
    ~Heading,        ~Font,   ~FontSize, ~FontColour,
    "Level 1",       "Arial", 30,        "#004785",
    "Level 2",       "Arial", 20,        "444040"
  )
  return(fonts_lookup)
}

CoralTrends_make_strip_text_grob__ <- function(region) {
  headings_lookup <- CoralTrends_get_headings_lookup__()
  strip_text_grob <- textGrob(  #tG
    headings_lookup |>
      filter(Region == region) |>
      pull(Heading),
    x = unit(0, "npc"),
    y = unit(1, "npc"),
    vjust = 0,
    hjust = 0,
    gp = gpar(fontsize = 30,
      fontfamily = "Arial",
      fontface = "bold",
      col = "#004785")
  )
  return(strip_text_grob)
}
CoralTrends_make_strip_subtext_grob__ <- function(region) {
  headings_lookup <- CoralTrends_get_headings_lookup__()
  strip_subtext_grob <- textGrob(
    headings_lookup |>
      filter(Region == region) |>
      pull(Subheading),
    x = unit(0, "npc"),
    y = unit(1, "npc"),
    vjust = 0,
    hjust = 0,
    gp = gpar(fontsize = 20,
      fontfamily = "Arial",
      fontface = "bold",
      col = "#989898")
  )
  return(strip_subtext_grob)
}
CoralTrends_logo_to_raster_grob__ <- function(aims_logo) {
  rG <- rasterGrob(aims_logo,
    x = unit(0, "npc"),
    y = unit(0, "npc"),
    hjust = 0,
    vjust = 0,
    width = unit(1, "npc"))
  return(rG)
}

## The following function is the workhorse that blends the single trend
## with the required addornment
CoralTrends_addorn_trend_plot__ <- function(trend_plot, fig_label, region) {
  png(file = tempfile(), width = 1, height = 1, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  ## need to temporarily move the cwd since ggplot_gtable insists on
  ## generating a Rplots.pdf and the current user does not have write access
  ## to this location
  tem <- getwd()
  setwd(tempdir())
  g <- ggplot_gtable(ggplot_build(trend_plot$plot[[1]]))
  setwd(tem)

  facets <- grep("strip-t-1-1", g$layout$name)

  if (trend_plot$theme_type[[1]] == "corporate") {
    strip_text_grob <- CoralTrends_make_strip_text_grob__(region)
    strip_subtext_grob <- CoralTrends_make_strip_subtext_grob__(region)
    ## Incorporate the strip text into the strip of the plot
    gg <- with(
      g$layout[facets,],
      gtable_add_grob(g,
        arrangeGrob(strip_text_grob, padding = unit(0, "line"),
          vp = viewport(x = 1, y = -0.4, width = 1,
            just = c("right", "bottom"))),
        t=t, l=l, b=b, r=l, name="pic_predator1")
    )
    ## Incorporate the strip subtext into the strip of the plot
    gg <- with(
      g$layout[facets,],
      gtable_add_grob(gg,
        arrangeGrob(strip_subtext_grob,
          padding = unit(0, "line"),
          vp = viewport(x = 1, y = -0.8,
            width = 1,
            just = c("right", "bottom"))),
        t=t, l=l, b=b, r=l,
        name="pic_predator2")
    )
    ## print(getwd())
    aims_logo <- CoralTrends_get_aims_logo(wch = "gov")
    aims_logo_rg <- CoralTrends_logo_to_raster_grob__(aims_logo)
    gg_with_logo <- with(
      g$layout[facets,],
      gtable_add_grob(
        gg,
        arrangeGrob(aims_logo_rg,
          padding = unit(0, "line"),
          vp = viewport(x = 1, y = 0.22,
            width = 0.3,
            just = c("right", "bottom"))),
        t=t, l=l, b=b, r=l,
        name="pic_predator")
    )
    plot <- gg_with_logo
  } else {
    minimap <- CoralTrends_generate_single_banner(region = domain)
    gg <- with(
      g$layout[facets,],
      gtable_add_grob(g,
        ggplotGrob(minimap),
        t=t, l=4, b=b, r=6, name="pic_predator")
    )
    aims_logo <- CoralTrends_get_aims_logo(wch = "AIMS")
    aims_logo_rg <- CoralTrends_logo_to_raster_grob__(aims_logo)
    gg_with_logo <- with(
      g$layout[facets,],
      gtable_add_grob(
        gg,
        arrangeGrob(aims_logo_rg,
          padding = unit(0, "line"),
          vp = viewport(x = 1, y = 0.15, #y = 0.22,
            width = 0.3,
            just = c("right", "bottom"))),
        t=t, l=l, b=b, r=l,
        name="pic_predator")
    )
    plot <- gg_with_logo
  }

  ## fname <- paste0(fig_label, trend_plot$file_str, ".png")
  fname <- paste0(fig_path, "gg_fig_", fig_label,
    trend_plot$file_str[[1]],
    ".png")

  ggsave(file = fname,
    plot,
    width = 801,
    height = 525,
    units = "px",
    dpi = 85
  )
  ggsave(file = gsub(".png", "hires.png", fname),
    plot,
    width = 801,
    height = 525,
    units = "px",
    dpi = 600
  )
  ggsave(file = gsub(".png", ".pdf", fname),
    plot,
    device = cairo_pdf,
    width = 801,
    height = 525,
    units = "px",
    dpi = 85
  )
  ## place marker for saving plot
  return(fname)
}

## The following function divies up the multiple rows to call out
## the get a single plot per row
CoralTrends_addorn_trend_plots__ <- function(trend_plot, fig_label, region) {
  plots <- seq_len(nrow(trend_plot)) |>
    map(function(i) {
     CoralTrends_addorn_trend_plot__(trend_plot[i, ], fig_label, region)
    })
  return(plots)
}

## In all likelihood, data here will only have a single row, so the
## map2 seems a little over the top.  That said, it means it can be
## expanded to accommodate multi-row datasets
CoralTrends_addorn_plots <- function(data, region) {
  data <- data |>
    mutate(trend_plot = map2(
      .x = trend_plot,
      .y = fig_label,
      .f =  ~ {
        trend_plot <- .x
        fig_label <- .y
        plots <- CoralTrends_addorn_trend_plots__(trend_plot, fig_label, region)
        return(trend_plot |> mutate(addorned_plot = plots))
      }))
  data <- data |>
    mutate(trend_plot = list(trend_plot))
  return(data)
}
