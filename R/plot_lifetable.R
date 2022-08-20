#' Plot a life table
#'
#' Create a log probability of death plot of the given life table, optionally
#' including multiple lines to indicate different input open age intervals.
#'
#' @param lt Life table to plot. Should have at least columns `loc_id`, `year`,
#'   `sex`, `age_start`, and `qx`.
#' @param open_age_var Optional column name indicating column with the input
#'   open age to use as the color aesthetic.
#'
#' @return a plot of the life table.
#' @export
#'
#' @examples
#' plot_lifetable(usa_2000_processed_lt, "input_open_age")
plot_lifetable <- function(lt, open_age_var = NULL) {

  gg <-
    ggplot(lt, aes(x = age_start, y = qx)) +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
    scale_y_log10(
      labels = scales::label_percent(drop0trailing = TRUE),
    ) +
    annotation_logticks(sides = "l") +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.minor.y = element_blank()
    ) +
    labs(
      title = paste(lt[1, loc_id], lt[1, year], lt[1, sex]),
      x = "Age Start",
      y = bquote(q[x])
    )

  if (!is.null(open_age_var)) {

    gg <- gg +
      geom_line(aes(color = factor(get(open_age_var)))) +
      labs(color = "Input open age")

  } else {

    gg <- gg + geom_line()

  }

  gg

}
