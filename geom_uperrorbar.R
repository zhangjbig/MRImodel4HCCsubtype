library(ggplot2)
library(magrittr)
#' @export
#' @rdname geom_linerange
geom_uperrorbar <- function(mapping = NULL, data = NULL,
                            stat = "identity", position = "identity",
                            ...,
                            na.rm = FALSE,
                            orientation = NA,
                            show.legend = NA,
                            inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomUperrorbar,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      orientation = orientation,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomUperrorbar <- ggproto("GeomUperrorbar", Geom,
                          default_aes = aes(colour = "black", size = 0.5, linetype = 1, width = 0.5,
                                            alpha = NA),
                          
                          draw_key = draw_key_path,
                          
                          required_aes = c("x|y", "ymin|xmin", "ymax|xmax"),
                          
                          setup_params = function(data, params) {
                            GeomLinerange$setup_params(data, params)
                          },
                          
                          extra_params = c("na.rm", "orientation"),
                          
                          setup_data = function(data, params) {
                            data$flipped_aes <- params$flipped_aes
                            data <- flip_data(data, params$flipped_aes)
                            data$width <- data$width %||%
                              params$width %||% (resolution(data$x, FALSE) * 0.9)
                            data <- transform(data,
                                              xmin = x - width / 2, xmax = x + width / 2, width = NULL
                            )
                            flip_data(data, params$flipped_aes)
                          },
                          
                          draw_panel = function(data, panel_params, coord, width = NULL, flipped_aes = FALSE) {
                            data <- flip_data(data, flipped_aes)
                            #x <- as.vector(rbind(data$xmin, data$xmax, NA, data$x,    data$x,    NA, data$xmin, data$xmax))
                            #y <- as.vector(rbind(data$ymax, data$ymax, NA, data$ymax, data$ymin, NA, data$ymin, data$ymin))
                            sel <- data$y < 0 
                            data$ymax[sel] <- data$ymin[sel]
                            x <- as.vector(rbind(data$xmin, data$xmax, NA, data$x,    data$x))
                            y <- as.vector(rbind(data$ymax, data$ymax, NA, data$ymax, data$y))
                            data <- new_data_frame(list(
                              x = x,
                              y = y,
                              colour = rep(data$colour, each = 5),
                              alpha = rep(data$alpha, each = 5),
                              size = rep(data$size, each = 5),
                              linetype = rep(data$linetype, each = 5),
                              group = rep(1:(nrow(data)), each = 5),
                              row.names = 1:(nrow(data) * 5)
                            ))
                            data <- flip_data(data, flipped_aes)
                            GeomPath$draw_panel(data, panel_params, coord)
                          }
)

new_data_frame <- function(x = list(), n = NULL) {
  if (length(x) != 0 && is.null(names(x))) {
    abort("Elements must be named")
  }
  lengths <- vapply(x, length, integer(1))
  if (is.null(n)) {
    n <- if (length(x) == 0 || min(lengths) == 0) 0 else max(lengths)
  }
  for (i in seq_along(x)) {
    if (lengths[i] == n) next
    if (lengths[i] != 1) {
      abort("Elements must equal the number of rows or 1")
    }
    x[[i]] <- rep(x[[i]], n)
  }
  
  class(x) <- "data.frame"
  
  attr(x, "row.names") <- .set_row_names(n)
  x
}