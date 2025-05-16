#' @title Custom Panel Function for Histograms in Lattice Plots
#' @description Utility functions used in lattice plots to show distribution overlays in bwplot comparisons.
#' @keywords internal
#' @import lattice
#' @export
panel.hanoi <- function(x, y, horizontal, breaks = "Sturges", ...) {
  if (horizontal) {
    condvar <- y
    datavar <- x
  } else {
    condvar <- x
    datavar <- y
  }

  conds <- sort(unique(condvar))

  for (i in seq_along(conds)) {
    h <- hist(datavar[condvar == conds[i]], plot = FALSE, breaks = breaks)
    brks.cnts <- stripOuterZeros(h$breaks, h$counts)
    brks <- brks.cnts[[1]]
    cnts <- brks.cnts[[2]]

    halfrelfs <- (cnts / sum(cnts)) / 2
    center <- i

    if (horizontal) {
      lattice::panel.rect(head(brks, -1), center - halfrelfs, tail(brks, -1), center + halfrelfs, ...)
    } else {
      lattice::panel.rect(center - halfrelfs, head(brks, -1), center + halfrelfs, tail(brks, -1), ...)
    }
  }
}

# Strip zeros on both sides of histogram bins
stripOuterZeros <- function(brks, cnts) {
  do.call("stripLeftZeros", stripRightZeros(brks, cnts))
}

stripLeftZeros <- function(brks, cnts) {
  if (cnts[1] == 0) {
    stripLeftZeros(brks[-1], cnts[-1])
  } else {
    list(brks, cnts)
  }
}

stripRightZeros <- function(brks, cnts) {
  len <- length(cnts)
  if (cnts[len] == 0) {
    stripRightZeros(brks[-(len + 1)], cnts[-len])
  } else {
    list(brks, cnts)
  }
}
