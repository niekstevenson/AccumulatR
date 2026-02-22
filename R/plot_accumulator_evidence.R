.pae_extract_tables <- function(model) {
  if (is_model_tables(model)) {
    return(model)
  }
  if (is.list(model) && !is.null(model$model_spec) && is.list(model$model_spec)) {
    model <- model$model_spec
  }
  model_to_tables(model)
}

.pae_extract_sources <- function(expr) {
  if (is.null(expr)) return(character(0))
  kind <- expr$kind %||% "event"
  if (identical(kind, "event")) {
    src <- expr$source %||% NA_character_
    return(if (!is.na(src) && nzchar(src)) src else character(0))
  }
  if (identical(kind, "guard")) {
    out <- c(
      .pae_extract_sources(expr$reference),
      .pae_extract_sources(expr$blocker)
    )
    unless_list <- expr$unless %||% list()
    if (length(unless_list) > 0L) {
      for (i in seq_along(unless_list)) {
        out <- c(out, .pae_extract_sources(unless_list[[i]]))
      }
    }
    return(unique(out))
  }
  if (identical(kind, "and") || identical(kind, "or")) {
    args <- expr$args %||% list()
    if (length(args) == 0L) return(character(0))
    out <- unlist(lapply(args, .pae_extract_sources), use.names = FALSE)
    return(unique(out))
  }
  if (identical(kind, "not")) {
    return(.pae_extract_sources(expr$arg))
  }
  character(0)
}

.pae_collect_and_source_sets <- function(expr) {
  sets <- list()
  walk <- function(node) {
    if (is.null(node)) return(invisible(NULL))
    kind <- node$kind %||% "event"
    if (identical(kind, "and")) {
      sets[[length(sets) + 1L]] <<- unique(.pae_extract_sources(node))
      args <- node$args %||% list()
      if (length(args) > 0L) lapply(args, walk)
      return(invisible(NULL))
    }
    if (identical(kind, "guard")) {
      walk(node$reference)
      walk(node$blocker)
      unless_list <- node$unless %||% list()
      if (length(unless_list) > 0L) lapply(unless_list, walk)
      return(invisible(NULL))
    }
    if (identical(kind, "or")) {
      args <- node$args %||% list()
      if (length(args) > 0L) lapply(args, walk)
      return(invisible(NULL))
    }
    if (identical(kind, "not")) {
      walk(node$arg)
      return(invisible(NULL))
    }
    invisible(NULL)
  }
  walk(expr)
  sets
}

.pae_collect_guard_source_links <- function(expr) {
  out <- list()
  walk <- function(node) {
    if (is.null(node)) return(invisible(NULL))
    kind <- node$kind %||% "event"
    if (identical(kind, "guard")) {
      refs <- unique(.pae_extract_sources(node$reference))
      blks <- unique(.pae_extract_sources(node$blocker))
      if (length(refs) > 0L && length(blks) > 0L) {
        out[[length(out) + 1L]] <<- list(blockers = blks, references = refs)
      }
      walk(node$reference)
      walk(node$blocker)
      unless_list <- node$unless %||% list()
      if (length(unless_list) > 0L) lapply(unless_list, walk)
      return(invisible(NULL))
    }
    if (identical(kind, "and") || identical(kind, "or")) {
      args <- node$args %||% list()
      if (length(args) > 0L) lapply(args, walk)
      return(invisible(NULL))
    }
    if (identical(kind, "not")) {
      walk(node$arg)
      return(invisible(NULL))
    }
    invisible(NULL)
  }
  walk(expr)
  out
}

.pae_pair_key <- function(a, b) {
  parts <- sort(c(a, b))
  paste(parts[[1]], parts[[2]], sep = "||")
}

.pae_directed_key <- function(from, to) {
  paste(from, to, sep = "->")
}

.pae_draw_curved_arrow <- function(x0, y0, x1, y1,
                                   col = "gray50", lty = 3, lwd = 1,
                                   curvature = 0.2, both = FALSE,
                                   length = 0.08, n = 80L) {
  dx <- x1 - x0
  dy <- y1 - y0
  dist <- sqrt(dx^2 + dy^2)
  if (!is.finite(dist) || dist <= .Machine$double.eps) return(invisible(NULL))

  mx <- (x0 + x1) / 2
  my <- (y0 + y1) / 2
  px <- -dy / dist
  py <- dx / dist
  offset <- curvature * dist
  cx <- mx + px * offset
  cy <- my + py * offset

  tt <- seq(0, 1, length.out = max(20L, as.integer(n)))
  bx <- (1 - tt)^2 * x0 + 2 * (1 - tt) * tt * cx + tt^2 * x1
  by <- (1 - tt)^2 * y0 + 2 * (1 - tt) * tt * cy + tt^2 * y1
  graphics::lines(bx, by, col = col, lty = lty, lwd = lwd)

  qbezier <- function(t) {
    c(
      (1 - t)^2 * x0 + 2 * (1 - t) * t * cx + t^2 * x1,
      (1 - t)^2 * y0 + 2 * (1 - t) * t * cy + t^2 * y1
    )
  }

  t_head <- 0.92
  e0 <- qbezier(t_head)
  e1 <- qbezier(1)
  if (both) {
    s0 <- qbezier(0.08)
    s1 <- qbezier(0)
    graphics::arrows(s0[[1]], s0[[2]], s1[[1]], s1[[2]],
      col = col, lty = lty, lwd = lwd, code = 2, length = length
    )
  }
  graphics::arrows(e0[[1]], e0[[2]], e1[[1]], e1[[2]],
    col = col, lty = lty, lwd = lwd, code = 2, length = length
  )
  invisible(NULL)
}

#' Plot accumulator trajectories with structural links
#'
#' Draws a model-wide accumulator schematic with time on the x-axis and evidence
#' on the y-axis. Each accumulator is drawn as an onset-aware diagonal trajectory.
#' Blocker links from guard relationships are overlaid as red dotted curves.
#'
#' @param model A race model (`race_spec`, finalized model structure, or `model_tables`).
#' @param xlim Optional x-axis limits. If `NULL`, computed from onsets and line length.
#' @param ylim Y-axis limits.
#' @param line_length Trajectory length in plot units.
#' @param angle_dodge Target angular spacing in degrees between neighbouring accumulators.
#'   Angles are centered at 45 degrees and constrained to the interval [10, 80].
#'   If there are too many accumulators sharing the same onset, spacing is reduced
#'   uniformly within that onset group.
#' @param accumulator_order Optional character vector of accumulator labels, ordered
#'   from top-to-bottom. This only affects accumulators sharing the same onset.
#' @param curve_strength Curvature multiplier for relationship arcs.
#' @param line_col Color for accumulator trajectories.
#' @param line_lwd Line width for accumulator trajectories.
#' @param axis_lwd Line width for x/y axis arrows.
#' @param link_lwd Reserved for reciprocal relationship arrows (currently disabled).
#' @param blocker_lwd Line width for blocker arrows.
#' @param show_labels If `TRUE`, draw accumulator labels near trajectory endpoints.
#' @param main Main title.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param ... Additional arguments passed to `graphics::plot`.
#'
#' @return Invisibly returns a list with plotting coordinates and inferred links.
#' @examples
#' spec <- race_spec() |>
#'   add_accumulator("go", "lognormal") |>
#'   add_accumulator("stop", "lognormal", onset = 0.15) |>
#'   add_pool("P", c("go", "stop")) |>
#'   add_outcome("R", "P") |>
#'   finalize_model()
#' plot_accumulators(spec)
#' @export
plot_accumulators <- function(model,
                              xlim = NULL,
                              ylim = c(0, 1),
                              line_length = 1,
                              angle_dodge = 10,
                              accumulator_order = NULL,
                              curve_strength = 0.18,
                              line_col = "gray20",
                              line_lwd = 3,
                              axis_lwd = 3,
                              link_lwd = 2.2,
                              blocker_lwd = 2.6,
                              show_labels = TRUE,
                              main = "Accumulator Schematic",
                              xlab = "Time",
                              ylab = "Evidence",
                              ...) {
  tables <- .pae_extract_tables(model)

  struct_nodes <- tables$struct_nodes
  acc_rows <- struct_nodes[struct_nodes$node_type == "accumulator", , drop = FALSE]
  pool_rows <- struct_nodes[struct_nodes$node_type == "pool", , drop = FALSE]
  outcome_rows <- struct_nodes[struct_nodes$node_type == "outcome", , drop = FALSE]
  if (nrow(acc_rows) == 0L) {
    stop("Model contains no accumulators", call. = FALSE)
  }

  acc_ids <- as.character(acc_rows$label)
  pool_members <- setNames(
    lapply(seq_len(nrow(pool_rows)), function(i) {
      payload <- pool_rows$payload[[i]] %||% list()
      as.character(payload$members %||% character(0))
    }),
    as.character(pool_rows$label)
  )

  expand_to_accumulators <- local({
    cache <- new.env(parent = emptyenv())
    rec <- function(source, stack = character(0)) {
      source <- as.character(source %||% "")
      if (!nzchar(source)) return(character(0))
      if (exists(source, envir = cache, inherits = FALSE)) {
        return(get(source, envir = cache, inherits = FALSE))
      }
      if (source %in% stack) return(character(0))
      if (source %in% acc_ids) {
        assign(source, source, envir = cache)
        return(source)
      }
      members <- pool_members[[source]] %||% character(0)
      if (length(members) == 0L) return(character(0))
      out <- unique(unlist(lapply(members, rec, stack = c(stack, source)), use.names = FALSE))
      assign(source, out, envir = cache)
      out
    }
    rec
  })

  onset_raw <- setNames(
    lapply(seq_len(nrow(acc_rows)), function(i) {
      payload <- acc_rows$payload[[i]] %||% list()
      payload$onset %||% 0
    }),
    acc_ids
  )

  onset_value <- local({
    cache <- new.env(parent = emptyenv())
    rec <- function(acc_id, stack = character(0)) {
      if (exists(acc_id, envir = cache, inherits = FALSE)) {
        return(get(acc_id, envir = cache, inherits = FALSE))
      }
      if (acc_id %in% stack) return(0)
      raw <- onset_raw[[acc_id]] %||% 0
      if (is.numeric(raw) && length(raw) > 0L && is.finite(raw[[1]])) {
        val <- as.numeric(raw[[1]])
        assign(acc_id, val, envir = cache)
        return(val)
      }
      if (is.list(raw) && identical(raw$kind %||% "", "after")) {
        src <- as.character(raw$source %||% "")
        lag <- suppressWarnings(as.numeric(raw$lag %||% 0)[1])
        if (!is.finite(lag)) lag <- 0
        src_accs <- expand_to_accumulators(src)
        if (length(src_accs) == 0L) {
          val <- lag
          assign(acc_id, val, envir = cache)
          return(val)
        }
        src_vals <- vapply(src_accs, rec, numeric(1), stack = c(stack, acc_id))
        val <- max(src_vals, na.rm = TRUE) + lag
        assign(acc_id, val, envir = cache)
        return(val)
      }
      assign(acc_id, 0, envir = cache)
      0
    }
    function(ids) setNames(vapply(ids, rec, numeric(1)), ids)
  })(acc_ids)

  ord <- order(onset_value, acc_ids)
  acc_ids <- acc_ids[ord]
  onset_value <- onset_value[acc_ids]

  n_acc <- length(acc_ids)
  yr <- diff(ylim)
  if (!is.finite(yr) || yr <= 0) stop("ylim must define a positive range", call. = FALSE)
  if (!is.numeric(line_length) || length(line_length) != 1L || !is.finite(line_length) || line_length <= 0) {
    stop("line_length must be a positive finite numeric scalar", call. = FALSE)
  }
  if (!is.numeric(angle_dodge) || length(angle_dodge) != 1L || !is.finite(angle_dodge) || angle_dodge < 0) {
    stop("angle_dodge must be a non-negative finite numeric scalar", call. = FALSE)
  }
  if (!is.null(accumulator_order)) {
    if (!is.character(accumulator_order)) {
      stop("accumulator_order must be NULL or a character vector of accumulator labels", call. = FALSE)
    }
    if (anyNA(accumulator_order) || any(!nzchar(accumulator_order))) {
      stop("accumulator_order cannot contain NA or empty labels", call. = FALSE)
    }
    if (anyDuplicated(accumulator_order)) {
      stop("accumulator_order cannot contain duplicate labels", call. = FALSE)
    }
    unknown <- setdiff(accumulator_order, acc_ids)
    if (length(unknown) > 0L) {
      stop(
        sprintf(
          "Unknown accumulator label(s) in accumulator_order: %s",
          paste(unknown, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  if (!is.numeric(axis_lwd) || length(axis_lwd) != 1L || !is.finite(axis_lwd) || axis_lwd <= 0) {
    stop("axis_lwd must be a positive finite numeric scalar", call. = FALSE)
  }
  if (!is.numeric(link_lwd) || length(link_lwd) != 1L || !is.finite(link_lwd) || link_lwd <= 0) {
    stop("link_lwd must be a positive finite numeric scalar", call. = FALSE)
  }
  if (!is.numeric(blocker_lwd) || length(blocker_lwd) != 1L || !is.finite(blocker_lwd) || blocker_lwd <= 0) {
    stop("blocker_lwd must be a positive finite numeric scalar", call. = FALSE)
  }

  mean_angle_deg <- 55
  angle_bounds_deg <- c(10, 80)
  angle_deg <- rep(mean_angle_deg, n_acc)
  rank_top_to_bottom <- setNames(seq_len(n_acc), acc_ids)
  if (!is.null(accumulator_order) && length(accumulator_order) > 0L) {
    remainder <- acc_ids[!acc_ids %in% accumulator_order]
    ordered_full <- c(accumulator_order, remainder)
    rank_top_to_bottom <- setNames(seq_along(ordered_full), ordered_full)
  }
  onset_group <- split(seq_len(n_acc), sprintf("%.8f", onset_value))
  for (grp in onset_group) {
    m <- length(grp)
    if (m <= 1L) {
      angle_deg[grp] <- mean_angle_deg
      next
    }
    max_step <- diff(angle_bounds_deg) / (m - 1L)
    step_deg <- min(angle_dodge, max_step)
    half_span <- step_deg * (m - 1L) / 2
    vals <- seq(mean_angle_deg - half_span, mean_angle_deg + half_span, length.out = m)
    vals <- pmin(angle_bounds_deg[[2]], pmax(angle_bounds_deg[[1]], vals))
    if (!is.null(accumulator_order) && length(accumulator_order) > 0L) {
      grp_ids <- acc_ids[grp]
      grp_sorted <- grp[order(rank_top_to_bottom[grp_ids], grp_ids)]
      angle_deg[grp_sorted] <- rev(vals)
    } else {
      angle_deg[grp] <- vals
    }
  }
  angle_rad <- angle_deg * pi / 180

  x0 <- onset_value
  x1 <- x0 + line_length * cos(angle_rad)
  y0 <- rep(ylim[[1]], n_acc)
  y1 <- y0 + line_length * sin(angle_rad)
  y1 <- pmin(ylim[[2]] - 0.02 * yr, pmax(y0 + 0.08 * yr, y1))

  if (is.null(xlim)) {
    xmin <- min(0, min(x0, na.rm = TRUE)) - 0.05 * line_length
    xmax <- max(x1, na.rm = TRUE) + 0.20 * line_length
    xlim <- c(xmin, xmax)
  }

  graphics::plot(
    NA, NA,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    xlab = "",
    ylab = "",
    axes = FALSE,
    main = main,
    ...
  )
  xr <- diff(xlim)
  x_origin <- min(0, min(x0, na.rm = TRUE))
  y_origin <- ylim[[1]]
  graphics::arrows(x_origin, y_origin, xlim[[2]], y_origin, code = 2, length = 0.08, xpd = NA, lwd = axis_lwd)
  graphics::arrows(x_origin, y_origin, x_origin, ylim[[2]], code = 2, length = 0.08, xpd = NA, lwd = axis_lwd)
  graphics::text(xlim[[2]], y_origin - 0.04 * yr, labels = xlab, pos = 2, xpd = NA)
  graphics::text(x_origin - 0.03 * xr, ylim[[2]], labels = ylab, srt = 90, pos = 2, xpd = NA)

  idx_map <- setNames(seq_along(acc_ids), acc_ids)

  reciprocal_pairs <- matrix(character(0), ncol = 2L)
  guard_links <- list()
  if (nrow(outcome_rows) > 0L) {
    for (i in seq_len(nrow(outcome_rows))) {
      payload <- outcome_rows$payload[[i]] %||% list()
      expr <- payload$expr
      if (is.null(expr)) next
      guard_links <- c(guard_links, .pae_collect_guard_source_links(expr))
    }
  }

  blocker_keys <- character(0)
  blocker_pairs <- matrix(character(0), ncol = 2L)
  add_block_pair <- function(from, to) {
    if (identical(from, to)) return(invisible(NULL))
    key <- .pae_directed_key(from, to)
    if (key %in% blocker_keys) return(invisible(NULL))
    blocker_keys <<- c(blocker_keys, key)
    blocker_pairs <<- rbind(blocker_pairs, c(from, to))
    invisible(NULL)
  }

  if (length(guard_links) > 0L) {
    for (lnk in guard_links) {
      blk_acc <- unique(unlist(lapply(lnk$blockers, expand_to_accumulators), use.names = FALSE))
      ref_acc <- unique(unlist(lapply(lnk$references, expand_to_accumulators), use.names = FALSE))
      blk_acc <- intersect(blk_acc, acc_ids)
      ref_acc <- intersect(ref_acc, acc_ids)
      if (length(blk_acc) == 0L || length(ref_acc) == 0L) next
      for (b in blk_acc) {
        for (r in ref_acc) add_block_pair(b, r)
      }
    }
  }

  blocker_pair_fracs <- function(n_pairs) {
    if (n_pairs <= 1L) return(0.62)
    center <- 0.62
    span <- 0.18
    base <- seq(center - span / 2, center + span / 2, length.out = n_pairs)
    idx <- order(abs(base - center), base)
    as.numeric(base[idx])
  }
  pair_anchor_fracs <- blocker_pair_fracs(nrow(blocker_pairs))
  point_on_line <- function(i, frac) {
    c(
      x0[[i]] + frac * (x1[[i]] - x0[[i]]),
      y0[[i]] + frac * (y1[[i]] - y0[[i]])
    )
  }

  if (nrow(blocker_pairs) > 0L) {
    for (k in seq_len(nrow(blocker_pairs))) {
      from <- blocker_pairs[[k, 1]]
      to <- blocker_pairs[[k, 2]]
      i_from <- idx_map[[from]]
      i_to <- idx_map[[to]]
      if (is.null(i_from) || is.null(i_to)) next
      pair_frac <- pair_anchor_fracs[[k]]
      p_from <- point_on_line(i_from, pair_frac)
      p_to <- point_on_line(i_to, pair_frac)
      sign <- if (i_from < i_to) -1 else 1
      curv <- curve_strength * sign * (0.95 + 0.10 * ((k - 1L) %% 3L))
      .pae_draw_curved_arrow(
        p_from[[1]], p_from[[2]],
        p_to[[1]], p_to[[2]],
        col = "firebrick3", lty = "11", lwd = blocker_lwd,
        curvature = curv, both = FALSE, length = 0.11
      )
    }
  }

  graphics::arrows(x0, y0, x1, y1, col = line_col, lwd = line_lwd, code = 2, length = 0.11)
  if (isTRUE(show_labels)) {
    graphics::text(x1, y1, labels = acc_ids, pos = 4, cex = 0.85, col = line_col, xpd = NA)
  }

  invisible(list(
    accumulators = data.frame(
      accumulator = acc_ids,
      onset = as.numeric(onset_value),
      angle_deg = as.numeric(angle_deg),
      x0 = as.numeric(x0),
      y0 = as.numeric(y0),
      x1 = as.numeric(x1),
      y1 = as.numeric(y1),
      stringsAsFactors = FALSE
    ),
    reciprocal_pairs = reciprocal_pairs,
    blocker_pairs = blocker_pairs
  ))
}

plot_accumulator_evidence <- function(...) {
  plot_accumulators(...)
}
