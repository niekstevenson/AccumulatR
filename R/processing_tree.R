.pt_node_definition <- function(source_id, model_view) {
  source_id <- as.character(source_id %||% "")
  if (length(source_id) != 1L || is.na(source_id) || !nzchar(source_id)) {
    source_id <- "?"
  }
  if (!is.null(model_view$accumulators[[source_id]])) {
    lines <- c(source_id, "accumulator")
    return(list(type = "accumulator", label = paste(lines, collapse = "\\n")))
  }

  pool <- model_view$pools[[source_id]]
  if (!is.null(pool)) {
    members <- pool$members %||% character(0)
    k <- suppressWarnings(as.integer(pool$k %||% NA_integer_)[1])
    n <- length(members)
    if (!is.na(k) && n > 0) {
      pool_line <- sprintf("pool: %d-of-%d", k, n)
    } else if (n > 0) {
      pool_line <- sprintf("pool: 1-of-%d", n)
    } else {
      pool_line <- "pool"
    }
    lines <- c(source_id, pool_line)
    return(list(type = "pool", label = paste(lines, collapse = "\\n")))
  }

  list(type = "event", label = paste(c(source_id, "event"), collapse = "\\n"))
}

.pt_to_dot <- function(nodes, edges) {
  if (nrow(nodes) == 0) stop("No nodes to render")

  esc <- function(x) gsub("\"", "\\\\\"", as.character(x %||% ""))

  node_shape <- function(type) {
    switch(type,
      outcome = "box",
      accumulator = "ellipse",
      pool = "box",
      operator = "diamond",
      event = "ellipse",
      "ellipse"
    )
  }

  edge_style <- function(role) {
    switch(role %||% "",
      contributes = list(style = "solid", penwidth = "1.3", dir = "back"),
      reference = list(style = "solid", penwidth = "1.8", dir = "back"),
      blocker = list(style = "bold", color = "#B42318", fontcolor = "#8A1C13", penwidth = "2.0", dir = "back"),
      unless = list(style = "dotted", color = "#B54708", fontcolor = "#8A3A06", dir = "back"),
      member = list(style = "dashed", color = "#4B5563", fontcolor = "#374151", dir = "back"),
      list(style = "solid", dir = "back")
    )
  }

  format_attrs <- function(attrs) {
    parts <- character(0)
    for (nm in names(attrs)) {
      value <- attrs[[nm]]
      if (is.null(value) || is.na(value) || !nzchar(as.character(value))) next
      if (nm %in% c("penwidth", "constraint")) {
        parts <- c(parts, sprintf("%s=%s", nm, value))
      } else {
        parts <- c(parts, sprintf('%s="%s"', nm, esc(value)))
      }
    }
    paste(parts, collapse = ", ")
  }

  node_lines <- apply(nodes, 1, function(row) {
    sprintf('  %s [label="%s", shape=%s];',
      row[["id"]], esc(row[["label"]]), node_shape(row[["type"]])
    )
  })

  edge_lines <- apply(edges, 1, function(row) {
    attrs <- edge_style(row[["role"]])
    label <- row[["label"]]
    if (!is.null(label) && !is.na(label) && nzchar(label)) {
      attrs$label <- label
    }
    attr_txt <- format_attrs(attrs)
    if (nzchar(attr_txt)) {
      sprintf("  %s -> %s [%s];", row[["from"]], row[["to"]], attr_txt)
    } else {
      sprintf("  %s -> %s;", row[["from"]], row[["to"]])
    }
  })

  rank_lines <- character(0)
  if (any(nodes$type == "outcome")) {
    rank_lines <- c(rank_lines, sprintf("  { rank=min; %s; }", paste(nodes$id[nodes$type == "outcome"], collapse = "; ")))
  }
  if (any(nodes$type == "pool")) {
    rank_lines <- c(rank_lines, sprintf("  { rank=same; %s; }", paste(nodes$id[nodes$type == "pool"], collapse = "; ")))
  }
  if (any(nodes$type %in% c("accumulator", "event"))) {
    rank_lines <- c(rank_lines, sprintf("  { rank=max; %s; }", paste(nodes$id[nodes$type %in% c("accumulator", "event")], collapse = "; ")))
  }

  paste0(
    "digraph processing_tree {\n",
    "  rankdir=TB;\n",
    "  graph [fontname=\"Helvetica\", bgcolor=\"white\"];\n",
    "  node [fontname=\"Helvetica\"];\n",
    "  edge [fontname=\"Helvetica\"];\n",
    paste(node_lines, collapse = "\n"), "\n",
    if (length(rank_lines) > 0) paste(rank_lines, collapse = "\n") else "",
    if (length(rank_lines) > 0) "\n" else "",
    paste(edge_lines, collapse = "\n"), "\n",
    "}"
  )
}

#' Draw a processing tree for the responses in a model
#'
#' This is a compact visual summary of the response rules in a model. It avoids
#' equations and instead shows the observed responses, the accumulators or pools
#' that feed them, and any blocking relationships.
#'
#' @param model A race model or finalized model structure.
#' @param outcome_label Optional response label. If supplied, only that response
#'   is shown.
#' @param return_dot If `TRUE`, return the Graphviz DOT string instead of a plot.
#'
#' @return If `DiagrammeR` is available and `return_dot = FALSE`, a
#'   `DiagrammeR` graph. Otherwise, a list with `dot`, `nodes`, and `edges`.
#' @export
processing_tree <- function(model, outcome_label = NULL, return_dot = FALSE) {
  view <- .model_view(model)
  outcomes <- view$outcomes
  if (length(outcomes) == 0L) stop("No outcomes found in model", call. = FALSE)

  labels <- vapply(outcomes, function(out) out$label %||% "", character(1))
  missing_labels <- is.na(labels) | !nzchar(labels)
  labels[missing_labels] <- names(outcomes)[missing_labels]
  if (!is.null(outcome_label)) {
    keep <- labels == outcome_label
    outcomes <- outcomes[keep]
    labels <- labels[keep]
    if (length(outcomes) == 0L) {
      stop(sprintf("Outcome '%s' not found in model", outcome_label), call. = FALSE)
    }
  }

  nodes <- data.frame(id = character(0), label = character(0), type = character(0), stringsAsFactors = FALSE)
  edges <- data.frame(
    from = character(0),
    to = character(0),
    label = character(0),
    role = character(0),
    stringsAsFactors = FALSE
  )
  node_counter <- 0L
  event_nodes <- list()
  pool_expanded <- list()

  next_id <- function(prefix = "n") {
    node_counter <<- node_counter + 1L
    sprintf("%s%d", prefix, node_counter)
  }

  add_node <- function(id, label, type) {
    nodes <<- rbind(nodes, data.frame(id = id, label = label, type = type, stringsAsFactors = FALSE))
    id
  }

  add_edge <- function(from, to, label = "", role = "contributes") {
    edges <<- rbind(
      edges,
      data.frame(from = from, to = to, label = label, role = role, stringsAsFactors = FALSE)
    )
    invisible(NULL)
  }

  ensure_event_node <- function(source_id) {
    source_id <- as.character(source_id %||% "?")
    if (length(source_id) != 1L || is.na(source_id) || !nzchar(source_id)) {
      source_id <- "?"
    }
    key <- sprintf("event:%s", source_id)
    existing <- event_nodes[[key]] %||% NA_character_
    if (!is.na(existing) && nzchar(existing)) return(existing)
    node_id <- next_id("e")
    def <- .pt_node_definition(source_id, view)
    add_node(node_id, def$label, def$type)
    event_nodes[[key]] <<- node_id
    node_id
  }

  expand_pool <- function(source_id) {
    source_id <- as.character(source_id %||% "")
    if (length(source_id) != 1L || is.na(source_id) || !nzchar(source_id)) {
      return(invisible(NULL))
    }
    pool <- view$pools[[source_id]]
    if (is.null(pool)) return(invisible(NULL))
    key <- sprintf("pool:%s", source_id)
    if (isTRUE(pool_expanded[[key]])) return(invisible(NULL))

    pool_id <- ensure_event_node(source_id)
    members <- pool$members %||% character(0)
    if (length(members) > 0) {
      for (i in seq_along(members)) {
        member <- members[[i]]
        member_id <- ensure_event_node(member)
        add_edge(pool_id, member_id, sprintf("member %d", i), role = "member")
        expand_pool(member)
      }
    }
    pool_expanded[[key]] <<- TRUE
    invisible(NULL)
  }

  build_expr <- function(expr, parent_id = NULL, edge_label = "", parent_role = "reference") {
    if (is.null(expr)) return(NA_character_)
    kind <- expr$kind %||% "event"

    if (identical(kind, "event")) {
      source_id <- expr$source %||% "?"
      node_id <- ensure_event_node(source_id)
      expand_pool(source_id)
      if (!is.null(parent_id)) add_edge(parent_id, node_id, edge_label, role = parent_role)
      return(node_id)
    }

    if (identical(kind, "guard")) {
      ref_id <- build_expr(expr$reference, parent_id = parent_id, edge_label = edge_label, parent_role = parent_role)
      block_id <- build_expr(expr$blocker, parent_id = NULL, edge_label = "")
      if (!is.na(block_id) && !is.na(ref_id)) {
        add_edge(ref_id, block_id, "blocks", role = "blocker")
      }
      unless_list <- expr$unless %||% list()
      if (length(unless_list) > 0) {
        for (i in seq_along(unless_list)) {
          unless_id <- build_expr(unless_list[[i]], parent_id = NULL, edge_label = "")
          if (!is.na(unless_id) && !is.na(block_id)) {
            add_edge(block_id, unless_id, sprintf("unless %d", i), role = "unless")
          }
        }
      }
      return(ref_id)
    }

    node_id <- next_id("o")
    op_label <- switch(kind,
      and = "AND",
      or = "OR",
      not = "NOT",
      toupper(kind)
    )
    add_node(node_id, op_label, "operator")
    if (!is.null(parent_id)) add_edge(parent_id, node_id, edge_label, role = parent_role)

    if (identical(kind, "and") || identical(kind, "or")) {
      args <- expr$args %||% list()
      if (length(args) > 0) {
        for (i in seq_along(args)) {
          build_expr(args[[i]], node_id, sprintf("part %d", i), parent_role = "contributes")
        }
      }
    } else if (identical(kind, "not")) {
      build_expr(expr$arg, node_id, "part", parent_role = "contributes")
    }
    node_id
  }

  for (i in seq_along(outcomes)) {
    outcome <- outcomes[[i]]
    expr <- outcome$expr
    outcome_id <- next_id("out")
    add_node(outcome_id, paste(c(labels[[i]], "outcome"), collapse = "\\n"), "outcome")
    if (!is.null(expr)) {
      build_expr(expr, parent_id = outcome_id, edge_label = "defined by", parent_role = "reference")
    }
  }

  dot <- .pt_to_dot(nodes, edges)
  if (!return_dot && requireNamespace("DiagrammeR", quietly = TRUE)) {
    return(DiagrammeR::grViz(dot))
  }
  list(dot = dot, nodes = nodes, edges = edges)
}
