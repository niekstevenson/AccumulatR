`%||%` <- function(lhs, rhs) if (is.null(lhs) || length(lhs) == 0) rhs else lhs

#' Build a processing tree representation for a model outcome
#'
#' @param model A race model (`race_spec`, `race_model_spec`) or `model_tables`.
#' @param outcome_label Outcome label to visualise.
#' @param component Optional mixture component (unused in structure, retained for API symmetry).
#' @param return_dot If `TRUE`, always return the DOT string. If `FALSE` (default) and
#'   DiagrammeR is available, returns a DiagrammeR graph.
#' @return DiagrammeR graph object (if available) or list with `dot`, `nodes`, `edges`.
processing_tree <- function(model, outcome_label, component = NULL, return_dot = FALSE) {
  if (is_model_tables(model)) {
    tables <- model
  } else {
    tables <- model_to_tables(model)
  }
  outcome_struct_id <- sprintf("outcome:%s", outcome_label)
  struct_row <- tables$struct_nodes[tables$struct_nodes$struct_id == outcome_struct_id, , drop = FALSE]
  if (nrow(struct_row) == 0) stop(sprintf("Outcome '%s' not found in model tables", outcome_label))
  outcome_payload <- struct_row$payload[[1]] %||% list()
  expr <- outcome_payload$expr
  if (is.null(expr)) stop(sprintf("Outcome '%s' is missing expression payload", outcome_label))

  struct_lookup <- build_struct_lookup(tables$struct_nodes)

  tree <- build_expr_tree(expr, struct_lookup)
  # Prepend outcome node
  root_id <- tree$new_node_id()
  outcome_node <- data.frame(
    id = root_id,
    label = sprintf("Outcome: %s", outcome_label),
    type = "outcome",
    stringsAsFactors = FALSE
  )
  nodes <- rbind(outcome_node, tree$nodes)
  edges <- rbind(
    data.frame(from = root_id, to = tree$root_id, label = "expr", stringsAsFactors = FALSE),
    tree$edges
  )

  dot <- to_dot(nodes, edges)
  if (!return_dot && requireNamespace("DiagrammeR", quietly = TRUE)) {
    return(DiagrammeR::grViz(dot))
  }
  list(dot = dot, nodes = nodes, edges = edges)
}

build_struct_lookup <- function(struct_nodes) {
  lookup <- list()
  for (i in seq_len(nrow(struct_nodes))) {
    row <- struct_nodes[i, ]
    if (is.null(lookup[[row$label]])) {
      lookup[[row$label]] <- list(
        struct_id = row$struct_id,
        node_type = row$node_type,
        label = row$label,
        dist_family = row$dist_family %||% NA_character_,
        payload = row$payload[[1]] %||% list()
      )
    }
  }
  lookup
}

build_expr_tree <- function(expr, struct_lookup) {
  nodes <- data.frame(id = character(0), label = character(0), type = character(0), stringsAsFactors = FALSE)
  edges <- data.frame(from = character(0), to = character(0), label = character(0), stringsAsFactors = FALSE)
  counter <- 0L
  new_node_id <- function() {
    counter <<- counter + 1L
    sprintf("n%d", counter)
  }

  build_rec <- function(node, parent_id = NULL, edge_label = NULL) {
    node_id <- new_node_id()
    kind <- node$kind %||% "event"
    label <- kind_label(node, struct_lookup)
    nodes <<- rbind(nodes, data.frame(id = node_id, label = label, type = kind, stringsAsFactors = FALSE))
    if (!is.null(parent_id)) {
      edges <<- rbind(edges, data.frame(from = parent_id, to = node_id,
                                        label = edge_label %||% "", stringsAsFactors = FALSE))
    }

    if (identical(kind, "guard")) {
      build_rec(node$reference, node_id, "reference")
      build_rec(node$blocker, node_id, "blocker")
      unless_list <- node$unless %||% list()
      if (length(unless_list) > 0) {
        for (i in seq_along(unless_list)) {
          build_rec(unless_list[[i]], node_id, sprintf("unless[%d]", i))
        }
      }
    } else if (identical(kind, "and") || identical(kind, "or")) {
      args <- node$args %||% list()
      for (i in seq_along(args)) {
        build_rec(args[[i]], node_id, sprintf("arg[%d]", i))
      }
    } else if (identical(kind, "not")) {
      build_rec(node$arg, node_id, "arg")
    }
    node_id
  }

  root_id <- build_rec(expr)
  list(nodes = nodes, edges = edges, root_id = root_id, new_node_id = new_node_id)
}

kind_label <- function(node, struct_lookup) {
  kind <- node$kind %||% "event"
  if (identical(kind, "event")) {
    src <- node$source %||% "?"
    info <- struct_lookup[[src]] %||% list(node_type = "event")
    sprintf("Event: %s\n(%s)", src, info$node_type)
  } else if (identical(kind, "guard")) {
    sprintf("Guard\nref:%s\nblock:%s", label_for_event(node$reference, struct_lookup),
            label_for_event(node$blocker, struct_lookup))
  } else if (identical(kind, "and")) {
    "AND"
  } else if (identical(kind, "or")) {
    "OR"
  } else if (identical(kind, "not")) {
    "NOT"
  } else {
    toupper(kind)
  }
}

label_for_event <- function(node, struct_lookup) {
  if (is.null(node)) return("?")
  src <- node$source %||% "?"
  info <- struct_lookup[[src]] %||% list(node_type = "event")
  sprintf("%s (%s)", src, info$node_type)
}

.simple_node_definition <- function(source_id, struct_lookup) {
  info <- struct_lookup[[source_id]] %||% list(node_type = "event", payload = list())
  node_type <- info$node_type %||% "event"
  payload <- info$payload %||% list()

  if (identical(node_type, "accumulator")) {
    lines <- c(source_id, "accumulator")
    onset <- payload$onset %||% 0
    onset_num <- suppressWarnings(as.numeric(onset)[1])
    if (!is.na(onset_num) && onset_num != 0) {
      lines <- c(lines, sprintf("onset: %.3f", onset_num))
    }
    return(list(type = "accumulator", label = paste(lines, collapse = "\\n")))
  }

  if (identical(node_type, "pool")) {
    members <- payload$members %||% character(0)
    rule <- payload$rule %||% list()
    k <- suppressWarnings(as.integer(rule$k %||% NA_integer_)[1])
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

  list(type = "event", label = paste(c(source_id, node_type), collapse = "\\n"))
}

to_dot_simple <- function(nodes, edges) {
  if (nrow(nodes) == 0) stop("No nodes to render")

  esc <- function(x) gsub("\"", "\\\\\"", as.character(x %||% ""))

  node_shape <- function(type) {
    switch(type,
           outcome = "box",
           accumulator = "ellipse",
           pool = "box",
           operator = "diamond",
           event = "ellipse",
           "ellipse")
  }

  node_lines <- apply(nodes, 1, function(row) {
    sprintf('  %s [label="%s", shape=%s];',
            row[["id"]], esc(row[["label"]]), node_shape(row[["type"]]))
  })

  edge_style <- function(role) {
    switch(role %||% "",
           contributes = list(style = "solid", penwidth = "1.3", dir = "back"),
           reference = list(style = "solid", penwidth = "1.8", dir = "back"),
           blocker = list(style = "bold", color = "#B42318", fontcolor = "#8A1C13", penwidth = "2.0", dir = "back"),
           unless = list(style = "dotted", color = "#B54708", fontcolor = "#8A3A06", dir = "back"),
           member = list(style = "dashed", color = "#4B5563", fontcolor = "#374151", dir = "back"),
           list(style = "solid", dir = "back"))
  }

  format_attrs <- function(attrs) {
    parts <- character(0)
    for (nm in names(attrs)) {
      v <- attrs[[nm]]
      if (is.null(v) || is.na(v) || !nzchar(as.character(v))) next
      if (nm %in% c("penwidth", "constraint")) {
        parts <- c(parts, sprintf("%s=%s", nm, v))
      } else {
        parts <- c(parts, sprintf('%s="%s"', nm, esc(v)))
      }
    }
    paste(parts, collapse = ", ")
  }

  edge_lines <- apply(edges, 1, function(row) {
    label <- row[["label"]]
    role <- if ("role" %in% names(edges)) row[["role"]] else ""
    attrs <- edge_style(role)
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

  paste0(
    "digraph processing_tree_simple {\n",
    "  rankdir=TB;\n",
    "  graph [fontname=\"Helvetica\", bgcolor=\"white\"];\n",
    "  node [fontname=\"Helvetica\"];\n",
    "  edge [fontname=\"Helvetica\"];\n",
    paste(node_lines, collapse = "\n"), "\n",
    if (any(nodes$type == "outcome")) {
      sprintf("  { rank=min; %s; }\n", paste(nodes$id[nodes$type == "outcome"], collapse = "; "))
    } else "",
    paste(edge_lines, collapse = "\n"), "\n",
    "}"
  )
}

#' Build a minimal processing tree without mathematical formulas.
#'
#' @param model Race model (spec/tables).
#' @param outcome_label Optional single outcome label. If NULL, include all outcomes.
#' @param component Optional component (reserved for future use).
#' @param return_dot If TRUE, return DOT string.
#' @return DiagrammeR graph object (if available) or list with DOT/nodes/edges.
processing_simple_tree <- function(model, outcome_label = NULL, component = NULL, return_dot = FALSE) {
  if (is_model_tables(model)) {
    tables <- model
  } else {
    tables <- model_to_tables(model)
  }

  struct_lookup <- build_struct_lookup(tables$struct_nodes)
  outcome_rows <- tables$struct_nodes[tables$struct_nodes$node_type == "outcome", , drop = FALSE]
  if (nrow(outcome_rows) == 0) stop("No outcomes found in model tables")

  if (!is.null(outcome_label)) {
    outcome_rows <- outcome_rows[outcome_rows$label == outcome_label, , drop = FALSE]
    if (nrow(outcome_rows) == 0) {
      stop(sprintf("Outcome '%s' not found in model tables", outcome_label))
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
    key <- sprintf("event:%s", source_id)
    existing <- event_nodes[[key]] %||% NA_character_
    if (!is.na(existing) && nzchar(existing)) return(existing)
    node_id <- next_id("e")
    def <- .simple_node_definition(source_id, struct_lookup)
    add_node(node_id, def$label, def$type)
    event_nodes[[key]] <<- node_id
    node_id
  }

  expand_pool <- function(source_id) {
    info <- struct_lookup[[source_id]] %||% list(node_type = "event", payload = list())
    if (!identical(info$node_type %||% "", "pool")) return(invisible(NULL))
    key <- sprintf("pool:%s", source_id)
    if (isTRUE(pool_expanded[[key]])) return(invisible(NULL))

    pool_id <- ensure_event_node(source_id)
    members <- (info$payload %||% list())$members %||% character(0)
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

  build_expr <- function(expr, parent_id = NULL, edge_label = "") {
    if (is.null(expr)) return(NA_character_)
    kind <- expr$kind %||% "event"

    if (identical(kind, "event")) {
      source_id <- expr$source %||% "?"
      node_id <- ensure_event_node(source_id)
      expand_pool(source_id)
      if (!is.null(parent_id)) add_edge(parent_id, node_id, edge_label, role = "contributes")
      return(node_id)
    }

    if (identical(kind, "guard")) {
      ref_id <- build_expr(expr$reference, parent_id = parent_id, edge_label = edge_label)
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
                       toupper(kind))
    add_node(node_id, op_label, "operator")
    if (!is.null(parent_id)) add_edge(parent_id, node_id, edge_label, role = "reference")

    if (identical(kind, "and") || identical(kind, "or")) {
      args <- expr$args %||% list()
      if (length(args) > 0) {
        for (i in seq_along(args)) {
          build_expr(args[[i]], node_id, sprintf("part %d", i))
        }
      }
    } else if (identical(kind, "not")) {
      build_expr(expr$arg, node_id, "part")
    }
    node_id
  }

  for (i in seq_len(nrow(outcome_rows))) {
    row <- outcome_rows[i, ]
    payload <- row$payload[[1]] %||% list()
    expr <- payload$expr
    outcome_id <- next_id("out")
    add_node(outcome_id, paste(c(row$label, "outcome"), collapse = "\\n"), "outcome")
    if (!is.null(expr)) {
      build_expr(expr, parent_id = outcome_id, edge_label = "defined by")
    }
  }

  dot <- to_dot_simple(nodes, edges)
  if (!return_dot && requireNamespace("DiagrammeR", quietly = TRUE)) {
    return(DiagrammeR::grViz(dot))
  }
  list(dot = dot, nodes = nodes, edges = edges)
}

to_dot <- function(nodes, edges) {
  if (nrow(nodes) == 0) stop("No nodes to render")

  dot_escape <- function(x) {
    gsub("\"", "\\\\\"", as.character(x %||% ""))
  }

  type_shape <- function(type) {
    switch(type,
           outcome = "box",
           guard = "octagon",
           and = "diamond",
           or = "diamond",
           not = "diamond",
           event = "ellipse",
           likelihood = "doublecircle",
           operator = "diamond",
           expr = "box",
           factor = "plaintext",
           competitor = "ellipse",
           definitions = "note",
           definitions_general = "note",
           "ellipse")
  }

  type_style <- function(type) {
    switch(type,
           likelihood = list(fillcolor = "#16324F", color = "#16324F", fontcolor = "white", style = "filled,bold", penwidth = "1.8"),
           outcome = list(fillcolor = "#DDF2E3", color = "#2C7A4B", fontcolor = "#153D2A", style = "filled"),
           event = list(fillcolor = "#E8F1F8", color = "#2A5778", fontcolor = "#1B364B", style = "filled"),
           operator = list(fillcolor = "#FFF4E6", color = "#C77D1F", fontcolor = "#6C3D00", style = "filled"),
           definitions = list(fillcolor = "#F4F5F7", color = "#6B7280", fontcolor = "#374151", style = "filled"),
           definitions_general = list(fillcolor = "#F4F5F7", color = "#6B7280", fontcolor = "#374151", style = "filled"),
           expr = list(fillcolor = "#F5F7FA", color = "#5A6570", fontcolor = "#1F2933", style = "filled"),
           list(fillcolor = "#F5F7FA", color = "#5A6570", fontcolor = "#1F2933", style = "filled"))
  }

  edge_style <- function(edge_type) {
    switch(edge_type %||% "",
           likelihood_main = list(color = "#1F4E79", fontcolor = "#1F4E79", penwidth = "1.8"),
           likelihood_competitor = list(color = "#8B2F2F", fontcolor = "#8B2F2F", style = "dashed"),
           arg = list(color = "#6B7280", fontcolor = "#4B5563"),
           reference = list(color = "#2C7A4B", fontcolor = "#1F5B38", penwidth = "1.5"),
           blocker = list(color = "#B42318", fontcolor = "#8A1C13", penwidth = "1.5"),
           unless = list(color = "#B54708", fontcolor = "#8A3A06", style = "dashed"),
           legend = list(color = "#6B7280", fontcolor = "#6B7280", style = "dotted"),
           layout = list(style = "invis"),
           list(color = "#6B7280", fontcolor = "#4B5563"))
  }

  cluster_label <- function(cluster_id) {
    switch(cluster_id,
           target = "Target Outcome Path",
           competitor = "Competitor Paths",
           helper = "Helper / Logic Nodes",
           cluster_id)
  }

  cluster_style <- function(cluster_id) {
    switch(cluster_id,
           target = list(color = "#9DBDD5", style = "rounded,dashed"),
           competitor = list(color = "#D9A3A3", style = "rounded,dashed"),
           helper = list(color = "#C2C7CF", style = "rounded,dashed"),
           list(color = "#C2C7CF", style = "rounded,dashed"))
  }

  format_attrs <- function(attrs) {
    vals <- character(0)
    for (nm in names(attrs)) {
      value <- attrs[[nm]]
      if (is.null(value) || is.na(value) || !nzchar(as.character(value))) next
      if (nm %in% c("shape", "penwidth")) {
        vals <- c(vals, sprintf("%s=%s", nm, value))
      } else {
        vals <- c(vals, sprintf('%s="%s"', nm, dot_escape(value)))
      }
    }
    paste(vals, collapse = ", ")
  }

  node_lines <- apply(nodes, 1, function(row) {
    type <- row[['type']] %||% "expr"
    attrs <- c(list(
      label = row[['label']],
      shape = type_shape(type)
    ), type_style(type))
    if ("fillcolor" %in% names(nodes) && !is.na(row[['fillcolor']]) && row[['fillcolor']] != "") {
      attrs$fillcolor <- row[['fillcolor']]
    }
    if ("color" %in% names(nodes) && !is.na(row[['color']]) && row[['color']] != "") {
      attrs$color <- row[['color']]
    }
    if ("fontcolor" %in% names(nodes) && !is.na(row[['fontcolor']]) && row[['fontcolor']] != "") {
      attrs$fontcolor <- row[['fontcolor']]
    }
    sprintf("  %s [%s];", row[['id']], format_attrs(attrs))
  })
  edge_lines <- apply(edges, 1, function(row) {
    lbl <- row[['label']]
    style <- if ("style" %in% names(edges)) row[['style']] else ""
    edge_type <- if ("edge_type" %in% names(edges)) row[['edge_type']] else ""
    attrs <- edge_style(edge_type)
    if (!is.null(style) && !is.na(style) && style != "") {
      attrs$style <- style
    }
    if (!is.null(lbl) && !is.na(lbl) && lbl != "") {
      attrs$label <- lbl
    }
    if (length(attrs) > 0) {
      sprintf("  %s -> %s [%s];", row[['from']], row[['to']], format_attrs(attrs))
    } else {
      sprintf("  %s -> %s;", row[['from']], row[['to']])
    }
  })

  rank_lines <- character(0)
  if ("rank" %in% names(nodes)) {
    ranks <- unique(nodes$rank[!is.na(nodes$rank) & nodes$rank != ""])
    for (rk in ranks) {
      idx <- which(!is.na(nodes$rank) & nodes$rank == rk)
      ids <- nodes$id[idx]
      if (length(ids) > 0) {
        rank_lines <- c(rank_lines, sprintf("  { rank=%s; %s; }", rk, paste(ids, collapse = "; ")))
      }
    }
  }

  cluster_lines <- character(0)
  if ("cluster" %in% names(nodes)) {
    cluster_ids <- unique(nodes$cluster[!is.na(nodes$cluster) & nodes$cluster != ""])
    for (cid in cluster_ids) {
      ids <- nodes$id[!is.na(nodes$cluster) & nodes$cluster == cid]
      if (length(ids) == 0L) next
      st <- cluster_style(cid)
      attrs <- c(
        sprintf('label="%s"', dot_escape(cluster_label(cid))),
        sprintf('color="%s"', dot_escape(st$color %||% "#C2C7CF")),
        sprintf('style="%s"', dot_escape(st$style %||% "rounded,dashed"))
      )
      cluster_lines <- c(
        cluster_lines,
        sprintf("  subgraph cluster_%s {", .math_sanitize_token(cid)),
        sprintf("    graph [%s];", paste(attrs, collapse = ", ")),
        sprintf("    %s;", paste(ids, collapse = "; ")),
        "  }"
      )
    }
  }

  paste0(
    "digraph processing_tree {\n",
    "  rankdir=TB;\n",
    "  graph [fontname=\"Helvetica\", bgcolor=\"white\"];\n",
    "  node [fontname=\"Helvetica\", margin=\"0.08,0.04\"];\n",
    "  edge [fontname=\"Helvetica\"];\n",
    paste(node_lines, collapse = "\n"), "\n",
    if (length(cluster_lines) > 0) paste(cluster_lines, collapse = "\n") else "",
    if (length(cluster_lines) > 0) "\n" else "",
    if (length(rank_lines) > 0) paste(rank_lines, collapse = "\n") else "",
    if (length(rank_lines) > 0) "\n" else "",
    paste(edge_lines, collapse = "\n"), "\n",
    "}"
  )
}

.is_trivial_formula <- function(line) {
  if (!grepl("=", line, fixed = TRUE)) return(FALSE)
  rhs <- trimws(sub("^[^=]+=", "", line))
  if (rhs == "") return(TRUE)
  if (grepl("integral", rhs, fixed = TRUE) || grepl("\\*", rhs) ||
      grepl("\\+", rhs) || grepl("sum", rhs, fixed = TRUE) ||
      grepl("prod", rhs, fixed = TRUE) || grepl("∑", rhs) ||
      grepl("∏", rhs)) return(FALSE)
  if (grepl("^1 - F_", rhs)) return(TRUE)
  if (grepl("^F_", rhs)) return(TRUE)
  if (grepl("^S_.*\\(t\\)$", rhs)) return(TRUE)
  FALSE
}

._filter_formulas <- function(formulas) {
  formulas[!vapply(formulas, .is_trivial_formula, logical(1))]
}

.describe_expr_short <- function(expr) {
  if (is.null(expr)) return("?")
  kind <- expr$kind %||% "event"
  if (identical(kind, "event")) {
    return(expr$source %||% "?")
  }
  if (identical(kind, "guard")) {
    ref_txt <- .describe_expr_short(expr$reference)
    block_txt <- .describe_expr_short(expr$blocker)
    unless_txt <- expr$unless %||% list()
    if (length(unless_txt) > 0) {
      unl_lab <- paste(vapply(unless_txt, .describe_expr_short, character(1)), collapse = ", ")
      return(sprintf("guard(%s | !%s, unless %s)", ref_txt, block_txt, unl_lab))
    }
    return(sprintf("guard(%s | !%s)", ref_txt, block_txt))
  }
  if (identical(kind, "and") || identical(kind, "or")) {
    op <- if (identical(kind, "and")) " & " else " | "
    args <- expr$args %||% list()
    if (length(args) == 0) return(sprintf("%s()", kind))
    return(paste(vapply(args, .describe_expr_short, character(1)), collapse = op))
  }
  if (identical(kind, "not")) {
    return(sprintf("not %s", .describe_expr_short(expr$arg)))
  }
  kind
}

.math_sanitize_token <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) x <- "node"
  if (grepl("^[0-9]", x)) x <- paste0("N_", x)
  x
}

.math_truncate <- function(x, max_chars = 64L) {
  x <- as.character(x %||% "")
  if (nchar(x, type = "chars") <= max_chars) return(x)
  paste0(substr(x, 1L, max_chars - 3L), "...")
}

.math_unique_formula_name <- function(ctx, candidate) {
  candidate <- .math_sanitize_token(candidate)
  counts <- ctx$formula_name_counts
  n <- counts[[candidate]] %||% 0L
  n <- as.integer(n) + 1L
  counts[[candidate]] <- n
  ctx$formula_name_counts <- counts
  if (n <= 1L) candidate else sprintf("%s_%d", candidate, n)
}

.math_formula_name <- function(ctx, symbol) {
  info <- ctx$node_info[[symbol]] %||% list()
  info$formula_name %||% symbol
}

.math_assign_outcome_alias <- function(ctx, symbol, outcome_label) {
  alias <- .math_unique_formula_name(ctx, paste0("O_", outcome_label))
  ctx$outcome_alias[[symbol]] <- alias
  alias
}

.math_symbol_ref <- function(ctx, symbol, prefer_outcome_alias = FALSE) {
  if (prefer_outcome_alias) {
    alias <- ctx$outcome_alias[[symbol]] %||% NA_character_
    if (!is.na(alias) && nzchar(alias)) return(alias)
  }
  .math_formula_name(ctx, symbol)
}

.math_label_from_info <- function(ctx, info) {
  mode <- ctx$render_mode %||% "compact"
  formulas <- info$formulas %||% character(0)
  name <- info$formula_name %||% info$symbol %||% "node"
  if (identical(mode, "verbose")) {
    if (length(formulas) == 0) {
      formulas <- c(name, info$description %||% info$display %||% info$label %||% "")
      formulas <- formulas[nzchar(formulas)]
    }
    return(paste(formulas, collapse = "\\n"))
  }
  detail <- info$description %||% info$display %||% info$label %||% ""
  detail <- gsub("[\r\n]+", " ", detail)
  detail <- .math_truncate(detail, max_chars = 70L)
  if (!nzchar(detail) || identical(detail, name)) {
    return(name)
  }
  paste(name, detail, sep = "\\n")
}

.math_update_node_label <- function(ctx, node_id, label) {
  idx <- which(ctx$nodes$id == node_id)
  if (length(idx) > 0) {
    ctx$nodes$label[idx[1]] <- label
  }
}

.math_set_node_cluster <- function(ctx, node_ids, cluster) {
  if (length(node_ids) == 0L) return(invisible(NULL))
  idx <- which(ctx$nodes$id %in% node_ids)
  if (length(idx) > 0L) {
    ctx$nodes$cluster[idx] <- cluster
  }
  invisible(NULL)
}

.collect_descendants <- function(edges, root_ids) {
  roots <- unique(root_ids[!is.na(root_ids) & nzchar(root_ids)])
  if (length(roots) == 0L || nrow(edges) == 0L) return(character(0))
  out <- character(0)
  queue <- roots
  while (length(queue) > 0L) {
    current <- queue[[1]]
    queue <- queue[-1]
    if (current %in% out) next
    out <- c(out, current)
    kids <- edges$to[edges$from == current]
    kids <- kids[!is.na(kids) & nzchar(kids)]
    queue <- c(queue, kids)
  }
  unique(out)
}

.math_new_ctx <- function(struct_lookup, mode = c("compact", "verbose")) {
  mode <- match.arg(mode)
  ctx <- new.env(parent = emptyenv())
  ctx$struct_lookup <- struct_lookup
  ctx$id_counter <- 0L
  ctx$symbol_counter <- 0L
  ctx$render_mode <- mode
  ctx$formula_name_counts <- list()
  ctx$symbol_lookup <- list()
  ctx$node_lookup <- list()
  ctx$node_info <- list()
  ctx$node_to_symbol <- list()
  ctx$symbol_sequence <- character(0)
  ctx$outcome_map <- list()
  ctx$outcome_order <- character(0)
  ctx$outcome_alias <- list()
  ctx$name_to_symbol <- list()
  ctx$nodes <- data.frame(
    id = character(0),
    label = character(0),
    type = character(0),
    rank = character(0),
    cluster = character(0),
    stringsAsFactors = FALSE
  )
  ctx$edges <- data.frame(
    from = character(0),
    to = character(0),
    label = character(0),
    style = character(0),
    edge_type = character(0),
    stringsAsFactors = FALSE
  )
  ctx
}

.math_next_id <- function(ctx, prefix = "m") {
  ctx$id_counter <- ctx$id_counter + 1L
  sprintf("%s%d", prefix, ctx$id_counter)
}

.math_next_symbol <- function(ctx, prefix = "N", hint = NULL) {
  if (!is.null(hint) && nzchar(hint)) {
    candidate <- sprintf("%s_%s", prefix, .math_sanitize_token(hint))
    return(.math_unique_formula_name(ctx, candidate))
  }
  ctx$symbol_counter <- ctx$symbol_counter + 1L
  .math_unique_formula_name(ctx, sprintf("%s_%d", prefix, ctx$symbol_counter))
}

.math_add_node <- function(ctx, id, label, type, rank = NA_character_, cluster = NA_character_) {
  ctx$nodes <- rbind(
    ctx$nodes,
    data.frame(id = id, label = label, type = type,
               rank = rank %||% NA_character_,
               cluster = cluster %||% NA_character_,
               stringsAsFactors = FALSE)
  )
}

.math_add_edge <- function(ctx, from, to, label = "", style = "", edge_type = NA_character_) {
  ctx$edges <- rbind(
    ctx$edges,
    data.frame(from = from, to = to, label = label, style = style, edge_type = edge_type,
               stringsAsFactors = FALSE)
  )
}

.math_register_symbol <- function(ctx, symbol, info) {
  info$symbol <- symbol
  if (is.null(info$description) && !is.null(info$label)) {
    info$description <- info$label
  }
  if (is.null(info$formula_name)) {
    info$formula_name <- symbol
  }
  if (is.null(info$display_name)) {
    info$display_name <- info$description %||% symbol
  }
  ctx$node_info[[symbol]] <- info
  ctx$name_to_symbol[[info$formula_name]] <- symbol
  node_id <- info$node_id %||% NA_character_
  if (!is.na(node_id) && nzchar(node_id)) {
    ctx$node_to_symbol[[node_id]] <- symbol
    .math_update_node_label(ctx, node_id, .math_label_from_info(ctx, info))
  }
  if (!(symbol %in% ctx$symbol_sequence)) {
    ctx$symbol_sequence <- c(ctx$symbol_sequence, symbol)
  }
}

.math_explicit_sum <- function(child_syms, base_prefix, other_prefix) {
  if (length(child_syms) == 0) return("0")
  terms <- vapply(seq_along(child_syms), function(i) {
    base <- sprintf("%s_%s(t)", base_prefix, child_syms[[i]])
    others <- child_syms[-i]
    if (length(others) == 0) return(base)
    paste0(base, " * ", paste(sprintf("%s_%s(t)", other_prefix, others), collapse = " * "))
  }, character(1))
  paste(terms, collapse = " + ")
}

.math_and_density_line <- function(symbol, child_syms) {
  if (length(child_syms) == 0) return(sprintf("f_%s(t) = 0", symbol))
  if (length(child_syms) <= 4) {
    rhs <- .math_explicit_sum(child_syms, "f", "F")
    return(sprintf("f_%s(t) = %s", symbol, rhs))
  }
  sprintf("f_%s(t) = sum_i f_{child_i}(t) * prod_{j!=i} F_{child_j}(t)", symbol)
}

.math_and_cdf_line <- function(symbol, child_syms) {
  if (length(child_syms) == 0) return(sprintf("F_%s(t) = 1", symbol))
  if (length(child_syms) <= 4) {
    rhs <- paste(sprintf("F_%s(t)", child_syms), collapse = " * ")
    return(sprintf("F_%s(t) = %s", symbol, rhs))
  }
  sprintf("F_%s(t) = prod_i F_{child_i}(t)", symbol)
}

.math_or_density_line <- function(symbol, child_syms) {
  if (length(child_syms) == 0) return(sprintf("f_%s(t) = 0", symbol))
  if (length(child_syms) <= 4) {
    rhs <- .math_explicit_sum(child_syms, "f", "S")
    return(sprintf("f_%s(t) = %s", symbol, rhs))
  }
  sprintf("f_%s(t) = sum_i f_{child_i}(t) * prod_{j!=i} S_{child_j}(t)", symbol)
}

.math_or_cdf_line <- function(symbol, child_syms) {
  if (length(child_syms) == 0) return(sprintf("F_%s(t) = 0", symbol))
  if (length(child_syms) <= 4) {
    rhs <- paste(sprintf("S_%s(t)", child_syms), collapse = " * ")
    return(sprintf("F_%s(t) = 1 - %s", symbol, rhs))
  }
  sprintf("F_%s(t) = 1 - prod_i S_{child_i}(t)", symbol)
}

.math_build_expr <- function(ctx, expr) {
  kind <- expr$kind %||% "event"

  if (identical(kind, "event")) {
    source_id <- expr$source %||% "?"
    key <- sprintf("event:%s", source_id)
    if (!is.null(ctx$symbol_lookup[[key]])) {
      symbol <- ctx$symbol_lookup[[key]]
      node_id <- ctx$node_lookup[[key]]
      info <- ctx$node_info[[symbol]] %||% list()
      display <- info$display %||% info$label %||% source_id
      return(list(node_id = node_id, symbol = symbol, display = display))
    }
    symbol <- .math_next_symbol(ctx, prefix = "E", hint = source_id)
    node_id <- .math_next_id(ctx, "m")
    info <- ctx$struct_lookup[[source_id]] %||% list(node_type = "event", payload = list())
    display <- source_id
    type_txt <- info$node_type %||% "event"
    payload <- info$payload %||% list()
    members <- payload$members %||% character(0)
    rule <- payload$rule %||% list()
    event_name <- symbol
    base_formulas <- c(
      sprintf("f_%s(t) = f_%s(t)", event_name, event_name),
      sprintf("F_%s(t) = ∫ f_%s(u) du", event_name, event_name),
      sprintf("S_%s(t) = 1 - F_%s(t)", event_name, event_name)
    )
    formulas <- ._filter_formulas(base_formulas)
    if (length(formulas) == 0) formulas <- base_formulas[1]
    .math_add_node(ctx, node_id, event_name, "event")
    ctx$symbol_lookup[[key]] <- symbol
    ctx$node_lookup[[key]] <- node_id
    .math_register_symbol(ctx, symbol, list(
      kind = "event",
      label = source_id,
      description = sprintf("%s (%s)", source_id, type_txt),
      display = display,
      display_name = source_id,
      type = type_txt,
      members = members,
      rule = list(
        k = if (!is.null(rule$k)) rule$k else if (length(members) > 0) 1L else NA_integer_,
        n = length(members)
      ),
      formulas = formulas,
      formula_name = event_name,
      node_id = node_id
    ))
    return(list(node_id = node_id, symbol = symbol, display = display))
  }

  if (identical(kind, "and")) {
    node_id <- .math_next_id(ctx, "m")
    args <- expr$args %||% list()
    child_results <- lapply(seq_along(args), function(i) .math_build_expr(ctx, args[[i]]))
    child_syms <- if (length(child_results) > 0) vapply(child_results, `[[`, character(1), "symbol") else character(0)
    child_names <- if (length(child_syms) > 0) vapply(child_syms, function(s) .math_formula_name(ctx, s), character(1)) else character(0)
    child_disp <- if (length(child_results) > 0) vapply(child_results, `[[`, character(1), "display") else character(0)
    display <- if (length(child_disp) > 0) paste(child_disp, collapse = " & ") else "AND()"
    desc_text <- if (length(child_disp) > 0) sprintf("AND(%s)", paste(child_disp, collapse = ", ")) else "AND()"
    symbol <- .math_next_symbol(
      ctx,
      prefix = "A",
      hint = if (length(child_names) > 0) paste(child_names, collapse = "_and_") else "empty"
    )
    self_name <- symbol
    formulas <- c(.math_and_density_line(self_name, child_names), .math_and_cdf_line(self_name, child_names),
                  sprintf("S_%s(t) = 1 - F_%s(t)", self_name, self_name))
    formulas <- ._filter_formulas(formulas)
    if (length(formulas) == 0) formulas <- .math_and_density_line(self_name, child_names)
    .math_add_node(ctx, node_id, self_name, "operator")
    .math_register_symbol(ctx, symbol, list(
      kind = "and",
      label = "AND",
      display = display,
      description = desc_text,
      formulas = formulas,
      formula_name = self_name,
      node_id = node_id
    ))
    if (length(child_results) > 0) {
      for (i in seq_along(child_results)) {
        .math_add_edge(ctx, node_id, child_results[[i]]$node_id, sprintf("arg[%d]", i), edge_type = "arg")
      }
    }
    return(list(node_id = node_id, symbol = symbol, display = display))
  }

  if (identical(kind, "or")) {
    node_id <- .math_next_id(ctx, "m")
    args <- expr$args %||% list()
    child_results <- lapply(seq_along(args), function(i) .math_build_expr(ctx, args[[i]]))
    child_syms <- if (length(child_results) > 0) vapply(child_results, `[[`, character(1), "symbol") else character(0)
    child_names <- if (length(child_syms) > 0) vapply(child_syms, function(s) .math_formula_name(ctx, s), character(1)) else character(0)
    child_disp <- if (length(child_results) > 0) vapply(child_results, `[[`, character(1), "display") else character(0)
    display <- if (length(child_disp) > 0) paste(child_disp, collapse = " | ") else "OR()"
    desc_text <- if (length(child_disp) > 0) sprintf("OR(%s)", paste(child_disp, collapse = ", ")) else "OR()"
    symbol <- .math_next_symbol(
      ctx,
      prefix = "R",
      hint = if (length(child_names) > 0) paste(child_names, collapse = "_or_") else "empty"
    )
    self_name <- symbol
    formulas <- c(.math_or_density_line(self_name, child_names), .math_or_cdf_line(self_name, child_names),
                  sprintf("S_%s(t) = 1 - F_%s(t)", self_name, self_name))
    formulas <- ._filter_formulas(formulas)
    if (length(formulas) == 0) formulas <- .math_or_density_line(self_name, child_names)
    .math_add_node(ctx, node_id, self_name, "operator")
    .math_register_symbol(ctx, symbol, list(
      kind = "or",
      label = "OR",
      display = display,
      description = desc_text,
      formulas = formulas,
      formula_name = self_name,
      node_id = node_id
    ))
    if (length(child_results) > 0) {
      for (i in seq_along(child_results)) {
        .math_add_edge(ctx, node_id, child_results[[i]]$node_id, sprintf("arg[%d]", i), edge_type = "arg")
      }
    }
    return(list(node_id = node_id, symbol = symbol, display = display))
  }

  if (identical(kind, "guard")) {
    node_id <- .math_next_id(ctx, "m")
    reference <- expr$reference
    blocker <- expr$blocker
    unless_list <- expr$unless %||% list()
    ref_res <- .math_build_expr(ctx, reference)
    block_res <- .math_build_expr(ctx, blocker)
    unless_res <- lapply(seq_along(unless_list), function(i) .math_build_expr(ctx, unless_list[[i]]))
    unless_syms <- if (length(unless_res) > 0) vapply(unless_res, `[[`, character(1), "symbol") else character(0)
    display <- sprintf("%s guarded by %s", ref_res$display, block_res$display)
    ref_name <- .math_formula_name(ctx, ref_res$symbol)
    block_name <- .math_formula_name(ctx, block_res$symbol)
    unless_names <- if (length(unless_syms) > 0) {
      vapply(unless_syms, function(s) .math_formula_name(ctx, s), character(1))
    } else character(0)
    symbol <- .math_next_symbol(
      ctx,
      prefix = "G",
      hint = paste(c(ref_name, "not", block_name, unless_names), collapse = "_")
    )
    self_name <- symbol
    formulas <- c(
      sprintf("f_%s(t) = f_%s(t) * S_%s^{eff}(t)", self_name, ref_name, block_name),
      sprintf("F_%s(t) = ∫ f_%s(u) du", self_name, self_name),
      sprintf("S_%s(t) = 1 - F_%s(t)", self_name, self_name)
    )
    if (length(unless_syms) == 0) {
      formulas <- c(formulas, sprintf("S_%s^{eff}(t) = S_%s(t)", block_name, block_name))
    } else {
      protector_terms <- paste(sprintf("S_%s(u)", unless_names), collapse = " * ")
      formulas <- c(formulas, sprintf("S_%s^{eff}(t) = 1 - ∫ f_%s(u) * %s du",
                                      block_name, block_name, protector_terms))
    }
    formulas <- ._filter_formulas(formulas)
    if (length(formulas) == 0) formulas <- sprintf("f_%s(t) = f_%s(t)", self_name, ref_name)
    .math_add_node(ctx, node_id, self_name, "operator")
    .math_register_symbol(ctx, symbol, list(
      kind = "guard",
      label = "GUARD",
      display = display,
      description = sprintf("guard(%s | !%s)", ref_res$display, block_res$display),
      formulas = formulas,
      formula_name = self_name,
      node_id = node_id
    ))
    .math_add_edge(ctx, node_id, ref_res$node_id, "reference", edge_type = "reference")
    .math_add_edge(ctx, node_id, block_res$node_id, "blocker", edge_type = "blocker")
    if (length(unless_res) > 0) {
      for (i in seq_along(unless_res)) {
        .math_add_edge(ctx, node_id, unless_res[[i]]$node_id, sprintf("unless[%d]", i), edge_type = "unless")
      }
    }
    return(list(node_id = node_id, symbol = symbol, display = display))
  }

  if (identical(kind, "not")) {
    node_id <- .math_next_id(ctx, "m")
    child <- expr$arg
    child_res <- .math_build_expr(ctx, child)
    display <- sprintf("not %s", child_res$display)
    child_name <- .math_formula_name(ctx, child_res$symbol)
    symbol <- .math_next_symbol(ctx, prefix = "N", hint = paste("not", child_name, sep = "_"))
    self_name <- symbol
    formulas <- c(
      sprintf("f_%s(t) = 0", self_name),
      sprintf("F_%s(t) = 1 - F_%s(t)", self_name, child_name),
      sprintf("S_%s(t) = 1 - F_%s(t)", self_name, self_name)
    )
    formulas <- ._filter_formulas(formulas)
    if (length(formulas) == 0) formulas <- sprintf("f_%s(t) = 0", self_name)
    .math_add_node(ctx, node_id, self_name, "operator")
    .math_register_symbol(ctx, symbol, list(
      kind = "not",
      label = "NOT",
      display = display,
      description = display,
      formulas = formulas,
      formula_name = self_name,
      node_id = node_id
    ))
    .math_add_edge(ctx, node_id, child_res$node_id, "arg", edge_type = "arg")
    return(list(node_id = node_id, symbol = symbol, display = display))
  }

  symbol <- .math_next_symbol(ctx, prefix = "X", hint = kind %||% "expr")
  node_id <- .math_next_id(ctx, "m")
  display <- .describe_expr_short(expr)
  label_lines <- c(symbol, toupper(kind %||% "expr"))
  .math_add_node(ctx, node_id, paste(label_lines, collapse = "\\n"), "operator")
  .math_register_symbol(ctx, symbol, list(
    kind = kind,
    label = toupper(kind %||% "expr"),
    display = display,
    description = display,
    formula_name = symbol,
    node_id = node_id
  ))
  list(node_id = node_id, symbol = symbol, display = display)
}

#' Build a mathematical processing tree describing likelihood contributions.
#'
#' @param model Race model (spec/tables).
#' @param outcome_label Outcome label of interest.
#' @param component Optional component (reserved for future use).
#' @param return_dot If TRUE, return DOT string.
#' @param mode Rendering mode: `"compact"` (short readable labels) or `"verbose"` (full formulas).
#' @return DiagrammeR graph object (if available) or list with DOT/nodes/edges.
processing_math_tree <- function(model, outcome_label, component = NULL, return_dot = FALSE,
                                 mode = c("compact", "verbose")) {
  mode <- match.arg(mode)
  if (is_model_tables(model)) {
    tables <- model
  } else {
    tables <- model_to_tables(model)
  }

  outcome_struct_id <- sprintf("outcome:%s", outcome_label)
  struct_row <- tables$struct_nodes[tables$struct_nodes$struct_id == outcome_struct_id, , drop = FALSE]
  if (nrow(struct_row) == 0) {
    stop(sprintf("Outcome '%s' not found in model tables", outcome_label))
  }
  outcome_expr <- (struct_row$payload[[1]] %||% list())$expr
  if (is.null(outcome_expr)) {
    stop(sprintf("Outcome '%s' is missing expression payload", outcome_label))
  }

  struct_lookup <- build_struct_lookup(tables$struct_nodes)
  ctx <- .math_new_ctx(struct_lookup, mode = mode)
  outcome_res <- .math_build_expr(ctx, outcome_expr)
  .math_assign_outcome_alias(ctx, outcome_res$symbol, outcome_label)
  ctx$outcome_map[[outcome_res$symbol]] <- outcome_label
  if (!(outcome_res$symbol %in% ctx$outcome_order)) {
    ctx$outcome_order <- c(ctx$outcome_order, outcome_res$symbol)
  }
  if (!is.null(ctx$node_info[[outcome_res$symbol]])) {
    info <- ctx$node_info[[outcome_res$symbol]]
    if (is.null(info$description)) info$description <- info$label
    info$label <- outcome_label
    info$display <- outcome_label
    ctx$node_info[[outcome_res$symbol]] <- info
  }

  outcome_rows <- tables$struct_nodes[tables$struct_nodes$node_type == "outcome", , drop = FALSE]
  competitor_rows <- outcome_rows[outcome_rows$label != outcome_label, , drop = FALSE]
  competitor_results <- list()
  if (nrow(competitor_rows) > 0) {
    competitor_results <- lapply(seq_len(nrow(competitor_rows)), function(i) {
      row <- competitor_rows[i, ]
      payload <- row$payload[[1]] %||% list()
      expr <- payload$expr
      if (is.null(expr)) return(NULL)
      res <- .math_build_expr(ctx, expr)
      res$outcome_label <- row$label %||% sprintf("outcome[%d]", i)
      .math_assign_outcome_alias(ctx, res$symbol, res$outcome_label)
      ctx$outcome_map[[res$symbol]] <- res$outcome_label
      if (!(res$symbol %in% ctx$outcome_order)) {
        ctx$outcome_order <- c(ctx$outcome_order, res$symbol)
      }
      if (!is.null(ctx$node_info[[res$symbol]])) {
        info <- ctx$node_info[[res$symbol]]
        if (is.null(info$description)) info$description <- info$label
        info$label <- res$outcome_label
        info$display <- res$outcome_label
        ctx$node_info[[res$symbol]] <- info
      }
      res
    })
    competitor_results <- competitor_results[!vapply(competitor_results, is.null, logical(1))]
  }

  root_id <- .math_next_id(ctx, "L")
  outcome_ref <- .math_symbol_ref(ctx, outcome_res$symbol, prefer_outcome_alias = TRUE)
  comp_refs <- if (length(competitor_results) > 0) {
    vapply(competitor_results, function(x) .math_symbol_ref(ctx, x$symbol, prefer_outcome_alias = TRUE), character(1))
  } else character(0)

  product_term <- if (length(comp_refs) > 0) {
    paste(sprintf("[1 - F_%s(t)]", comp_refs), collapse = " * ")
  } else ""
  root_symbol <- .math_unique_formula_name(ctx, paste0("L_", outcome_label))
  root_lines <- c(
    sprintf("%s(t)", root_symbol),
    if (nzchar(product_term)) {
      sprintf("= f_%s(t) * %s", outcome_ref, product_term)
    } else {
      sprintf("= f_%s(t)", outcome_ref)
    }
  )
  .math_add_node(ctx, root_id, paste(root_lines, collapse = "\\n"), "likelihood")
  .math_add_edge(
    ctx, root_id, outcome_res$node_id, sprintf("f_%s(t)", outcome_ref),
    edge_type = "likelihood_main"
  )
  if (length(competitor_results) > 0) {
    for (i in seq_along(competitor_results)) {
      comp_ref <- .math_symbol_ref(ctx, competitor_results[[i]]$symbol, prefer_outcome_alias = TRUE)
      edge_lab <- sprintf("1 - F_%s(t)", comp_ref)
      .math_add_edge(ctx, root_id, competitor_results[[i]]$node_id, edge_lab, edge_type = "likelihood_competitor")
    }
  }
  general_lines <- if (identical(mode, "verbose")) {
    c(
      "General:",
      "  f_E(t) : density of event E",
      "  F_E(t) : cumulative probability of E by t",
      "  S_E(t) : survival of E past t"
    )
  } else {
    c(
      "General",
      "f_E(t): density",
      "F_E(t): cumulative by t",
      "S_E(t): survival past t"
    )
  }
  general_id <- .math_next_id(ctx, "note")
  .math_add_node(ctx, general_id, paste(general_lines, collapse = "\\n"),
                 "definitions_general")

  legend_lines <- c("Definitions")
  if (length(ctx$outcome_order) > 0) {
    legend_lines <- c(legend_lines, "Outcomes:")
    for (sym in ctx$outcome_order) {
      sym_txt <- .math_symbol_ref(ctx, sym, prefer_outcome_alias = TRUE)
      label_txt <- ctx$outcome_map[[sym]] %||% "[unknown]"
      info <- ctx$node_info[[sym]]
      desc <- info$description %||% info$label %||% label_txt
      line <- if (!identical(label_txt, desc)) {
        sprintf("  %s = %s [%s]", sym_txt, desc, label_txt)
      } else {
        sprintf("  %s = %s", sym_txt, desc)
      }
      legend_lines <- c(legend_lines, line)
    }
  }
  other_syms <- setdiff(ctx$symbol_sequence, ctx$outcome_order)
  if (length(other_syms) > 0) {
    other_lines <- character(0)
    for (sym in other_syms) {
      info <- ctx$node_info[[sym]]
      if (is.null(info)) next
      if (identical(info$kind, "event")) {
        extras <- character(0)
        if (!is.null(info$type)) extras <- c(extras, info$type)
        members <- info$members %||% character(0)
        if (length(members) > 0) {
          extras <- c(extras, sprintf("members: %s", paste(members, collapse = ", ")))
        }
        rule <- info$rule %||% list()
        if (!is.null(rule$k) && !is.na(rule$k) && !is.null(rule$n) && rule$n > 0) {
          extras <- c(extras, sprintf("rule: %d-of-%d", rule$k, rule$n))
        }
        if (length(extras) > 0) {
          desc <- paste(extras, collapse = "; ")
        } else {
          desc <- info$description %||% info$label %||% info$display %||% "[unknown]"
        }
      } else {
        desc <- info$description %||% info$label %||% info$display %||% "[unknown]"
      }
      display_name <- .math_formula_name(ctx, sym)
      other_lines <- c(other_lines, sprintf("  %s = %s", display_name, desc))
    }
    if (length(other_lines) > 0) {
      legend_lines <- c(legend_lines, "Other nodes:", other_lines)
    }
  }
  legend_id <- .math_next_id(ctx, "note")
  .math_add_node(ctx, legend_id, paste(legend_lines, collapse = "\\n"), "definitions")
  .math_add_edge(ctx, root_id, legend_id, "legend", style = "dotted", edge_type = "legend")
  .math_add_edge(ctx, general_id, root_id, "", style = "invis", edge_type = "layout")

  target_nodes <- .collect_descendants(ctx$edges, outcome_res$node_id)
  competitor_root_ids <- if (length(competitor_results) > 0) {
    vapply(competitor_results, `[[`, character(1), "node_id")
  } else character(0)
  competitor_nodes <- .collect_descendants(ctx$edges, competitor_root_ids)
  shared_nodes <- intersect(target_nodes, competitor_nodes)

  operator_node_ids <- character(0)
  if (length(ctx$symbol_sequence) > 0) {
    for (sym in ctx$symbol_sequence) {
      info <- ctx$node_info[[sym]] %||% list()
      if ((info$kind %||% "") %in% c("and", "or", "guard", "not")) {
        nid <- info$node_id %||% NA_character_
        if (!is.na(nid) && nzchar(nid)) operator_node_ids <- c(operator_node_ids, nid)
      }
    }
  }
  helper_ids <- unique(c(shared_nodes, operator_node_ids, general_id, legend_id))
  target_only <- setdiff(target_nodes, c(helper_ids, competitor_nodes))
  competitor_only <- setdiff(competitor_nodes, c(helper_ids, target_nodes))
  .math_set_node_cluster(ctx, target_only, "target")
  .math_set_node_cluster(ctx, competitor_only, "competitor")
  .math_set_node_cluster(ctx, helper_ids, "helper")
  .math_set_node_cluster(ctx, root_id, "target")

  nodes <- ctx$nodes
  edges <- ctx$edges
  dot <- to_dot(nodes, edges)
  if (!return_dot && requireNamespace("DiagrammeR", quietly = TRUE)) {
    return(DiagrammeR::grViz(dot))
  }
  list(dot = dot, nodes = nodes, edges = edges)
}
