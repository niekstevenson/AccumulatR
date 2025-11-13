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

to_dot <- function(nodes, edges) {
  if (nrow(nodes) == 0) stop("No nodes to render")

  type_shape <- function(type) {
    switch(type,
           outcome = "box",
           guard = "octagon",
           and = "diamond",
           or = "diamond",
           not = "diamond",
           event = "ellipse",
           likelihood = "doublecircle",
           expr = "box",
           factor = "plaintext",
           competitor = "ellipse",
           definitions = "note",
           definitions_general = "note",
           "ellipse")
  }
  node_lines <- apply(nodes, 1, function(row) {
    sprintf('  %s [label="%s", shape=%s];',
            row[['id']], row[['label']], type_shape(row[['type']]))
  })
  edge_lines <- apply(edges, 1, function(row) {
    lbl <- row[['label']]
    style <- if ("style" %in% names(edges)) row[['style']] else ""
    attrs <- character(0)
    if (!is.null(lbl) && lbl != "") {
      attrs <- c(attrs, sprintf('label="%s"', lbl))
    }
    if (!is.null(style) && style != "") {
      attrs <- c(attrs, sprintf('style="%s"', style))
    }
    if (length(attrs) > 0) {
      sprintf("  %s -> %s [%s];", row[['from']], row[['to']], paste(attrs, collapse = ", "))
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

  paste0(
    "digraph processing_tree {\n",
    "  rankdir=TB;\n",
    "  node [fontname=\"Helvetica\"];\n",
    paste(node_lines, collapse = "\n"), "\n",
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

.math_new_ctx <- function(struct_lookup) {
  ctx <- new.env(parent = emptyenv())
  ctx$struct_lookup <- struct_lookup
  ctx$id_counter <- 0L
  ctx$symbol_counter <- 0L
  ctx$symbol_lookup <- list()
  ctx$node_lookup <- list()
  ctx$node_info <- list()
  ctx$symbol_sequence <- character(0)
  ctx$outcome_map <- list()
  ctx$outcome_order <- character(0)
  ctx$name_to_symbol <- list()
  ctx$nodes <- data.frame(
    id = character(0),
    label = character(0),
    type = character(0),
    rank = character(0),
    stringsAsFactors = FALSE
  )
  ctx$edges <- data.frame(
    from = character(0),
    to = character(0),
    label = character(0),
    stringsAsFactors = FALSE
  )
  ctx
}

.math_next_id <- function(ctx, prefix = "m") {
  ctx$id_counter <- ctx$id_counter + 1L
  sprintf("%s%d", prefix, ctx$id_counter)
}

.math_next_symbol <- function(ctx) {
  ctx$symbol_counter <- ctx$symbol_counter + 1L
  sprintf("E%d", ctx$symbol_counter)
}

.math_add_node <- function(ctx, id, label, type, rank = NA_character_) {
  ctx$nodes <- rbind(
    ctx$nodes,
    data.frame(id = id, label = label, type = type,
               rank = rank %||% NA_character_, stringsAsFactors = FALSE)
  )
}

.math_add_edge <- function(ctx, from, to, label = "", style = "") {
  ctx$edges <- rbind(
    ctx$edges,
    data.frame(from = from, to = to, label = label, style = style, stringsAsFactors = FALSE)
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
    symbol <- .math_next_symbol(ctx)
    node_id <- .math_next_id(ctx, "m")
    info <- ctx$struct_lookup[[source_id]] %||% list(node_type = "event", payload = list())
    display <- source_id
    type_txt <- info$node_type %||% "event"
    payload <- info$payload %||% list()
    members <- payload$members %||% character(0)
    rule <- payload$rule %||% list()
    event_label <- source_id
    base_formulas <- c(
      sprintf("f_%s(t) = f_%s(t)", event_label, event_label),
      sprintf("F_%s(t) = ∫ f_%s(u) du", event_label, event_label),
      sprintf("S_%s(t) = 1 - F_%s(t)", event_label, event_label)
    )
    formulas <- ._filter_formulas(base_formulas)
    if (length(formulas) == 0) formulas <- base_formulas[1]
    .math_add_node(ctx, node_id, event_label, "expr")
    ctx$symbol_lookup[[key]] <- symbol
    ctx$node_lookup[[key]] <- node_id
    .math_register_symbol(ctx, symbol, list(
      kind = "event",
      label = source_id,
      description = source_id,
      display = display,
      display_name = event_label,
      type = type_txt,
      members = members,
      rule = list(
        k = if (!is.null(rule$k)) rule$k else if (length(members) > 0) 1L else NA_integer_,
        n = length(members)
      ),
      formulas = formulas,
      formula_name = event_label,
      node_id = node_id
    ))
    return(list(node_id = node_id, symbol = symbol, display = display))
  }

  if (identical(kind, "and")) {
    symbol <- .math_next_symbol(ctx)
    node_id <- .math_next_id(ctx, "m")
    args <- expr$args %||% list()
    child_results <- lapply(seq_along(args), function(i) .math_build_expr(ctx, args[[i]]))
    child_syms <- if (length(child_results) > 0) vapply(child_results, `[[`, character(1), "symbol") else character(0)
    child_names <- if (length(child_syms) > 0) vapply(child_syms, function(s) ctx$node_info[[s]]$formula_name, character(1)) else character(0)
    child_disp <- if (length(child_results) > 0) vapply(child_results, `[[`, character(1), "display") else character(0)
    display <- if (length(child_disp) > 0) paste(child_disp, collapse = " & ") else "AND()"
    desc_text <- if (length(child_disp) > 0) sprintf("AND(%s)", paste(child_disp, collapse = ", ")) else "AND()"
    self_name <- ctx$node_info[[symbol]]$formula_name %||% symbol
    formulas <- c(.math_and_density_line(self_name, child_names), .math_and_cdf_line(self_name, child_names),
                  sprintf("S_%s(t) = 1 - F_%s(t)", self_name, self_name))
    formulas <- ._filter_formulas(formulas)
    if (length(formulas) == 0) formulas <- .math_and_density_line(self_name, child_names)
    .math_add_node(ctx, node_id, paste(formulas, collapse = "\\n"), "expr")
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
        .math_add_edge(ctx, node_id, child_results[[i]]$node_id, sprintf("arg[%d]", i))
      }
    }
    return(list(node_id = node_id, symbol = symbol, display = display))
  }

  if (identical(kind, "or")) {
    symbol <- .math_next_symbol(ctx)
    node_id <- .math_next_id(ctx, "m")
    args <- expr$args %||% list()
    child_results <- lapply(seq_along(args), function(i) .math_build_expr(ctx, args[[i]]))
    child_syms <- if (length(child_results) > 0) vapply(child_results, `[[`, character(1), "symbol") else character(0)
    child_names <- if (length(child_syms) > 0) vapply(child_syms, function(s) ctx$node_info[[s]]$formula_name %||% s, character(1)) else character(0)
    child_disp <- if (length(child_results) > 0) vapply(child_results, `[[`, character(1), "display") else character(0)
    display <- if (length(child_disp) > 0) paste(child_disp, collapse = " | ") else "OR()"
    desc_text <- if (length(child_disp) > 0) sprintf("OR(%s)", paste(child_disp, collapse = ", ")) else "OR()"
    self_name <- ctx$node_info[[symbol]]$formula_name %||% symbol
    formulas <- c(.math_or_density_line(self_name, child_names), .math_or_cdf_line(self_name, child_names),
                  sprintf("S_%s(t) = 1 - F_%s(t)", self_name, self_name))
    formulas <- ._filter_formulas(formulas)
    if (length(formulas) == 0) formulas <- .math_or_density_line(self_name, child_names)
    .math_add_node(ctx, node_id, paste(formulas, collapse = "\\n"), "expr")
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
        .math_add_edge(ctx, node_id, child_results[[i]]$node_id, sprintf("arg[%d]", i))
      }
    }
    return(list(node_id = node_id, symbol = symbol, display = display))
  }

  if (identical(kind, "guard")) {
    symbol <- .math_next_symbol(ctx)
    node_id <- .math_next_id(ctx, "m")
    reference <- expr$reference
    blocker <- expr$blocker
    unless_list <- expr$unless %||% list()
    ref_res <- .math_build_expr(ctx, reference)
    block_res <- .math_build_expr(ctx, blocker)
    unless_res <- lapply(seq_along(unless_list), function(i) .math_build_expr(ctx, unless_list[[i]]))
    unless_syms <- if (length(unless_res) > 0) vapply(unless_res, `[[`, character(1), "symbol") else character(0)
    display <- sprintf("%s guarded by %s", ref_res$display, block_res$display)
    self_name <- ctx$node_info[[symbol]]$formula_name %||% symbol
    ref_name <- ctx$node_info[[ref_res$symbol]]$formula_name %||% ref_res$symbol
    block_name <- ctx$node_info[[block_res$symbol]]$formula_name %||% block_res$symbol
    formulas <- c(
      sprintf("f_%s(t) = f_%s(t) * S_%s^{eff}(t)", self_name, ref_name, block_name),
      sprintf("F_%s(t) = ∫ f_%s(u) du", self_name, self_name),
      sprintf("S_%s(t) = 1 - F_%s(t)", self_name, self_name)
    )
    if (length(unless_syms) == 0) {
      formulas <- c(formulas, sprintf("S_%s^{eff}(t) = S_%s(t)", block_name, block_name))
    } else {
      protector_names <- vapply(unless_syms, function(s) ctx$node_info[[s]]$formula_name %||% s, character(1))
      protector_terms <- paste(sprintf("S_%s(u)", protector_names), collapse = " * ")
      formulas <- c(formulas, sprintf("S_%s^{eff}(t) = 1 - ∫ f_%s(u) * %s du",
                                      block_name, block_name, protector_terms))
    }
    formulas <- ._filter_formulas(formulas)
    if (length(formulas) == 0) formulas <- sprintf("f_%s(t) = f_%s(t)", self_name, ref_name)
    .math_add_node(ctx, node_id, paste(formulas, collapse = "\\n"), "expr")
    .math_register_symbol(ctx, symbol, list(
      kind = "guard",
      label = "GUARD",
      display = display,
      description = sprintf("guard(%s | !%s)", ref_res$display, block_res$display),
      formulas = formulas,
      formula_name = self_name,
      node_id = node_id
    ))
    .math_add_edge(ctx, node_id, ref_res$node_id, "reference")
    .math_add_edge(ctx, node_id, block_res$node_id, "blocker")
    if (length(unless_res) > 0) {
      for (i in seq_along(unless_res)) {
        .math_add_edge(ctx, node_id, unless_res[[i]]$node_id, sprintf("unless[%d]", i))
      }
    }
    return(list(node_id = node_id, symbol = symbol, display = display))
  }

  if (identical(kind, "not")) {
    symbol <- .math_next_symbol(ctx)
    node_id <- .math_next_id(ctx, "m")
    child <- expr$arg
    child_res <- .math_build_expr(ctx, child)
    display <- sprintf("not %s", child_res$display)
    formulas <- c(
      sprintf("f_%s(t) = 0", symbol),
      sprintf("F_%s(t) = 1 - F_%s(t)", symbol, child_res$symbol),
      sprintf("S_%s(t) = 1 - F_%s(t)", symbol, symbol)
    )
    formulas <- ._filter_formulas(formulas)
    if (length(formulas) == 0) formulas <- sprintf("f_%s(t) = 0", symbol)
    .math_add_node(ctx, node_id, paste(formulas, collapse = "\\n"), "expr")
    .math_register_symbol(ctx, symbol, list(
      kind = "not",
      label = "NOT",
      display = display,
      description = display,
      formulas = formulas
    ))
    .math_add_edge(ctx, node_id, child_res$node_id, "arg")
    return(list(node_id = node_id, symbol = symbol, display = display))
  }

  symbol <- .math_next_symbol(ctx)
  node_id <- .math_next_id(ctx, "m")
  display <- .describe_expr_short(expr)
  label_lines <- c(symbol, toupper(kind %||% "expr"))
  .math_add_node(ctx, node_id, paste(label_lines, collapse = "\\n"), "expr")
  .math_register_symbol(ctx, symbol, list(
    kind = kind,
    label = toupper(kind %||% "expr"),
    display = display,
    description = display
  ))
  list(node_id = node_id, symbol = symbol, display = display)
}

#' Build a mathematical processing tree describing likelihood contributions.
#'
#' @param model Race model (spec/tables).
#' @param outcome_label Outcome label of interest.
#' @param component Optional component (reserved for future use).
#' @param return_dot If TRUE, return DOT string.
#' @return DiagrammeR graph object (if available) or list with DOT/nodes/edges.
processing_math_tree <- function(model, outcome_label, component = NULL, return_dot = FALSE) {
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
  ctx <- .math_new_ctx(struct_lookup)
  outcome_res <- .math_build_expr(ctx, outcome_expr)
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
  comp_syms <- if (length(competitor_results) > 0) vapply(competitor_results, `[[`, character(1), "symbol") else character(0)

  product_term <- if (length(comp_syms) > 0) {
    paste(sprintf("[1 - F_%s(t)]", comp_syms), collapse = " * ")
  } else ""
  root_lines <- c(
    sprintf("ℒ_%s(t)", outcome_label),
    if (nzchar(product_term)) {
      sprintf("= f_%s(t) * %s", outcome_res$symbol, product_term)
    } else {
      sprintf("= f_%s(t)", outcome_res$symbol)
    }
  )
  .math_add_node(ctx, root_id, paste(root_lines, collapse = "\\n"), "likelihood")
  .math_add_edge(ctx, root_id, outcome_res$node_id, sprintf("f_%s(t)", outcome_res$symbol))
  if (length(competitor_results) > 0) {
    for (i in seq_along(competitor_results)) {
      edge_lab <- sprintf("1 - F_%s(t)", competitor_results[[i]]$symbol)
      .math_add_edge(ctx, root_id, competitor_results[[i]]$node_id, edge_lab)
    }
  }
  general_lines <- c(
    "General:",
    "  f_E(t) : density of event E",
    "  F_E(t) : cumulative probability of E by t",
    "  S_E(t) : survival of E past t"
  )
  general_id <- .math_next_id(ctx, "note")
  .math_add_node(ctx, general_id, paste(general_lines, collapse = "\\n"),
                 "definitions_general")

  legend_lines <- c("Definitions")
  if (length(ctx$outcome_order) > 0) {
    legend_lines <- c(legend_lines, "Outcomes:")
    for (sym in ctx$outcome_order) {
      label_txt <- ctx$outcome_map[[sym]] %||% "[unknown]"
      info <- ctx$node_info[[sym]]
      desc <- info$description %||% info$label %||% label_txt
      line <- if (!identical(label_txt, desc)) {
        sprintf("  %s = %s [%s]", sym, desc, label_txt)
      } else {
        sprintf("  %s = %s", sym, desc)
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
      display_name <- if (identical(info$kind, "event")) info$label %||% sym else sym
      other_lines <- c(other_lines, sprintf("  %s = %s", display_name, desc))
    }
    if (length(other_lines) > 0) {
      legend_lines <- c(legend_lines, "Other nodes:", other_lines)
    }
  }
  legend_id <- .math_next_id(ctx, "note")
  .math_add_node(ctx, legend_id, paste(legend_lines, collapse = "\\n"), "definitions")
  .math_add_edge(ctx, root_id, legend_id, "legend", style = "dotted")
  .math_add_edge(ctx, general_id, root_id, "", style = "invis")

  nodes <- ctx$nodes
  edges <- ctx$edges
  dot <- to_dot(nodes, edges)
  if (!return_dot && requireNamespace("DiagrammeR", quietly = TRUE)) {
    return(DiagrammeR::grViz(dot))
  }
  list(dot = dot, nodes = nodes, edges = edges)
}
