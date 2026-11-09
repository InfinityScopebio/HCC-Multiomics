# Vendored from CoVarNet v0.3.0:
# https://github.com/QiangShiPKU/CoVarNet
# commit d4c52a4b2ea821dc4b75861c3772e21232d1ee51
# The implementation is namespaced and lightly rewritten to remove pipe and
# non-standard-evaluation dependencies while preserving numerical behavior.

vendor_freq_calculate <- function(meta) {
  required <- c("sampleID", "majorCluster", "subCluster")
  missing <- setdiff(required, colnames(meta))
  if (length(missing) > 0L) {
    stop("CoVarNet metadata missing: ", paste(missing, collapse = ", "))
  }
  state_annotation <- meta[
    !duplicated(meta$subCluster),
    c("majorCluster", "subCluster"),
    drop = FALSE
  ]
  state_annotation <- state_annotation[
    order(state_annotation$subCluster),
    ,
    drop = FALSE
  ]
  rownames(state_annotation) <- state_annotation$subCluster
  contingency <- table(meta$subCluster, meta$sampleID)
  counts <- as.data.frame(
    matrix(
      as.numeric(contingency),
      nrow = nrow(contingency),
      ncol = ncol(contingency),
      dimnames = dimnames(contingency)
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  state_annotation <- state_annotation[rownames(counts), , drop = FALSE]
  frequencies <- counts
  for (major in unique(state_annotation$majorCluster)) {
    state_rows <- which(state_annotation$majorCluster == major)
    denominator <- colSums(counts[state_rows, , drop = FALSE])
    frequencies[state_rows, ] <- sweep(
      counts[state_rows, , drop = FALSE],
      2L,
      denominator,
      "/"
    )
  }
  frequencies[is.na(frequencies)] <- 0
  data.frame(
    majorCluster = factor(state_annotation$majorCluster),
    frequencies,
    check.names = FALSE,
    row.names = rownames(frequencies)
  )
}

vendor_freq_normalize <- function(mat_fq, normalize = "minmax") {
  values <- as.matrix(mat_fq[, -1L, drop = FALSE])
  storage.mode(values) <- "numeric"
  if (normalize == "minmax") {
    values <- values - apply(values, 1L, min)
    values <- values / apply(values, 1L, max)
  } else if (normalize == "zscore") {
    values <- values - apply(values, 1L, mean)
    values <- values / apply(values, 1L, sd)
    values[values < 0] <- 0
  } else {
    stop("Unknown CoVarNet normalization: ", normalize)
  }
  values
}

vendor_pair_correlation <- function(mat_fq, method = "pearson") {
  state_annotation <- data.frame(
    majorCluster = as.character(mat_fq[[1L]]),
    subCluster = rownames(mat_fq),
    stringsAsFactors = FALSE
  )
  values <- as.matrix(mat_fq[, -1L, drop = FALSE])
  storage.mode(values) <- "numeric"
  correlation <- psych::corr.test(t(values), method = method)

  correlation_matrix <- correlation$r
  correlation_matrix[lower.tri(correlation_matrix, diag = TRUE)] <- NA
  correlation_table <- reshape2::melt(
    correlation_matrix,
    na.rm = TRUE,
    value.name = "correlation"
  )
  colnames(correlation_table)[1:2] <- c("subCluster1", "subCluster2")

  pvalue_matrix <- correlation$p
  pvalue_matrix[lower.tri(pvalue_matrix, diag = TRUE)] <- NA
  pvalue_table <- reshape2::melt(
    pvalue_matrix,
    na.rm = TRUE,
    value.name = "pval"
  )
  colnames(pvalue_table)[1:2] <- c("subCluster1", "subCluster2")
  pvalue_table$pval_fdr <- p.adjust(pvalue_table$pval, method = "fdr")

  full_correlation <- correlation$r
  specificity_matrix <- matrix(
    NA_real_,
    nrow(full_correlation),
    ncol(full_correlation),
    dimnames = dimnames(full_correlation)
  )
  state_count <- nrow(specificity_matrix)
  for (left_index in seq_len(state_count - 1L)) {
    for (right_index in (left_index + 1L):state_count) {
      background <- c(
        full_correlation[left_index, -left_index],
        full_correlation[-c(left_index, right_index), right_index]
      )
      specificity_matrix[left_index, right_index] <- mean(
        background <= full_correlation[left_index, right_index]
      )
    }
  }
  specificity_table <- reshape2::melt(
    specificity_matrix,
    na.rm = TRUE,
    value.name = "spe"
  )
  colnames(specificity_table)[1:2] <- c("subCluster1", "subCluster2")

  output <- Reduce(
    merge,
    list(correlation_table, pvalue_table, specificity_table)
  )
  output$majorCluster1 <- state_annotation$majorCluster[
    match(output$subCluster1, state_annotation$subCluster)
  ]
  output$majorCluster2 <- state_annotation$majorCluster[
    match(output$subCluster2, state_annotation$subCluster)
  ]
  output
}

vendor_cm_network <- function(
    nmf_result,
    correlation_result,
    top_n = 10L,
    corr = -2,
    pval_fdr = 0.05) {
  annotation <- rbind(
    setNames(
      correlation_result[, c("subCluster1", "majorCluster1")],
      c("subCluster", "majorCluster")
    ),
    setNames(
      correlation_result[, c("subCluster2", "majorCluster2")],
      c("subCluster", "majorCluster")
    )
  )
  annotation <- annotation[!duplicated(annotation$subCluster), , drop = FALSE]
  rownames(annotation) <- annotation$subCluster

  weights <- NMF::basis(nmf_result)
  weights <- weights[, order(colnames(weights)), drop = FALSE]
  weight_table <- reshape2::melt(weights)
  colnames(weight_table) <- c("subCluster", "cm", "weight")
  weight_table <- data.table::as.data.table(weight_table)
  data.table::setorder(weight_table, cm, -weight)
  weight_table <- weight_table[
    ,
    head(.SD, top_n),
    by = cm
  ]

  passes <- correlation_result$correlation > corr &
    correlation_result$pval_fdr <= pval_fdr
  edge_table <- correlation_result[passes, , drop = FALSE]
  state_count <- nrow(annotation)
  specificity_cutoff <- 1 -
    ((top_n - 1) * 2 - 1) / ((state_count - 1) * 2 - 1)
  edge_table <- edge_table[
    edge_table$spe >= specificity_cutoff,
    ,
    drop = FALSE
  ]
  graph <- igraph::graph_from_data_frame(edge_table, directed = FALSE)
  global_nodes <- data.frame(
    subCluster = igraph::V(graph)$name,
    stringsAsFactors = FALSE
  )
  global_nodes$majorCluster <- annotation$majorCluster[
    match(global_nodes$subCluster, annotation$subCluster)
  ]
  global_nodes <- global_nodes[order(global_nodes$subCluster), , drop = FALSE]

  node_each <- data.frame(
    cm = character(),
    subCluster = character(),
    stringsAsFactors = FALSE
  )
  edge_each <- data.frame()
  for (module in unique(weight_table$cm)) {
    selected_nodes <- intersect(
      igraph::V(graph)$name,
      weight_table[cm == module, subCluster]
    )
    sub_graph <- igraph::induced_subgraph(graph, vids = selected_nodes)
    isolates <- igraph::V(sub_graph)[igraph::degree(sub_graph) == 0]
    if (length(isolates) > 0L) {
      sub_graph <- igraph::delete_vertices(sub_graph, isolates)
    }
    if (igraph::vcount(sub_graph) == 0L) {
      next
    }
    module_nodes <- data.frame(
      cm = as.character(module),
      subCluster = igraph::V(sub_graph)$name,
      stringsAsFactors = FALSE
    )
    node_each <- rbind(node_each, module_nodes)
    module_edges <- edge_table[
      edge_table$subCluster1 %in% module_nodes$subCluster &
        edge_table$subCluster2 %in% module_nodes$subCluster,
      ,
      drop = FALSE
    ]
    edge_each <- rbind(
      edge_each,
      data.frame(cm = as.character(module), module_edges)
    )
  }
  node_each$majorCluster <- annotation$majorCluster[
    match(node_each$subCluster, annotation$subCluster)
  ]
  filtered_weights <- merge(node_each, weight_table)
  filtered_weights <- filtered_weights[
    order(filtered_weights$cm, -filtered_weights$weight),
    ,
    drop = FALSE
  ]
  rownames(filtered_weights) <- NULL
  rownames(edge_each) <- NULL

  list(
    global = list(node = global_nodes, edge = edge_table),
    each = list(node = node_each, edge = edge_each),
    raw = nmf_result,
    filter = filtered_weights,
    ann = annotation
  )
}
