## Extracted from GP_four_model_illustrtation.R. See docs/SOURCE_MAP.csv.
## Original mathematical operations and their order are retained.

choose_q_modularity_eigengap <- function(
    values,
    tol = MOD_EIG_TOL,
    max_q = MOD_MAX_Q
) {
  positive_values <- sort(
    as.numeric(values[values > tol]),
    decreasing = TRUE
  )

  n_positive <- length(positive_values)

  if (n_positive == 0L) {
    return(list(
      q = 1L,
      positive_values = numeric(0),
      gaps = numeric(0),
      max_gap = NA_real_
    ))
  }

  if (n_positive == 1L) {
    return(list(
      q = 1L,
      positive_values = positive_values,
      gaps = numeric(0),
      max_gap = NA_real_
    ))
  }

  gaps <- positive_values[-n_positive] - positive_values[-1L]
  q <- which.max(gaps)

  if (is.finite(max_q)) {
    q <- min(q, as.integer(max_q))
  }

  q <- max(1L, min(q, n_positive))

  list(
    q = as.integer(q),
    positive_values = positive_values,
    gaps = gaps,
    max_gap = gaps[which.max(gaps)]
  )
}

orient_eigenvectors_deterministically <- function(V) {
  V <- as.matrix(V)
  for (column_index in seq_len(ncol(V))) {
    largest_loading_index <- which.max(abs(V[, column_index]))
    if (V[largest_loading_index, column_index] < 0) {
      V[, column_index] <- -V[, column_index]
    }
  }
  V
}

modularity_features_single_network <- function(
    A,
    seed,
    eig_tol = MOD_EIG_TOL,
    max_q = MOD_MAX_Q,
    k_max = MOD_K_MAX
) {
  set.seed(seed)
  A <- sanitize_network_adjacency(A)
  n <- nrow(A)
  degree <- rowSums(A)
  n_edges <- sum(degree) / 2

  if (!is.finite(n_edges) || n_edges <= 0) {
    stop("Cannot construct modularity features from an edgeless network.")
  }

  modularity_operator <- function(x, args) {
    as.vector(
      A %*% x -
        degree * sum(degree * x) / (2 * n_edges)
    )
  }

  k_use <- min(as.integer(k_max), n - 2L)

  eig <- RSpectra::eigs_sym(
    A = modularity_operator,
    k = k_use,
    n = n,
    which = "LA",
    opts = list(retvec = TRUE)
  )

  eigenvalues <- Re(eig$values)
  eigenvectors <- Re(eig$vectors)
  ordering <- order(eigenvalues, decreasing = TRUE)
  eigenvalues <- eigenvalues[ordering]
  eigenvectors <- eigenvectors[, ordering, drop = FALSE]
  eigenvectors <- orient_eigenvectors_deterministically(eigenvectors)

  positive_index <- which(eigenvalues > eig_tol)

  if (length(positive_index) == 0L) {
    keep <- 1L
  } else {
    gap_information <- choose_q_modularity_eigengap(
      values = eigenvalues,
      tol = eig_tol,
      max_q = max_q
    )
    keep <- positive_index[seq_len(gap_information$q)]
  }

  Z_raw <- eigenvectors[, keep, drop = FALSE]
  colnames(Z_raw) <- paste0("network_mod", seq_len(ncol(Z_raw)))

  standardized <- standardize_columns(Z_raw)
  Z <- standardized$scaled
  colnames(Z) <- colnames(Z_raw)

  list(
    Z = Z,
    raw = Z_raw,
    q = ncol(Z),
    selected_eigenvalues = eigenvalues[keep],
    all_eigenvalues = eigenvalues,
    center = standardized$center,
    scale = standardized$scale
  )
}

prepare_network_components_for_blocks <- function(
    A,
    network_name = "network"
) {
  A <- sanitize_network_adjacency(A)
  n <- nrow(A)

  g <- igraph::graph_from_adjacency_matrix(
    A,
    mode = "undirected",
    diag = FALSE
  )

  comp <- igraph::components(g)$membership
  comp_ids <- sort(unique(comp))
  components <- vector("list", length(comp_ids))

  for (cc_pos in seq_along(comp_ids)) {
    cc <- comp_ids[cc_pos]
    idx <- which(comp == cc)
    s <- length(idx)

    object <- list(idx = idx, size = s)

    if (s == 1L) {
      object$type <- "isolate"
      components[[cc_pos]] <- object
      next
    }

    A_sub <- A[idx, idx, drop = FALSE]
    edge_total <- sum(A_sub) / 2
    is_clique <- abs(edge_total - s * (s - 1) / 2) < 1e-8

    if (is_clique) {
      object$type <- "clique"
      object$lambda_orth <- s / (s - 1)
      components[[cc_pos]] <- object
      next
    }

    degree <- rowSums(A_sub)
    inverse_sqrt_degree <- rep(0, s)
    positive_degree <- which(degree > 0)
    inverse_sqrt_degree[positive_degree] <-
      1 / sqrt(degree[positive_degree])

    D_inverse <- diag(inverse_sqrt_degree, s)
    L_normalized <- diag(1, s) -
      D_inverse %*% A_sub %*% D_inverse
    L_normalized <- symmetrize(L_normalized)

    eig <- eigen(L_normalized, symmetric = TRUE)

    object$type <- "general"
    object$U <- Re(eig$vectors)
    object$U2 <- object$U^2
    object$lambda <- pmax(Re(eig$values), 0)
    components[[cc_pos]] <- object
  }

  component_sizes <- vapply(components, `[[`, integer(1), "size")
  component_types <- vapply(components, `[[`, character(1), "type")

  list(
    name = network_name,
    A = A,
    components = components,
    component_membership = comp,
    component_summary = data.frame(
      network = network_name,
      n_components = length(components),
      max_component_size = max(component_sizes),
      n_isolates = sum(component_types == "isolate"),
      n_cliques = sum(component_types == "clique"),
      n_general = sum(component_types == "general"),
      stringsAsFactors = FALSE
    )
  )
}

network_kernel_from_geom <- function(geom, ell, sigma) {
  n <- nrow(geom$A)
  K <- matrix(0, n, n)
  sigma2 <- sigma^2

  for (component in geom$components) {
    idx <- component$idx
    s <- component$size

    if (component$type == "isolate") {
      q <- exp(-1 / (2 * ell^2))
      K[idx, idx] <- K[idx, idx] + sigma2 * q
      next
    }

    if (component$type == "clique") {
      q <- exp(-component$lambda_orth / (2 * ell^2))
      off_diagonal <- sigma2 * (1 - q) / s
      K[idx, idx] <- K[idx, idx] + off_diagonal
      K[cbind(idx, idx)] <- K[cbind(idx, idx)] + sigma2 * q
      next
    }

    gamma <- sigma2 * exp(-component$lambda / (2 * ell^2))
    U_scaled <- sweep(component$U, 2L, sqrt(gamma), "*")
    K_component <- symmetrize(tcrossprod(U_scaled))
    K[idx, idx] <- K[idx, idx] + K_component
  }

  symmetrize(K)
}

