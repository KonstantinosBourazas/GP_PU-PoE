## GCN is the graph convolutional network of Kipf and Welling (2017). The original
## competitor code and the archived files call this method GNN.
scar_registry <- function() {
  data.frame(
    Method = c(
      "Lin-Cov",
      "Lin-PU-Cov",
      "Lin-Cov+Net",
      "Lin-PU-Cov+Net",
      "GP-Cov",
      "GP-PU-Cov",
      "GP-PoE",
      "GP-PU-PoE",
      "GP-PU-PoE_beta0",
      "Ada",
      "Ada+Net",
      "nnPU",
      "nnPU+Net",
      "WLR",
      "WLR+Net",
      "INLA",
      "INLA+Net",
      "GCN"
    ),
    source_method = c(
      "Bayesian Linear Cov",
      "Bayesian Linear Cov + PU",
      "Bayesian Linear Cov + Net",
      "Bayesian Linear Cov + Net + PU",
      "GP-Cov",
      "GP-Cov + PU",
      "PoE",
      "PU-PoE",
      "PU-PoE mu=0",
      "Ada",
      "Ada + Net",
      "nnPU",
      "nnPU + Net",
      "WLR",
      "WLR + Net",
      "INLA",
      "INLA + Net",
      "GNN"
    ),
    group = c(
      rep("Bayesian_linear", 4),
      "GP_Cov",
      "GP_Cov_PU",
      "GP_PoE",
      "GP_PU_PoE",
      "GP_PU_PoE_beta0",
      "Competitors",
      "Competitors",
      "nnPU",
      "nnPU",
      rep("Competitors", 5)
    ),
    stringsAsFactors = FALSE
  )
}

scar_groups <- function() {
  c(
    "Competitors",
    "Bayesian_linear",
    "GP_Cov",
    "GP_Cov_PU",
    "GP_PoE",
    "GP_PU_PoE",
    "GP_PU_PoE_beta0",
    "nnPU"
  )
}

scar_group_spec <- function() {
  data.frame(
    group = scar_groups(),
    prefix = c(
      "competing_methods",
      "Bayesian_linear",
      "GP_Cov",
      "GP_Cov_PU",
      "PoE",
      "PU_PoE",
      "PU_PoE_mu0",
      "nnpu"
    ),
    checkpoint_dir = c(
      "replicate_checkpoints",
      rep("production_checkpoints", 6),
      "replicate_checkpoints"
    ),
    pattern = c(
      "replicate_%03d.rds",
      "replicate_%03d.rds",
      "GP_Cov_replicate_%03d.rds",
      "GP_Cov_PU_replicate_%03d.rds",
      "PoE_replicate_%03d.rds",
      "PU_PoE_replicate_%03d.rds",
      "PU_PoE_mu0_replicate_%03d.rds",
      "replicate_%03d.rds"
    ),
    stringsAsFactors = FALSE
  )
}

scar_metrics <- function() {
  c(
    "Overall_AUC",
    "Hidden_AUC",
    "Overall_Recall_at_45",
    "Hidden_Recall_at_15",
    "Overall_LogLoss",
    "Hidden_LogLoss",
    "Overall_Coverage_80",
    "Hidden_Coverage_80",
    "Overall_Interval_Score_80",
    "Hidden_Interval_Score_80",
    "Overall_Brier",
    "Hidden_Brier"
  )
}

scar_canonical <- function(method) {
  r <- scar_registry()
  k <- match(as.character(method), r$source_method)
  scar_assert(
    !anyNA(k),
    paste(
      "Unrecognized original method label:",
      paste(unique(method[is.na(k)]), collapse = ", ")
    )
  )
  r$Method[k]
}
