`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

SCENARIOS <- c("M0", "M1", "M2", "M3", "M4", "M5")
ALLELES <- matrix(c("A", "G", "C", "T", "G", "A", "T", "C"), ncol = 2, byrow = TRUE)

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

load_config <- function(path) {
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

benchmark_scenarios <- function(config, benchmark) {
  bench_cfg <- config[[benchmark]] %||% list()
  scenarios <- bench_cfg[["scenarios"]]
  if (is.null(scenarios) || length(scenarios) == 0) {
    return(SCENARIOS)
  }
  as.character(unlist(scenarios, use.names = FALSE))
}

cell_scenarios <- function(config, benchmark, cell) {
  scenarios <- cell[["scenarios"]]
  if (is.null(scenarios) || length(scenarios) == 0) {
    return(benchmark_scenarios(config, benchmark))
  }
  as.character(unlist(scenarios, use.names = FALSE))
}

scenario_alias_map <- function(config, benchmark) {
  bench_cfg <- config[[benchmark]] %||% list()
  aliases <- bench_cfg[["scenario_aliases"]]
  if (is.null(aliases)) {
    return(character())
  }
  stats::setNames(as.character(unname(unlist(aliases))), names(aliases))
}

effective_scenario <- function(name, alias_map) {
  if (length(name) == 0 || is.na(name[[1]])) {
    return(NA_character_)
  }
  if (!(name %in% names(alias_map))) {
    return(name)
  }
  mapped <- alias_map[[name]]
  if (is.null(mapped) || identical(mapped, character(0))) name else mapped
}

norm_sf <- function(z) {
  stats::pnorm(abs(z), lower.tail = FALSE)
}

choose_alleles <- function(n) {
  ALLELES[sample.int(nrow(ALLELES), size = n, replace = TRUE), , drop = FALSE]
}

write_table <- function(path, rows, fieldnames, delimiter) {
  if (is.null(rows) || nrow(rows) == 0) {
    rows <- as.data.frame(setNames(replicate(length(fieldnames), character(0), simplify = FALSE), fieldnames))
  }
  rows <- rows[, fieldnames, drop = FALSE]
  for (nm in names(rows)) {
    if (is.list(rows[[nm]])) {
      rows[[nm]] <- vapply(
        rows[[nm]],
        function(x) {
          if (length(x) == 0 || all(is.na(x))) {
            ""
          } else if (length(x) == 1) {
            as.character(x)
          } else {
            paste(as.character(x), collapse = ",")
          }
        },
        character(1)
      )
    }
  }
  utils::write.table(rows, file = path, sep = delimiter, quote = FALSE, row.names = FALSE, col.names = TRUE)
}

read_table <- function(path, delimiter) {
  df <- utils::read.table(
    path,
    header = TRUE,
    sep = delimiter,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (nm in names(df)) {
    df[[nm]] <- utils::type.convert(df[[nm]], as.is = TRUE)
  }
  df
}

scenario_parameters <- function(cell, scenario) {
  params <- list(
    beta1 = if (scenario %in% c("M1", "M2", "M5")) as.numeric(cell[["beta1"]]) else 0.0,
    beta2 = if (scenario %in% c("M1", "M4")) as.numeric(cell[["beta2"]]) else 0.0,
    beta3 = if (scenario %in% c("M1", "M3", "M5")) as.numeric(cell[["beta3"]] %||% 0.0) else 0.0,
    delta_sd = as.numeric(cell[["delta_sd"]]),
    psi_sd = as.numeric(cell[["psi_sd"]]),
    phi_sd = as.numeric(cell[["phi_sd"]]),
    rho_m5 = as.numeric(cell[["rho_m5"]] %||% cell[["rho_m3"]] %||% 0.0),
    direction_flip_rate = as.numeric(cell[["direction_flip_rate"]] %||% 0.0)
  )
  overrides <- cell[["scenario_overrides"]][[scenario]]
  if (!is.null(overrides)) {
    for (nm in names(overrides)) {
      params[[nm]] <- as.numeric(overrides[[nm]])
    }
  }
  params
}

sample_set_counts <- function(benchmark, cell, scenario) {
  scenario_counts <- cell[["scenario_set_counts"]][[scenario]]
  if (!is.null(scenario_counts)) {
    return(c(A = as.integer(scenario_counts[["A"]]), B = as.integer(scenario_counts[["B"]]), C = as.integer(scenario_counts[["C"]])))
  }
  scenario_dist <- cell[["scenario_set_count_distribution"]][[scenario]]
  if (!is.null(scenario_dist)) {
    return(c(
      A = sample(as.integer(unlist(scenario_dist[["A"]])), 1),
      B = sample(as.integer(unlist(scenario_dist[["B"]])), 1),
      C = sample(as.integer(unlist(scenario_dist[["C"]])), 1)
    ))
  }
  if (identical(benchmark, "classification")) {
    sets <- cell[["set_counts"]]
    return(c(A = as.integer(sets[["A"]]), B = as.integer(sets[["B"]]), C = as.integer(sets[["C"]])))
  }
  dist <- cell[["set_count_distribution"]]
  c(
    A = sample(as.integer(unlist(dist[["A"]])), 1),
    B = sample(as.integer(unlist(dist[["B"]])), 1),
    C = sample(as.integer(unlist(dist[["C"]])), 1)
  )
}

rmvnorm_simple <- function(n, sigma) {
  if (n <= 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 2))
  }
  z <- matrix(stats::rnorm(n * 2), ncol = 2)
  z %*% chol(sigma)
}

simulate_effects <- function(scenario, cell, n_a, n_b, n_c) {
  params <- scenario_parameters(cell, scenario)
  beta1 <- params[["beta1"]]
  beta2 <- params[["beta2"]]
  beta3 <- params[["beta3"]]
  delta_sd <- params[["delta_sd"]]
  psi_sd <- params[["psi_sd"]]
  phi_sd <- params[["phi_sd"]]
  rho <- params[["rho_m5"]]
  direction_flipped <- identical(scenario, "M1") && stats::runif(1) < params[["direction_flip_rate"]]
  if (direction_flipped) {
    beta2 <- -beta2
  }
  gamma_range <- as.numeric(unlist(cell[["gamma_range"]]))
  cis_range <- as.numeric(unlist(cell[["cis_effect_range"]]))

  gamma_a <- if (n_a > 0) stats::runif(n_a, gamma_range[1], gamma_range[2]) else numeric()
  gamma_c <- if (n_c > 0) stats::runif(n_c, gamma_range[1], gamma_range[2]) else numeric()
  cis_b <- if (n_b > 0) stats::runif(n_b, cis_range[1], cis_range[2]) else numeric()
  cis_c <- if (n_c > 0) stats::runif(n_c, cis_range[1], cis_range[2]) else numeric()

  if (identical(scenario, "M5")) {
    cov_mat <- matrix(
      c(delta_sd^2, rho * delta_sd * psi_sd, rho * delta_sd * psi_sd, psi_sd^2),
      nrow = 2,
      byrow = TRUE
    )
    joint <- rmvnorm_simple(n_a + n_c, cov_mat)
    delta_ac <- if (nrow(joint) > 0) joint[, 1] else numeric()
    psi_ac <- if (nrow(joint) > 0) joint[, 2] else numeric()
  } else {
    delta_ac <- if ((n_a + n_c) > 0) stats::rnorm(n_a + n_c, sd = delta_sd) else numeric()
    psi_ac <- if ((n_a + n_c) > 0) stats::rnorm(n_a + n_c, sd = psi_sd) else numeric()
  }

  phi_b <- if (n_b > 0) stats::rnorm(n_b, sd = phi_sd) else numeric()
  delta_a <- if (n_a > 0) delta_ac[seq_len(n_a)] else numeric()
  delta_c <- if (n_c > 0) delta_ac[seq.int(n_a + 1, n_a + n_c)] else numeric()
  psi_a <- if (n_a > 0) psi_ac[seq_len(n_a)] else numeric()
  psi_c <- if (n_c > 0) psi_ac[seq.int(n_a + 1, n_a + n_c)] else numeric()

  alpha_a <- beta1 * gamma_a + delta_a
  alpha_b <- cis_b
  alpha_c <- beta1 * gamma_c + cis_c + delta_c
  gamma_b <- if (n_b > 0) stats::rnorm(n_b, sd = 0.005) else numeric()

  gamma_out_a <- (beta2 * beta1 + beta3) * gamma_a + beta2 * delta_a + psi_a
  gamma_out_b <- beta2 * cis_b + phi_b
  gamma_out_c <- (beta2 * beta1 + beta3) * gamma_c + beta2 * (cis_c + delta_c) + psi_c

  if (identical(scenario, "M0")) {
    alpha_a[] <- 0.0
    alpha_b[] <- 0.0
    alpha_c[] <- 0.0
    gamma_out_a[] <- 0.0
    gamma_out_b[] <- 0.0
    gamma_out_c[] <- 0.0
  } else if (identical(scenario, "M2")) {
    gamma_out_a <- psi_a
    gamma_out_b <- phi_b
    gamma_out_c <- psi_c
  } else if (identical(scenario, "M3")) {
    alpha_a[] <- 0.0
    alpha_b[] <- 0.0
    alpha_c[] <- 0.0
    gamma_out_a <- beta3 * gamma_a + psi_a
    gamma_out_b <- phi_b
    gamma_out_c <- beta3 * gamma_c + psi_c
  } else if (identical(scenario, "M4")) {
    alpha_a[] <- 0.0
    alpha_c[] <- 0.0
    gamma_out_a <- psi_a
    gamma_out_b <- beta2 * cis_b + phi_b
    gamma_out_c <- beta2 * cis_c + psi_c
  } else if (identical(scenario, "M5")) {
    gamma_out_a <- beta3 * gamma_a + psi_a
    gamma_out_b <- phi_b
    gamma_out_c <- beta3 * gamma_c + psi_c
  }

  list(
    true_beta1 = beta1,
    true_beta2 = beta2,
    true_beta3 = beta3,
    direction_flipped = direction_flipped,
    gamma_a = gamma_a,
    gamma_b = gamma_b,
    gamma_c = gamma_c,
    alpha_a = alpha_a,
    alpha_b = alpha_b,
    alpha_c = alpha_c,
    Gamma_a = gamma_out_a,
    Gamma_b = gamma_out_b,
    Gamma_c = gamma_out_c
  )
}

generate_proteins <- function(benchmark, cell, scenario_sequence) {
  if (isTRUE(cell[["shared_rf_panel"]] %||% FALSE)) {
    return(generate_shared_panel_proteins(benchmark, cell, scenario_sequence))
  }

  proteins <- vector("list", length(scenario_sequence))
  for (idx in seq_along(scenario_sequence)) {
    scenario <- scenario_sequence[[idx]]
    set_counts <- sample_set_counts(benchmark, cell, scenario)
    n_a <- unname(set_counts[["A"]])
    n_b <- unname(set_counts[["B"]])
    n_c <- unname(set_counts[["C"]])
    effects <- simulate_effects(scenario, cell, n_a, n_b, n_c)
    chr_num <- sample.int(22, 1)
    base <- 10000000 + idx * 50000
    gene_start <- base
    gene_end <- base + 20000
    allele_mat <- choose_alleles(n_a + n_b + n_c)
    freqs <- stats::runif(n_a + n_b + n_c, 0.08, 0.45)
    rf_se <- as.numeric(cell[["rf_se"]])
    pqtl_se <- as.numeric(cell[["pqtl_se"]])
    outcome_se <- as.numeric(cell[["outcome_se"]])
    cursor <- 1L
    snps <- list(A = list(), B = list(), C = list())

    for (set_name in c("A", "B", "C")) {
      count <- switch(set_name, A = n_a, B = n_b, C = n_c)
      if (count == 0) {
        next
      }
      set_rows <- vector("list", count)
      for (j in seq_len(count)) {
        a1 <- allele_mat[cursor, 1]
        a2 <- allele_mat[cursor, 2]
        freq <- round(as.numeric(freqs[cursor]), 4)
        if (identical(set_name, "A")) {
          rf_beta <- effects[["gamma_a"]][j]
          pqtl_beta <- effects[["alpha_a"]][j]
          outcome_beta <- effects[["Gamma_a"]][j]
          bp <- gene_start - 2000000 - j * 2500
        } else if (identical(set_name, "B")) {
          rf_beta <- effects[["gamma_b"]][j]
          pqtl_beta <- effects[["alpha_b"]][j]
          outcome_beta <- effects[["Gamma_b"]][j]
          bp <- gene_start + 2000 + j * 1200
        } else {
          rf_beta <- effects[["gamma_c"]][j]
          pqtl_beta <- effects[["alpha_c"]][j]
          outcome_beta <- effects[["Gamma_c"]][j]
          bp <- gene_start + 9000 + j * 1400
        }
        set_rows[[j]] <- data.frame(
          rsid = sprintf("rs%s_%d_%s%d", tolower(scenario), idx, tolower(set_name), j),
          a1 = a1,
          a2 = a2,
          freq = freq,
          rf_beta = rf_beta,
          pqtl_beta = pqtl_beta,
          outcome_beta = outcome_beta,
          rf_p = 2.0 * norm_sf(rf_beta / rf_se),
          pqtl_p = 2.0 * norm_sf(pqtl_beta / pqtl_se),
          outcome_p = 2.0 * norm_sf(outcome_beta / outcome_se),
          bp = as.integer(bp),
          stringsAsFactors = FALSE
        )
        cursor <- cursor + 1L
      }
      snps[[set_name]] <- do.call(rbind, set_rows)
    }

    params <- scenario_parameters(cell, scenario)
    proteins[[idx]] <- list(
      protein_id = sprintf("%s_P%04d", scenario, idx),
      gene_name = sprintf("GENE_%s_%04d", scenario, idx),
      scenario = scenario,
      chr = chr_num,
      gene_start = gene_start,
      gene_end = gene_end,
      rf_se = rf_se,
      pqtl_se = pqtl_se,
      outcome_se = outcome_se,
      true_beta1 = effects[["true_beta1"]],
      true_beta2 = effects[["true_beta2"]],
      true_beta3 = effects[["true_beta3"]],
      direction_flipped = effects[["direction_flipped"]],
      set_counts = set_counts,
      snps = snps
    )
  }
  proteins
}

generate_shared_panel_proteins <- function(benchmark, cell, scenario_sequence) {
  n_shared <- as.integer(cell[["shared_rf_count"]] %||% 24L)
  gamma_range <- as.numeric(unlist(cell[["gamma_range"]]))
  shared_gamma <- stats::runif(n_shared, gamma_range[1], gamma_range[2])
  shared_rf_outcome_beta <- as.numeric(cell[["shared_rf_outcome_beta"]] %||% 0.0)
  shared_outcome_sd <- as.numeric(cell[["shared_rf_outcome_sd"]] %||% 0.0)
  shared_outcome <- shared_rf_outcome_beta * shared_gamma +
    if (n_shared > 0) stats::rnorm(n_shared, sd = shared_outcome_sd) else numeric()
  shared_alleles <- choose_alleles(n_shared)
  shared_freqs <- stats::runif(n_shared, 0.08, 0.45)
  rf_se <- as.numeric(cell[["rf_se"]])
  pqtl_se <- as.numeric(cell[["pqtl_se"]])
  outcome_se <- as.numeric(cell[["outcome_se"]])
  shared_rows <- vector("list", n_shared)
  for (i in seq_len(n_shared)) {
    shared_rows[[i]] <- data.frame(
      rsid = sprintf("rsrf_shared_%04d", i),
      a1 = shared_alleles[i, 1],
      a2 = shared_alleles[i, 2],
      freq = round(as.numeric(shared_freqs[i]), 4),
      rf_beta = shared_gamma[i],
      pqtl_beta = NA_real_,
      outcome_beta = shared_outcome[i],
      rf_p = 2.0 * norm_sf(shared_gamma[i] / rf_se),
      pqtl_p = NA_real_,
      outcome_p = 2.0 * norm_sf(shared_outcome[i] / outcome_se),
      chr = 22L,
      bp = as.integer(1000000 + i * 2500000),
      stringsAsFactors = FALSE
    )
  }
  shared_template <- do.call(rbind, shared_rows)

  proteins <- vector("list", length(scenario_sequence))
  for (idx in seq_along(scenario_sequence)) {
    scenario <- scenario_sequence[[idx]]
    set_counts <- sample_set_counts(benchmark, cell, scenario)
    n_b <- unname(set_counts[["B"]])
    n_c <- unname(set_counts[["C"]])
    params <- scenario_parameters(cell, scenario)
    beta1 <- params[["beta1"]]
    beta2 <- params[["beta2"]]
    beta3 <- params[["beta3"]]
    direction_flipped <- identical(scenario, "M1") && stats::runif(1) < params[["direction_flip_rate"]]
    if (direction_flipped) {
      beta2 <- -beta2
    }
    delta_sd <- params[["delta_sd"]]
    psi_sd <- params[["psi_sd"]]
    phi_sd <- params[["phi_sd"]]
    rho <- params[["rho_m5"]]
    cis_range <- as.numeric(unlist(cell[["cis_effect_range"]]))
    chr_num <- 1L + ((idx - 1L) %% 21L)
    base <- 10000000 + idx * 50000
    gene_start <- base
    gene_end <- base + 20000
    snps <- list(A = list(), B = list(), C = list())

    total_ac <- n_shared + n_c
    if (identical(scenario, "M5") && total_ac > 0) {
      cov_mat <- matrix(
        c(delta_sd^2, rho * delta_sd * psi_sd, rho * delta_sd * psi_sd, psi_sd^2),
        nrow = 2,
        byrow = TRUE
      )
      joint <- rmvnorm_simple(total_ac, cov_mat)
      delta_ac <- joint[, 1]
      psi_ac <- joint[, 2]
    } else {
      delta_ac <- if (total_ac > 0) stats::rnorm(total_ac, sd = delta_sd) else numeric()
      psi_ac <- if (total_ac > 0) stats::rnorm(total_ac, sd = psi_sd) else numeric()
    }
    delta_shared <- if (n_shared > 0) delta_ac[seq_len(n_shared)] else numeric()
    psi_shared <- if (n_shared > 0) psi_ac[seq_len(n_shared)] else numeric()
    delta_c <- if (n_c > 0) delta_ac[seq.int(n_shared + 1L, n_shared + n_c)] else numeric()
    psi_c <- if (n_c > 0) psi_ac[seq.int(n_shared + 1L, n_shared + n_c)] else numeric()

    shared_for_protein <- shared_template
    shared_for_protein$pqtl_beta <- 0.0
    shared_for_protein$pqtl_p <- 1.0
    if (scenario %in% c("M1", "M2", "M5")) {
      shared_for_protein$pqtl_beta <- beta1 * shared_gamma + delta_shared
      shared_for_protein$pqtl_p <- 2.0 * norm_sf(shared_for_protein$pqtl_beta / pqtl_se)
    }
    snps$A <- shared_for_protein

    if (n_b > 0) {
      allele_b <- choose_alleles(n_b)
      cis_b <- stats::runif(n_b, cis_range[1], cis_range[2])
      phi_b <- stats::rnorm(n_b, sd = phi_sd)
      b_rows <- vector("list", n_b)
      for (j in seq_len(n_b)) {
        pqtl_beta <- if (scenario %in% c("M1", "M2", "M4", "M5")) cis_b[j] else 0.0
        outcome_beta <- if (scenario %in% c("M1", "M4")) beta2 * cis_b[j] + phi_b[j] else phi_b[j]
        b_rows[[j]] <- data.frame(
          rsid = sprintf("rs%s_%d_b%d", tolower(scenario), idx, j),
          a1 = allele_b[j, 1],
          a2 = allele_b[j, 2],
          freq = round(stats::runif(1, 0.08, 0.45), 4),
          rf_beta = 0.0,
          pqtl_beta = pqtl_beta,
          outcome_beta = outcome_beta,
          rf_p = 1.0,
          pqtl_p = 2.0 * norm_sf(pqtl_beta / pqtl_se),
          outcome_p = 2.0 * norm_sf(outcome_beta / outcome_se),
          chr = chr_num,
          bp = as.integer(gene_start + 2000 + j * 1200),
          stringsAsFactors = FALSE
        )
      }
      snps$B <- do.call(rbind, b_rows)
    }

    if (n_c > 0) {
      allele_c <- choose_alleles(n_c)
      gamma_c <- stats::runif(n_c, gamma_range[1], gamma_range[2])
      cis_c <- stats::runif(n_c, cis_range[1], cis_range[2])
      c_rows <- vector("list", n_c)
      for (j in seq_len(n_c)) {
        if (scenario %in% c("M1", "M2", "M5")) {
          pqtl_beta <- beta1 * gamma_c[j] + cis_c[j] + delta_c[j]
        } else if (identical(scenario, "M4")) {
          pqtl_beta <- cis_c[j]
        } else {
          pqtl_beta <- 0.0
        }
        if (identical(scenario, "M1")) {
          outcome_beta <- (beta2 * beta1 + beta3) * gamma_c[j] + beta2 * (cis_c[j] + delta_c[j]) + psi_c[j]
        } else if (identical(scenario, "M3")) {
          outcome_beta <- beta3 * gamma_c[j] + psi_c[j]
        } else if (identical(scenario, "M4")) {
          outcome_beta <- beta2 * cis_c[j] + psi_c[j]
        } else if (identical(scenario, "M5")) {
          outcome_beta <- beta3 * gamma_c[j] + psi_c[j]
        } else {
          outcome_beta <- psi_c[j]
        }
        c_rows[[j]] <- data.frame(
          rsid = sprintf("rs%s_%d_c%d", tolower(scenario), idx, j),
          a1 = allele_c[j, 1],
          a2 = allele_c[j, 2],
          freq = round(stats::runif(1, 0.08, 0.45), 4),
          rf_beta = gamma_c[j],
          pqtl_beta = pqtl_beta,
          outcome_beta = outcome_beta,
          rf_p = 2.0 * norm_sf(gamma_c[j] / rf_se),
          pqtl_p = 2.0 * norm_sf(pqtl_beta / pqtl_se),
          outcome_p = 2.0 * norm_sf(outcome_beta / outcome_se),
          chr = chr_num,
          bp = as.integer(gene_start + 9000 + j * 1400),
          stringsAsFactors = FALSE
        )
      }
      snps$C <- do.call(rbind, c_rows)
    }

    proteins[[idx]] <- list(
      protein_id = sprintf("%s_P%04d", scenario, idx),
      gene_name = sprintf("GENE_%s_%04d", scenario, idx),
      scenario = scenario,
      chr = chr_num,
      gene_start = gene_start,
      gene_end = gene_end,
      rf_se = rf_se,
      pqtl_se = pqtl_se,
      outcome_se = outcome_se,
      true_beta1 = beta1,
      true_beta2 = beta2,
      true_beta3 = beta3,
      direction_flipped = direction_flipped,
      set_counts = set_counts,
      snps = snps
    )
  }
  proteins
}

row_value <- function(row, name, default) {
  if (name %in% names(row)) row[[name]][[1]] else default
}

expected_instrument_counts <- function(rf_rows, pqtl_rows, info_rows, truth_rows, global_cfg = NULL) {
  global_cfg <- global_cfg %||% list()
  rf_threshold <- as.numeric(global_cfg[["rf_p_threshold"]] %||% 5e-6)
  cis_threshold <- as.numeric(global_cfg[["cis_p_threshold"]] %||% 5e-6)
  cis_window_bp <- as.integer(global_cfg[["cis_window_bp"]] %||% 1000000L)
  rf_inst <- unique(rf_rows$SNP[as.numeric(rf_rows$P) < rf_threshold])

  for (i in seq_len(nrow(truth_rows))) {
    pid <- truth_rows$protein_id[[i]]
    info <- info_rows[info_rows$PROTEIN == pid, , drop = FALSE]
    if (nrow(info) != 1L) {
      next
    }
    pq <- pqtl_rows[pqtl_rows$PROTEIN == pid, , drop = FALSE]
    pq_p <- suppressWarnings(as.numeric(pq$P))
    cis <- pq[
      !is.na(pq_p) &
        pq_p < cis_threshold &
        as.integer(pq$CHR) == as.integer(info$CHR[[1]]) &
        as.integer(pq$BP) >= as.integer(info$START[[1]]) - cis_window_bp &
        as.integer(pq$BP) <= as.integer(info$END[[1]]) + cis_window_bp,
      ,
      drop = FALSE
    ]
    cis_snps <- unique(cis$SNP)
    truth_rows$nA_true[[i]] <- sum(!(rf_inst %in% cis_snps))
    truth_rows$nB_true[[i]] <- sum(!(cis_snps %in% rf_inst))
    truth_rows$nC_true[[i]] <- sum(cis_snps %in% rf_inst)
  }
  truth_rows
}

write_dataset <- function(outdir, benchmark, cell_name, replicate_id, proteins, global_cfg = NULL) {
  rep_dir <- ensure_dir(file.path(outdir, benchmark, cell_name, sprintf("rep_%04d", replicate_id)))
  rf_rows <- list()
  pqtl_rows <- list()
  cancer_rows <- list()
  info_rows <- list()
  truth_rows <- list()
  rf_idx <- pqtl_idx <- cancer_idx <- info_idx <- truth_idx <- 1L

  for (protein in proteins) {
    pid <- protein[["protein_id"]]
    gene <- protein[["gene_name"]]
    chr_num <- protein[["chr"]]
    gene_start <- protein[["gene_start"]]
    gene_end <- protein[["gene_end"]]
    set_counts <- protein[["set_counts"]]
    rf_se <- protein[["rf_se"]]
    pqtl_se <- protein[["pqtl_se"]]
    outcome_se <- protein[["outcome_se"]]

    info_rows[[info_idx]] <- data.frame(PROTEIN = pid, GENE = gene, CHR = chr_num, START = gene_start, END = gene_end, stringsAsFactors = FALSE)
    info_idx <- info_idx + 1L
    truth_rows[[truth_idx]] <- data.frame(
      protein_id = pid,
      gene_name = gene,
      true_scenario = protein[["scenario"]],
      nA_true = as.integer(set_counts[["A"]]),
      nB_true = as.integer(set_counts[["B"]]),
      nC_true = as.integer(set_counts[["C"]]),
      true_beta1 = protein[["true_beta1"]],
      true_beta2 = protein[["true_beta2"]],
      true_beta3 = protein[["true_beta3"]],
      true_mediated_effect = protein[["true_beta1"]] * protein[["true_beta2"]],
      direction_flipped = as.integer(isTRUE(protein[["direction_flipped"]])),
      stringsAsFactors = FALSE
    )
    truth_idx <- truth_idx + 1L

    for (set_name in c("A", "B", "C")) {
      snp_df <- protein[["snps"]][[set_name]]
      if (is.null(snp_df) || !is.data.frame(snp_df) || nrow(snp_df) == 0) {
        next
      }
      for (j in seq_len(nrow(snp_df))) {
        rec <- snp_df[j, , drop = FALSE]
        rec_chr <- as.integer(row_value(rec, "chr", chr_num))
        rec_bp <- as.integer(row_value(rec, "bp", rec$bp))
        rf_rows[[rf_idx]] <- data.frame(
          SNP = rec$rsid, A1 = rec$a1, A2 = rec$a2, FREQ = rec$freq, BETA = rec$rf_beta, SE = rf_se,
          P = rec$rf_p, CHR = rec_chr, BP = rec_bp, stringsAsFactors = FALSE
        )
        rf_idx <- rf_idx + 1L
        pqtl_rows[[pqtl_idx]] <- data.frame(
          PROTEIN = pid, SNP = rec$rsid, A1 = rec$a1, A2 = rec$a2, FREQ = rec$freq, BETA = rec$pqtl_beta, SE = pqtl_se,
          P = rec$pqtl_p, CHR = rec_chr, BP = rec_bp, stringsAsFactors = FALSE
        )
        pqtl_idx <- pqtl_idx + 1L
        cancer_rows[[cancer_idx]] <- data.frame(
          SNP = rec$rsid, A1 = rec$a1, A2 = rec$a2, FREQ = rec$freq, BETA = rec$outcome_beta, SE = outcome_se,
          P = rec$outcome_p, CHR = rec_chr, BP = rec_bp, stringsAsFactors = FALSE
        )
        cancer_idx <- cancer_idx + 1L
      }
    }
  }

  rf_df <- do.call(rbind, rf_rows)
  pqtl_df <- do.call(rbind, pqtl_rows)
  cancer_df <- do.call(rbind, cancer_rows)
  info_df <- do.call(rbind, info_rows)
  truth_df <- do.call(rbind, truth_rows)
  rf_df <- rf_df[!duplicated(rf_df$SNP), , drop = FALSE]
  cancer_df <- cancer_df[!duplicated(cancer_df$SNP), , drop = FALSE]
  truth_df <- expected_instrument_counts(rf_df, pqtl_df, info_df, truth_df, global_cfg = global_cfg)

  rf <- file.path(rep_dir, "rf_sumstat.txt")
  pqtl <- file.path(rep_dir, "pqtl_sumstat.txt")
  cancer <- file.path(rep_dir, "cancer_sumstat.txt")
  protein_info <- file.path(rep_dir, "protein_info.txt")
  truth <- file.path(rep_dir, "truth.tsv")

  write_table(rf, rf_df, c("SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "CHR", "BP"), " ")
  write_table(pqtl, pqtl_df, c("PROTEIN", "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "CHR", "BP"), " ")
  write_table(cancer, cancer_df, c("SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "CHR", "BP"), " ")
  write_table(protein_info, info_df, c("PROTEIN", "GENE", "CHR", "START", "END"), " ")
  write_table(
    truth,
    truth_df,
    c("protein_id", "gene_name", "true_scenario", "nA_true", "nB_true", "nC_true", "true_beta1", "true_beta2", "true_beta3", "true_mediated_effect", "direction_flipped"),
    "\t"
  )
  list(rf = rf, pqtl = pqtl, cancer = cancer, protein_info = protein_info, truth = truth)
}

classification_sequence <- function(proteins_per_scenario, scenarios = SCENARIOS) {
  rep(scenarios, each = proteins_per_scenario)
}

calibration_sequence <- function(proteins_per_replicate, scenario_mix, scenarios = SCENARIOS) {
  probs <- vapply(scenarios, function(s) as.numeric(scenario_mix[[s]]), numeric(1))
  probs <- probs / sum(probs)
  sample(scenarios, size = proteins_per_replicate, replace = TRUE, prob = probs)
}

global_arg <- function(global_cfg, key, flag, transform = identity) {
  value <- global_cfg[[key]]
  if (is.null(value)) {
    return(character())
  }
  c(flag, as.character(transform(value)))
}

config_extra_args <- function(global_cfg) {
  args <- c(
    global_arg(global_cfg, "cis_window_bp", "--cis-window", function(x) as.integer(as.numeric(x) / 1000)),
    global_arg(global_cfg, "eb_tol", "--eb-tol"),
    global_arg(global_cfg, "elbo_tol", "--elbo-tol"),
    global_arg(global_cfg, "prior_p0", "--prior-p0"),
    global_arg(global_cfg, "prior_p1", "--prior-p1"),
    global_arg(global_cfg, "prior_p2", "--prior-p2"),
    global_arg(global_cfg, "prior_p3", "--prior-p3"),
    global_arg(global_cfg, "prior_p4", "--prior-p4"),
    global_arg(global_cfg, "prior_p5", "--prior-p5"),
    global_arg(global_cfg, "eb_prior_strength", "--eb-prior-strength"),
    global_arg(global_cfg, "m1_min_cis_only", "--m1-min-cis-only"),
    global_arg(global_cfg, "m1_min_first_stage_z", "--m1-min-first-stage-z"),
    global_arg(global_cfg, "m1_min_second_stage_z", "--m1-min-second-stage-z"),
    global_arg(global_cfg, "m1_resid_corr_threshold", "--m1-resid-corr-threshold"),
    global_arg(global_cfg, "m1_resid_corr_penalty", "--m1-resid-corr-penalty"),
    global_arg(global_cfg, "sigma2_beta1", "--sigma2-beta1"),
    global_arg(global_cfg, "sigma2_beta2", "--sigma2-beta2"),
    global_arg(global_cfg, "sigma2_beta3", "--sigma2-beta3"),
    global_arg(global_cfg, "direction_mode", "--direction-mode"),
    global_arg(global_cfg, "direction_weight", "--direction-weight"),
    global_arg(global_cfg, "direction_min_prob", "--direction-min-prob")
  )
  if (isTRUE(global_cfg[["fixed_priors"]])) {
    args <- c(args, "--fixed-priors")
  }
  extra <- global_cfg[["binary_options"]]
  if (!is.null(extra)) {
    if (is.list(extra) && is.null(names(extra))) {
      args <- c(args, as.character(unlist(extra, use.names = FALSE)))
    } else {
      for (nm in names(extra)) {
        value <- extra[[nm]]
        flag <- if (startsWith(nm, "--")) nm else paste0("--", gsub("_", "-", nm))
        if (isTRUE(value)) {
          args <- c(args, flag)
        } else if (!identical(value, FALSE) && !is.null(value)) {
          args <- c(args, flag, as.character(value))
        }
      }
    }
  }
  args
}

analysis_config <- function(global_cfg, cell) {
  cfg <- global_cfg
  overrides <- cell[["analysis_overrides"]]
  if (is.null(overrides)) {
    return(cfg)
  }
  for (nm in names(overrides)) {
    cfg[[nm]] <- overrides[[nm]]
  }
  cfg
}

run_bmediator <- function(binary, files, out_prefix, global_cfg) {
  args <- c(
    "--rf-sumstat", files[["rf"]],
    "--pqtl-sumstat", files[["pqtl"]],
    "--cancer-sumstat", files[["cancer"]],
    "--protein-info", files[["protein_info"]],
    "--out", out_prefix,
    "--p-thresh-rf", as.character(global_cfg[["rf_p_threshold"]]),
    "--p-thresh-cis", as.character(global_cfg[["cis_p_threshold"]]),
    "--max-eb-iter", as.character(global_cfg[["max_eb_iter"]]),
    "--max-cavi-iter", as.character(global_cfg[["max_cavi_iter"]]),
    "--threads", as.character(global_cfg[["threads"]]),
    config_extra_args(global_cfg)
  )
  status <- system2(binary, args = args)
  if (!identical(status, 0L)) {
    stop(sprintf("bmediator failed with status %s", status), call. = FALSE)
  }
  list(mediation = sprintf("%s.mediation", out_prefix), hyp = sprintf("%s.hyp", out_prefix))
}

read_truth <- function(path) {
  read_table(path, "\t")
}

read_results <- function(path) {
  read_table(path, "\t")
}

resolve_column <- function(df, candidates) {
  hits <- candidates[candidates %in% names(df)]
  if (length(hits) == 0) {
    stop(sprintf("could not resolve any of columns: %s", paste(candidates, collapse = ", ")), call. = FALSE)
  }
  hits[[1]]
}

attach_truth <- function(results, truth) {
  protein_col <- resolve_column(results, c("Protein", "protein_id"))
  results <- results[!duplicated(results[[protein_col]]), , drop = FALSE]
  match_idx <- match(truth$protein_id, results[[protein_col]])
  extra_cols <- setdiff(names(results), protein_col)
  extra <- results[match_idx, extra_cols, drop = FALSE]
  merged <- cbind(truth, extra, stringsAsFactors = FALSE)
  rownames(merged) <- NULL
  merged
}

compute_benchmark_metrics <- function(rows, alias_map = character()) {
  prob_cols <- c(
    M0 = resolve_column(rows, c("P_M0", "prob_M0")),
    M1 = resolve_column(rows, c("P_M1", "prob_M1")),
    M2 = resolve_column(rows, c("P_M2", "prob_M2")),
    M3 = resolve_column(rows, c("P_M3", "prob_M3")),
    M4 = resolve_column(rows, c("P_M4", "prob_M4")),
    M5 = resolve_column(rows, c("P_M5", "prob_M5"))
  )

  prob_mat <- as.matrix(rows[, prob_cols, drop = FALSE])
  storage.mode(prob_mat) <- "numeric"
  raw_pred <- apply(prob_mat, 1, function(x) {
    if (all(is.na(x))) {
      return(NA_character_)
    }
    names(prob_cols)[which.max(replace(x, is.na(x), -Inf))]
  })
  true_scenario <- vapply(rows$true_scenario, effective_scenario, character(1), alias_map = alias_map)
  pred_scenario <- vapply(raw_pred, effective_scenario, character(1), alias_map = alias_map)

  rows$raw_true_scenario <- rows$true_scenario
  rows$raw_pred_scenario <- raw_pred
  rows$true_scenario <- true_scenario
  rows$pred_scenario <- pred_scenario
  rows$analyzed <- as.integer(!is.na(rows$raw_pred_scenario))
  rows$correct_class <- ifelse(rows$analyzed == 1L, as.integer(rows$true_scenario == rows$pred_scenario), NA_integer_)
  rows$p_m1 <- as.numeric(rows[[prob_cols[["M1"]]]])
  rows$true_state_prob <- vapply(
    seq_len(nrow(rows)),
    function(i) {
      raw_state <- rows$raw_true_scenario[[i]]
      if (raw_state %in% names(prob_cols)) as.numeric(rows[[prob_cols[[raw_state]]]][[i]]) else NA_real_
    },
    numeric(1)
  )
  rows$is_true_m1 <- as.integer(rows$true_scenario == "M1")

  rows <- rows[order(-rows$p_m1, na.last = TRUE), , drop = FALSE]
  rows$rank <- NA_integer_
  rows$estimated_bfdr <- NA_real_
  analyzed_idx <- which(!is.na(rows$p_m1))
  rows$rank[analyzed_idx] <- seq_along(analyzed_idx)
  rows$estimated_bfdr[analyzed_idx] <- cumsum(1.0 - rows$p_m1[analyzed_idx]) / seq_along(analyzed_idx)
  rownames(rows) <- NULL
  rows
}

summarize_fdr_power <- function(rows, thresholds = c(0.01, 0.05, 0.10)) {
  rows <- rows[!is.na(rows$estimated_bfdr), , drop = FALSE]
  total_true_m1 <- sum(rows$is_true_m1 == 1L, na.rm = TRUE)
  out <- vector("list", length(thresholds))
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[[i]]
    selected <- rows[rows$estimated_bfdr <= threshold, , drop = FALSE]
    true_positives <- sum(selected$is_true_m1 == 1L, na.rm = TRUE)
    false_positives <- nrow(selected) - true_positives
    out[[i]] <- data.frame(
      bfdr_threshold = threshold,
      n_selected = nrow(selected),
      true_positives = true_positives,
      false_positives = false_positives,
      empirical_fdr = if (nrow(selected) > 0) false_positives / nrow(selected) else NA_real_,
      power = if (total_true_m1 > 0) true_positives / total_true_m1 else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

summarize_classification <- function(rows, scenarios = NULL) {
  if (is.null(scenarios)) {
    scenarios <- unique(rows$true_scenario)
  }
  by_rows <- list()
  conf_rows <- list()
  bi <- ci <- 1L
  for (true_scenario in scenarios) {
    subset <- rows[rows$true_scenario == true_scenario, , drop = FALSE]
    if (nrow(subset) > 0) {
      analyzed <- subset[subset$analyzed == 1L, , drop = FALSE]
      by_rows[[bi]] <- data.frame(
        true_scenario = true_scenario,
        n = nrow(subset),
        n_analyzed = nrow(analyzed),
        accuracy = mean(analyzed$correct_class, na.rm = TRUE),
        mean_p_m1 = mean(analyzed$p_m1, na.rm = TRUE),
        mean_true_state_prob = mean(analyzed$true_state_prob, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      bi <- bi + 1L
    }
    for (pred in scenarios) {
      conf_rows[[ci]] <- data.frame(
        true_scenario = true_scenario,
        pred_scenario = pred,
        count = sum(rows$true_scenario == true_scenario & rows$pred_scenario == pred, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      ci <- ci + 1L
    }
  }
  list(by_scenario = do.call(rbind, by_rows), confusion = do.call(rbind, conf_rows))
}

summarize_calibration <- function(rows, n_bins = 10L) {
  valid <- !is.na(rows$p_m1) & !is.na(rows$is_true_m1)
  rows <- rows[valid, , drop = FALSE]
  if (nrow(rows) == 0) {
    return(list(calibration = data.frame(), bfdr = data.frame()))
  }

  bins <- seq(0, 1, length.out = n_bins + 1L)
  cal_rows <- list()
  ci <- 1L
  for (idx in seq_len(n_bins)) {
    lo <- bins[idx]
    hi <- bins[idx + 1L]
    if (idx == n_bins) {
      subset <- rows[rows$p_m1 >= lo & rows$p_m1 <= hi, , drop = FALSE]
    } else {
      subset <- rows[rows$p_m1 >= lo & rows$p_m1 < hi, , drop = FALSE]
    }
    if (nrow(subset) == 0) {
      next
    }
    cal_rows[[ci]] <- data.frame(
      calibration_bin = idx - 1L,
      n = nrow(subset),
      mean_pred = mean(subset$p_m1, na.rm = TRUE),
      observed_m1 = mean(subset$is_true_m1, na.rm = TRUE),
      bin_mid = 0.5 * (lo + hi),
      stringsAsFactors = FALSE
    )
    ci <- ci + 1L
  }

  rank_breaks <- unique(c(10L, 25L, 50L, 100L, nrow(rows)))
  start <- 1L
  bfdr_rows <- list()
  bi <- 1L
  for (stop in rank_breaks) {
    if (stop < start) {
      next
    }
    subset <- rows[start:min(stop, nrow(rows)), , drop = FALSE]
    if (nrow(subset) == 0) {
      next
    }
    bfdr_rows[[bi]] <- data.frame(
      rank_bin = sprintf("%d-%d", start, min(stop, nrow(rows))),
      n = nrow(subset),
      mean_estimated_bfdr = mean(subset$estimated_bfdr, na.rm = TRUE),
      empirical_fdr = 1.0 - mean(subset$is_true_m1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    bi <- bi + 1L
    start <- stop + 1L
  }

  list(
    calibration = if (length(cal_rows)) do.call(rbind, cal_rows) else data.frame(),
    bfdr = if (length(bfdr_rows)) do.call(rbind, bfdr_rows) else data.frame()
  )
}
