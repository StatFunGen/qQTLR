#' Quantile TWAS Pipeline
#'
#' Performs quantile TWAS analysis for multiple contexts. This is the quantile-specific
#' counterpart to pecotmr's \code{twas_pipeline}, reusing its harmonization and
#' TWAS z-score computation but handling weight extraction, model selection, and
#' output formatting differently for quantile methods.
#'
#' @param twas_weights_data List of list of twas weights output from generate_twas_db function.
#' @param ld_meta_file_path Path to LD reference: either a PLINK2/PLINK1 prefix, or a tab-delimited
#'   metadata file with columns "#chrom", "start", "end", "path" (auto-detected).
#' @param gwas_meta_file A file path for a dataframe table with column of "study_id", "chrom" (integer), "file_path",
#'   "column_mapping_file".
#' @param region_block A string with LD region information of chromosome number, start and end position
#'   of LD block connected with "_".
#' @param ld_reference_sample_size Sample size of the LD reference panel (integer). Required.
#' @param output_twas_data Logical, whether to output the formatted TWAS data. Default is FALSE.
#' @param event_filters Optional event filters.
#' @param column_file_path Optional path to column mapping file.
#' @param comment_string Comment string for file parsing. Default is "#".
#' @return A list containing twas_result, twas_data, and mr_result.
#' @importFrom stats rnorm
#' @export
quantile_twas_pipeline <- function(twas_weights_data,
                                   ld_meta_file_path,
                                   gwas_meta_file,
                                   region_block,
                                   ld_reference_sample_size,
                                   output_twas_data = FALSE,
                                   event_filters = NULL,
                                   column_file_path = NULL,
                                   comment_string = "#") {
  # Internal function to format TWAS output for quantile mode
  format_quantile_twas_data <- function(post_qc_twas_data, twas_table) {
    weights_list <- do.call(c, lapply(names(post_qc_twas_data), function(molecular_id) {
      contexts <- names(post_qc_twas_data[[molecular_id]][["weights_qced"]])
      chrom <- post_qc_twas_data[[molecular_id]][["chrom"]]
      do.call(c, lapply(contexts, function(context) {
        weight <- list()
        data_type <- post_qc_twas_data[[molecular_id]][["data_type"]][[context]]
        postqc_scaled_weight <- list()
        gwas_studies <- names(post_qc_twas_data[[molecular_id]][["weights_qced"]][[context]])

        # For quantile TWAS: extract all available methods from weight matrix columns
        if (length(gwas_studies) > 0) {
          sample_weight_matrix <- post_qc_twas_data[[molecular_id]][["weights_qced"]][[context]][[gwas_studies[1]]][["scaled_weights"]]
          all_columns <- colnames(sample_weight_matrix)
          potential_methods <- unique(gsub("_weights$", "", all_columns))
          methods_with_suffix <- potential_methods[paste0(potential_methods, "_weights") %in% all_columns]

          if (length(methods_with_suffix) == 0) {
            available_methods <- all_columns
            use_direct_columns <- TRUE
          } else {
            available_methods <- methods_with_suffix
            use_direct_columns <- FALSE
          }

          for (method in available_methods) {
            for (study in gwas_studies) {
              weight_matrix <- post_qc_twas_data[[molecular_id]][["weights_qced"]][[context]][[study]][["scaled_weights"]]

              if (use_direct_columns) {
                selected_col <- method
              } else {
                weight_col_candidates <- c(paste0(method, "_weights"), method, "weight")
                selected_col <- NULL
                for (col_candidate in weight_col_candidates) {
                  if (!is.null(col_candidate) && col_candidate %in% colnames(weight_matrix)) {
                    selected_col <- col_candidate
                    break
                  }
                }
              }

              if (!is.null(selected_col) && selected_col %in% colnames(weight_matrix)) {
                postqc_scaled_weight[[study]] <- weight_matrix[, selected_col, drop = FALSE]
                colnames(postqc_scaled_weight[[study]]) <- "weight"
                context_variants <- rownames(weight_matrix)
                context_range <- as.integer(sapply(context_variants, function(variant) {
                  parts <- strsplit(variant, ":")[[1]]
                  if (length(parts) >= 2) as.integer(parts[2]) else NA
                }))
                context_range <- context_range[!is.na(context_range)]

                if (length(context_range) > 0) {
                  quantile_info <- twas_table[twas_table$molecular_id == molecular_id &
                                                twas_table$context == context &
                                                twas_table$method == method, ]

                  if (nrow(quantile_info) > 0) {
                    quantile_start <- if ("quantile_start" %in% colnames(quantile_info)) quantile_info$quantile_start[1] else NA
                    quantile_end <- if ("quantile_end" %in% colnames(quantile_info)) quantile_info$quantile_end[1] else NA
                  } else {
                    quantile_start <- NA
                    quantile_end <- NA
                  }

                  quantile_range <- if (!is.na(quantile_start) && !is.na(quantile_end)) {
                    paste0(quantile_start, "-", quantile_end)
                  } else {
                    NA
                  }

                  safe_data_type <- if (is.null(data_type) || any(is.na(data_type)) || length(data_type) == 0) {
                    "unknown"
                  } else {
                    data_type[1]
                  }

                  weight_id <- paste0(molecular_id, "|", safe_data_type, "_", context, "_", method)
                  if (is.null(weight[[weight_id]])) {
                    weight[[weight_id]] <- list()
                  }

                  weight[[weight_id]][[study]] <- list(
                    chrom = chrom,
                    p0 = min(context_range),
                    p1 = max(context_range),
                    wgt = postqc_scaled_weight[[study]],
                    molecular_id = molecular_id,
                    weight_name = paste0(safe_data_type, "_", context, "_", method),
                    type = safe_data_type,
                    context = context,
                    method = method,
                    quantile_start = quantile_start,
                    quantile_end = quantile_end,
                    quantile_range = quantile_range,
                    n_wgt = length(context_variants)
                  )
                }
              }
            }
          }
        }
        return(weight)
      }))
    }))
    weights <- weights_list[!sapply(weights_list, is.null)]

    # gene_z table — quantile TWAS includes all methods
    if (nrow(twas_table) > 0) {
      twas_table$id <- paste0(twas_table$molecular_id, "|",
                              ifelse(is.null(twas_table$type) | any(is.na(twas_table$type)), "unknown", twas_table$type),
                              "_", twas_table$context, "_", twas_table$method)
      twas_table$group <- paste0(twas_table$context, "|",
                                 ifelse(is.null(twas_table$type) | any(is.na(twas_table$type)), "unknown", twas_table$type),
                                 "|", twas_table$method)

      twas_table$z <- twas_table$twas_z

      output_columns <- c("id", "z", "type", "context", "group", "gwas_study", "method")
      if ("quantile_start" %in% colnames(twas_table)) output_columns <- c(output_columns, "quantile_start", "quantile_end")
      if ("pseudo_R2_avg" %in% colnames(twas_table)) output_columns <- c(output_columns, "pseudo_R2_avg")
      twas_table <- twas_table[, intersect(output_columns, colnames(twas_table)), drop = FALSE]
      studies <- unique(twas_table$gwas_study)
      z_gene_list <- list()
      for (study in studies) {
        z_gene_list[[study]] <- twas_table[twas_table$gwas_study == study, , drop = FALSE]
      }
      result <- list(weights = weights, z_gene = z_gene_list)
    } else {
      result <- list(weights = weights, z_gene = list())
    }
    return(result)
  }

  # Step 1: Harmonize TWAS data
  twas_data_qced_result <- pecotmr::harmonize_twas(twas_weights_data, ld_meta_file_path, gwas_meta_file,
                                                    ld_reference_sample_size = ld_reference_sample_size,
                                                    column_file_path = column_file_path, comment_string = comment_string)
  twas_results_db <- lapply(names(twas_weights_data), function(weight_db) {
    twas_weights_data[[weight_db]][["molecular_id"]] <- weight_db
    twas_data_qced <- twas_data_qced_result$twas_data_qced
    if (length(twas_data_qced[[weight_db]]) == 0 | is.null(twas_data_qced[[weight_db]])) {
      warning(paste0("No data harmonized for ", weight_db, ". Returning NULL for TWAS result for this region."))
      return(NULL)
    }
    # Quantile TWAS: skip model selection
    message("Quantile TWAS: skipping best model selection.")
    twas_data_qced[[weight_db]][["model_selection"]] <- setNames(
      rep(NA, length(names(twas_weights_data[[weight_db]]$weights))),
      names(twas_weights_data[[weight_db]]$weights)
    )
    if (!"data_type" %in% names(twas_weights_data[[weight_db]])) {
      twas_data_qced[[weight_db]][["data_type"]] <- setNames(rep(
        list(NA),
        length(names(twas_weights_data[[weight_db]]$weights))
      ), names(twas_weights_data[[weight_db]]$weights))
    }
    if (length(weight_db) < 1) stop(paste0("No data harmonized for ", weight_db, ". "))
    contexts <- names(twas_data_qced[[weight_db]][["weights_qced"]])
    gwas_studies <- names(twas_data_qced[[weight_db]][["gwas_qced"]])

    # TWAS analysis (same as standard TWAS)
    twas_gene_results <- lapply(contexts, function(context) {
      study_results <- lapply(gwas_studies, function(study) {
        twas_variants <- Reduce(intersect, list(
          rownames(twas_data_qced[[weight_db]][["weights_qced"]][[context]][[study]][["weights"]]),
          twas_data_qced[[weight_db]][["variant_names"]][[context]][[study]],
          twas_data_qced[[weight_db]][["gwas_qced"]][[study]]$variant_id
        ))
        if (length(twas_variants) == 0) {
          return(list(twas_rs_df = data.frame()))
        }
        twas_rs <- pecotmr::twas_analysis(
          twas_data_qced[[weight_db]][["weights_qced"]][[context]][[study]][["weights"]],
          twas_data_qced[[weight_db]][["gwas_qced"]][[study]],
          twas_data_qced[[weight_db]][["LD"]], twas_variants
        )
        if (is.null(twas_rs)) {
          return(list(twas_rs_df = data.frame()))
        }
        twas_rs_df <- data.frame(
          gwas_study = study, method = sub("_[^_]+$", "", names(twas_rs)),
          twas_z = pecotmr::find_data(twas_rs, c(2, "z")),
          twas_pval = pecotmr::find_data(twas_rs, c(2, "pval")),
          context = context, molecular_id = weight_db
        )
        return(list(twas_rs_df = twas_rs_df))
      })
      twas_context_table <- do.call(rbind, lapply(study_results, function(x) x$twas_rs_df))
      return(list(twas_context_table = twas_context_table))
    })
    twas_data_qced[[weight_db]][["LD"]] <- NULL
    twas_weights_data[[weight_db]] <- NULL
    twas_gene_table <- do.call(rbind, lapply(twas_gene_results, function(x) x$twas_context_table))
    return(list(twas_table = twas_gene_table, twas_data_qced = twas_data_qced[weight_db]))
  })
  rm(twas_data_qced_result)
  gc()
  twas_results_db <- twas_results_db[!sapply(twas_results_db, function(x) is.null(x) || (is.list(x) && all(sapply(x, is.null))))]
  if (length(twas_results_db) == 0) {
    return(list(NULL))
  }
  twas_results_table <- do.call(rbind, lapply(twas_results_db, function(x) x$twas_table))
  twas_data <- do.call(c, lapply(twas_results_db, function(x) x$twas_data_qced))
  rm(twas_results_db)
  gc()

  # Step 2: Summarize and merge quantile CV results
  twas_table <- do.call(rbind, lapply(names(twas_data), function(molecular_id) {
    contexts <- names(twas_weights_data[[molecular_id]]$weights)
    gene_table <- do.call(rbind, lapply(contexts, function(context) {
      methods <- sub("_[^_]+$", "", names(twas_weights_data[[molecular_id]]$twas_cv_performance[[context]]))
      cv_performance <- twas_weights_data[[molecular_id]]$twas_cv_performance[[context]]
      if (length(methods) == 0) {
        return(data.frame())
      }
      method_results <- list()
      for (method in methods) {
        if (!is.null(cv_performance[[paste0(method, "_performance")]])) {
          performance_data <- cv_performance[[paste0(method, "_performance")]]
          method_results[[method]] <- data.frame(
            context = context,
            method = method,
            quantile_start = performance_data[, "quantile_start"],
            quantile_end = performance_data[, "quantile_end"],
            pseudo_R2_avg = performance_data[, "pseudo_R2_avg"],
            type = twas_weights_data[[molecular_id]][["data_type"]][[context]]
          )
        }
      }
      if (length(method_results) > 0) {
        do.call(rbind, method_results)
      } else {
        data.frame()
      }
    }))
    gene_table$molecular_id <- molecular_id
    return(gene_table)
  }))
  twas_table$chr <- as.integer(sub("^chr", "", gsub("\\_.*", "", region_block)))
  twas_table$block <- region_block

  # Step 3: Merge and output
  colname_ordered <- c("chr", "molecular_id", "context", "gwas_study", "method",
                       "quantile_start", "quantile_end", "pseudo_R2_avg",
                       "twas_z", "twas_pval", "type", "block")
  if (nrow(twas_results_table) == 0) {
    return(list(twas_result = NULL, twas_data = NULL, mr_result = NULL))
  }
  twas_table <- merge(twas_table, twas_results_table, by = c("molecular_id", "context", "method"))
  if (output_twas_data & nrow(twas_table) > 0) {
    twas_data_subset <- format_quantile_twas_data(twas_data, twas_table)
  } else {
    twas_data_subset <- NULL
  }
  return(list(twas_result = twas_table[, colname_ordered], twas_data = twas_data_subset, mr_result = NULL))
}
