#' Plot prior vs posterior for a parameter in a PEPI object
#'
#' @param pepi A PEPI object with `posterior` and `prior` extracted.
#' @param parameter Name of the parameter to plot.
#' @param log_scale Logical, whether to plot x-axis on log scale.
#' @return A ggplot object.
#' @export
plot_inference <- function(pepi, parameter, log_scale = FALSE) {
  if (is.null(pepi$posterior)) stop("Posterior not found. Run get_posterior() first.")
  if (is.null(pepi$prior)) stop("Prior not found. Run get_prior() first.")
  
  draws <- bind_rows(
    pepi$posterior %>% filter(variable == parameter),
    pepi$prior %>% filter(variable == parameter)
  )
  
  p <- ggplot(draws, aes(x = value, fill = type, alpha = type)) +
    geom_histogram(position = "identity", bins = 100) +
    theme_minimal() +
    scale_alpha_manual(values = c(posterior = 1, prior = 0.5)) +
    labs(x = parameter, y = "Count") +
    theme(legend.position = "top")
  
  if (log_scale) p <- p + scale_x_log10()
  
  return(p)
}


#' Plot posterior predictive for a parameter
#'
#' @param pepi A PEPI object with `posterior` extracted.
#' @param parameter Name of the posterior predictive variable (starts with "pred").
#' @param obs Optional data frame with columns: variable and exp (observed value).
#' @param ncol Number of columns in facet_wrap if assembling multiple variables.
#' @param log_scale Logical, whether to plot x-axis on log scale.
#' @return A ggplot object.
#' @export
plot_posterior_predictive <- function(pepi, parameter, obs = NULL, ncol = 4, log_scale = TRUE) {
  if (is.null(pepi$posterior)) stop("Posterior not found. Run get_posterior() first.")
  
  draws <- pepi$posterior %>% filter(variable == parameter)
  if (nrow(draws) == 0) stop("Parameter not found in posterior draws.")
  
  p <- ggplot(draws, aes(x = value)) +
    geom_histogram(fill = "steelblue", bins = 100, position = "stack") +
    theme_minimal() +
    labs(x = parameter, y = "Count")
  
  if (!is.null(obs)) {
    p <- p + geom_vline(data = obs %>% filter(variable == parameter), 
                        aes(xintercept = exp), linetype = "dashed", color = "red")
  }
  
  if (log_scale) p <- p + scale_x_log10()
  
  return(p)
}


#' Plot a forest of cells
#'
#' @param forest A data.frame representing the forest. Must contain columns:
#'   \code{cell_id}, \code{ancestor}, \code{mutant}, \code{epistate}, \code{sample}, \code{birth_time}.
#'
#' @return A ggplot2 object representing the forest as a tree.
#' @export
plot_forest <- function(forest) {
  
  required_cols <- c("cell_id", "ancestor", "mutant", "epistate", "sample", "birth_time")
  missing_cols <- setdiff(required_cols, colnames(forest))
  if (length(missing_cols) > 0) {
    stop("The forest object is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Create species vector
  species <- forest %>%
    dplyr::select(mutant, epistate) %>%
    dplyr::distinct() %>%
    dplyr::mutate(species = paste0(mutant, epistate)) %>%
    dplyr::pull(species)
  
  # Assign a ggsci color palette
  species_colors <- ggsci::pal_igv()(length(species))
  names(species_colors) <- species
  
  nodes <- forest
  
  if (nrow(nodes) == 0) {
    warning("The forest does not contain any node")
    return(ggplot2::ggplot())
  }
  
  forest_data <- nodes
  forest_data[nrow(forest_data) + 1, ] <- c(NA, NA, NA, NA, NA, 0)
  
  forest_data <- forest_data %>%
    dplyr::as_tibble() %>%
    dplyr::rename(from = .data$ancestor, to = .data$cell_id) %>%
    dplyr::mutate(
      from = ifelse(is.na(.data$from), "WT", .data$from),
      to = ifelse(is.na(.data$to), "WT", .data$to),
      species = paste0(.data$mutant, .data$epistate),
      sample = ifelse(is.na(.data$sample), "N/A", .data$sample),
      highlight = FALSE
    )
  
  edges <- forest_data %>% dplyr::select("from", "to", "highlight")
  
  graph <- tidygraph::as_tbl_graph(edges, directed = TRUE) %>%
    tidygraph::activate("nodes") %>%
    dplyr::left_join(
      forest_data %>%
        dplyr::rename(name = .data$to) %>%
        dplyr::mutate(name = as.character(.data$name)),
      by = "name"
    )
  
  layout <- ggraph::create_layout(graph, layout = "tree", root = "WT")
  max_Y <- max(layout$birth_time, na.rm = TRUE)
  layout$reversed_btime <- max_Y - layout$birth_time
  layout$y <- layout$reversed_btime
  
  nsamples <- forest %>% filter(!is.na(sample)) %>% pull(sample) %>% unique() %>% length()
  point_size <- c(.5, rep(1, nsamples))
  names(point_size) <- c("N/A", forest %>% filter(!is.na(sample)) %>% pull(sample) %>% unique())
  labels_every <- max_Y / 10
  
  graph_plot <- ggraph::ggraph(layout, "tree") +
    ggraph::geom_edge_link(edge_width = .1,
                           ggplot2::aes(edge_color = ifelse(highlight, "indianred3", "black"))) +
    ggraph::geom_node_point(ggplot2::aes(color = .data$species,
                                         shape = ifelse(is.na(.data$sample), "N/A", .data$sample),
                                         size = .data$sample)) +
    ggplot2::scale_shape_manual(values = c(0:nsamples + 1)) +
    ggplot2::scale_color_manual(values = species_colors) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(shape = "Sample", x = NULL, y = "Cell division") +
    ggplot2::guides(size = "none",
                    shape = ggplot2::guide_legend("Sample"),
                    fill = ggplot2::guide_legend("Species")) +
    ggplot2::scale_size_manual(values = point_size) +
    ggplot2::scale_y_continuous(labels = seq(0, max_Y, labels_every) %>% round() %>% rev,
                                breaks = seq(0, max_Y, labels_every) %>% round()) +
    ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
  
  return(graph_plot)
}


#' Plot VAF distributions and pairwise scatterplots
#'
#' @param muts A data.frame containing mutation VAFs. Column names must be like vaf_1_n, vaf_1_p, etc.
#' @param variables Character vector of columns to plot. Must exist in muts.
#' @param vaf_cut Minimum VAF to include in plots (default 0.02).
#'
#' @return A list of ggplot objects: marginals and, if more than 1 variable, pairwise scatterplots.
#' @export
plot_vaf <- function(muts, variables = NULL, vaf_cut = 0.02) {
  
  # Check label column
  has_label <- "label" %in% colnames(muts)
  
  # Default: all VAF columns
  vaf_cols <- colnames(muts)[grepl("^vaf_\\d+_[np]$", colnames(muts))]
  
  if (!is.null(variables)) {
    missing_vars <- setdiff(variables, colnames(muts))
    if (length(missing_vars) > 0) {
      stop("These variables are missing in muts: ", paste(missing_vars, collapse = ", "))
    }
    vaf_cols <- variables
  }
  
  # Check all columns are of correct format
  if (!all(grepl("^vaf_\\d+_[np]$", vaf_cols))) {
    stop("All VAF variables must be of the form vaf_X_n or vaf_X_p")
  }
  
  # Color palette for clusters
  if (has_label) {
    clusters <- unique(muts$label)
    cls <- ggsci::pal_simpsons()(length(clusters))
    names(cls) <- clusters
  }
  
  # Marginal histograms
  marginals <- lapply(vaf_cols, function(col) {
    y <- muts
    colnames(y)[colnames(y) == col] <- "VAF"
    
    p <- ggplot(y %>% filter(VAF > vaf_cut & VAF < 0.7)) +
      geom_histogram(aes(x = VAF), binwidth = 0.01, fill = "steelblue") +
      labs(x = col)
    
    if (has_label) {
      p <- ggplot(y %>% filter(VAF > vaf_cut & VAF < 0.7)) +
        geom_histogram(aes(x = VAF, fill = label), binwidth = 0.01) +
        scale_fill_manual(values = cls) +
        labs(x = col)
    }
    
    p
  })
  
  # Pairwise scatterplots if more than 1 variable
  multi_plots <- NULL
  if (length(vaf_cols) > 1) {
    g <- expand.grid(l1 = 1:length(vaf_cols), l2 = 1:length(vaf_cols)) %>%
      dplyr::filter(l2 > l1) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(lab1 = vaf_cols[l1], lab2 = vaf_cols[l2]) %>%
      dplyr::ungroup()
    
    multi_plots <- lapply(1:nrow(g), function(i) {
      y <- muts
      colnames(y)[colnames(y) %in% c(g$lab1[i], g$lab2[i])] <- c("VAF_1", "VAF_2")
      
      p <- ggplot(y %>% filter(VAF_1 < 0.7 & VAF_2 < 0.7)) +
        geom_point(aes(x = VAF_1, y = VAF_2), color = "steelblue") +
        labs(x = g$lab1[i], y = g$lab2[i])
      
      if (has_label) {
        p <- ggplot(y %>% filter(VAF_1 < 0.7 & VAF_2 < 0.7)) +
          geom_point(aes(x = VAF_1, y = VAF_2, color = label)) +
          scale_color_manual(values = cls) +
          labs(x = g$lab1[i], y = g$lab2[i])
      }
      
      p
    })
  }
  
  if (!is.null(multi_plots)) {
    return(list(marginals = marginals, scatterplots = multi_plots))
  } else {
    return(list(marginals = marginals))
  }
}

