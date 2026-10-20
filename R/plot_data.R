#' Get the branch of a cell in a phylogenetic forest.
#'
#' Walks the forest from `cell_id` back through its ancestors, stopping
#' before `stop_id`.
#'
#' @param forest Table of forest nodes with columns `cell_id` and `ancestor`
#'   (as returned by ProCESS' \code{forest$get_nodes()}).
#' @param cell_id Id of the cell the branch starts from.
#' @param stop_id Id of the ancestor where the walk stops (excluded).
#' @return Vector of cell ids, from `cell_id` up to the child of `stop_id`.
#' @examples
#' data("forest")
#' driver_id = min(forest$cell_id[forest$mutant == "B"])
#' PEPI:::get_branch(forest, cell_id = driver_id, stop_id = 0)
#' @keywords internal

get_branch = function(forest, cell_id, stop_id){

  s = cell_id
  branch = s

  while(s > stop_id){

    s = forest %>% dplyr::filter(.data$cell_id == s) %>% dplyr::pull(ancestor)
    if(s > stop_id){branch = c(branch, s)}
  }

  return(branch)

}


#' Plot a phylogenetic forest coloured by species.
#'
#' Cells are placed by birth time, coloured by species (mutant + epistate)
#' and shaped by sample.
#'
#' @param forest Table of forest nodes with columns `cell_id`, `ancestor`,
#'   `mutant`, `epistate`, `sample` and `birth_time`.
#' @param ids Optional vector of cell ids to plot (e.g. the branches of a
#'   subset of sampled cells). If NULL, the whole forest is plotted.
#' @param colors Named vector of species colors.
#' @return A ggplot object.
#' @examples
#' data("forest")
#' set.seed(31)
#' cells = sample(forest$cell_id[!is.na(forest$sample)], 20)
#' ids = unique(c(0, unlist(lapply(cells, function(id) PEPI:::get_branch(forest, id, 0)))))
#' PEPI:::plot_forest(forest, ids)
#' @keywords internal

plot_forest = function(forest, ids = NULL,
                       colors = c("wt-" = "#82CCE0", "wt+" = "#CE6D6E",
                                  "B-" = "steelblue", "B+" = "darkorange3")){

  if(!is.null(ids)){
    nodes = forest %>% dplyr::filter(cell_id %in% ids)
  }else{
    nodes = forest
  }

  if(nrow(nodes) == 0){
    warning("The forest does not contain any node")
    return(ggplot2::ggplot())
  }

  forest_data = nodes %>%
    dplyr::as_tibble() %>%
    dplyr::select(cell_id, ancestor, mutant, epistate, sample, birth_time) %>%
    dplyr::bind_rows(tibble::tibble(birth_time = 0)) %>%
    dplyr::rename(from = ancestor, to = cell_id) %>%
    dplyr::mutate(
      from = ifelse(is.na(from), "WT", from),
      to = ifelse(is.na(to), "WT", to),
      species = paste0(mutant, epistate),
      sample = ifelse(is.na(sample), "N/A", sample),
      highlight = FALSE
    )

  edges = forest_data %>% dplyr::select(from, to, highlight)

  graph = tidygraph::as_tbl_graph(edges, directed = TRUE)

  graph = graph %>%
    tidygraph::activate("nodes") %>%
    dplyr::left_join(
      forest_data %>%
        dplyr::rename(name = to) %>%
        dplyr::mutate(name = as.character(name)),
      by = "name"
    )

  # tree layout, root on top
  layout = ggraph::create_layout(graph, layout = "tree", root = "WT")

  max_Y = max(layout$birth_time, na.rm = TRUE)
  layout$y = max_Y - layout$birth_time

  samples = forest %>% dplyr::filter(!is.na(sample)) %>% dplyr::pull(sample) %>% unique()
  nsamples = length(samples)

  labels_every = max_Y/10

  point_size = c(.5, rep(1, nsamples))
  names(point_size) = c("N/A", samples)

  ggraph::ggraph(layout, "tree") +
    ggraph::geom_edge_link(edge_width = .1,
                           ggplot2::aes(edge_color = ifelse(highlight, "indianred3", "black"))) +
    ggraph::geom_node_point(ggplot2::aes(color = species,
                                         shape = ifelse(is.na(sample), "N/A", sample),
                                         size = sample)) +
    ggplot2::scale_shape_manual(values = c(0:nsamples + 1)) +
    ggplot2::scale_color_manual(values = colors) +
    ggraph::scale_edge_color_identity() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(shape = "Sample", color = "Species", x = NULL, y = "Cell division") +
    ggplot2::guides(size = "none",
                    shape = ggplot2::guide_legend("Sample"),
                    fill = ggplot2::guide_legend("Species")) +
    ggplot2::scale_size_manual(values = point_size) +
    ggplot2::scale_y_continuous(labels = seq(0, max_Y, labels_every) %>% round %>% rev,
                                breaks = seq(0, max_Y, labels_every) %>% round) +
    ggplot2::theme(
      axis.line.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

}


#' Plot a phylogenetic forest with nodes labelled by mutation class.
#'
#' Nodes carrying a label (e.g. "driver", "truncal") are highlighted;
#' all other nodes are labelled "subclonal".
#'
#' @param forest Table of forest nodes with columns `cell_id`, `ancestor`,
#'   `mutant`, `epistate`, `sample` and `birth_time`.
#' @param labels Table with columns `cell_id` and `label`.
#' @param cls Named vector of label colors. If NULL a palette is generated.
#' @return A ggplot object.
#' @examples
#' data("forest")
#' set.seed(31)
#' cells = sample(forest$cell_id[!is.na(forest$sample)], 20)
#' ids = unique(c(0, unlist(lapply(cells, function(id) PEPI:::get_branch(forest, id, 0)))))
#' driver_id = min(forest$cell_id[forest$mutant == "B"])
#' labels = data.frame(cell_id = PEPI:::get_branch(forest, driver_id, 0), label = "driver")
#' PEPI:::plot_sticks(forest[forest$cell_id %in% ids, ], labels,
#'   cls = c(subclonal = "gainsboro", driver = "forestgreen"))
#' @keywords internal

plot_sticks = function(forest, labels, cls = NULL){

  forest_data = forest %>% dplyr::as_tibble() %>%
    dplyr::rename(from = ancestor, to = cell_id) %>%
    dplyr::select(from, to, mutant, epistate, sample, birth_time) %>%
    dplyr::mutate(
      from = ifelse(is.na(from), "WT", from),
      species = paste0(mutant, epistate),
      sample = ifelse(is.na(sample), "N/A", sample),
      highlight = FALSE
    ) %>%
    dplyr::full_join(labels %>%
                       dplyr::rename(to = cell_id) %>%
                       dplyr::select(dplyr::all_of(c("to", "label"))),
                     by = "to") %>%
    dplyr::mutate(label = ifelse(is.na(label), "subclonal", label))

  edges = forest_data %>% dplyr::select(from, to, highlight)

  graph = tidygraph::as_tbl_graph(edges, directed = TRUE)
  graph = graph %>% tidygraph::activate("nodes") %>%
    dplyr::left_join(forest_data %>% dplyr::rename(name = to) %>%
                       dplyr::mutate(name = as.character(name)),
                     by = "name")

  layout = ggraph::create_layout(graph, layout = "tree", root = "WT")
  max_Y = max(layout$birth_time, na.rm = TRUE)
  layout$y = max_Y - layout$birth_time

  if(is.null(cls)){
    cls = ggsci::pal_simpsons()(forest_data %>% dplyr::pull(label) %>% unique() %>% length())
    names(cls) = forest_data %>% dplyr::pull(label) %>% unique()
    cls["subclonal"] = "gainsboro"
  }

  ncells = graph %>% tidygraph::activate("nodes") %>% dplyr::pull(name) %>% length()
  ncells_sampled = graph %>% tidygraph::activate("nodes") %>%
    dplyr::filter(sample != "N/A") %>% dplyr::pull(name) %>% length()
  nsamples = forest %>% dplyr::filter(!is.na(sample)) %>% dplyr::pull(sample) %>% unique() %>% length()

  labels_every = max_Y/10

  highlighted = unique(forest_data$label[forest_data$label != "subclonal"])
  point_size = c(0, rep(1, length(highlighted)))
  names(point_size) = c("subclonal", highlighted)

  ggraph::ggraph(layout, "tree") +
    ggraph::geom_edge_link(edge_width = 0.5,
                           ggplot2::aes(edge_color = ifelse(highlight, "indianred3", "gainsboro"))) +
    ggraph::geom_node_point(ggplot2::aes(color = label, size = label)) +
    ggplot2::scale_color_manual(values = cls, na.translate = FALSE) +
    ggraph::scale_edge_color_identity() +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      title = "Phylogenetic tree",
      subtitle = paste0(ncells, " cells (", ncells_sampled, " sampled from ", nsamples, " samples)"),
      color = NULL,
      x = NULL,
      y = "Birth time"
    ) +
    ggplot2::guides(size = "none") +
    ggplot2::scale_size_manual(values = point_size) +
    ggplot2::scale_y_continuous(labels = seq(0, max_Y, labels_every) %>% round %>% rev,
                                breaks = seq(0, max_Y, labels_every) %>% round) +
    ggplot2::theme(
      axis.line.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

}


#' Plot marginal and pairwise VAF distributions of the sorted samples.
#'
#' Every column whose name contains "VAF" is treated as one sample
#' (columns starting with "normal" are dropped). Columns named as returned
#' by ProCESS' \code{simulate_seq()} with an epistate labelling
#' (e.g. `Sampling_1_N.VAF`, `Sampling_1_P.VAF`) are renamed to
#' `VAF 1-`, `VAF 1+`, ...
#'
#' @param muts Mutation table with VAF columns and, optionally, a `label` column.
#' @param colors Named vector of label colors.
#' @param vaf_cut Minimum VAF shown in the marginal histograms.
#' @param size Point size of the pairwise scatter plots.
#' @param assemble If TRUE, return one figure with the pairwise scatter plots
#'   on top and the marginal histograms below; otherwise a list with elements
#'   `marginals` and `pairs`.
#' @return A ggplot object, or a list of ggplot objects.
#' @examples
#' data("muts")
#' PEPI:::plot_vaf(muts, vaf_cut = 0.03, size = 0.1)
#' @keywords internal

plot_vaf = function(muts, colors = NULL, vaf_cut = 0.02, size = 0.01, assemble = TRUE){

  muts = muts %>% dplyr::select(-dplyr::starts_with("normal"))

  colnames(muts) = sub("^Sampling_(.*)_N\\.VAF$", "VAF \\1-", colnames(muts))
  colnames(muts) = sub("^Sampling_(.*)_P\\.VAF$", "VAF \\1+", colnames(muts))

  has_label = "label" %in% colnames(muts)

  if(has_label){
    if(!is.null(colors)){
      cls = colors
    }else{
      cls = ggsci::pal_simpsons()(muts %>% dplyr::pull(label) %>% unique() %>% length())
      names(cls) = muts %>% dplyr::pull(label) %>% unique()
      cls["subclonal"] = "gainsboro"
    }
  }

  lab = colnames(muts)[grepl(x = colnames(muts), pattern = "VAF")]

  marginals = lapply(lab, function(l){

    y = muts
    colnames(y)[colnames(y) == l] = "VAF"

    p = ggplot2::ggplot(y %>% dplyr::filter(VAF > vaf_cut)) +
      ggplot2::geom_histogram(ggplot2::aes(x = VAF), binwidth = 0.01) +
      get_pepi_theme() +
      ggplot2::labs(x = l)

    if(has_label){
      p = p + ggplot2::geom_histogram(ggplot2::aes(x = VAF, fill = label), binwidth = 0.01) +
        ggplot2::scale_fill_manual(values = cls)
    }
    p

  })

  if(length(lab) < 2){
    return(marginals)
  }

  g = expand.grid(l1 = seq_along(lab), l2 = seq_along(lab)) %>% dplyr::filter(l2 > l1)

  pairs = lapply(seq_len(nrow(g)), function(i){

    y = muts
    y$VAF_1 = y[[lab[g$l1[i]]]]
    y$VAF_2 = y[[lab[g$l2[i]]]]

    p = ggplot2::ggplot(y) +
      ggplot2::geom_point(ggplot2::aes(x = VAF_1, y = VAF_2), size = size) +
      get_pepi_theme() +
      ggplot2::labs(x = lab[g$l1[i]], y = lab[g$l2[i]])

    if(has_label){
      p = p + ggplot2::geom_point(ggplot2::aes(x = VAF_1, y = VAF_2, color = label), size = size) +
        ggplot2::scale_color_manual(values = cls)
    }
    p

  })

  if(!assemble){
    return(list(marginals = marginals, pairs = pairs))
  }

  nm = ceiling(sqrt(length(marginals)))

  ggpubr::ggarrange(
    ggpubr::ggarrange(plotlist = pairs, ncol = 3, nrow = ceiling(length(pairs)/3),
                      common.legend = TRUE, legend = "none"),
    ggpubr::ggarrange(plotlist = marginals, ncol = nm, nrow = ceiling(length(marginals)/nm),
                      common.legend = TRUE, legend = "none"),
    ncol = 2, nrow = 1, widths = c(3, 2)
  )

}


#' Plot the simulated data: phylogenetic forest and VAF spectra.
#'
#' Samples `n_cells` cells from each sample and restricts the forest to their
#' branches. Top: the forest coloured by species (\code{plot_forest()}) and
#' with the driver branch and the truncal node highlighted
#' (\code{plot_sticks()}). Bottom: pairwise and marginal VAF spectra of the
#' sorted samples (\code{plot_vaf()}), with mutations labelled as truncal,
#' driver (acquired on the branch from the root to the first cell of
#' `driver_mutant`) or subclonal.
#'
#' @param forest Table of forest nodes (see \code{plot_forest()}).
#' @param muts Mutation table with VAF columns and `first_cell_id`.
#' @param n_cells Number of cells drawn from each sample.
#' @param seed Optional seed for the cell sampling.
#' @param driver_mutant Name of the driver mutant.
#' @param vaf_cut Minimum VAF shown in the marginal histograms.
#' @param cls Named vector of colors for the labels `subclonal`, `driver`, `truncal`.
#' @return A ggplot object.
#' @examples
#' data("forest")
#' data("muts")
#' PEPI:::plot_data(forest, muts, seed = 31)
#' @keywords internal

plot_data = function(forest, muts, n_cells = 50, seed = NULL, driver_mutant = "B", vaf_cut = 0.03,
                     cls = c("subclonal" = "gainsboro", "driver" = "forestgreen", "truncal" = "darkgoldenrod3")){

  if(!is.null(seed)) set.seed(seed)

  samples = forest %>% dplyr::filter(!is.na(sample)) %>% dplyr::pull(sample) %>% unique() %>% sort()

  ids = unlist(lapply(samples, function(s){
    cells = forest %>% dplyr::filter(!is.na(sample), sample == s) %>% dplyr::pull(cell_id)
    cells[sample(x = seq_along(cells), size = min(n_cells, length(cells)), replace = FALSE)]
  }))

  all_ids = c(0)
  for(id in ids){
    all_ids = unique(c(all_ids, get_branch(forest = forest, cell_id = id, stop_id = 0)))
  }

  forest = forest %>% dplyr::mutate(sample = gsub(x = sample, pattern = "Sampling_", replacement = "T"))

  driver_id = forest %>% dplyr::filter(mutant == driver_mutant) %>% dplyr::pull(cell_id) %>% min()
  branch_driver = get_branch(forest = forest, cell_id = driver_id, stop_id = 0)

  muts = muts %>%
    dplyr::mutate(label = ifelse(first_cell_id %in% branch_driver, "driver",
                                 ifelse(first_cell_id == 0, "truncal", "subclonal")))

  p_forest = plot_forest(forest, all_ids)

  p_sticks = plot_sticks(
    forest = forest %>% dplyr::filter(cell_id %in% all_ids),
    labels = muts %>% dplyr::filter(first_cell_id %in% all_ids) %>%
      dplyr::select(cell_id = first_cell_id, label) %>% unique(),
    cls = cls
  )

  p_vaf = plot_vaf(muts, vaf_cut = vaf_cut, colors = cls, size = 0.1, assemble = FALSE)

  nm = ceiling(sqrt(length(p_vaf$marginals)))

  ggpubr::ggarrange(
    p_forest,
    p_sticks,
    ggpubr::ggarrange(plotlist = p_vaf$pairs, ncol = 3, nrow = ceiling(length(p_vaf$pairs)/3), legend = "none"),
    ggpubr::ggarrange(plotlist = p_vaf$marginals, ncol = nm, nrow = ceiling(length(p_vaf$marginals)/nm), legend = "none"),
    nrow = 2, ncol = 2
  )

}
