#' plotting functions

#' Plot multivariate VAF distributions with cluster associated to tree nodes.
#'
#' A multivariate plot is generated from a labelled dataset
#'
#' @param spectrum VAF spectrum with cluster labels
#' @return A multivariate plot
#' @examples
#' plot_multivariate(spectrum)
#' @export

plot_multivariate = function(spectrum){


if(!"node" %in% colnames(spectrum)){
  
  stop("no cluster labels") 
  
}
  
max_level = spectrum %>% pull(node) %>% nchar() %>% max() - 1
cls = get_colors(max_level)

ggplot(spectrum %>% mutate(vaf_x = Nx/DPx, vaf_y = Ny/DPy)) + geom_point(aes(x = vaf_x, y = vaf_y,color = node)) +
  get_pepi_theme() + scale_colour_manual(values = cls) + 
  labs(title = "Multivariate Spectrum", x = "VAF -", y = "VAF +")
  
}

#' Plot marginal VAF distributions with cluster associated to tree nodes.
#'
#' Two marginal histograms are generated from a labelled dataset.
#'
#' @param spectrum VAF spectrum with cluster labels
#' @return Two marginal histograms.
#' @examples
#' plot_marginal(spectrum)
#' @export

plot_marginal = function(spectrum){
  
  if(!"node" %in% colnames(spectrum)){
    
    stop("no cluster labels") 
    
  }
  
  max_level = spectrum %>% pull(node) %>% nchar() %>% max() - 1
  cls = get_colors(max_level)
  
px =   ggplot(spectrum %>% mutate(vaf_x = Nx/DPx) %>% filter(vaf_x > 0)) + 
  geom_histogram(aes(x = vaf_x, fill = node), bins = 50) +
    get_pepi_theme() + scale_fill_manual(values = cls) + labs(title = "Marginal -", x = "VAF")

py =   ggplot(spectrum %>% mutate(vaf_y = Ny/DPy) %>% filter(vaf_y > 0)) + 
  geom_histogram(aes(x = vaf_y, fill= node), bins = 50) +
  get_pepi_theme() + scale_fill_manual(values = cls) + labs(title = "Marginal +", x = "VAF")

 ggarrange(plotlist = list(px,py),ncol = 2, nrow = 1)
  
 }

#' Plot sample tree with branches associated to nodes.
#'
#' A ggtree plot of the sample tree is generated from the inferred tree.
#'
#' @param x Tree object
#' @return Plot of the tree.
#' @examples
#' plot_tree(tree)
#' @export

plot_tree = function(tree){
 
  max_level = tree %>% pull(level) %>% max 
  
  tree = tree %>% mutate(level = ifelse(! paste0(node,"-") %in% node,max_level,level))
  
  tree = tree %>% mutate(text = ifelse(level == max_level, 
                                       paste0(substr(node,start = nchar(node),stop = nchar(node)),
                                                    ":",round(m)),paste0(":",round(m)))) %>% 
    mutate(text = paste0(text,"[&&NHX:S=",node,"]"))
  
  for( lev in 1:max_level){
    
    nodes =  tree %>% filter(level == max_level - lev) %>% pull(node)
    
     for (nod in nodes){
       
        t1 = tree %>% filter(node == paste0(nod,"+")) %>% pull(text)
        t2 = tree %>% filter(node == paste0(nod,"-")) %>% pull(text)
        tree = tree %>% mutate(text = ifelse(node == nod,paste0("(",t1,",",t2,")",text),text))
     }
    
  }
  
  text = paste0("(",tree %>% filter(level == 0) %>% pull(text),");")
  
  phylo <- read.nhx(textConnection(text))
  
  tree_plot = ggtree(phylo) + geom_tiplab() + 
    geom_label(aes(x=branch, label=S, fill = S)) 
  
  cls = get_colors(max_level)
  
  tree_plot = tree_plot +
    scale_fill_manual(values = cls) + theme(legend.position = "none")

 tree_plot + ggtitle(
   label = "Switching Tree",
 )
 
}


#' Plot predicted vs input counts.
#'
#' A plot of predicted and input counts is generated.
#'
#' @param x Pepi object
#' @return Plot of the counts.
#' @examples
#' plot_counts(x)
#' @export


plot_counts = function(x){
  
  if(is.null(x$predicted_counts)){
    
    stop("no predicted counts")
  }
  
  obj = rbind(x$counts %>% mutate(type = "data") %>% 
                dplyr::select(-genotype),
        x$predicted_counts %>% mutate(type = "prediction"))
  
  cls = c("dodgerblue","black")
  names(cls) = c("prediction","data")
    
  plot = ggplot(obj) + geom_point(aes(x = time,y = counts,color = type)) + 
          scale_color_manual(values = cls) + facet_grid(~epistate) + 
          get_pepi_theme()
  
   return(plot)
}

#' Plot posterior and prior distributions.
#'
#' Posterior and prior draws histograms are plotted for any required parameter.
#'
#' @param post Posterior draws
#' @param params A list of parameters 
#' @return A plot with posterior and prior distributions
#' @examples
#' plot_inference(x,params = c("rn","rp"))
#' @export

plot_inference = function(x,params = NULL){
  
  prior = rbind(as_tibble(x$prior$tree) %>% mutate(class = "tree") %>% reshape2::melt(),
                as_tibble(x$prior$fitness) %>% mutate(class = "fitness") %>% reshape2::melt(),
                as_tibble(x$prior$counts) %>% mutate(class = "counts") %>% reshape2::melt()) %>% 
    mutate(type = "prior") 
  
  post = rbind(as_tibble(x$posterior$tree) %>% mutate(class = "tree") %>% reshape2::melt(),
               as_tibble(x$posterior$fitness) %>% mutate(class = "fitness") %>% reshape2::melt(),
               as_tibble(x$posterior$counts) %>% mutate(class = "counts") %>% reshape2::melt()) %>% 
    mutate(type = "post")  
  
  sampling = rbind(prior,post)
  
  if(!is.null(params)){
    
    sampling = sampling %>% filter(variable %in% params) 
  }
  
  if (is.null(sampling)) {
    stop("required parameters are not present")
    
  }
  
  nr = round(length(sampling$variable %>% unique())/4 + 1)
  nc = min(length(sampling$variable %>% unique()) + 1,4) 
  cls = c("indianred","steelblue","darkgreen")
  names(cls) = c("tree","fitness","counts")
  
  ggplot(sampling) + geom_density(aes(x = value, alpha = type,fill = class)) + 
    facet_wrap(~variable,scales = "free",nrow = nr, ncol = nc)  + 
    scale_alpha_manual(values = c("prior" = 0.4, "posterior" = 1)) + 
    scale_fill_manual(values = cls) + 
    get_pepi_theme() + theme(legend.position="none")
  
}



