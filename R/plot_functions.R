# plotting functions

# Plot multivariate VAF distributions with cluster associated to tree nodes.
#
# A multivariate plot is generated from a labelled dataset
#
# @param x Pepi object with labelled spectrum
# @return A multivariate plot
# @examples
# plot_multivariate(spectrum)
# @export

plot_multivariate = function(x){

spectrum = x$VAF

if(!"node" %in% colnames(spectrum)){
  
  stop("no cluster labels") 
  
}
  
max_depth = spectrum %>% pull(node) %>% unique() %>% nchar() %>% max() - 1
cls = get_colors(max_depth)

ggplot(spectrum %>% mutate(vaf_x = Nx/DPx, vaf_y = Ny/DPy)) + geom_point(aes(x = vaf_x, y = vaf_y,color = node)) +
  CNAqc:::my_ggplot_theme() + scale_colour_manual(values = cls) + 
  labs(title = "Multivariate Spectrum", x = "VAF -", y = "VAF +")
  
}

# Plot marginal VAF distributions with cluster associated to tree nodes.
#
# Two marginal histograms are generated from a labelled dataset.
#
# @param spectrum Dataframe with number of variants,depth and node label for any mutation
# @return Two marginal histograms.
# @examples
# plot_marginal(spectrum)
# @export

plot_marginal = function(x){
  
  spectrum = x$VAF
  
  max_depth = spectrum %>% pull(node) %>% unique() %>% nchar() %>% max() - 1
  cls = get_colors(max_depth)
  
  if(!"node" %in% colnames(spectrum)){
    
    stop("no cluster labels") 
    
  }
  
px =   ggplot(spectrum %>% mutate(vaf_x = Nx/DPx) %>% filter(vaf_x > 0)) + 
  geom_histogram(aes(x = vaf_x, fill = node), bins = 50) +
    CNAqc:::my_ggplot_theme() + scale_fill_manual(values = cls) + labs(title = "Marginal -", x = "VAF")

py =   ggplot(spectrum %>% mutate(vaf_y = Ny/DPy) %>% filter(vaf_y > 0)) + 
  geom_histogram(aes(x = vaf_y, fill= node), bins = 50) +
  CNAqc:::my_ggplot_theme() + scale_fill_manual(values = cls) + labs(title = "Marginal +", x = "VAF")

 ggarrange(plotlist = list(px,py),ncol = 2, nrow = 1)
  
 }

# Plot sample tree with branches associated to nodes.
#
# A ggtree plot of the sample tree is generated from the inferred tree.
#
# @param x Pepi object
# @return Plot of the tree.
# @examples
# plot_tree(tree)
# @export

plot_tree = function(x){
  
  if(! "inferred_tree" %in% names(x)){
    
    stop("no inferred tree") 
    
  }
  
  tree = x$inferred_tree
  
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


# to do: plot predicted counts vs data


# Associate colors to nodes.
#
# A list of colors labelled by nodes is generated.
#
# @param max_depth Maximal number of levels of the tree
# @return Named list of colors.
# @examples
# get_colors(max_depth = 2)
# @export

get_colors = function(max_depth){
  
  tree = data.frame(node = "-",level = 0) 
  
  for(l in 1:max_depth){
    
    epsilon = tree %>% filter(level == l-1) 
    node = epsilon %>% pull(node)
    
    new =  lapply(1:length(node), function(i){
      
     new_nodes = tibble(node = c(paste0(node[i],"-"),paste0(node[i],"+")), level = l)
      
       }) %>% bind_rows()
    
    tree = rbind(tree,new)
    
  }
  
   nodes = tree %>% pull(node)
   cls = ggsci::pal_igv()(nodes %>% length())
   names(cls) = nodes

  return(cls)
  
}


