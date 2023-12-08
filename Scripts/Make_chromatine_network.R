Do_net_features <- function(chromatine_net, Features, filtered = F){
  fragment_of_interest <- Features$name
  net <- chromatine_net
  igraph <- igraph::graph_from_data_frame(net, directed = F)
  eigen_centrality_result <- igraph::eigen_centrality(igraph, directed = F)$vector %>% as.data.frame()
  page_rank_result <- igraph::page_rank(igraph, directed = F)$vector %>% as.data.frame()
  degree_result <- igraph::degree(igraph) %>% as.data.frame()
  closenessness_result <- igraph::closeness(igraph) %>% as.data.frame()
  betweenness_result <- igraph::betweenness(igraph, directed = F) %>% as.data.frame()
  res <- merge(Features, eigen_centrality_result, by.x = "name", by.y = 0, all.x = T, all.y = T)
  colnames(res)[4] <- "eigenvalue"
  res <- merge(res, page_rank_result, by.x = "name", by.y = 0, all.x = T, all.y = T)
  colnames(res)[5] <- "Page_Rank"
  res <- merge(res, degree_result, by.x = "name", by.y = 0, all.x = T, all.y = T)
  colnames(res)[6] <- "degree"
  res <- merge(res, closenessness_result, by.x = "name", by.y = 0, all.x = T, all.y = T)
  colnames(res)[7] <- "closeness"
  res <- merge(res, betweenness_result, by.x = "name", by.y = 0, all.x = T, all.y = T)
  colnames(res)[8] <- "betweenness"
  min_eigen <- res$eigenvalue %>% .[.!= 0] %>% na.omit(.) %>% min() %>% log()
  res$eigenvalue <- log(res$eigenvalue)
  
  res$eigenvalue <- sapply(res$eigenvalue, function(eigen){
    test <- (is.na(eigen) | eigen == -Inf)
    if (test){
      min_eigen
    }else{
      eigen
    }
  })
  res$Page_Rank <- log(res$Page_Rank)
  deal_with_missing_values(res, "logFC", 0)
  deal_with_missing_values(res, "P.Value", 1)
  list("features" = res, "net" = net)
}

Do_cool_scatterplot <- function(Feature, title, filtered_plot = T){
  min_eigen <- min(Feature$eigenvalue)
  DCpGs <- sapply(1:nrow(Feature), function(n){
    res <- "NoSign"
    if (Feature[n,"P.Value"] < 0.05 & Feature[n,"eigenvalue"] > min_eigen){
      if (Feature[n, "logFC"] > 0){
        res <- "UP"
      }else{
        res <- "DOWN"
      }
    }
    res
  })
  color <-  c("#00FF00", "#888888", "#0000FF")
  if(filtered_plot){
    Feature <- Feature[DCpGs != "NoSign",]
    DCpGs <- DCpGs[DCpGs != "NoSign"]
    color =  c("#00FF00", "#FF0000")
  }
  ggplot(Feature, aes(x = Page_Rank, y = eigenvalue, label = name, colour = DCpGs))+
    geom_text(check_overlap = T, size = 4, nudge_x = 0.0005, hjust = 0, outlier.size = 0)+
    geom_point(size = 0.5)+
    labs(title = paste0("Network-based node prioritization ", title))+
    xlab("Page Rank (log)")+
    ylab("Eigen Centrality (log)")+
    scale_colour_manual(values= color)
}

deal_with_missing_values <- function(df, column, final_value){
  for (i in column){
    data.table::set(df, which(is.na(df[[i]])), i, final_value)
  }
}

From_network_to_TF_activity <- function(chromatine_net, Features, output_folder, filtered = F, title){
  features <- Do_net_features(chromatine_net, Features, filtered)
  params <- paste0("_params", "_filtered"[filtered])
  scat <- Do_cool_scatterplot(features[["features"]], title = title)
  folder_output <- output_folder
  dir.create(folder_output, showWarnings = F)
  
  png(paste0(folder_output, "/Chromatine_Network_based_nodes.png"), width = 1280, height = 720)
  plot(scat)
  dev.off()
  
  min_eigen <- min(na.omit(features[["features"]]$eigenvalue))
  deal_with_missing_values(features[["features"]], c("eigenvalue"), min_eigen)
  deal_with_missing_values(features[["features"]], c("P.Value"), 1)
  res <- list( "features" = features[["features"]], "scat" = scat, "net" = features[["net"]])
}

Make_Cytoscape_network <- function(net, feat, title, collection){
  Genes_n_frag <- feat %>% 
    dplyr::filter(abs(logFC) > 0.3 & P.Value < 0.05) %>% 
    unlist(.$name) %>% 
    unique()
  filtered <- net %>% dplyr::filter(source %in% Genes_n_frag | target %in% Genes_n_frag)
  filtered_feat <- feat %>% 
    dplyr::filter(name %in% filtered$source | name %in% filtered$target)
  colnames(filtered_feat)[1] <- c("id")
  not_featured <- filtered$source
  filtered_feat$Promoter <- ifelse(stringr::str_detect(filtered_feat$id, "[:alpha:]"), "Promoter", "Other_end")
  createNetworkFromDataFrames(nodes = filtered_feat, edges = filtered, title = title, collection = collection)
  Modify_Cytoscape_network(filtered, filtered_feat, title, collection)
  return(list(net = filtered, nodes = filtered_feat))
}

Modify_Cytoscape_network <- function(e, n, title, collection){
  defaults <- list(NODE_SHAPE="diamond",
                   NODE_SIZE=10,
                   EDGE_TRANSPARENCY=200,
                   NODE_LABEL_POSITION="c,c,c,0.00,0.00")
  nodeLabels <- mapVisualProperty(visual.prop = "node label", table.column = 'id', mapping.type = 'p')
  createVisualStyle(title, defaults, list(nodeLabels))
  setVisualStyle(title)
  
  min_logFC <- min(n$logFC)
  max_logFC <- max(n$logFC)
  setNodeShapeMapping(table.column = "Promoter",
                      table.column.values = c("Promoter", "Other_end"),
                      shapes = c('diamond', 'ELLIPSE'),
                      style.name = title)
  
  setNodeColorMapping(table.column = 'logFC',
                      table.column.values = c(min_logFC, 0.0, max_logFC),
                      colors = c('#0000FF', '#FFFFFF', '#FF0000'),
                      style.name = title)
  
  setNodeComboOpacityMapping(table.column = 'Promoter', 
                             table.column.values = c("Other_end", "Promoter"), 
                             opacities = c(20, 255), 
                             style.name = title, 
                             mapping.type = 'd')
  
  setNodeBorderOpacityMapping(table.column = 'P.Value',
                              table.column.values = c(0, 0.1, 1),
                              opacities = c(255, 200, 20),
                              style.name = title)
  
  createColumnFilter(filter.name = "fragment_name", column = "Promoter", criterion = "Other_end", type = "nodes", predicate = "IS")
  setNodeBorderOpacityBypass(getSelectedNodes(), new.values = 1)

  max_degree <- max(n$degree)
  
  setNodeSizeMapping (table.column = 'degree',
                      table.column.values = c(0, max_degree),
                      sizes = c(10, 200),
                      style.name = title)
  
  setNodeFontSizeMapping(table.column = 'degree',
                         table.column.values = c(0, max_degree),
                         sizes = c(15, 75),
                         style.name = title)
  
  layoutNetwork()
  clearSelection()
  exportImage(paste0("../Results/Network_analysis/", collection, "_", title, ".svg"), 'SVG', zoom=200)
}
