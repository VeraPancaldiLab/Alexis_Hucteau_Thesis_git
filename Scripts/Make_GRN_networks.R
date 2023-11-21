Do_net_features <- function(vip, DEG){
  genes_of_interest <- dplyr::filter(DEG, abs(logFC) > 1.5 & P.Value < 0.1) %>%
    .$ID
  TF_of_interest <- dplyr::filter(vip$mrs_table, pval < 0.1) %>% 
    .$TF
  
  element_of_interest <- c(genes_of_interest, TF_of_interest) %>% unique()
  lower_mor <- summary(vip$regulons$mor[vip$regulons$mor < 0])[2]
  upper_mor <- summary(vip$regulons$mor[vip$regulons$mor > 0])[5]
  reg <- vip$regulons %>% 
    dplyr::filter((tf %in% element_of_interest | target %in% element_of_interest) & (mor < lower_mor | mor > upper_mor))
  igraph <- igraph::graph_from_data_frame(reg, directed = T)
  igraph_4_eigenvalue <- igraph::reverse_edges(igraph)
  eigen_centrality_result <- igraph::eigen_centrality(igraph, directed = F, scale = T)$vector %>% as.data.frame()
  page_rank_result <- igraph::page_rank(igraph_4_eigenvalue, directed = T)$vector %>% as.data.frame()
  degree_result <- igraph::degree(igraph) %>% as.data.frame()
  closenessness_result <- igraph::closeness(igraph) %>% as.data.frame()
  betweenness_result <- igraph::betweenness(igraph, directed = F) %>% as.data.frame()
  res <- merge(DEG, eigen_centrality_result, by.x = "ID", by.y = 0, all.x = T, all.y = T)
  colnames(res)[8] <- "eigenvalue"
  res <- merge(res, page_rank_result, by.x = "ID", by.y = 0, all.x = T, all.y = T)
  colnames(res)[9] <- "Page_Rank"
  res <- merge(res, degree_result, by.x = "ID", by.y = 0, all.x = T, all.y = T)
  colnames(res)[10] <- "degree"
  res <- merge(res, closenessness_result, by.x = "ID", by.y = 0, all.x = T, all.y = T)
  colnames(res)[11] <- "closeness"
  res <- merge(res, betweenness_result, by.x = "ID", by.y = 0, all.x = T, all.y = T)
  colnames(res)[12] <- "betweenness"
  min_eigen <- min(res$eigenvalue) %>% log()
  res$eigenvalue <- log(res$eigenvalue)
  
  res <- merge(res, vip$mrs_table, by.x = "ID", by.y = "TF", all.x = T)
  res$eigenvalue <- sapply(res$eigenvalue, function(eigen){
    test <- (is.na(eigen) | eigen == -Inf)
    if (test){
      min_eigen
    }else{
      eigen
    }
  })
  res$Page_Rank <- log(res$Page_Rank)
  res <- dplyr::filter(res, ID %in% c(reg$tf, reg$target))
  deal_with_missing_values(res, c("logFC", "AveExpr", "t", "B", "nes", "size"), 0)
  deal_with_missing_values(res, c("P.Value", "pval"), 1)
  list("features" = res, "net" = reg)
}

Do_cool_scatterplot <- function(Feature, title, filtered_plot = T){
  min_eigen <- min(na.omit(Feature$eigenvalue))
  min_Page_rank <- min(Feature$Page_Rank)
  Feature <- dplyr::filter(Feature, eigenvalue > min_eigen & Page_Rank > min_Page_rank )#& ((P.Value < 0.05 & abs(logFC) > 1.5) | pval < 0.05))
  DEG <- sapply(1:nrow(Feature), function(n){
    res <- "NoSign"
    if (Feature[n,"P.Value"] < 0.1){
      if (Feature[n, "logFC"] > 0){
        res <- "UP"
      }else{
        res <- "DOWN"
      }
    }
    if(!is.na(Feature[n, "pval"])){
      if(Feature[n, "pval"] < 0.1){
        if (Feature[n, "nes"] > 0){
          res <- "TF_UP"
        }else{
          res <- "TF_DOWN"
        }
      }
    }
    res
  })
  color <-  c("#00FF00", "#888888", "#0000FF", "#FF0000", "#FF00FF")
  if(filtered_plot){
    Feature <- Feature[DEG != "NoSign",]
    DEG <- DEG[DEG != "NoSign"]
    color =  c("#00FF00", "#0000FF", "#FF0000", "#FF00FF")
  }
  ggplot(Feature, aes(x = Page_Rank, y = eigenvalue, label = ID, colour = DEG))+
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

From_network_to_TF_activity <- function(ms_vip, DEG_analysis, output_folder, filtered = F, title){
  features <- Do_net_features(ms_vip, DEG_analysis)
  scat <- Do_cool_scatterplot(features[["features"]], title = title, filtered)
  folder_output <- output_folder
  dir.create(folder_output, showWarnings = F)
  
  png(paste0(folder_output, "/GRN_Network_based_nodes.png"), width = 1280, height = 720)
  plot(scat)
  dev.off()
  
  ms_vip_table <- ms_vip$mrs_table %>%
    dplyr::filter(TF %in% c(features[["net"]]$tf, features[["net"]]$target))
  min_eigen <- min(na.omit(features[["features"]]$eigenvalue))
  deal_with_missing_values(features[["features"]], c("eigenvalue"), min_eigen)
  deal_with_missing_values(features[["features"]], c("nes", "size", "closeness"), 0)
  deal_with_missing_values(features[["features"]], c("pval", "pval.fdr"), 1)
  res <- list("ms_vip" = ms_vip_table, "features" = features[["features"]], "scat" = scat, "mrs" = ms_vip$mrs, "net" = features[["net"]])
}

Make_Cytoscape_network <- function(net, feat, title, collection){
  Genes_n_TFs <- feat %>% 
    dplyr::filter((abs(logFC) > 1.5 & P.Value < 0.05) | pval < 0.1) %>% 
    unlist(.$ID) %>% 
    unique()
  filtered <- net %>% dplyr::filter(tf %in% Genes_n_TFs | target %in% Genes_n_TFs)
  filtered_feat <- feat %>% 
    dplyr::filter(ID %in% filtered$tf | ID %in% filtered$target)
  colnames(filtered)[1:2] <- c("source", "target")
  colnames(filtered_feat)[1] <- c("id")
  filtered <- filtered %>%
    dplyr::filter(source %in% filtered_feat$id & target %in% filtered_feat$id)
  filtered_feat$TF <- ifelse(filtered_feat$nes != 0, "TF", "Gene")
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

  setNodeShapeMapping(table.column = "TF",
                      table.column.values = c("Gene", "TF"),
                      shapes = c('diamond', 'ELLIPSE'),
                      style.name = title)
  
  setNodeColorMapping(table.column = 'logFC',
                      table.column.values = c(min_logFC, 0, max_logFC),
                      colors = c('#0000FF', '#FFFFFF', '#FF0000'),
                      style.name = title)
  
  setNodeFillOpacityMapping(table.column = 'P.Value',
                            table.column.values = c(0, 0.1, 1),
                            opacities = c(255, 200, 20),
                            style.name = title)
  
  createColumnFilter(filter.name = "TF_diff", column = "pval", criterion = 0.1, type = "nodes", predicate = "LESS_THAN")
  setNodeOpacityBypass(getSelectedNodes(), new.values = 255)
  
  setNodeBorderOpacityMapping(table.column = 'P.Value',
                              table.column.values = c(0, 0.1, 1),
                              opacities = c(255, 200, 20),
                              style.name = title)
  
  setNodeLabelOpacityMapping(table.column = 'P.Value',
                             table.column.values = c(0, 0.1, 1),
                             opacities = c(255, 200, 20),
                             style.name = title)
  
  setEdgeTargetArrowShapeDefault(style.name = title, 
                                 new.shape = "ARROW")
  
  setEdgeColorMapping(table.column = 'mor', 
                      table.column.values = c(-1, 0, 1), 
                      colors = c("#FF0000", "#FFFFFF", "#0000FF"), 
                      style.name = title)
  
  min_eigen <- min(n$eigenvalue)
  max_eigen <- max(n$eigenvalue)
  
  setNodeSizeMapping (table.column = 'eigenvalue',
                      table.column.values = c(min_eigen, max_eigen),
                      sizes = c(10, 200),
                      style.name = title)
  
  setNodeFontSizeMapping(table.column = 'eigenvalue',
                         table.column.values = c(min_eigen, max_eigen),
                         sizes = c(15, 75),
                         style.name = title)
  
  layoutNetwork()
  clearSelection()
  exportImage(paste0("~/GitHub/Thesis_paper/Results/Network_analysis/", collection, "_", title), 'SVG', zoom=200)
}
