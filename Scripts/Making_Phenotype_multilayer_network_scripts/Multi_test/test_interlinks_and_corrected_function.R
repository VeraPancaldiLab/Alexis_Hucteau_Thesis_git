rm(DmIALite, DOLite, DOLiteTerm, HsIALite, MmIALite, RnIALite)

library(dplyr)
library(RCy3)
library(igraph)
library(data.table)
library(stringr)
library(clusterProfiler)
library(ggplot2)
library(muxViz)
library(grid)
library(ggraph)
library(RColorBrewer)
library(viridis)

buildMultilayerNetworkFromMuxvizFiles_corrected <-
  function(config.file,
           isDirected,
           isWeighted,
           MultisliceType,
           LayerCouplingStrength = 1,
           format = "muxviz edge-colored",
           verbose = T) {
    if (format == "muxviz edge-colored") {
      # Expected format: layer_file;layer_label;layout_file
      df.config <- utils::read.table(config.file, sep = ";", header = F)
      colnames(df.config) <-
        c("layer.file", "layer.label", "layout.file")
      Layers <- nrow(df.config)
      layerTensor <-
        BuildLayersTensor(
          Layers = Layers,
          OmegaParameter = LayerCouplingStrength,
          MultisliceType = MultisliceType
        )
      layerLabels <- df.config$layer.label
      
      if (verbose)
        cat(paste("Found", Layers, "layers...\n"))
      
      layout.file <- unique(as.character(df.config$layout.file))
      if (length(layout.file) > 1)
        stop("More than one layout file specified.")
      
      #Expected format: nodeID nodeLabel (optional)
      df.layout <- utils::read.table(layout.file, sep = " ", header = T)
      
      nodeIDs <- df.layout$nodeID
      nodeLabels <- df.layout$nodeLabel
      nodeX <- df.layout$nodeX
      nodeY <- df.layout$nodeY
      
      Nodes <- length(nodeIDs)
      if (verbose)
        cat(paste("Found", Nodes, "nodes\n"))
      
      #Read the edge-colored network
      nodeTensor <- list()
      g.list <- list()
      l <- 1
      for (input.file in as.character(df.config$layer.file)) {
        if (verbose)
          cat(paste("  Reading layer from file", input.file, "...\n"))
        edges <- utils::read.table(input.file, header = T)
        if (ncol(edges) == 3) {
          colnames(edges) <- c("from", "to", "weight")
          if (!isWeighted) {
            cat(
              paste(
                "  WARNING! You asked for an unweighted network but weights are found. Assigning 1 by default.\n"
              )
            )
            edges$weight <- 1
          }
        } else if (ncol(edges) == 2) {
          colnames(edges) <- c("from", "to")
          if (isWeighted) {
            cat(
              paste(
                "  WARNING! You assume a weighted network but no weights are found. Assigning 1 by default.\n"
              )
            )
            edges$weight <- 1
          }
        }
        
        g.list[[l]] <-
          igraph::graph.data.frame(edges, directed = isDirected, vertices = nodeIDs)
        #g.list[[l]] <- simplify(g.list[[l]])
        nodeTensor[[l]] <-
          igraph::as_adjacency_matrix(g.list[[l]], attr = "weight")
        l <- l + 1
      }
      
      #Build the multilayer adjacency tensor
      M <-
        BuildSupraAdjacencyMatrixFromEdgeColoredMatrices(nodeTensor, layerTensor, Layers, Nodes)
      
      #Build the aggregate matrix and network
      aggregateTensor <- GetAggregateMatrixFromNetworkList(g.list)
      g.agg <- GetAggregateNetworkFromNetworkList(g.list)
      
      return(
        list(
          Nodes = Nodes,
          nodeIDs = nodeIDs,
          nodeLabels = nodeLabels,
          nodeX = nodeX,
          nodeY = nodeY,
          nodeTensor = nodeTensor,
          Layers = Layers,
          layerIDs = 1:Layers,
          layerLabels = layerLabels,
          layerTensor = layerTensor,
          adjacencyTensor = M,
          g.list = g.list,
          isDirected = isDirected,
          isWeighted = isWeighted,
          aggregateTensor = aggregateTensor,
          g.agg = g.agg
        )
      )
    } else if (format == "muxviz general") {
      #Expected format: layers_file;layer_label;layout_file
      df.config <- utils::read.table(config.file, sep = ";", header = F)
      colnames(df.config) <-
        c("layers.file", "layer.label.file", "layout.file")
      
      layerLabels <-
        utils::read.table(as.character(df.config$layer.label.file), header = T)
      
      #Expected format: node layer node layer weight
      mEdges <-
        utils::read.table(as.character(df.config$layers.file), header = F)
      
      inter.edges <- mEdges[mEdges[, 2] != mEdges[, 4], ]
      intra.edges <- mEdges[mEdges[, 2] == mEdges[, 4], ]
      
      Layers <- max(max(mEdges[, 2]), max(mEdges[, 4]))
      Nodes <- max(max(mEdges[, 1]), max(mEdges[, 3]))
      
      if (verbose)
        cat(paste("Found", Layers, "layers...\n"))
      
      if (nrow(inter.edges) == 0) {
        cat("Warning: no inter-layer links found, input network is edge-colored.\n")
        cat(
          paste(
            "Applying",
            MultisliceType,
            "coupling with intensity",
            LayerCouplingStrength,
            ".\n"
          )
        )
        layerTensor <-
          BuildLayersTensor(
            Layers = Layers,
            OmegaParameter = LayerCouplingStrength,
            MultisliceType = MultisliceType
          )
      }
      layerTensor <-
        BuildLayersTensor(
          Layers = Layers,
          OmegaParameter = LayerCouplingStrength,
          MultisliceType = MultisliceType
        )
      
      layout.file <- unique(as.character(df.config$layout.file))
      if (length(layout.file) > 1)
        stop("Error: More than one layout file specified.")
      
      #Expected format: nodeID nodeLabel (optional)
      df.layout <- utils::read.table(layout.file, sep = " ", header = T)
      
      nodeIDs <- df.layout$nodeID
      nodeLabels <- df.layout$nodeLabel
      nodeX <- df.layout$nodeX
      nodeY <- df.layout$nodeY
      
      Nodes2 <- length(nodeIDs)
      
      if (Nodes != Nodes2) {
        stop("Error: Nodes specified in the layout do not match nodes used in the edges list.\n")
      }
      
      if (verbose)
        cat(paste("Found", Nodes, "nodes\n"))
      
      if (verbose) {
        cat(paste("  Inter-links:", nrow(inter.edges), "\n"))
        cat(paste("  Intra-links:", nrow(intra.edges), "\n"))
      }
      
      #Read the layers
      nodeTensor <- list()
      g.list <- list()
      layerEdges <- list()
      
      for (l in 1:Layers) {
        if (verbose)
          cat(paste("  Reading layer from file", df.config$layers.file, "...\n"))
        
        if (ncol(mEdges) == 5) {
          layerEdges[[l]] <- mEdges[mEdges[, 2] == l &
                                      mEdges[, 4] == l, c(1, 3, 5)]
        } else {
          layerEdges[[l]] <- mEdges[mEdges[, 2] == l & mEdges[, 4] == l, c(1, 3)]
        }
        
        if (ncol(layerEdges[[l]]) == 3) {
          colnames(layerEdges[[l]]) <- c("from", "to", "weight")
          if (!isWeighted) {
            cat(
              paste(
                "  WARNING! You asked for an unweighted network but weights are found. Assigning 1 by default.\n"
              )
            )
            layerEdges[[l]]$weight <- 1
          }
        } else if (ncol(layerEdges[[l]]) == 2) {
          colnames(layerEdges[[l]]) <- c("from", "to")
          if (isWeighted) {
            cat(
              paste(
                "  WARNING! You assume a weighted network but no weights are found. Assigning 1 by default.\n"
              )
            )
            layerEdges[[l]]$weight <- 1
          }
        }
        
        g.list[[l]] <-
          igraph::graph.data.frame(layerEdges[[l]], directed = isDirected, vertices =
                                     nodeIDs)
        #g.list[[l]] <- simplify(g.list[[l]])
        nodeTensor[[l]] <-
          igraph::as_adjacency_matrix(g.list[[l]], attr = "weight")
      }
      
      #Build the multilayer adjacency tensor
      M <-
        BuildSupraAdjacencyMatrixFromExtendedEdgelist(mEdges, Layers, Nodes, isDirected)
      
      #Build the aggregate matrix and network
      aggregateTensor <- GetAggregateMatrixFromNetworkList(g.list)
      g.agg <- GetAggregateNetworkFromNetworkList(g.list)
      
      return(
        list(
          Nodes = Nodes,
          nodeIDs = nodeIDs,
          nodeLabels = nodeLabels,
          nodeX = nodeX,
          nodeY = nodeY,
          nodeTensor = nodeTensor,
          Layers = Layers,
          layerIDs = 1:Layers,
          layerLabels = layerLabels,
          layerTensor = layerTensor,
          adjacencyTensor = M,
          g.list = g.list,
          isDirected = isDirected,
          isWeighted = isWeighted,
          aggregateTensor = aggregateTensor,
          g.agg = g.agg
        )
      )
    } else {
      stop("Format not recognized.")
    }
  }


Multiplex_test <- buildMultilayerNetworkFromMuxvizFiles_corrected(
  config.file = "GitHub/Thesis_paper/Scripts/Making Phenotype_multilayer_network_scripts/Multi_test/config.txt", 
  MultisliceType = "ordinal", isWeighted = T, isDirected = F, verbose = T, format = "muxviz general")

layer.colors <- brewer.pal(8, "Set2")

lay <- layoutMultiplex(Multiplex_test$g.list, layout="fr", ggplot.format=F, box=T)

plot_multiplex3D(Multiplex_test$g.list, layer.layout=lay, layer.colors=layer.colors,
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels="auto", layer.labels.cex=1.5,
                 node.size.values="auto", node.size.scale=0.8,
                 show.aggregate=T)

rgl::snapshot3d("GitHub/Thesis_paper/Scripts/Making Phenotype_multilayer_network_scripts/Multi_test/test_visu_interlinks.png", fmt="png", width = 1920, height = 1080)
