# Set up R error handling to go to stderr
options(show.error.messages = FALSE,
        error = function() {
        cat(geterrmessage(), file = stderr())
        q("no", 1, FALSE)})
# Avoid crashing Galaxy with an UTF8 error on German LC settings
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
# Import required libraries and data
suppressPackageStartupMessages({
  library(netboxr)
  library(igraph)
  library(RColorBrewer)
})

data(netbox2010)
args <- commandArgs(TRUE)
# Vars
gene_list <- scan(args[2], what = character(), sep = "\n")
cutoff <- args[4]
community <- args[6]
global_model <- args[8]
global_iterations <- args[10]
global_number <- args[12]
local_model <- args[14]
local_iterations <- args[16]

network_plot <- args[18]
plot_width <- args[20]
output_sif <- args[22]
neighbor_list <- args[24]
modmem <- args[26]
nt <- args[28]

sink("metadata.txt")
sink(stdout(), type = "message")
# Network analysis as described in netboxr vignette
sif_network <- netbox2010$network
graph_reduced <- networkSimplify(sif_network, directed = FALSE)
threshold <- cutoff
results <- print(geneConnector(geneList = gene_list, networkGraph = graph_reduced,
                         directed = FALSE, pValueAdj = "BH", pValueCutoff = threshold,
                         communityMethod = community, keepIsolatedNodes = FALSE))

# Check the p-value of the selected linker
linker_df <- results$neighborData
linker_df[linker_df$pValueFDR < threshold, ]
graph_layout <- layout_with_fr(results$netboxGraph)

# Global Network Null Model
if (global_model) {
  global_test <- globalNullModel(netboxGraph = results$netboxGraph, networkGraph = graph_reduced,
                                iterations = global_iterations, numOfGenes = global_number)
  global_test
}

# Local Network Null Model
if (local_model) {
  local_test <- localNullModel(netboxGraph = results$netboxGraph, iterations = local_iterations)
  local_test
}

## Output
# Plot the edge annotated graph
if (network_plot) {

  edges <- results$netboxOutput
  interaction_type <- unique(edges[, 2])
  interaction_type_color <- brewer.pal(length(interaction_type), name = "Spectral")
  edge_colors <- data.frame(interaction_type, interaction_type_color, stringsAsFactors = FALSE)
  colnames(edge_colors) <- c("INTERACTION_TYPE", "COLOR")
  netbox_graph_annotated <- annotateGraph(netboxResults = results, edgeColors =
                                          edge_colors, directed = FALSE, linker = TRUE)
  pdf("network_plot.pdf", width = plot_width)
  plot(results$netboxCommunity, netbox_graph_annotated, layout = graph_layout,
       vertex.size = 10, vertex.shape = V(netbox_graph_annotated)$shape, edge.color
       = E(netbox_graph_annotated)$interactionColor, edge.width = 3)

  # Add interaction type annotations
  legend(x = -1.8, y = -1, legend = interaction_type, col =
           interaction_type_color, lty = 1, lwd = 2, bty = "n", cex = 1)
  dev.off()
}

# Local Network Null Model
if (local_model) {
  pdf("localModel_histogram.pdf")
  h <- hist(local_test$randomModularityScore, breaks = 35, plot = FALSE)
  h$density <- h$counts / sum(h$counts)
  plot(h, freq = FALSE, ylim = c(0, 0.1), xlim = c(0.1, 0.6), col = "lightblue")
  abline(v = local_test$modularityScoreObs, col = "red")
  dev.off()
}

# NetBox algorithm output in SIF format.
if (output_sif) {
  write.table(results$netboxOutput, file = "network.sif", sep = "\t", quote = FALSE,
              col.names = FALSE, row.names = FALSE)
}

# Save neighbor data
if (neighbor_list) {
  write.table(results$neighborData, file = "neighbor_data.txt", sep = "\t",
              quote = FALSE, col.names = TRUE, row.names = FALSE)
}

#Save identified pathway module numbers
if (modmem) {
 write.table(results$moduleMembership, file = "community.membership.txt", sep = "\t",
              quote = FALSE, col.names = FALSE, row.names = FALSE)
}

# Save file that indicates whether the node is a 'linker' or 'candidate'
if (nt) {
  write.table(results$nodeType, file = "nodeType.txt", sep = "\t", quote = FALSE, col.names = FALSE,
              row.names = FALSE)
}
sink(NULL)
