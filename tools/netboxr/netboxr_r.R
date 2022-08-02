# Set up R error handling to go to stderr
options(show.error.messages = FALSE,
        error = function() {
        cat(geterrmessage(), file = stderr())
        q("no", 1, FALSE)})
# Avoid crashing Galaxy with an UTF8 error on German LC settings
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
# Import required libraries and data
suppressPackageStartupMessages({
  library(optparse)
  library(netboxr)
  library(igraph)
  library(RColorBrewer)
})

data(netbox2010)
option_list <- list(
  make_option("--geneList", type = "character", help = "Tab-delimited list of genes of interest"),
  make_option("--cutoff", type = "double", help = "p-value cutoff value"),
  make_option("--community", type = "character", help = "community detection method"),
  #make_option("--resolutionParam", type = "integer", help = "community size"),
  #make_option("--networkType", type = "character", help = "edge weights"),
  #make_option("--weightsInput", type = "", help = "edge weights"),
  make_option("--globalModel", type = "logical", help = "Used to assess the global connectivity
              (number of nodes and edges) of the largest module in the identified network
              compared with the same number but randomly selected gene list"),
  make_option("--globalIterations", type = "integer", help = "Global model iterations"),
  make_option("--globalNumber", type = "integer", help = "Global model number of genes"),
  make_option("--localModel", type = "logical", help = "Used to assess the network modularity in
              the identified network compared with random re-wired network"),
  make_option("--localIterations", type = "integer", help = "Local model iterations"),

  make_option("--networkPlot", type = "logical", help = "Plot of edge-annotated netboxr graph"),
  make_option("--plotWidth", type = "integer", help = "Plot width"),
  make_option("--outputSIF", type = "logical", help = "NetBox algorithm output in SIF format"),
  make_option("--neighborList", type = "logical", help = "Contains information of all neighbor nodes"),
  make_option("--modmem", type = "logical", help = "Identified pathway module numbers"),
  make_option("--nt", type = "logical", help = "Indicates whether node is linker or candidate")
)
parser <- OptionParser(usage = "%prog [options] file", option_list =
                         option_list)
args <- parse_args(parser)
# Vars
gene_list <- scan(args$geneList, what = character(), sep = "\n")
cutoff <- args$cutoff
community <- args$community
global_model <- args$globalModel
global_iterations <- args$globalIterations
global_number <- args$globalNumber
local_model <- args$localModel
local_iterations <- args$localIterations

network_plot <- args$networkPlot
plot_width <- args$plotWidth
output_sif <- args$outputSIF
neighbor_list <- args$neighborList
modmem <- args$modmem
nt <- args$nt

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
