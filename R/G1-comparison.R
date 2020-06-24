library("phyloseq"); packageVersion("phyloseq")
library(qiime2R)
library("ggplot2"); packageVersion("ggplot2")

setwd("/home/ben/Documents/gradients-amplicon/data/g1")


physeq<-qza_to_phyloseq(
  features="pr2-dada2-final.qza",
  taxonomy="g1-taxonomy-pr2-dada2.qza",
  metadata = "dada2-metadata.txt"
)

# Check out the top 20 across samples.
theme_set(theme_bw())

top20 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, fill="Class", facet_grid=~size.fraction)


# For fun lets make a heatmap and see how it looks.
gpt <- prune_taxa(names(sort(taxa_sums(physeq),TRUE)[1:100]), physeq)
plot_heatmap(gpt, "NMDS", "bray", "latitude", "Family" )

# Hopefully the samples should cluster by size and lat.
ps.ord <- ordinate(physeq, "PCoA", "bray")
p2 = plot_ordination(physeq, ps.ord, type="samples", color="latitude") 
# Well that's pretty messy.

# Let's try a network and see what happens.
top20 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:500]
ps.top20 <- prune_taxa(top20, ps.top20)
ig <- make_network(ps.top20, dist.fun = 'jaccard', max.dist=0.4, type = "taxa")
plot_network(ig, ps.top20, label = NULL)
# Nothing super interesting.

