# 01_Diversity_Calculation.R
# Taxonomic, functional, and phylogenetic diversity (Hill q = 1) and LCBD
# Packages: hillR, vegan, ape, cluster, adespatial, SYNCSA

library(hillR)
library(vegan)
library(ape)
library(cluster)
library(adespatial)
library(SYNCSA)

comm   <- read.table("abundance.txt",     header=TRUE, sep="\t", row.names=1, check.names=FALSE)
traits <- read.table("traits.txt",        header=TRUE, sep="\t", row.names=1, stringsAsFactors=TRUE)
tree   <- read.tree("phylogeny.nwk")
meta   <- read.table("site_metadata.txt", header=TRUE, sep="\t", row.names=1)

tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$edge.length[tree$edge.length <= 0] <- 1e-6

common_taxa <- intersect(intersect(colnames(comm), rownames(traits)), tree$tip.label)
comm   <- comm[, common_taxa]
traits <- traits[common_taxa, , drop=FALSE]
tree   <- keep.tip(tree, common_taxa)

# Alpha diversity (Hill q = 1)
comm_rel   <- decostand(comm, method="total", MARGIN=1)
gower_dist <- daisy(traits, metric="gower")

TD     <- hill_taxa(comm_rel, q=1)
FD_mat <- hill_func(comm_rel, gower_dist, traits_as_is=TRUE, q=1)
FD     <- FD_mat["FD_q", ]
PD     <- hill_phylo(comm_rel, tree, q=1)

alpha_result <- data.frame(
  Site = rownames(comm),
  TD   = as.numeric(TD[rownames(comm)]),
  FD   = as.numeric(FD[rownames(comm)]),
  PD   = as.numeric(PD[rownames(comm)]))
write.table(alpha_result, "alpha_diversity.txt", sep="\t", row.names=FALSE, quote=FALSE)

# LCBD
compute_lcbd <- function(mat) {
  mat <- mat[rowSums(mat) > 0, , drop=FALSE]
  hel     <- decostand(mat, method="hellinger")
  res_tax <- beta.div(hel, method="euclidean", nperm=999)
  sp_f    <- intersect(colnames(mat), rownames(traits))
  mx      <- matrix.x(comm=mat[, sp_f], traits=traits[sp_f, ],
                      scale=TRUE, notification=FALSE)$matrix.X
  res_func <- beta.div(mx, method="chord", nperm=999)
  sp_p    <- intersect(colnames(mat), tree$tip.label)
  tree_p  <- keep.tip(tree, sp_p)
  mp      <- matrix.p(comm=mat[, sp_p], phylodist=cophenetic.phylo(tree_p),
                      notification=FALSE)$matrix.P
  res_phylo <- beta.div(mp, method="chord", nperm=999)
  
  data.frame(
    Site      = rownames(mat),
    LCBDtax   = res_tax$LCBD   * 1000, P_tax    = res_tax$p.LCBD,
    LCBDfunc  = res_func$LCBD  * 1000, P_func   = res_func$p.LCBD,
    LCBDphylo = res_phylo$LCBD * 1000, P_phylo  = res_phylo$p.LCBD
  )
}

comm_mat <- as.matrix(comm)

# National
write.table(compute_lcbd(comm_mat), "LCBD_national.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Large / small rivers
for (rtype in c("large", "small")) {
  sites_sub <- intersect(rownames(meta)[meta$river_type == rtype], rownames(comm_mat))
  write.table(compute_lcbd(comm_mat[sites_sub, ]),
              paste0("LCBD_", rtype, "_river.txt"), sep="\t", row.names=FALSE, quote=FALSE)}

# Basin level
res_basin <- do.call(rbind, lapply(unique(meta$basin), function(b) {
  sites_b <- intersect(rownames(meta)[meta$basin == b], rownames(comm_mat))
  df <- compute_lcbd(comm_mat[sites_b, ]); df$Basin <- b; df
}))
write.table(res_basin, "LCBD_basin.txt", sep="\t", row.names=FALSE, quote=FALSE)