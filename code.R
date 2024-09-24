load('SO.1008cells_14celltypes.RData')


library(tidyverse)
library(Seurat)
library(celltalker)


DefaultAssay(SO) <- "RNA"
so_interactions <- celltalk(input_object=SO,
                               metadata_grouping="celltype",
                               ligand_receptor_pairs=ramilowski_pairs,
                               number_cells_required=10,
                               min_expression=1000,
                               max_expression=20000,
                               scramble_times=5)


## Identify top statistically significant interactions
top_stats <- so_interactions %>%
  mutate(fdr=p.adjust(p_val,method="fdr")) %>%
  filter(fdr<0.05) %>%
  group_by(cell_type1) %>%
  top_n(3,interact_ratio) %>%
  ungroup()

###
source("celltalker_R_utils.R")
# no adj

group_stats <- compare_group_interactions(seurat_object=SO,
                                          interaction_stats=top_stats,
                                          sample_replicates="Sample ID",
                                          sample_groups="Ancestry",
                                          metadata_grouping="celltype")



mod_p_vals <- do.call(rbind,
                      lapply(group_stats,function(x) {
                        if (class(x) != "lm") stop("Not an object of class 'lm' ")
                        f <- summary(x)$fstatistic
                        p <- pf(f[1],f[2],f[3],lower.tail=F)
                        attributes(p) <- NULL
                        return(p)
                      })
)

mod_p_vals <- data.frame(mod_p_vals,p_val_adj=p.adjust(mod_p_vals[,1],method="fdr"))

View(mod_p_vals)



# adj 
group_stats_adj <- compare_group_interactions_adj(seurat_object=SO, adj_var = c("Age_dx","BMI"),
                                          interaction_stats=top_stats,
                                          sample_replicates="Sample ID",
                                          sample_groups="Ancestry",
                                          metadata_grouping="celltype")



mod_p_vals_adj2 <- do.call(rbind,
                      lapply(group_stats_adj,function(x) {
                        if (class(x) != "lm") stop("Not an object of class 'lm' ")
                        f <- summary(x)$fstatistic
                        p <- pf(f[1],f[2],f[3],lower.tail=F)
                        attributes(p) <- NULL
                        return(p)
                      })
)

mod_p_vals_adj2 <- data.frame(mod_p_vals_adj2,p_val_adj=p.adjust(mod_p_vals[,1],method="fdr"))

View(mod_p_vals_adj2)





### https://arc85.github.io/celltalker/articles/consistent_interactions.html

