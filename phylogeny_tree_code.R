library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(tidytree)
library(ggstar)
library(ggforce)
library(RColorBrewer)

## Define input files
meta_file <- "results_29_new_seq/ksnp3/total_meta_data_phylogeny.csv"
tree_file <- "results_29_new_seq/ksnp3/tree.core.tre"

## Define output files
tree_plotfile <- "results_29_new_seq/ksnp3/tree.png"


# PLOT THE TREE ----------------------------------------------------------------
## Read and prep metadata
meta <- read_csv(meta_file, show_col_types = FALSE) |>
  janitor::clean_names() |> 
  rename(Serotype = serotype) |>
  mutate(year = factor(year),
         source = factor(source),
         farm=factor(farm),
         mlst=factor(mlst),
         serogroup=factor(serogroup))

## Read the tree file and prep the tree
tree <- read.tree(tree_file)
#tree$tip.label <- sub("AG21-0", "", tree$tip.label)
nseqs <- length(tree$tip.label)
tree$node.label <- round(as.numeric(tree$node.label) * 100)

## Define colors
#randomcoloR::distinctColorPalette(20)
mycols <- c("#a74bb4", "#62b54f", "#7064d3", "#b5b348", "#dd6fc5",
            "#4db18a", "#ce416e", "#45aecf", "#d55035", "#7784cb",
            "#cc8e3e", "#ac6090", "#647934", "#df837e", "#9c5736")
sero_cols <- mycols[1:length(unique(meta$Serotype))]
year_cols <- c("#595959", "#B30000", "#020099", '#1B9E77', "#7570b3", "#9c5736")
source_cols <- c("#a74bb4", "#62b54f", "#7064d3", "#b5b348")
farm_cols<-brewer.pal(n = 8, name = "Paired")
mlst_cols<- c("#a74bb4", "#62b54f", "#7064d3", "#b5b348", "#dd6fc5",
              "#4db18a", "#ce416e", "#45aecf", "#d55035", "#7784cb",
              "#cc8e3e", "#ac6090", "#647934", "#df837e", "#9c5736",
              "#595959", "#B30000", "#020099")
serogroup_cols<-brewer.pal(n = 7, name = "Dark2")

## Make the plot
plot<-ggtree(tree, layout = "circular") %<+% meta +
  geom_tiplab(aes(color = Serotype),
              size = 3, offset = .1, align = TRUE, fontface = 'bold') +
  scale_color_manual(values = sero_cols) +
  geom_fruit(geom = geom_tile,
             mapping = aes(y = genome, fill = year),
             offset = 1, pwidth = 0.05, width = 0.05) +
  scale_fill_manual(
    name = "Isolation Year",
    values = year_cols,
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 1, nrow=1)
  ) +
  new_scale_fill() +
  geom_fruit(geom = geom_tile,
             mapping = aes(y = genome, fill = source),
             offset = 0.3, pwidth = 0.05, width = 0.05) +
  scale_fill_manual(
    name = "Isolation Source",
    values= source_cols,
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 2, nrow=1))+
  new_scale_fill() +
  geom_fruit(geom = geom_tile,
             mapping = aes(y = genome, fill = farm),
             offset = 0.3, pwidth = 0.05, width = 0.05) +
  scale_fill_manual(
    name = "Farm",
    values= farm_cols,
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 3, nrow=2))+
theme(
  legend.title=element_text(size=9, face="bold"), 
  legend.text=element_text(size=8),
  legend.spacing.y = unit(0.02, "cm")) +
  new_scale_fill() +
  geom_fruit(geom = geom_tile,
             mapping = aes(y = genome, fill = mlst),
             offset = 0.3, pwidth = 0.05, width = 0.05) +
  scale_fill_manual(
    name = "MLST",
    values= mlst_cols,
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 4, nrow = 4))+
  theme(
    legend.title=element_text(size=9, face="bold"), 
    legend.text=element_text(size=8),
    legend.spacing.y = unit(0.02, "cm")) +
  new_scale_fill() +
  geom_fruit(geom = geom_tile,
             mapping = aes(y = genome, fill = serogroup),
             offset = 0.3, pwidth = 0.05, width = 0.05) +
  scale_fill_manual(
    name = "Serogroup",
    values= serogroup_cols,
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 5, nrow = 2))+
  theme(
    legend.title=element_text(size=9, face="bold"), 
    legend.text=element_text(size=8),
    legend.spacing.y = unit(0.02, "cm")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))


plot
