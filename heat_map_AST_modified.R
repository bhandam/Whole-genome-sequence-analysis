# SETUP ------------------------------------------------------------------------
## Load packages
library(tidyverse)
library(patchwork)    # For combining multiple plots into one

## Define input files
suscep_file <- "results/amrfinderplus/pheno_ast_data/heat_map/susceptibility_data.csv"
antibiotic_file <- "results/amrfinderplus/pheno_ast_data/heat_map/antibiotics.csv"

## Define output files
plotfile <- "results/amrfinderplus/pheno_ast_data/heat_map/heatmap.png"

## Set ggplot theme
#useful-you should not change the theme of difefrent plot
my_theme <- theme_bw(base_size = 13) +
  theme(axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_text(vjust = 0, size = 15),
      axis.title.y = element_text(vjust = 2, size = 15))
      
##override the default theme
theme_set(my_theme)

## Define susceptibility levels and color
suscep_levels_long <- c("Susceptible", "Intermediate", "Resistant")
suscep_levels <- c("S", "I", "R")
suscep_colors <- c("#e0ebdf", "#ffaa00", "darkred")
names(suscep_colors) <- suscep_levels
print(suscep_colors)

# PREP DATA --------------------------------------------------------------------
## Dataframe with antibiotic to antibiotic class mapping
anti_cats <- read.csv(antibiotic_file)

## Read susceptibility data
mock_data_raw <- read.csv(suscep_file, row.names = 1) %>%
  rownames_to_column("Isolate")

## Separate serotype data
sero <- mock_data_raw %>% select(Isolate, Serotype)

## Main dataframe
mock_data <- mock_data_raw %>%
  #drop serotype column
  select(-Serotype) %>% 
  pivot_longer(cols = -Isolate,
               names_to = "Antibiotic", values_to = "Susceptibility") %>%
  left_join(anti_cats, by = "Antibiotic") %>% 
  mutate(Susceptibility = factor(Susceptibility, levels = suscep_levels),
         Antibiotic = fct_inorder(Antibiotic),
         Isolate = sub("AG21_", "", Isolate),
         Class2 = Class)
## In the `Class` 2 column, repeated Classes are only listed once
#becuase we are putting class only once at the top of the plot
mock_data$Class2[which(duplicated(mock_data$Class))] <- ""

## Check antiobiotics classes:
mock_data %>% pull(Class) %>% unique() %>% sort()


# MAKE PLOTS -------------------------------------------------------------------
## Main plot
p_main <- ggplot(mock_data) +
  aes(x = Abbreviation, y = Isolate, fill = Susceptibility) +
  geom_tile(height = 0.9, width = 0.9, color = "black", size = 0.4) +
  scale_fill_manual(values = suscep_colors, labels = suscep_levels_long) +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme(plot.margin = margin(0, 0, 7, 7), 
        axis.text.x = element_text(angle = 50, hjust = 0, size = 12, face="bold", colour="black"),
        axis.text.y=element_text(face="bold", size=15, colour="black"))
  
        #axis.text.y = element_blank()) # Removes the isolate labels

## Antibiotics categories
p_cat <- mock_data %>%
  mutate(mock_column = "X") %>% 
  ggplot() +
  aes(x = Antibiotic, y = mock_column, fill = Class) +
  geom_tile() +
  geom_text(aes(label = substr(Class2, 1, 3)),
            size = 4, fontface = "plain", color = "grey20") +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  labs(x = "Antibiotic & Class", y = NULL) +
  guides(fill = "none") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(5.5, 5.5, 0, 5.5))
p_cat

## Serovars
p_sero <- sero %>%
  mutate(mock_column = "X") %>% 
  ggplot() +
  aes(x = mock_column, y = Isolate, fill = Serotype) +
  scale_fill_brewer(palette = "Dark2") +
  geom_tile() +
  coord_fixed() +
  labs(x = NULL, y = "Isolate/Serovar") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 0))

## Combine all the plots
## See https://patchwork.data-imaginist.com/articles/guides/layout.html
p_fill <- ggplot() + theme_void() # Empty filler plot
p_fill +  p_sero + p_cat + p_main + plot_layout(guides = "collect")
#p_cat 

## Save the final figure
ggsave(plotfile, width = 8, height = 6)
