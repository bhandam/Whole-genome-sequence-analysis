mut_files <- list.files("/fs/ess/PAS0471/jelmer/assist/2022-04_menuka/results/ncbi_genomes/amrfinderplus", pattern = "*_genomic.txt", full.names = TRUE)
a <-read_tsv(mut_files, show_col_types = FALSE)
##To select few columns
c<-a %>% select(1, 7:13)

#Join the column with the our sample with the NCBI sample
d<-data.frame(unique(c$`Gene symbol`))

my_files<-list.files("/fs/project/PAS0471/Menuka/WGS_Menuka/results/amrfinderplus",
                     pattern = "*[0-9].txt", full.names = TRUE)


my_AMR <-read_tsv(my_files, show_col_types = FALSE)
our_sample<-my_AMR %>% select(1, 7:13)
merge_our_sample_ncbi<-rbind(c, our_sample)
merge_our_sample_ncbi$Name<- sub("\\_[A-Z].*", "", merge_our_sample_ncbi$Name)
total_AMR_data<-merge_our_sample_ncbi %>% rename(Assembly = Name)
## Add metadata and merge with the AMR genes file
dir<-"/fs/ess/PAS0471/Menuka/WGS_Menuka"
setwd("/fs/ess/PAS0471/Menuka/WGS_Menuka")
metadata<-read.csv("metadata_for_PCA.csv")

## Data for PCA
PCA_data<-left_join(total_AMR_data, metadata, by="Assembly")
PCA_data<-PCA_data %>% rename(gene = "Gene symbol")
write_xlsx(PCA_data, "PCA_total_data.xlsx")
b<-PCA_data  %>% group_by('Assembly', 'gene')%>% count()

length(PCA_data$gene)

genes<- PCA_data %>% group_by(Assembly, gene)%>% summarise()

d<-pivot_wider(b, names_from = gene, values_from = n)
d[is.na(d)] <- 0
d[d ==2] <-1
d


#Codes from another R
getwd()
library(tidyr)
library(readxl)
library(writexl)
library(dplyr)
library(ggplot2)
a<-read_excel("PCA_total_data.xlsx")
meta<- read.csv("metadata_PCA.csv")
## Join two dataset to find genes present just in humans and environment
genes<-left_join(b, meta, by ="Assembly") 

##Split the dataframe to count number of genes
envi<-split(genes, genes$source)[[1]]
human<-split(genes, genes$source)[[2]]
our_sample<-split(genes, genes$source)[[3]]
## Unique genes in envi_and human sample
envi_unique<-data.frame(unique(envi$gene))|> rename(gene="unique.envi.gene.")
human_gene<-data.frame(unique(human$gene)) |> rename(gene="unique.human.gene.")
our_sample<- data.frame(unique(our_sample$gene)) |> rename(gene="unique.our_sample.gene.")
our_sample<- our_sample |> rename(gene="unique.our_sample.gene.")
## common genes between three sample
common_env_human<- intersect(envi_unique, human_gene)
common_our_human<-intersect(our_sample, human_gene)
common_our_env<-intersect(our_sample, envi_unique)

## save common genes for venn diagram
write_xlsx(common_env_human, "common_env_human.xlsx")
write_xlsx(envi_unique, "unique_env.xlsx")
write_xlsx(common_our_human, "common_our_human.xlsx")
write_xlsx(human_gene, "unique_human.xlsx")
write_xlsx(common_our_env, "common_our_env.xlsx")

## To filter out four samples which were porcine and missing
b<- a |>  group_by(Assembly, gene) |>  count()
d<-pivot_wider(b, names_from = gene, values_from = n)
##If any gene is unavailable/NA convert it to the 0
d[is.na(d)] <- 0
## If some samples have gene more than one then convert them to one

## Convert 1st column to row
e <- data.frame(d, row.names = 1)
pca <- prcomp(e, scale = TRUE)
summary(pca)
head(pca$x)
pca$x
m<-data.frame(pca$x)
n<-m |> tibble::rownames_to_column("Assembly")
o<-n |> select(1:3)
#write_xlsx(o, "total_samples_with_AMR_genes.xlsx")

## join PCA data to the category of environment_human and our sample
pca_scores <- left_join(o, meta, by ="Assembly")
#new_pca_scores<-pca_scores |>  select(1, 2:20, 108)
#new_pca_scores<-pca_scores |> rename(Source=Isolation.type) 
new_df <- subset(pca_scores, PC1> -15) 
new_df_2 <- filter(new_df, PC2<10, PC2>-10, PC1<20)

bad_samples <-  filter(new_df, PC2>10 | PC2 < -10 | PC1 > 20) %>%
  pull(Assembly)

score_plot <- ggplot(new_df_2) +
  geom_point(aes(x = PC1, y = PC2, color =source))+#,
             #shape = 1) +
  geom_point(data = filter(new_df_2, source == "Our_Sample"),
             aes(x = PC1, y = PC2), color = "blue") +
  theme_classic()+
  theme( axis.text.x = element_text(size = 15, face="bold", colour="black"),
         axis.text.y=element_text(face="bold", size=15, colour="black"),
         axis.title.x = element_text(face="bold", size=15, colour="black"),
         axis.title.y = element_text(face="bold", size=15, colour="black"),
         legend.title = element_text(face = "bold",size=15, colour="black" ),
         legend.text = element_text(face = "bold",size=10, colour="black"),
         legend.direction = "vertical")
score_plot
score_plot

## to find the genes that were common and different in the human, env and our sample

genes<-left_join(b, metadata, by ="Assembly")
genes<-genes |> rename(source= "Isolation.type")
envi<-split(genes, genes$source)[[1]]
human<-split(genes, genes$source)[[2]]
our_sample<-split(genes, genes$source)[[3]]

colnames(our_sample)
## remove n column
envi_gene<-data.frame(unique(envi$gene))|> rename(gene="unique.envi.gene.")
human_gene<-data.frame(unique(human$gene)) |> rename(gene="unique.human.gene.")
our_sample<- our_sample |> rename(gene="unique.our_sample.gene.")

## common genes between three sample
common_env_human<- intersect(envi_gene, human_gene)
common_our_human<-intersect(our_sample, human_gene)
common_our_env<-intersect(our_sample, envi_gene)

## save common genes for venn diagram
write_xlsx(common_env_human, "common_env_human.xlsx")
write_xlsx(envi_gene, "unique_env.xlsx")
write_xlsx(common_our_human, "common_our_human.xlsx")
write_xlsx(human_gene, "unique_human.xlsx")
write_xlsx(common_our_env, "common_our_env.xlsx")

##venn diagram in R, venn diagram was having scaling problem:
k<-list(envi_gene, human_gene, our_sample)
r<-venn.diagram(x=c(envi_gene, human_gene, our_sample), category.names = c("Environment" , "Human " , "our_sample"),
                fill=c("#E69F00", "#56B4E9", "#009E73"),
                filename = "venn.png", output=TRUE)
r
a
##https://biodash.github.io/codeclub/s03e03_pca/






## To match pattern
c$Name<- sub("\\_[A-Z].*", "", c$Name)
write_xlsx(c, "NCBI_AMR_finder_2751__genes.xlsx")

## Separate clinical samples from other environmental samples

clinical_samples_ID<-c("GCA_006650465.1",
"GCA_006822015.1", 
"GCA_002056945.1",
"GCA_002057025.1",
"GCA_009242365.1",
"GCA_002033735.1",
"GCA_010908165.1",
"GCA_011083625.1",
"GCA_010530695.1",
"GCA_016105225.1",
"GCA_000329285.1",
"GCA_000329505.1",
"GCA_016105065.1",
"GCA_000284775.2",
"GCA_000430105.1",
"GCA_009648795.1",
"GCA_006713205.1",
"GCA_006769545.1",
"GCA_006713235.1",
"GCA_006650495.1",
"GCA_010661665.1",
"GCA_011412755.1",
"GCA_025218965.1",
"GCA_002066895.1",
"GCA_000329365.2",
"GCA_004185475.1",
"GCA_006769725.1")

human_samples_genes<- c%>% filter(Name %in% clinical_samples_ID)
unique_human_genes<-data.frame(unique(human_samples_genes$`Gene symbol`))
environmental_samples_genes<- c%>% filter(!Name %in% clinical_samples_ID)
unique_environment_genes<-data.frame(unique(environmental_samples_genes$`Gene symbol`))
env_gene<-rename(unique_environment_genes, Gene_symbol=unique.environmental_samples_genes..Gene.symbol..)
human_genes<- rename(unique_human_genes, Gene_symbol=unique.human_samples_genes..Gene.symbol..)
common_genes<-inner_join(human_genes, env_gene)
write_xlsx(common_genes, "NCBI_AMR_finder_2751_common_genes.xlsx")
write_xlsx(human_genes, "NCBI_AMR_finder_human_genes.xlsx")
write_xlsx(env_gene, "NCBI_AMR_finder_env_genes.xlsx")
