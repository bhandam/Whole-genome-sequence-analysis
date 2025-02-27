##Trimmomatic to trim the adapters
for R1 in new_salmonella_29_seq/*_R1_*fastq.gz; do
    sbatch mcic-scripts/trim/trimmomatic.sh -i "$R1" -o results_29_new_seq/trim \
        -a mcic-scripts/trim/adapters.fa \
        -A "2:30:10" \
        -p "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36"
done

mkdir results_29_new_seq/multiqc_after_trim


##Run multiqc to see if the adapters have been removed
bash scripts/multiqc_WGS.sh results_29_new_seq/trim results_29_new_seq/multiqc_after_trim

##Assemble seqeunces using SPADES
cd mcic-scripts
git pull
# to update the github

for R1 in results_29_new_seq/trim/*_R1_*fastq.gz; do
    sampleID=$(basename "$R1"| sed 's/_.*//')
    sbatch mcic-scripts/assembly/spades.sh -i "$R1" -o results_29_new_seq/spades/"$sampleID" \
    -m isolate -k "21,33,55,77,99,127"
    done

## To look for the requested job
squeue -u menuka

## To run QUAST
for fasta in results_29_new_seq/spades/GRS-S-1/GRS-S-1.fasta
##run fro one sample
sbatch scripts/new_quast.sh results_29_new_seq/quast/double_check/ results_29_new_seq/spades/GRS-S-1/GRS-S-1.fasta 
##Run for all samples
## It didnt work, so i had to type all my sample number in the code itself
#for fasta_seq in results_29_new_seq/spades/GRS-S-*/GRS-S-*.fasta; do
#sampleID=$(basename "$fasta_seq"| sed 's/.fasta//') \

## Final code used for QUAST
sbatch scripts/new_quast.sh results_29_new_seq/quast/all_results

## determinae the coverage
## SAM file contains the seqeunce alignment data, BAM file is the binary version of SAM file
## Run BBmap to get output as Bam file, BAM file is an alignment file needed to get the coverage of the assembly
# Run BBmap
for R1 in results_29_new_seq/trim/*_R1_*.fastq.gz; do
    R2=${R1/_R1_/_R2_}
    sampleID=$(basename "$R1" | sed 's/_S.*//')
    assembly=results_29_new_seq/spades/"$sampleID"/"$sampleID".fasta
    outfile=results_29_new_seq/bbmap/"$sampleID".bam
    sbatch scripts/bbmap.sh "$R1" "$R2" "$assembly" "$outfile"
done

##Run samtools to find te coverage using the bam files
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/samtools 
for bam in results_29_new_seq/bbmap/*bam; do
    sbatch scripts/samtool2.sh "$bam" results_29_new_seq/samtools_cvg
done

## Run checkM
## change extension of all files at once
for a in results_29_new_seq/all_assemblies/*.fna; do 
mv "$a" "${a%.fna}.fasta"; 
done

## checkM
source activate /fs/project/PAS0471/jelmer/conda/checkm-1.2.0
# See https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow

sbatch mcic-scripts/assembly/checkm.sh -i results_29_new_seq/all_assemblies -o results_29_new_seq/checkm

##Run MLST to find the seqeunce type
module load python/3.6-conda5.2
source activate /fs/ess/PAS0471/jelmer/conda/mlst
mkdir results_29_new_seq/mlst_salmonella

for multi_locus in results_29_new_seq/all_assemblies/GRS*; do 
    sample_ID=$(basename "$multi_locus")
    mlst "$multi_locus" --legacy --scheme senterica_achtman_2 >results/mlst_salmonella/"$sample_ID".txt
    done

## Run SIStr for the old and new assembly
##make a common directory
mkdir results/common_assembly
cp results_29_new_seq/all_assemblies/*fasta results/common_assembly
cp results/contigs_SPADES results/common_assembly

## RUN SIStr to find the cgmLST

for cgmlst in results/common_assembly; do 
    sample_ID=$(basename "$cgmlst" | sed 's/.fasta//')
    sbatch mcic-scripts/bact/sistr.sh --indir "$cgmlst" --outdir results/cgmlst_salmonella/"$sample_ID" 
    done

## Scoary for serotype
module load python/3.6-conda5.2
source activate /fs/ess/PAS0471/jelmer/conda/scoary
scoary -g results/roary_1677068211/gene_presence_absence.csv -t results/roary_1677068211/scoary_file.csv
##scoary for isolation source
scoary -g results/roary_1677068211/gene_presence_absence.csv -t results/roary_1677068211/scoary_isolation_scource.csv

## Prokka for the annotation of the genome
for assembly in results_29_new_seq/all_assemblies/*fasta; do
  sbatch mcic-scripts/bact/prokka.sh -i "$assembly" -o results_29_new_seq/prokka --genus Salmonella --species enterica
done

##Assemble seqeunces using SPADES
cd mcic-scripts
git pull
##AMR finder plus
for assembly in results_29_new_seq/all_assemblies/*fasta; do
    aa_fasta=results_29_new_seq/prokka/$(basename "$assembly" .fasta).faa
    gff=results_29_new_seq/prokka/$(basename "$assembly" .fasta).gff
    sbatch mcic-scripts/bact/amrfinderplus.sh --nucleotide "$assembly" --protein "$aa_fasta" --gff "$gff" -o results_29_new_seq/amrfinderplus --organism Salmonella
done

##Virulence finder-abricate
source activate abricate
abricate-get_db --db vfdb --force
for contigs_fasta in results_29_new_seq/all_assemblies/*fasta; do 
    sample_ID=$(basename "$contigs_fasta")
    sbatch scripts/abricate.sh "$contigs_fasta" results_29_new_seq/Abricate_VGF/"$sample_ID".tab vfdb --minid 90
done
source deactivate abricate

##Plasmid find-mobsuite
for asm in results_29_new_seq/all_assemblies/*.fasta; do
    outdir=results_29_new_seq/mobsuite/$(basename "$asm" .fasta)
    sbatch mcic-scripts/bact/mob-suite.sh -i "$asm" -o "$outdir"
done

rm results_29_new_seq/prokka/GRS-S-13*
##Pan genome analysis
##Makedirectory of gff files
for gff_without_seq in results_29_new_seq/prokka/GRS-S-*.gff ; do
gff=results_29_new_seq/prokka/$(basename "$gff_without_seq")
cp $gff results/prokka/gff_file
done

cp results_29_new_seq/prokka/GRS-S-*_withseqs.gff results/prokka/gff_with_seq

## run roary
## Roary needs files with the nucleotide seqeunce at the end
sbatch scripts/roary.sh results/roary results/prokka/gff_with_seq
##Make figures with Roary
wget -O roary_plots.py "https://raw.githubusercontent.com/sanger-pathogens/Roary/master/contrib/roary_plots/roary_plots.py"
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/roary-3.13
python roary_plots.py results/roary_1677068211/accessory_binary_genes.fa.newick results/roary_1677068211/gene_presence_absence.csv
conda install matplotlib seaborn biopython

##Scoary
module load python/3.6-conda5.2
source activate /fs/ess/PAS0471/jelmer/conda/scoary
scoary -t results/roary_1677068211/scoary_1.csv -g results/roary_1677068211/gene_presence_absence.csv 

##Scoary with acessory genes
scoary -t results/roary_1677068211/scoary_1.csv -g results/roary_1677068211/accessory_gene_presence_absence.csv 

## Run kSNP3: to get the core SNPS for the phylogeny tree
#alignment of core SNPs
mkdir -p results_29_new_seq/ksnp3
# making files with the directory and the file name side by side
paste <(ls $PWD/results/contigs_SPADES/*) <(ls results/contigs_SPADES/ | sed 's/.fasta//') > results/ksnp3/assembly_list.txt
#Run the Ksnp3 from the MCIC scripts
sbatch mcic-scripts/trees/ksnp3.sh -i results/ksnp3/assembly_list.txt -o results/ksnp3

##RUN COG using one seqeunce of the prokka
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/cogclassifier
mkdir cog
COGclassifier -i results/prokka/AG21-0050.faa -o cog
##Chnage the default color of the COG
plot_cog_classifier_barchart -i cog/new_count.tsv -o cog/07_barchart_custom_change_color.html

## Run COG on the core genes
#get fasta file of core_gene_scoary.csv
faa=results/prokka/AG21-0050.gff
module load miniconda3
source activate /fs/ess/PAS0471/jelmer/conda/seqkit
mkdir cog/core_gene
seqkit grep -f results/roary_1677068211/core_gene_scoary.txt $faa > cog/core_gene/core_unique.faa
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/cogclassifier
COGclassifier -i cog/core_gene/core_unique.faa -o cog/core_gene
plot_cog_classifier_barchart -i cog/core_gene/classifier_count.tsv -o cog/core_gene/07_barchart_custom_change_color.html

#3Run COG in accessory genes
faa=cog/fasta_merged
module load miniconda3
source activate /fs/ess/PAS0471/jelmer/conda/seqkit
seqkit grep -f ID.txt $faa > cog/accessory_unique.faa

# Subset FASTA to unique serovar genes
##Agbeni
faa=results_29_new_seq/prokka/GRS-S-11.faa
grep -f cog/agbeni_unique.txt $faa
module load miniconda3
source activate /fs/ess/PAS0471/jelmer/conda/seqkit
gff=results/prokka/gff_with_seq/GRS-S-11_withseqs.gff
grep -f cog/agbeni/agbeni_unique.txt $gff | sed -E 's/.*ID=([^;]+);.*/\1/' > cog/agbeni/agbeni_unique_ids.txt
seqkit grep -f cog/agbeni/agbeni_unique_ids.txt $faa > cog/agbeni/agbeni_unique.faa
##Make COG for agbeni
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/cogclassifier
mkdir cog
COGclassifier -i cog/agbeni/agbeni_unique.faa -o cog/agbeni
##Chnage the default color of the COG
##change the color in the file called classifier_count.tsv and save the file and make bar plot
plot_cog_classifier_barchart -i cog/agbeni/new_counts_R.tsv -o cog/agbeni/07_barchart_custom_change_color.html

##Anatum:merge files for Anatum and look for ID AG21-0053, AG21-0054
mkdir cog/anatum
module load miniconda3
source activate /fs/ess/PAS0471/jelmer/conda/seqkit
##merge gff files
cat results/prokka/AG21-0053.gff results/prokka/AG21-0054.gff >cog/anatum/merge_file.txt
##merge fasta files to extract seqeunces
cat results/prokka/AG21-0053.faa results/prokka/AG21-0054.faa > cog/anatum/merge_aa.faa
grep -f cog/anatum/anatum_unique.txt cog/anatum/unique_gene_merge_gff_file.txt | sed -E 's/.*ID=([^;]+);.*/\1/' > cog/anatum/anatum_unique_ids.txt
seqkit grep -f cog/anatum/anatum_unique_ids.txt cog/anatum/merge_aa.faa > cog/anatum/anatum_unique.faa

##Make COG for Anatum
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/cogclassifier
COGclassifier -i cog/anatum/anatum_unique.faa -o cog/anatum
##Chnage the default color of the COG
##change the color in the file called classifier_count.tsv and save the file and make bar plot
plot_cog_classifier_barchart -i cog/anatum/new_counts_R.tsv -o cog/anatum/07_barchart_custom_change_color.html


##COG for accessory genes

seqkit grep -f cog/berta/berta_unique_ids.txt cog/berta/merge_aa.faa > cog/berta/berta_unique.faa

##Berta:merge files for Berta and look for ID AG21-0002 & AG21-0005
mkdir cog/berta
##merge gff files
cat results_29_new_seq/prokka/GRS-S-2.gff results_29_new_seq/prokka/GRS-S-5.gff >cog/berta/merge_gff_file.txt
##merge fasta files to extract seqeunces
cat results_29_new_seq/prokka/GRS-S-2.faa results_29_new_seq/prokka/GRS-S-5.faa > cog/berta/merge_aa.faa
grep -f cog/berta/berta_unique_genes.txt cog/berta/unique_gene_merge_gff_file.txt | sed -E 's/.*ID=([^;]+);.*/\1/' > cog/berta/berta_unique_ids.txt
seqkit grep -f cog/berta/berta_unique_ids.txt cog/berta/merge_aa.faa > cog/berta/berta_unique.faa
##Make COG for Berta
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/cogclassifier
##Make COG for Berta
COGclassifier -i cog/berta/berta_unique.faa -o cog/berta
##Chnage the default color of the COG
##change the color in the file called classifier_count.tsv and save the file and make bar plot
plot_cog_classifier_barchart -i cog/berta/new_counts_R.tsv -o cog/berta/07_barchart_custom_change_color.html

## brenderup:AG21-0061, AG21-0062, AG21-0063, AG21-0064, AG21-0065
mkdir cog/brenderup
##merge gff files 
cat results/prokka/gff_file/AG21-0061.gff results/prokka/gff_file/AG21-0062.gff results/prokka/gff_file/AG21-0063.gff results/prokka/gff_file/AG21-0064.gff results/prokka/gff_file/AG21-0065.gff >cog/brenderup/merge_gff_file.txt
##merge fasta files to extract seqeunces
cat results/prokka/AG21-0061.faa results/prokka/AG21-0062.faa results/prokka/AG21-0063.faa results/prokka/AG21-0064.faa results/prokka/AG21-0065.faa > cog/brenderup/merge_aa.faa
grep -f cog/brenderup/brenderup_unique_genes.txt cog/brenderup/unique_gene_merge_gff_file.txt | sed -E 's/.*ID=([^;]+);.*/\1/' > cog/brenderup/brenderup_unique_ids.txt
seqkit grep -f cog/brenderup/brenderup_unique_ids.txt cog/brenderup/merge_aa.faa > cog/brenderup/brenderup_unique.faa

##Make COG for brenderup
COGclassifier -i cog/brenderup/brenderup_unique.faa -o cog/brenderup
plot_cog_classifier_barchart -i cog/brenderup/new_counts_R.tsv -o cog/brenderup/07_barchart_custom_change_color.html

##Dublin- AG21-0051
mkdir cog/Dublin
##merge gff files
gff=results/prokka/AG21-0051_withseqs.gff
faa=results/prokka/AG21-0051.faa
grep -f cog/Dublin/dublin_unique_genes.txt $gff | sed -E 's/.*ID=([^;]+);.*/\1/' > cog/Dublin/dublin_unique_ids.txt
seqkit grep -f cog/Dublin/dublin_unique_ids.txt $faa > cog/Dublin/dublin_unique.faa
##Make COG for Dublin
COGclassifier -i cog/Dublin/dublin_unique.faa -o cog/Dublin
plot_cog_classifier_barchart -i cog/Dublin/new_counts_R.tsv -o cog/Dublin/07_barchart_custom_change_color.html


##Enteritidis:merge files for Anatum and look for ID AG21-0060, AG21-0066
mkdir cog/Enteritidis
##merge gff files
cat results/prokka/gff_file/AG21-0060.gff results/prokka/gff_file/AG21-0066.gff >cog/Enteritidis/merge_gff_file.txt
##merge fasta files to extract seqeunces
cat results/prokka/AG21-0060.faa results/prokka/AG21-0066.faa > cog/Enteritidis/merge_aa.faa
grep -f cog/Enteritidis/Enteritidis_unique_genes.txt cog/Enteritidis/unique_gene_merge_gff_file.txt | sed -E 's/.*ID=([^;]+);.*/\1/' > cog/Enteritidis/Enteritidis_unique_ids.txt
seqkit grep -f cog/Enteritidis/Enteritidis_unique_ids.txt cog/Enteritidis/merge_aa.faa > cog/Enteritidis/Enteritidis_unique.faa
##Make COG for Enteritidis
COGclassifier -i cog/Enteritidis/Enteritidis_unique.faa -o cog/Enteritidis
plot_cog_classifier_barchart -i cog/Enteritidis/new_counts_R.tsv -o cog/Enteritidis/07_barchart_custom_change_color.html


## Give: AG21-0564, AG21-0565, AG21-0566, AG21-0567, AG21-0568, AG21-0572
mkdir cog/Give
##merge gff files
cat results/prokka/gff_file/AG21-0564.gff results/prokka/gff_file/AG21-0565.gff results/prokka/gff_file/AG21-0566.gff results/prokka/gff_file/AG21-0567.gff results/prokka/gff_file/AG21-0568.gff results/prokka/gff_file/AG21-0572.gff >cog/Give/merge_gff_file.txt
##merge fasta files to extract seqeunces
cat results/prokka/AG21-0564.faa results/prokka/AG21-0565.faa results/prokka/AG21-0566.faa results/prokka/AG21-0567.faa results/prokka/AG21-0568.faa results/prokka/AG21-0572.faa > cog/Give/merge_aa.faa
grep -f cog/Give/Give_unique_genes.txt cog/Give/unique_gene_merge_gff_file.txt | sed -E 's/.*ID=([^;]+);.*/\1/' > cog/Give/Give_unique_ids.txt
seqkit grep -f cog/Give/Give_unique_ids.txt cog/Give/merge_aa.faa > cog/Give/Give_unique.faa
##Make COG for Give
COGclassifier -i cog/Give/Give_unique.faa -o cog/Give
plot_cog_classifier_barchart -i cog/Give/new_counts_R.tsv -o cog/Give/07_barchart_custom_change_color.html


## Hartford: GRS-S-24, GRS-S-25, GRS-S-26, GRS-S-27, GRS-S-28, GRS-S-6, GRS-S-7
mkdir cog/Hartford
##merge gff files
cat results_29_new_seq/prokka/gff/GRS-S-6.gff results_29_new_seq/prokka/gff/GRS-S-7.gff results_29_new_seq/prokka/gff/GRS-S-24.gff results_29_new_seq/prokka/gff/GRS-S-25.gff results_29_new_seq/prokka/gff/GRS-S-26.gff results_29_new_seq/prokka/gff/GRS-S-27.gff results_29_new_seq/prokka/gff/GRS-S-28.gff >cog/Hartford/merge_gff_file.txt
##merge fasta files to extract seqeunces
cat results_29_new_seq/prokka/GRS-S-6.faa results_29_new_seq/prokka/GRS-S-7.faa results_29_new_seq/prokka/GRS-S-24.faa results_29_new_seq/prokka/GRS-S-25.faa results_29_new_seq/prokka/GRS-S-26.faa results_29_new_seq/prokka/GRS-S-27.faa results_29_new_seq/prokka/GRS-S-28.faa > cog/Hartford/merge_aa.faa
grep -f cog/Hartford/hartford_unique_genes.txt cog/Hartford/unique_gene_merge_gff_file.txt | sed -E 's/.*ID=([^;]+);.*/\1/' > cog/Hartford/hartford_unique_ids.txt
seqkit grep -f cog/Hartford/hartford_unique_ids.txt cog/Hartford/merge_aa.faa > cog/Hartford/hartford_unique.faa
##Make COG for Hartford
COGclassifier -i cog/Hartford/hartford_unique.faa -o cog/Hartford
plot_cog_classifier_barchart -i cog/Hartford/new_counts_R.tsv -o cog/Hartford/07_barchart_custom_change_color.html


## Montevideo: AG21-0050, AG21-0057, AG21-0058, AG21-0059, AG21-0592
mkdir cog/Montevideo
##merge gff files
cat results/prokka/gff_file/AG21-0050.gff results/prokka/gff_file/AG21-0057.gff results/prokka/gff_file/AG21-0058.gff results/prokka/gff_file/AG21-0059.gff results/prokka/gff_file/AG21-0592.gff >cog/Montevideo/merge_gff_file.txt
##merge fasta files to extract seqeunces
cat results/prokka/AG21-0050.faa results/prokka/AG21-0057.faa results/prokka/AG21-0058.faa results/prokka/AG21-0059.faa results/prokka/AG21-0592.faa > cog/Montevideo/merge_aa.faa
grep -f cog/Montevideo/montevideo_unique_genes.txt cog/Montevideo/unique_gene_merge_gff_file.txt | sed -E 's/.*ID=([^;]+);.*/\1/' > cog/Montevideo/Montevideo_unique_ids.txt
seqkit grep -f cog/Montevideo/Montevideo_unique_ids.txt cog/Montevideo/merge_aa.faa > cog/Montevideo/Montevideo_unique.faa
##Make COG for Montevideo
COGclassifier -i cog/Montevideo/Montevideo_unique.faa -o cog/Montevideo
plot_cog_classifier_barchart -i cog/Montevideo/new_counts_R.tsv -o cog/Montevideo/07_barchart_custom_change_color.html

## Muenchen: AG21-0589, AG21-0590, AG21-0591
mkdir cog/Muenchen
##merge gff files
cat results/prokka/gff_file/AG21-0589.gff results/prokka/gff_file/AG21-0590.gff results/prokka/gff_file/AG21-0591.gff >cog/Muenchen/merge_gff_file.txt
##merge fasta files to extract seqeunces
cat results/prokka/AG21-0589.faa results/prokka/AG21-0590.faa results/prokka/AG21-0591.faa > cog/Muenchen/merge_aa.faa
gff=results/prokka/gff_file/AG21-0589.gff
faa=results/prokka/AG21-0589.faa
grep -f cog/Muenchen/muenchen_unique_genes.txt $gff | sed -E 's/.*ID=([^;]+);.*/\1/' > cog/Muenchen/muenchen_unique_ids.txt
seqkit grep -f cog/Muenchen/muenchen_unique_ids.txt $faa > cog/Muenchen/Muenchen_unique.faa
COGclassifier -i cog/Muenchen/Muenchen_unique.faa -o cog/Muenchen
plot_cog_classifier_barchart -i cog/Muenchen/new_counts_R.tsv -o cog/Muenchen/07_barchart_custom_change_color.html


## Newport: AG21-0583, AG21-0584, AG21-0585, AG21-0586, AG21-0587, AG21-0588, AG21-0012, AG21-0014, AG21-0004
mkdir cog/Newport
##merge gff files
cat results/prokka/gff_file/AG21-0583.gff \
 results/prokka/gff_file/AG21-0584.gff \
 results/prokka/gff_file/AG21-0585.gff \
 results/prokka/gff_file/AG21-0586.gff \
 results/prokka/gff_file/AG21-0587.gff \
 results/prokka/gff_file/AG21-0588.gff \
 results_29_new_seq/prokka/gff/GRS-S-12.gff \
 results_29_new_seq/prokka/gff/GRS-S-14.gff \
 results_29_new_seq/prokka/gff/GRS-S-4.gff >cog/Newport/merge_gff_file.txt
##merge fasta files to extract seqeunces

COGclassifier -i cog/Newport_4/newport_4_unique.faa -o cog/Newport_4
plot_cog_classifier_barchart -i cog/Newport_4/new_counts_R.tsv -o cog/Newport_4/07_barchart_custom_change_color.html

COGclassifier -i cog/Newport_5/newport_5_unique.faa -o cog/Newport_5
plot_cog_classifier_barchart -i cog/Newport_5/new_counts_R.tsv -o cog/Newport_5/07_barchart_custom_change_color.html

## paratyphi-AG21-0052
mkdir cog/paratyphi
##merge gff files
gff=results/prokka/AG21-0052_withseqs.gff
faa=results/prokka/AG21-0052.faa
grep -f cog/paratyphi/paratyphi_unique_genes.txt $gff | sed -E 's/.*ID=([^;]+);.*/\1/' > cog/paratyphi/paratyphi_unique_ids.txt
seqkit grep -f cog/paratyphi/paratyphi_unique_ids.txt $faa > cog/paratyphi/paratyphi_unique.faa

COGclassifier -i cog/paratyphi/paratyphi_unique.faa -o cog/paratyphi
plot_cog_classifier_barchart -i cog/paratyphi/new_counts_R.tsv -o cog/paratyphi/07_barchart_custom_change_color.html


##Litchfiled-AG21-0055
mkdir cog/Litchfiled
##merge gff files
gff=results/prokka/AG21-0055_withseqs.gff
faa=results/prokka/AG21-0055.faa
grep -f cog/Litchfiled/litchfield_unique_genes.txt $gff | sed -E 's/.*ID=([^;]+);.*/\1/' > cog/Litchfiled/litchfield_unique_ids.txt
seqkit grep -f cog/Litchfiled/litchfield_unique_ids.txt $faa > cog/Litchfiled/litchfield_unique.faa
COGclassifier -i cog/Litchfiled/litchfield_unique.faa -o cog/Litchfiled
plot_cog_classifier_barchart -i cog/Litchfiled/new_counts_R.tsv -o cog/Litchfiled/07_barchart_custom_change_color.html


##There were no genes 

## change the name of the files and run knsp3 again so that phylogeny tree would have same label
for file in results_29_new_seq/all_assemblies/GRS*; do mv "$file" "${file/GRS-S-/AG21-00}"; done
##why this is not working
rename 's/0*([12345678])$\.fasta)/000$1\.fasta/' results_29_new_seq/all_assemblies/AG21-00[1-8].fasta

##rerun KSNP3
mkdir -p results_29_new_seq/ksnp3
# making files with the directory and the file name side by side
paste <(ls $PWD/results_29_new_seq/all_assemblies/*) <(ls results_29_new_seq/all_assemblies | sed 's/.fasta//') > results_29_new_seq/ksnp3/assembly_list.txt
#Run the Ksnp3 from the MCIC scripts
sbatch mcic-scripts/trees/ksnp3.sh -i results_29_new_seq/ksnp3/assembly_list.txt -o results_29_new_seq/ksnp3

##Ran ANI from MCIC
## Run sourmash to get the ANI
indir=results_29_new_seq/all_assemblies
sbatch mcic-scripts/misc/sourmash_compare.sh -i "$indir" -o results_29_new_seq/sourmash/ani --ani

##ANI for the NCBI seqeunces
## To download the genome seqeunce data from NCBI, seqeunce will be in the current directory in zip file
module load python/3.6-conda5.2
source activate /fs/ess/PAS0471/jelmer/conda/ncbi-datasets
datasets download genome accession --inputfile results/2689_assembly_ohio.csv
## Unzip file
mkdir results/ncbi_2989_seq
unzip ncbi_dataset.zip -d results/ncbi_2989_seq
mkdir results/ncbi_2989_seq_fasta
find results/ncbi_2989_seq/ncbi_dataset -name "*genomic.fna" -exec mv -v {} . \; mv -v *genomic.fna results/ncbi_2989_seq_fasta/

##Run ANI for all isolates from Ohio
cp results_29_new_seq/all_assemblies/*.fasta results/ncbi_2989_seq/fasta/
f=results/ncbi_2989_seq/fasta/GCA_000256625.1_ASM25662v1_genomic.fna
for f in results/ncbi_2989_seq/fasta/*.fna; do mv -- "$f" "${f%.fna}.fasta"; done
indir=results/ncbi_2989_seq/fasta_files
sbatch mcic-scripts/misc/sourmash_compare.sh -i "$indir" -o results/sourmash/ani --ani

## Run snippy
mcic-scripts/bact/snippy-multi.sh
find WGS_Salmonella -type f -name "*.fastq.gz" | awk '{print $0}' > new_file.txt
find WGS_Salmonella -type f -name "*.fastq.gz" | awk -F/ '{print $NF"\t"}' > new_file.txt

##Run snippy
sbatch mcic-scripts/bact/snippy-multi.sh -i snippy_data.tsv -r snippy/GCA_000006945.2_ASM694v2_genomic.fna -o results/snippy

##Run snippy
sbatch mcic-scripts/bact/snippy-multi.sh -i snippy_data.tsv -r snippy/GCA_000006945.2_ASM694v2_genomic.fna -o results/snippy
## the gene was located at severeal psoition so I plan on extracting the AA seqeunce and aligning them
grep -r 'ZraR' results_29_new_seq/prokka/faa_files/*faa > zraR.txt
awk -F'[>:]' '/ZraR/{print $3}' zraR.txt > new_file.txt
sed -i 's/Transcriptional regulatory protein ZraR//g' new_file.txt 

for file in results_29_new_seq/prokka/faa_files/*faa; do 
seqkit grep -n -r -f new_file.txt $file >> neew.faa; 
done
sed -i 's/Transcriptional regulatory protein ZraR//g' neew.faa

##align the sequences
mafft new_file.txt > output.aln

rm results_29_new_seq/prokka/faa_files/*extracted.fasta

## Do for the mgrB genes
grep -r 'MgrB' results_29_new_seq/prokka/faa_files/*faa > results/colistin_resiatnce_genes/mgrB.txt
awk -F'[>:]' '/MgrB/{print $3}' results/colistin_resiatnce_genes/mgrB_final.txt > results/colistin_resiatnce_genes/new_file.txt
module load miniconda3
source activate /fs/ess/PAS0471/jelmer/conda/seqkit
for file in results_29_new_seq/prokka/faa_files/*faa; do 
seqkit grep -n -r -f results/colistin_resiatnce_genes/new_file.txt $file >> results/colistin_resiatnce_genes/neew.faa; 
done


## Do for phop gene
grep -r 'Virulence transcriptional regulatory protein PhoP' results_29_new_seq/prokka/faa_files/*faa > results/colistin_resiatnce_genes/phop.txt
results/colistin_resiatnce_genes/phop_aa_final.txt
for file in results_29_new_seq/prokka/faa_files/*faa; do 
seqkit grep -n -r -f results/colistin_resiatnce_genes/phop_aa_final.txt $file >> results/colistin_resiatnce_genes/phop.faa; 
done
sed -i 's/Virulence transcriptional regulatory protein PhoP//g' results/colistin_resiatnce_genes/phop.faa
