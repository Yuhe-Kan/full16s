##：reference：https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2%E2%80%902024.5)

# Start Terminal
wsl -d Ubuntu-24.04
#____________________________________________________________________

cd /mnt/e/
# set path
path_root=/mnt/e/b_myproject/kan_2024_full16s
path_data=$path_root/03_data/analysis-individual/kan_2024
path_script=$path_root/04_scripts/analysis-individual/kan_2024/
#____________________________________________________________________

# Check if the fq file was transferred completely.
cd $path_data/00_raw_fastq
md5sum -c data_md5.txt
#____________________________________________________________________

# Test read quality
# conda create -n qc
conda activate qc
# sudo apt update
# sudo apt install fastqc
# sudo apt install multiqc
fastqc $path_data/00_raw_fastq/*.fq.gz -o $path_data/02_fastqc_out
multiqc $path_data/02_fastqc_out --filename $path_data/02_fastqc_out/multiqc_report.html
#____________________________________________________________________

#import read
conda activate qiime2-amplicon-2024.5

awk 'NR==1{print "sample-id\tabsolute-filepath"} \
	   NR>1{print $1"\t/mnt/e/b_myproject/kan_2024_full16s/03_data/analysis-individual/kan_2024/00_raw_fastq/"$1"_clean.fq.gz"}' \
	   $path_data/01_Metadata/metadata.txt > $path_data/03_qiime2/manifest_kan
head -n10 $path_data/03_qiime2/manifest_kan
#____________________________________________________________________

qiime tools import \
    --type SampleData[SequencesWithQuality] \
    --input-path $path_data/03_qiime2/manifest_kan \
    --output-path $path_data/03_qiime2/demux_seqs_kan.qza \
    --input-format SingleEndFastqManifestPhred33V2

qiime demux summarize \
--i-data $path_data/03_qiime2/demux_seqs_kan.qza \
--o-visualization $path_data/03_qiime2/demux_seqs_kan.qzv
#____________________________________________________________________

# filter processing
qiime dada2 denoise-ccs \
   --i-demultiplexed-seqs $path_data/03_qiime2/demux_seqs_kan.qza \
   --p-min-len 1000 --p-max-len 1600 \
   --p-front AGRGTTYGATYMTGGCTCAG --p-adapter RGYTACCTTGTTACGACTT \
   --p-n-threads 10 \
   --o-table $path_data/03_qiime2/otutab.qza \
   --o-representative-sequences $path_data/03_qiime2/repseq.qza \
   --o-denoising-stats $path_data/03_qiime2/dada2_stats.qza \
   --verbose

#View noise reduction results
  qiime metadata tabulate \
    --m-input-file $path_data/03_qiime2/dada2_stats.qza \
    --o-visualization $path_data/03_qiime2/dada2_stats.qzv

# Statistical analysis for otutab
	qiime feature-table summarize \
	 --i-table $path_data/03_qiime2/otutab.qza \
	 --m-sample-metadata-file $path_data/01_Metadata/metadata.txt \
	 --o-visualization $path_data/03_qiime2/otutab.qzv
	 
#Assess subset of taxonomic assignments with BLAST
  qiime feature-table tabulate-seqs \
  --i-data $path_data/03_qiime2/repseq.qza \
  --o-visualization $path_data/03_qiime2/repseq.qzv
#____________________________________________________________________

#Assign taxonomy to ASVs
cluster=$path_root/03_data/analysis-individual/CLUSTER/taxonomy/classifier_GSR

#download database: GSR full-length 16S database: GSR-DB_full-16S.tar.gz，from：https://manichanh.vhir.org/gsrdb/download_db_links2.php

#The error message "conda install killed" indicates that channels are only set to default.
conda activate qiime2-amplicon-2023.5

qiime feature-classifier classify-sklearn \
--i-reads $path_data/03_qiime2/repseq.qza \
--i-classifier $cluster/classifier_GSR-DB_full-16S.qza \
--o-classification $path_data/03_qiime2/taxnomy.qza

qiime metadata tabulate \
  --m-input-file $path_data/03_qiime2/taxnomy.qza \
  --o-visualization $path_data/03_qiime2/taxonomy.qzv
#____________________________________________________________________


#Build tree with SEPP QIIME 2 plugin (16S data)
# method 1
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $path_data/03_qiime2/repseq.qza \
  --o-alignment $path_data/03_qiime2/tree_fasttree/aligned-rep-seqs.qza \
  --o-masked-alignment $path_data/03_qiime2/tree_fasttree/masked-aligned-rep-seqs.qza \
  --o-tree $path_data/03_qiime2/tree_fasttree/unrooted-tree.qza \
  --o-rooted-tree $path_data/03_qiime2/tree_fasttree/rooted-tree.qza

# method 2
# wget https://data.qiime2.org/2019.10/common/sepp-refs-gg-13-8.qza
qiime fragment-insertion sepp \
   --i-representative-sequences $path_data/03_qiime2/repseq.qza \
   --i-reference-database $path_root/03_data/analysis-individual/CLUSTER/phylotree/input/sepp-refs-gg-13-8.qza \
   --o-tree $path_data/03_qiime2/tree.qza \
   --o-placements $path_data/03_qiime2/insertion-placements.qza \
   --p-threads 10


qiime tools export \
  --input-path $path_data/03_qiime2/tree.qza \
  --output-path $path_data/03_qiime2/exported-tree
#____________________________________________________________________
