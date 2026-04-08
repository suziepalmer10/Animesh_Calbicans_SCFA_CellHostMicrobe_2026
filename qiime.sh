#!/usr/bin/env bash

threads=16
classifier=silva-138-99-515-806-nb-classifier.qza

[ -f metadata.txt ] && perl NGS_scripts/fixNewline.pl metadata.txt
[ -f metadataColumns.txt ] && perl NGS_scripts/fixNewline.pl metadataColumns.txt
perl prepare_input.pl

# Space-separated metadata columns
[ "$metadataColumns" == "" ] && [ -f metadataColumns.txt ] && metadataColumns=$(cat metadataColumns.txt)


. "/opt/anaconda3/etc/profile.d/conda.sh"
conda activate qiime2-2022.2

# Import FASTQ files using manifest file
time qiime tools import \
	--type 'SampleData[PairedEndSequencesWithQuality]' \
	--input-path manifest.txt \
	--output-path paired-end-demux.qza \
	--input-format PairedEndFastqManifestPhred33V2

# Check primer trimming
time qiime cutadapt trim-paired \
	--i-demultiplexed-sequences paired-end-demux.qza \
	--p-cores $threads \
	--p-front-f GTGCCAGCMGCCGCGGTAA \
	--p-front-r GGACTACHVGGGTWTCTAAT \
	--o-trimmed-sequences paired-end-demux-trimmed.qza \
	--verbose


# Denoise and dereplicate paired-end sequences
time qiime dada2 denoise-paired \
	--i-demultiplexed-seqs paired-end-demux.qza \
	--p-trunc-len-f 0 \
	--p-trunc-len-r 0 \
	--p-trim-left-f 0 \
	--p-trim-left-r 0 \
	--p-max-ee-f 2 \
	--p-max-ee-r 2 \
	--p-n-threads $threads \
	--o-table table.qza \
	--o-representative-sequences rep-seqs.qza \
	--o-denoising-stats stats.qza

time qiime tools export --input-path table.qza --output-path .
time qiime tools export --input-path rep-seqs.qza --output-path .
time qiime tools export --input-path stats.qza --output-path .

time biom convert --input-fp feature-table.biom --output-fp table.tsv --to-tsv
minimumSampleDepth=$(awk '(NR > 2)' table.tsv | perl NGS_scripts/table.extendLines.pl - `head -n2 table.tsv | tail -n1 | awk -F'\t' '{print "1.."NF - 1}'` `head -n2 table.tsv | tail -n1 | cut -f2-` | cut -f2,3 | perl NGS_scripts/table.mergeLines.pl -f sum - 0 | cut -f2 | sort -n | head -n1)

time qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
time qiime metadata tabulate --m-input-file stats.qza --o-visualization stats.qzv


# Pre-fitted sklearn-based taxonomy classifier
time qiime feature-classifier classify-sklearn \
	--i-reads rep-seqs.qza \
	--i-classifier $classifier \
	--p-n-jobs $threads \
	--o-classification taxonomy.qza

time qiime tools export --input-path taxonomy.qza --output-path .

time qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv


# Build a phylogenetic tree using fasttree and mafft alignment
time qiime phylogeny align-to-tree-mafft-fasttree \
	--i-sequences rep-seqs.qza \
	--o-alignment aligned-rep-seqs.qza \
	--o-masked-alignment masked-aligned-rep-seqs.qza \
	--o-tree unrooted-tree.qza \
	--o-rooted-tree rooted-tree.qza

time qiime tools export --input-path unrooted-tree.qza --output-path unrooted-tree
time qiime tools export --input-path rooted-tree.qza --output-path rooted-tree

time perl NGS_scripts/table.rearrangeColumns.pl -c taxonomy.tsv 'Feature ID' 'Taxon' | awk '(NR > 1)' | sed 's/; .__unidentified.*$//' | sed 's/; .__uncultured.*$//' | sed 's/; .__metagenome.*$//' | sed 's/; .__gut_metagenome.*$//' | sed 's/\t.*; /\t/' | awk -vOFS='\t' '{print $o, "black"}' > feature.name.color.txt

time sed "s/'//g" unrooted-tree/tree.nwk | Rscript NGS_scripts/tree.R stdin feature.name.color.txt unrooted-tree/tree.pdf 60 80
time sed "s/'//g" unrooted-tree/tree.nwk | Rscript NGS_scripts/tree_cladogram.R stdin feature.name.color.txt unrooted-tree/tree_cladogram.pdf 60 60

time Rscript NGS_scripts/tree.R rooted-tree/tree.nwk feature.name.color.txt rooted-tree/tree.pdf 60 80
time Rscript NGS_scripts/tree_cladogram.R rooted-tree/tree.nwk feature.name.color.txt rooted-tree/tree_cladogram.pdf 60 60


# Following commands require metadata file
if [ -f metadata.txt ]
then

# Summarize dada2 result table
time qiime feature-table summarize \
	--i-table table.qza \
	--m-sample-metadata-file metadata.txt \
	--o-visualization table.qzv

# Visualize taxonomy with an interactive bar plot
time qiime taxa barplot \
	--i-table table.qza \
	--i-taxonomy taxonomy.qza \
	--m-metadata-file metadata.txt \
	--o-visualization barplot.qzv

# Core diversity metrics (phylogenetic and non-phylogenetic)
time qiime diversity core-metrics-phylogenetic \
	--i-table table.qza \
	--i-phylogeny rooted-tree.qza \
	--p-sampling-depth $minimumSampleDepth \
	--m-metadata-file metadata.txt \
	--p-n-jobs-or-threads $threads \
	--output-dir core_metrics_results

# Alpha diversity comparisons
time for prefix in $(ls core_metrics_results/*_vector.qza | sed 's/_vector\.qza$//'); do
	qiime diversity alpha-group-significance \
		--i-alpha-diversity ${prefix}_vector.qza \
		--m-metadata-file metadata.txt \
		--o-visualization ${prefix}_significance.qzv
done

# Beta diversity group significance
time for prefix in $(ls core_metrics_results/*_distance_matrix.qza | sed 's/_distance_matrix\.qza$//'); do for metadataColumn in $metadataColumns; do
	qiime diversity beta-group-significance \
		--i-distance-matrix ${prefix}_distance_matrix.qza \
		--m-metadata-file metadata.txt \
		--m-metadata-column $metadataColumn \
		--o-visualization ${prefix}_${metadataColumn}_significance.qzv \
		--p-pairwise
done; done

time qiime tools export --input-path core_metrics_results/weighted_unifrac_distance_matrix.qza --output-path .

time for metadataColumn in $metadataColumns; do
	if [ $(perl NGS_scripts/table.rearrangeColumns.pl -c metadata.txt $metadataColumn | awk '(NR > 1)' | sort -u | wc -l) -le $(cat NGS_scripts/color.txt | wc -l) ]
	then
		bash -c "paste <(perl NGS_scripts/table.rearrangeColumns.pl -c metadata.txt $metadataColumn | awk '(NR > 1)' | sort -u) NGS_scripts/color.txt" | awk -F'\t' '($1 != "")' | perl NGS_scripts/table.addColumns.pl "perl NGS_scripts/table.rearrangeColumns.pl -c metadata.txt sample-id $metadataColumn | awk '(NR > 1)' |" 1 - 0 1 | cut -f1,3 > sample.${metadataColumn}_color.txt
		Rscript NGS_scripts/hclust_fan.color.R distance-matrix.tsv sample.${metadataColumn}_color.txt sample.${metadataColumn}_color.hclust_fan.pdf 30 30
		awk -F'\t' -vOFS='\t' '{print $1, $1, $2}' sample.${metadataColumn}_color.txt | Rscript NGS_scripts/distance_matrix.mst_plot.R distance-matrix.tsv stdin sample.${metadataColumn}_color.mst_plot.pdf 30 30 0 15000
	fi
done

fi

conda deactivate


time awk '(NR > 2)' table.tsv | perl NGS_scripts/table.extendLines.pl - `awk '(NR > 1)' table.tsv | head -n1 | awk -F'\t' '{print "1.."NF - 1}'` `awk '(NR > 1)' table.tsv | head -n1 | cut -f2-` | perl NGS_scripts/table.addColumns.pl -m - 0 taxonomy.tsv 0 1 | awk -F'\t' -vOFS='\t' '($3 > 0) {print $4, $2, $3}' | perl NGS_scripts/table.mergeLines.pl -f sum - 0,1 > taxonomy.sample.abundance.txt
time cut -f2,3 taxonomy.sample.abundance.txt | perl NGS_scripts/table.mergeLines.pl -f sum - 0 | perl NGS_scripts/table.addColumns.pl taxonomy.sample.abundance.txt 1 - 0 1 | awk -F'\t' -vOFS='\t' '{print $1, $2, $3 / $4}' > taxonomy.sample.fraction.txt
time perl NGS_scripts/table.mergeLines.pl -f sum taxonomy.sample.fraction.txt 0 1 2 `awk '(NR > 1)' table.tsv | head -n1 | cut -f2-` | perl NGS_scripts/table.substitute_value.pl - '' 0 | bash -c "cat <(awk '(NR > 1)' table.tsv | head -n1 | cut -f2- | sed 's/^/taxonomy\t/') -" > taxonomy.fraction.table.tsv
time awk '(NR > 1)' taxonomy.fraction.table.tsv | perl extend_taxonomy.pl - | perl NGS_scripts/table.mergeLines.pl -f sum - 0 | bash -c "cat <(head -n1 taxonomy.fraction.table.tsv) -" > taxonomy.fraction.table.extended.tsv


#time for metadataColumn in $metadataColumns; do bash qiime.lefse.sh $metadataColumn 2 0.05; done


mkdir result
for file in stats.tsv barplot.qzv taxonomy.fraction.table.extended.tsv; do (cd result; ln -sf ../$file $file); done
(cd result; ln -sf ../core_metrics_results/weighted_unifrac_emperor.qzv weighted_unifrac_emperor.qzv)
(cd result; ln -sf ../rooted-tree/tree.pdf rooted-tree.pdf)
(cd result; ln -sf ../unrooted-tree/tree.pdf unrooted-tree.pdf)
for metadataColumn in $metadataColumns; do for file in sample.${metadataColumn}_color.hclust_fan.pdf sample.${metadataColumn}_color.mst_plot.pdf; do (cd result; ln -sf ../$file $file); done; done


rm -rf fastq
rm paired-end-demux.qza paired-end-demux-trimmed.qza
