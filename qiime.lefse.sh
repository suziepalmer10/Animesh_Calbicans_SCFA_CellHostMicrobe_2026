#!/usr/bin/env bash

#Code by Jiwoong Kim
metadataColumn=$1 # group
lda=$2 # 2
alpha=$3 # 0.05

time head -n1 taxonomy.fraction.table.extended.tsv | sed 's/\t/\n/g' | awk '(NR > 1)' | perl NGS_scripts/table.search.pl -h - 0 "perl NGS_scripts/table.rearrangeColumns.pl -c metadata.txt $metadataColumn sample-id |" 1 | awk -F'\t' '($1 != "")' > $metadataColumn.sample.txt
time perl NGS_scripts/table.rearrangeColumns.pl -c taxonomy.fraction.table.extended.tsv taxonomy `cut -f2 $metadataColumn.sample.txt | awk '(NR > 1)'` | awk '(NR > 1)' | sed 's/; /|/g' | sed 's/^d__/k__/' | bash -c "cat <(perl NGS_scripts/table.transpose.pl $metadataColumn.sample.txt) -" > $metadataColumn.taxonomy.fraction.table.all_taxonomy.tsv
time awk '(NR > 2)' $metadataColumn.taxonomy.fraction.table.all_taxonomy.tsv | perl NGS_scripts/table.extendLines.pl - `head -n1 $metadataColumn.taxonomy.fraction.table.all_taxonomy.tsv | awk -F'\t' '{print "1.."NF - 1}'` | perl NGS_scripts/table.mergeLines.pl -f max - 0 | awk -F'\t' '($2 > 0)' | cut -f1 | bash -c "cat <(head -n2 $metadataColumn.taxonomy.fraction.table.all_taxonomy.tsv | cut -f1) -" | perl NGS_scripts/table.search.pl - 0 $metadataColumn.taxonomy.fraction.table.all_taxonomy.tsv 0 > $metadataColumn.taxonomy.fraction.table.tsv

. "/opt/anaconda3/etc/profile.d/conda.sh"
conda activate lefse

time lefse_format_input.py $metadataColumn.taxonomy.fraction.table.tsv $metadataColumn.lefse_format_input -f r -c 1 -s -1 -u 2 -o 1000000.0
time lefse_run.py $metadataColumn.lefse_format_input $metadataColumn.lefse.txt -l $lda -a $alpha -w $alpha -e 0 -y 1 -f 0.9
awk -F'\t' '($3 != "")' $metadataColumn.lefse.txt > $metadataColumn.lefse.differential_feature.txt

[ -e $metadataColumn.lefse.png ] && rm $metadataColumn.lefse.png
if [ $(cut -f3 $metadataColumn.lefse.differential_feature.txt | sort -u | wc -l) -eq 2 ]
then
	time lefse_plot_res.py $metadataColumn.lefse.txt $metadataColumn.lefse.png --title '' --subclades 1 --max_feature_len 60 --title_font_size 14 --feature_font_size 10 --class_legend_font_size 11 --width 10.0 --left_space 0.1 --right_space 0.1 --background_color w --format png --dpi 300 --all_feats `head -n1 $metadataColumn.taxonomy.fraction.table.tsv | sed 's/\t/\n/g' | awk '(NR > 1)' | sort -u | tr '\n' ':' | sed 's/:$/\n/'`
else
	time lefse_plot_res.py $metadataColumn.lefse.txt $metadataColumn.lefse.png --title '' --subclades 1 --max_feature_len 60 --title_font_size 14 --feature_font_size 10 --class_legend_font_size 11 --width 10.0 --left_space 0.5 --right_space 0.1 --background_color w --format png --dpi 300 --all_feats `head -n1 $metadataColumn.taxonomy.fraction.table.tsv | sed 's/\t/\n/g' | awk '(NR > 1)' | sort -u | tr '\n' ':' | sed 's/:$/\n/'`
fi

[ -e $metadataColumn.lefse.cladogram.png ] && rm $metadataColumn.lefse.cladogram.png
time lefse_plot_cladogram.py $metadataColumn.lefse.txt $metadataColumn.lefse.cladogram.png --sub_clade '' --expand_void_lev 0 --max_lev 6 --title '' --title_font_size 14 --label_font_size 8 --class_legend_font_size 10 --labeled_start_lev 2 --labeled_stop_lev 6 --abrv_start_lev 4 --abrv_stop_lev 6 --radial_start_lev 1 --max_point_size 7.0 --min_point_size 1.5 --point_edge_width 0.25 --siblings_connector_width 2.0 --parents_connector_width 0.8 --alpha 0.2 --left_space_prop 0.01 --right_space_prop 0.25 --colored_label 0 --background_color w --clade_sep 1.5 --format png --dpi 300 --all_feats `head -n1 $metadataColumn.taxonomy.fraction.table.tsv | sed 's/\t/\n/g' | awk '(NR > 1)' | sort -u | tr '\n' ':' | sed 's/:$/\n/'`

conda deactivate

#. "/archive/PCDC/shared/jkim23/miniconda3/etc/profile.d/conda.sh"
#conda activate python2-lefse
#
#if [ $(awk -F'\t' '($3 != "")' $metadataColumn.lefse.txt | cut -f3 | sort -u | wc -l) = 1 ]
#then
#	time python /home2/jkim23/src/lefse/home/ubuntu/lefse_to_export/plot_res.modified.py $metadataColumn.lefse.txt $metadataColumn.lefse.png --title '' --subclades 1 --max_feature_len 60 --title_font_size 14 --feature_font_size 10 --class_legend_font_size 11 --width 10.0 --left_space 0.1 --right_space 0.1 --background_color w --format png --dpi 300 --all_feats `head -n1 $metadataColumn.taxonomy.fraction.table.tsv | sed 's/\t/\n/g' | awk '(NR > 1)' | sort -u | tr '\n' ':' | sed 's/:$/\n/'`
#else
#	time python /home2/jkim23/src/lefse/home/ubuntu/lefse_to_export/plot_res.modified.py $metadataColumn.lefse.txt $metadataColumn.lefse.png --title '' --subclades 1 --max_feature_len 60 --title_font_size 14 --feature_font_size 10 --class_legend_font_size 11 --width 10.0 --left_space 0.5 --right_space 0.1 --background_color w --format png --dpi 300 --all_feats `head -n1 $metadataColumn.taxonomy.fraction.table.tsv | sed 's/\t/\n/g' | awk '(NR > 1)' | sort -u | tr '\n' ':' | sed 's/:$/\n/'`
#fi
#time python /home2/jkim23/src/lefse/home/ubuntu/lefse_to_export/plot_cladogram.modified.py $metadataColumn.lefse.txt $metadataColumn.lefse.cladogram.png --sub_clade '' --expand_void_lev 0 --max_lev 6 --title '' --title_font_size 14 --label_font_size 8 --class_legend_font_size 10 --labeled_start_lev 2 --labeled_stop_lev 6 --abrv_start_lev 4 --abrv_stop_lev 6 --radial_start_lev 1 --max_point_size 7.0 --min_point_size 1.5 --point_edge_width 0.25 --siblings_connector_width 2.0 --parents_connector_width 0.8 --alpha 0.2 --left_space_prop 0.01 --right_space_prop 0.25 --colored_label 0 --background_color w --clade_sep 1.5 --format png --dpi 300 --all_feats `head -n1 $metadataColumn.taxonomy.fraction.table.tsv | sed 's/\t/\n/g' | awk '(NR > 1)' | sort -u | tr '\n' ':' | sed 's/:$/\n/'`
#
#conda deactivate
