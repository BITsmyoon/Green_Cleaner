#!/bin/bash -e

path_data=$1

path_output=${path_data}/../output/qiime

mkdir ${path_output}

# raw asv count
sed 's/asv_id/Feature ID/g' ${path_data}/input_taxa_table.txt > ${path_output}/taxa_table_tmp.txt

qiime tools import \
        --type 'FeatureData[Taxonomy]' \
        --input-path ${path_output}/taxa_table_tmp.txt \
        --output-path ${path_output}/taxa_table.qza

sed '1s/^/# Constructed from biom file\n/' ${path_data}/input_raw_asv_count.txt > ${path_output}/add_header_asv_count.txt

biom convert -i ${path_output}/add_header_asv_count.txt -o ${path_output}/asv_count.biom --table-type="OTU table" --to-hdf5

qiime tools import --input-path ${path_output}/asv_count.biom --type 'FeatureTable[Frequency]' --output-path ${path_output}/asv_count.qza

qiime taxa barplot --i-table ${path_output}/asv_count.qza --i-taxonomy ${path_output}/taxa_table.qza --m-metadata-file ${path_data}/input_qiime_meta.txt --o-visualization ${path_output}/qiime_barplot.qzv

qiime tools export --input-path ${path_output}/qiime_barplot.qzv --output-path ${path_output}/qiime_barplot

# scrub asv count
sed '1s/^/# Constructed from biom file\n/' ${path_data}/input_scrub_asv_count.txt > ${path_output}/add_header_scrub_asv_count.txt

biom convert -i ${path_output}/add_header_scrub_asv_count.txt -o ${path_output}/scrub_asv_count.biom --table-type="OTU table" --to-hdf5

qiime tools import --input-path ${path_output}/scrub_asv_count.biom --type 'FeatureTable[Frequency]' --output-path ${path_output}/scrub_asv_count.qza

qiime taxa barplot --i-table ${path_output}/scrub_asv_count.qza --i-taxonomy ${path_output}/taxa_table.qza --m-metadata-file ${path_data}/input_scrub_meta_data.txt --o-visualization ${path_output}/scrub_qiime_barplot.qzv

qiime tools export --input-path ${path_output}/scrub_qiime_barplot.qzv --output-path ${path_output}/scrub_qiime_barplot

