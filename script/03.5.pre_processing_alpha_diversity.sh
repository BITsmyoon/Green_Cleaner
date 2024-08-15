#!/bin/bash -e

path_data=$1

path_output=${path_data}/../output/qiime/alpha_diversity
path_green_cleaner=${path_data}/../output/dataset_output/asv_count
path_dataset_qiime_output=${path_data}/../output/qiime/dataset

mkdir ${path_output}
mkdir ${path_dataset_qiime_output}

# raw alpha
qiime diversity alpha \
        --i-table ${path_output}/../asv_count.qza \
        --p-metric 'chao1_ci' \
        --o-alpha-diversity ${path_output}/raw_chao1_ci

qiime tools export \
        --input-path ${path_output}/raw_chao1_ci.qza \
        --output-path ${path_output}/raw_chao1_ci

# scrub alpha
sed 's/asv_id/#OTU ID/g' ${path_data}/input_scrub_asv_count.txt > ${path_output}/scrub_tmp_asv_count.txt
sed -i '1s/^/# Constructed from biom file\n/' ${path_output}/scrub_tmp_asv_count.txt

biom convert -i  ${path_output}/scrub_tmp_asv_count.txt -o ${path_output}/scrub_asv_count.biom --table-type="OTU table" --to-hdf5
qiime tools import --input-path ${path_output}/scrub_asv_count.biom --type 'FeatureTable[Frequency]' --output-path ${path_output}/scrub_asv_count.qza

qiime diversity alpha \
        --i-table ${path_output}/scrub_asv_count.qza \
        --p-metric 'chao1_ci' \
        --o-alpha-diversity ${path_output}/scrub_chao1_ci

qiime tools export \
        --input-path ${path_output}/scrub_chao1_ci.qza \
        --output-path ${path_output}/scrub_chao1_ci

# green cleaner
for value_num in {1..10}
do

	value_tmp_dataset="dataset_${value_num}"

	sed 's/asv_id/#OTU ID/g' ${path_green_cleaner}/${value_tmp_dataset}_merged_asv_count.txt > ${path_output}/${value_tmp_dataset}_gc_tmp_asv_count.txt
	sed -i '1s/^/# Constructed from biom file\n/' ${path_output}/${value_tmp_dataset}_gc_tmp_asv_count.txt
	
	biom convert -i  ${path_output}/${value_tmp_dataset}_gc_tmp_asv_count.txt -o  ${path_output}/${value_tmp_dataset}_gc_asv_count.biom --table-type="OTU table" --to-hdf5
	qiime tools import --input-path ${path_output}/${value_tmp_dataset}_gc_asv_count.biom --type 'FeatureTable[Frequency]' --output-path ${path_output}/${value_tmp_dataset}_gc_asv_count.qza

	qiime diversity alpha \
		--i-table ${path_output}/${value_tmp_dataset}_gc_asv_count.qza \
		--p-metric 'chao1_ci' \
		--o-alpha-diversity ${path_output}/${value_tmp_dataset}_gc_chao1_ci

	qiime tools export \
		--input-path ${path_output}/${value_tmp_dataset}_gc_chao1_ci.qza \
		--output-path ${path_output}/${value_tmp_dataset}_gc_chao1_ci

	qiime taxa barplot --i-table ${path_output}/${value_tmp_dataset}_gc_asv_count.qza --i-taxonomy ${path_dataset_qiime_output}/../taxa_table.qza --m-metadata-file ${path_green_cleaner}/../meta_data_for_qiime/${value_tmp_dataset}_metadata.txt --o-visualization ${path_dataset_qiime_output}/${value_tmp_dataset}_qiime_barplot.qzv
	
	qiime tools export --input-path ${path_dataset_qiime_output}/${value_tmp_dataset}_qiime_barplot.qzv --output-path ${path_dataset_qiime_output}/${value_tmp_dataset}_qiime_barplot


done
