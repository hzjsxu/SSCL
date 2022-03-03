#!/bin/bash
# File: complot_anchor_TE.sh
# Author: jsxu
# Contact: <hzaujsxu@163.com>
# Usage: TE analysis within anchor
#----

TE_file_dir="./TE_bed_files"
wholeGenomeBin_dir="./bin_dir"
chrom_size_file="Xv10.chrom.size"
out_dir="./bw_files"

### 1. calculate TE density.
for wholeGenomeBin_file in ${wholeGenomeBin_dir}/WholeGenome*bed;do
	bin=$(basename $wholeGenomeBin_file k.bed)
	bin=${bin##*.}
    genome_bin=$(basename $wholeGenomeBin_file .bed)
	for TE_file in ${TE_file_dir}/xv10*bed; do
    	TE_class=$(basename $TE_file .bed)
        echo $genome_bin, $TE_class
        bedtools intersect -a ${wholeGenomeBin_file} -b ${TE_file} -C -nonamecheck | \
            awk -v OFS='\t' -v bin=${bin}  '{print $1,$2,$3,$4/bin}' | \
			bedtools sort -i - >${out_dir}/${genome_bin}_${TE_class}.count.bdg
		bedGraphToBigWig ${out_dir}/${genome_bin}_${TE_class}.count.bdg $chrom_size_file ${out_dir}/${genome_bin}_${TE_class}.count.bw
    done
done


### 2. plot TE profiles and hetmap using deeptools.
computeMatrix reference-point -R SuperLoop.Anchor.25kb.txt \
 			-S LTR_Copia.bw LTR_Gypsy.bw LTR_ERV.bw LTR_Pao.bw LTR_DIRS.bw LINE_Penelope.bw SINE.bw LINE_L1.bw LINE_L2.bw DNA_hAT-Charlie.bw DNA_TcMar.bw DNA_Maverick.bw DNA_PiggyBac.bw Helitron.bw \
			--outFileName wholeGenome.5k.anchor.TE.matrix.gz \
			--sortRegions keep \
			--binSize 5000 \
			--referencePoint center -a 250000 -b 250000 \
			--missingDataAsZero -p 10

plotHeatmap -m wholeGenome.5k.anchor.TE.matrix.gz -o wholeGenome.5k.anchor.TE.plot.pdf \
			--sortRegions descend \
			--outFileSortedRegions wholeGenome.5k.anchor.TE.plot.SortedRegions.bed \
			--outFileNameMatrix wholeGenome.5k.anchor.TE.plot.Peak.matrix.gz \
			--missingDataColor white \
			--colorMap Blues \
			--whatToShow "plot, heatmap and colorbar" \
			--refPointLabel Anchor