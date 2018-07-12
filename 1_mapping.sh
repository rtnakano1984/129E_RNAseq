


raw="/your_raw_data_dir"
qf="/your_quality_filtered_dir"
out="/your_mappping_output_dir"

gind="/your_index_dir/TAIR10_g/TAIR10"
tind="/your_index_dir/TAIR10_t/tair10"
gff="/your_indxx_dir/TAIR10_GFF3_exons_merged.gff"

# quality trimming
for lib in {A..Z} A{A..V}; do
	bsub -q normal -R "rusage[mem=8000]" -n 1 fastq_quality_filter -q 25 -p75 -Q 33 -i ${raw}/1106_${lib}.fastq -o ${qf}/1106_A.fastq > ${qf}/Lib${lib}_score
done

# read mapping
for lib in {A..Z} A{A..V}; do
	bsub -q multicore20 -R "rusage[mem=8000]" -n 20 tophat2 -p 20 -a 10 -g 10 --library-type fr-unstranded  --transcriptome-index=${tind} -o ${out}/Lib{lib} ${gind} ${qf}/1106_${lib}.fastq
done

# read counting
for lib in {A..Z} A{A..V}; do
	cd ${out}/Lib{lib}
	coverageBed -split -abam accepted_hits.bam -b ${gff} > accepted_hits.cov
done
