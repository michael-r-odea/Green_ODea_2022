# run FastQC to explore quality of FASTQ files
find $(pwd) -type f | xargs fastqc -o fastqc_results

# trim low quality reads and adaptor sequences 
find ./ -type f -name "*_1.fq.gz" | awk -F "_1.fq.gz" '{print $1}' | sort -u \
 | parallel trim_galore \
 --paired --nextera --output_dir trimming {}_1.fq.gz {}_2.fq.gz

# run FastQC on trimmed fastq files
fastqc *.fq.gz -o trimmed_fastqc_results -t 68

# run Salmon to map reads 
ls *_1_val_1.fq.gz | sort > R1_trimmed_list
ls *_2_val_2.fq.gz | sort > R2_trimmed_list

paste R1_trimmed_list R2_trimmed_list | while read R1 R2
do
  echo -e "\n"
  echo -e "--Input file(s) is/are:\t""${R1}"", ""${R2}"
  outdir=$(echo "${R2}" | cut -f1 -d"_")
  echo -e "--Output directory is:\toutput/""${outdir}"
  mkdir -p "../../quantification/output/""${outdir}"

  salmon quant -i salmon_indices/alias/mm10/salmon_sa_index/default -l A \
	-1 <(gunzip -c "${R1}") -2 <(gunzip -c "${R2}") \
	--validateMappings --seqBias --gcBias --output="../quantification/output/""${outdir}"
done
