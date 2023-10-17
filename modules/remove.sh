reads=$1
echo $reads
abc="" && [ $(echo "${reads}" | grep ".gz") ] && abc=" --readFilesCommand zcat "
echo "STAR --genomeDir ${index} $abc"
echo ${abc}