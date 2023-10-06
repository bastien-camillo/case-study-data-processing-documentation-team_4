```bash
conda activate day1
```

```bash
mkdir -p ~/course/wdir/day1
mkdir -p ~/course/wdir/day1/data

cd ~/course/wdir/day1

#if data is not empty remove everything in it
if [[ ! -z "$(ls -A ~/course/wdir/day1/data)" ]]
then
    rm -r ~/course/wdir/day1/data/*    
fi

cp -r ~/course/data/day1/* ~/course/wdir/day1/data
ls ~/course/wdir/day1/data/* .
du -sk ~/course/wdir/day1/data/*
```


```bash
# EXECRCIE 4
for file in ~/course/wdir/day1/data/*
do
    echo "File: $file"
    # get a variable with the name of the file without the extension
    # file_name=$(basename $file .fastq.gz)
    if [[ $file == *.fastq.gz ]]
    then
        file_name=$(basename $file .fastq.gz)
        echo $file_name

        zcat $file | wc -l
        zcat $file | awk '{s++}END{print s/4}'
        fastp -i $file -o $file_name.trimmed.fastq -l 30
    fi
done

for file in ~/course/wdir/day1/data/*
do
    if [[ $file == *.trimmed.fastq ]]
    then
        cat $file|awk '{s++}END{print s/4}'
    fi
done

# OR

# cat fastp.json|grep "total_reads" && cat fastp.json|grep "read1_mean_length" && cat fastp.json|grep "read1_adapter_counts"

# for file in ~/course/wdir/day1/data/*.gz
# do
#     echo ${file}
#     zcat ${file} | awk '{s++}END{print s/4}'
# done

# EXECRCIE 5

bwa index ~/course/wdir/day1/data/paeruginosa.fasta.gz

for file in ~/course/wdir/day1/data/*.trimmed.fastq
do
    file_name=$(basename $file .trimmed.fastq)
    bwa aln -t 5 ~/course/wdir/day1/data/paeruginosa.fasta.gz $file | bwa samse ~/course/wdir/day1/data/paeruginosa.fasta.gz - $file | samtools view -b - | samtools sort -o ~/course/wdir/day1/data/$file_name.trimmed.sorted.bam
done

# EXECRCIE 6

for file in ~/course/wdir/day1/data/*.trimmed.fastq
do
    file_name=$(basename $file .trimmed.fastq)
    samtools flagstat ~/course/wdir/day1/data/$file_name.trimmed.sorted.bam
    samtools view -b -F4 -q 30 ~/course/wdir/day1/data/$file_name.trimmed.sorted.bam | samtools rmdup -s - ~/course/wdir/day1/data/$file_name.filtered.rmdup.bam
    samtools view ~/course/wdir/day1/data/$file_name.filtered.rmdup.bam | awk '{print length($10)}' | datamash mean 1
done

# EXECRCIE 7

for file in ~/course/wdir/day1/data/*.trimmed.fastq
do
    file_name=$(basename $file .trimmed.fastq)
    samtools index ~/course/wdir/day1/data/$file_name.filtered.rmdup.bam
    mapDamage -i ~/course/wdir/day1/data/$file_name.filtered.rmdup.bam -r ~/course/wdir/day1/data/paeruginosa.fasta.gz --no-stats
done
```
