```bash
cd ~/course/wdir/
mkdir -p ~/course/wdir/day2
mkdir -p ~/course/wdir/day2/mapping
cd ~/course/wdir/day2/mapping
```

```bash
conda activate day2
```

```bash
for file in ~/course/data/day2/fastq/*
do
    if [[ $file == *.fq.gz ]]
    then
        ln -s $file .
    fi
done

for file in ~/course/wdir/day2/mapping/*
do
    if [[ $file == *.fq.gz ]]
    then
        file_name=$(basename $file .fq.gz)
        vsearch --fastx_uniques $file --fastqout $file_name.vs.fq --minseqlength 30 --strand both
        gzip $file_name.vs.fq
        zcat $file_name.vs.fq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $file_name.read_length.txt
    fi
done
```

```bash
conda activate day1
```

```bash
for file in ~/course/wdir/day2/mapping/*
do
    if [[ $file == *.vs.fq.gz ]]
    then
        file_name=$(basename $file .vs.fq.gz)
        bowtie2 --threads 5 -k 100 -x ~/course/data/shared/mapping/db/aegenomics.db -U $file --no-unal | samtools view -bS - > $file_name.bam
        samtools sort -n $file_name.bam -@ 5 > $file_name.sort.bam
    fi
done
```

```bash
conda activate metaDMG
```

```bash
metaDMG config *.sort.bam --names ~/course/data/shared/mapping/taxonomy/names.dmp --nodes ~/course/data/shared/mapping/taxonomy/nodes.dmp --acc2tax ~/course/data/shared/mapping/taxonomy/acc2taxid.map.gz -m /usr/local/bin/metaDMG-cpp
```

#Edit the config file according to wishes, study and data. set custom database = true, set number of threads to 4.
```bash
vim config.yaml
```

```bash
metaDMG compute config.yaml 
metaDMG convert --output metaDMGresults.csv --add-fit-predictions
metaDMG dashboard config.yaml
```




