# session 1 | setting the stage #
```bash
cd YOUR_DIRECTORY

mkdir -p metaDMG
cd metaDMG

ln -s ~/YOUR_DIRECTORY data
```

# session 2 | extension, dereplication and mapping #
```bash
conda activate day4
```

```bash
mamba install bbmap
pip install tabview
```

```bash
for file in ~/YOUR_DIRECTORY/data/*
do
    if [[ $file == *.fastq.gz ]]
    then
        file_name=$(basename $file .fastq.gz)
        seqkit stats -j 5 -T $file | csvtk -t pretty
    fi
done
```

## read extension of all fastq files from YOUR_DIRECTORY ##
```bash
for file in ~/YOUR_DIRECTORY/data/*
do
    if [[ $file == *.fastq.gz ]]
    then
        file_name=$(basename $file .fastq.gz)
        MEM=$(loglog.sh seed=1234 k=17 in=$file ignorebadquality 2> >(grep Cardinality) \
    | awk -vP=0.5 -vB=16 -vH=3 '{{print int( (((B*H)/8)*$2)/P )}}')

        tadpole.sh -Xmx10G \
            k=17 \
            in=$file \
            out=$file_name.extended.fastq.gz \
            mode=extend \
            ibb=f \
            prefilter=0 \
            el=100 er=100 \
            threads=5 \
            overwrite=true \
            trimends=9 \
            ecc=f ecco=f \
            filtermem="${MEM}" \
            conservative=t \
            ignorebadquality
    fi
done
```

## dereplication ##
```bash
for file in ~/YOUR_DIRECTORY/data/*
do
    if [[ $file == *.extended.fastq.gz ]]
    then
        file_name=$(basename $file .extended.fastq.gz)
        seqkit rmdup -j 5 -s -o $file_name.extended.derep.fastq.gz $file_name.extended.fastq.gz
    fi
done
```

## get the original reads ##
```bash
for file in ~/YOUR_DIRECTORY/data/*
do
    if [[ $file == *.extended.fastq.gz ]]
    then
        file_name=$(basename $file .extended.fastq.gz)
        filterbyname.sh in=$file_name.fastq.gz out=$file_name.mapping.fastq.gz names=$file_name.extended.derep.fastq.gz threads=5 overwrite=t include=t
    fi
done
```

## mapping ##
```bash
conda activate mapping
```

```bash
mkdir -d data/mapping_output
```

```bash
for file in ~/YOUR_DIRECTORY/data/*
do
    if [[ $file == *.mapping.fastq.gz ]]
    then
        file_name=$(basename $file .mapping.fastq.gz)
        bowtie2 -p 5 -k 100 -D 10 -R 2 \
            -N 0 -D 5 -R 1 -L 22 -i S,0,2.50 \
            --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.1" \
            -x data/db/aegenomics.db \ ## !!! DATA_BASE_PATH !!! ##
            -q $file_name.mapping.fastq.gz --no-unal \
            | samtools view -F 4 -b \
            | samtools sort -@ 32 -m 8G -o data/mapping_output/$file_name.sorted.bam
    fi
done
```

# session3 | filtering #
```bash
for file in ~/YOUR_DIRECTORY/data/mapping_output/*
do
    if [[ $file == *.sorted.bam ]]
    then
        file_name=$(basename $file .sorted.bam)
        sambamba markdup -r -t 5 -p $file_name.sorted.bam $file_name.sorted.rmdup.bam
        filterBAM --chunk-size 25 \
          --bam $file_name.sorted.rmdup.bam \
          -N \
          -r data/misc/aegenomics.db.len.map \  ## !!! DATA_BASE_PATH !!! ##
          -A 92 \
          -a 94 \
          -n 100 \
          -b 0.75 \
          -B 0.01 \
          -t 5 \
          --sort-memory 8G \
          --include-low-detection \
          --stats $file_name.stats.tsv.gz \
          --stats-filtered $file_name.tsv.gz \
          --bam-filtered $file_name.filtered.bam
    fi
done
```

## view the filtering result ##
```bash
conda deactivate
```
### chui pas convaincu de ce passage ###
```bash
for file in ~/YOUR_DIRECTORY/data/mapping_output/*
do
    if [[ $file == *.sorted.bam ]]
    then
        file_name=$(basename $file .sorted.bam)
        zcat $file_name.stats-filtered.tsv.gz \
          | csvtk cut -t -T -f "reference,n_reads,read_ani_mean,read_ani_std,coverage_mean,breadth,exp_breadth,breadth_exp_ratio,norm_entropy,norm_gini,cov_evenness,tax_abund_tad" \
          | csvtk grep -r -t -v -f reference -p _plas -p _mito \
          | csvtk sort -t -T -k "n_reads:Nr" \
          | tabview -
    fi
done

### ici si tu veux comparer faut voir avec la reférence, si tu veux le faire on en reparlera pck vu la quantité d'échantillons d'espèces différentes va falloir y réfléchir ultérieurement je pense ###
```bash
## je laisse le code du cours qu'on pourra adapter plus tard
#ecriture du fichier des réfs
printf "GCA_002781685.1\nIMGVR_UViG_3300027782_000260\nGCA_014380485.1" > ref-list.txt
#comparaison
getRPercId --bam PRI-TJPGK-CATN-160-162.sorted.rmdup.bam --reference-list ref-list.txt --threads 5 --sort-memory 8G

# faut activer l'env où y a bamcov (day1-2 je crois) ou l'installer
bamcov -w 0 -m GCA_002781685.1.bam 
```

# session 4 | damage #

```bash
conda activate metaDMG
```
```bash
metaDMG config --config-file config.local.yaml \
  --metaDMG-cpp /usr/local/bin/metaDMG-cpp \
  --parallel-samples 1 \
  --cores-per-sample 5 \
  --output-dir data/mapping_output/local \
  --max-position 35 \
  --min-similarity-score 0.92 \
  --damage-mode local \
  data/mapping_output/*.filtered.bam

metaDMG compute config.local.yaml

metaDMG convert --add-fit-predictions \
  --output outputs.csv.gz \
  --results data/mapping_output/local/results

metaDMG dashboard --results data/mapping_output/local/results/
```
