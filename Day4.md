
# session 1 | setting the stage #

```bash
cd ~/course/wdir
mkdir -p day4
cd day4

ln -s ~/course/data/day4/ data
```
```bash
conda activate day4
```

# session 2 | extension, dereplication and mapping #

```bash
mamba install bbmap
pip install tabview
```


seqkit stats -j 5 -T data/fastq/PRI-TJPGK-CATN-160-162.fq.gz | csvtk -t pretty



MEM=$(loglog.sh seed=1234 k=17 in=data/fastq/PRI-TJPGK-CATN-160-162.fq.gz ignorebadquality 2> >(grep Cardinality) \
    | awk -vP=0.5 -vB=16 -vH=3 '{{print int( (((B*H)/8)*$2)/P )}}')

tadpole.sh -Xmx10G \
    k=17 \
    in=data/fastq/PRI-TJPGK-CATN-160-162.fq.gz \
    out=PRI-TJPGK-CATN-160-162.extended.fq.gz \
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




seqkit rmdup -j 5 -s -o PRI-TJPGK-CATN-160-162.extended.derep.fq.gz PRI-TJPGK-CATN-160-162.extended.fq.gz


filterbyname.sh in=data/fastq/PRI-TJPGK-CATN-160-162.fq.gz out=PRI-TJPGK-CATN-160-162.mapping.fastq.gz names=PRI-TJPGK-CATN-160-162.extended.derep.fq.gz threads=5 overwrite=t include=t


conda activate mapping


bowtie2 -p 5 -k 100 -D 10 -R 2 \
    -N 0 -D 5 -R 1 -L 22 -i S,0,2.50 \
    --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.1" \
    -x data/db/aegenomics.db \
    -q PRI-TJPGK-CATN-160-162.mapping.fastq.gz --no-unal \
    | samtools view -F 4 -b \
    | samtools sort -@ 32 -m 8G -o PRI-TJPGK-CATN-160-162.sorted.bam






#session3 | filtering

sambamba markdup -r -t 5 -p PRI-TJPGK-CATN-160-162.sorted.bam PRI-TJPGK-CATN-160-162.sorted.rmdup.bam


filterBAM --chunk-size 25 \
  --bam PRI-TJPGK-CATN-160-162.sorted.rmdup.bam \
  -N \
  -r data/misc/aegenomics.db.len.map \
  -A 92 \
  -a 94 \
  -n 100 \
  -b 0.75 \
  -B 0.01 \
  -t 5 \
  --sort-memory 8G \
  --include-low-detection \
  --stats PRI-TJPGK-CATN-160-162.stats.tsv.gz \
  --stats-filtered PRI-TJPGK-CATN-160-162.stats-filtered.tsv.gz \
  --bam-filtered PRI-TJPGK-CATN-160-162.filtered.bam


conda deactivate

  zcat PRI-TJPGK-CATN-160-162.stats-filtered.tsv.gz \
  | csvtk cut -t -T -f "reference,n_reads,read_ani_mean,read_ani_std,coverage_mean,breadth,exp_breadth,breadth_exp_ratio,norm_entropy,norm_gini,cov_evenness,tax_abund_tad" \
  | csvtk grep -r -t -v -f reference -p _plas -p _mito \
  | csvtk sort -t -T -k "n_reads:Nr" \
  | tabview -


printf "GCA_002781685.1\nIMGVR_UViG_3300027782_000260\nGCA_014380485.1" > ref-list.txt
getRPercId --bam PRI-TJPGK-CATN-160-162.sorted.rmdup.bam --reference-list ref-list.txt --threads 5 --sort-memory 8G



/home/antonio/opt/conda/envs/day2/bin/bamcov -w 0 -m GCA_002781685.1.bam 
/home/antonio/opt/conda/envs/day2/bin/bamcov -w 0 -m GCA_014380485.1.bam
/home/antonio/opt/conda/envs/day2/bin/bamcov -w 0 -m IMGVR_UViG_3300027782_000260.bam 
/home/antonio/opt/conda/envs/day2/bin/bamcov -w 0 -m data/misc/example-4.bam





#session 4 | damage


conda activate metaDMG


metaDMG config --config-file PRI-TJPGK-CATN-160-162.local.yaml \
  --metaDMG-cpp /usr/local/bin/metaDMG-cpp \
  --parallel-samples 1 \
  --cores-per-sample 5 \
  --output-dir PRI-TJPGK-CATN-160-162.local \
  --max-position 35 \
  --min-similarity-score 0.92 \
  --damage-mode local \
  PRI-TJPGK-CATN-160-162.filtered.bam

metaDMG compute PRI-TJPGK-CATN-160-162.local.yaml

metaDMG convert --add-fit-predictions \
  --output PRI-TJPGK-CATN-160-162.csv.gz \
  --results PRI-TJPGK-CATN-160-162.local/results

metaDMG dashboard --results PRI-TJPGK-CATN-160-162.local/results/



###or
metaDMG config --config-file PRI-TJPGK-CATN-160-162.lca.yaml \
  --custom-database \
  --names data/taxonomy/names.dmp \
  --nodes data/taxonomy/nodes.dmp \
  --acc2tax data/taxonomy/acc2taxid.map.gz \
  --metaDMG-cpp /usr/local/bin/metaDMG-cpp \
  --parallel-samples 1 \
  --cores-per-sample 5 \
  --output-dir PRI-TJPGK-CATN-160-162.lca \
  --max-position 35 \
  --lca-rank '' \
  --min-similarity-score 0.92 \
  --damage-mode lca \
  --weight-type 1 \
  PRI-TJPGK-CATN-160-162.filtered.bam
