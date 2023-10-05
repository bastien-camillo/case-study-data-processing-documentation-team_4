# Setup, extension, deplication and mapping

We installed two things:

mamba install bbmap
pip install tabview

This is the suggested workflow for microbes
extension -> dereplication -> mapping -> filtering -> damage estimation

This tool gives you stats on the reads
seqkit stats -j 5 -T data/fastq/PRI-TJPGK-CATN-160-162.fq.gz | csvtk -t pretty

Then we extended the reads, with a necessary calculation of memory first

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

Then we deplicated the reads

seqkit rmdup -j 5 -s -o PRI-TJPGK-CATN-160-162.extended.derep.fq.gz PRI-TJPGK-CATN-160-162.extended.fq.gz

Then we calculated the different max read and average read lengths, as well as number of sequences in each to see the effects. You can use cat.

Then we switched back to the orginal reads for mapping.

filterbyname.sh in=data/fastq/PRI-TJPGK-CATN-160-162.fq.gz out=PRI-TJPGK-CATN-160-162.mapping.fastq.gz names=PRI-TJPGK-CATN-160-162.extended.derep.fq.gz threads=5 overwrite=t include=t

Then we activated mapping.
conda activate mapping

and performed the mapping using BOWTIE

bowtie2 -p 5 -k 100 -D 10 -R 2 \
    -N 0 -D 5 -R 1 -L 22 -i S,0,2.50 \
    --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.1" \
    -x data/db/aegenomics.db \
    -q PRI-TJPGK-CATN-160-162.mapping.fastq.gz --no-unal \
    | samtools view -F 4 -b \
    | samtools sort -@ 32 -m 8G -o PRI-TJPGK-CATN-160-162.sorted.bam

Paramater explanation:
There's a lot going on in this last command. First, with -k we ask for 100 alignments per read. This is because many of those reads are going to be mapping equally well to multiple references. We can adjust the % identity in --score-min "L,0,-0.1", where the -0.1 == 90, 0.05 == 95%
We will do a fast search by using -N 1 -D 5 -R 1 -L 22 -i S,0,2.50. Those are the preset parameters for a fast search in Bowtie2, but with the difference that we are using -N 1 to allow a mismatch in the seed to increase sensitivity and we will use.
The part -i S,0,2.5 sets the interval function f to f(x) = 1 + 2.5 * sqrt(x), where x is the read length. It means how many seed substrings will generate.


This is just an example for the tutorial, in real life you might want to use:


recommended, good compromise between sensitivity, specificity and speed: -D 15 -R 2 -N 1 -L 22 -i S,1,1.15 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.1"


faster: -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.1"


more sensitive, but ~10X slower: -D 15 -R 2 -N 1 -L 20 -i S,1,1.15  --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.1"


super sensitive, but ~20X slower (Reduce -L if you want to be more sensitive): -D 15 -R 3 -N 1 -L 20 -i S,1,0.5  --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.1"


# filtering
sambamba to remove duplicates (you can also use picard)
sambamba markdup -r -t 5 -p PRI-TJPGK-CATN-160-162.sorted.bam PRI-TJPGK-CATN-160-162.sorted.rmdup.bam

Then we filtered using bamfilter

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

  Paramters: There's a lot going on in this last command. We use -N to sort the filtered BAM files by name, so it can be used by metaDMG. The `-r` specifies where we can find the non-concatenated length of the references for the abundance estimation. With `-A` we will only keep those reads with at least 92% ANI and `-a`  will keep those references with at least 94% of average ANI. The `-n` will keep only those references with at least 100 reads mapping. 

  conda deactivate

  check out the filtered results

  zcat PRI-TJPGK-CATN-160-162.stats-filtered.tsv.gz \
  | csvtk cut -t -T -f "reference,n_reads,read_ani_mean,read_ani_std,coverage_mean,breadth,exp_breadth,breadth_exp_ratio,norm_entropy,norm_gini,cov_evenness,tax_abund_tad" \
  | csvtk grep -r -t -v -f reference -p _plas -p _mito \
  | csvtk sort -t -T -k "n_reads:Nr" \
  | tabview -

  You can use bamcov to check coverage patterns.

  bamcov -w 0 -m GCA_002781685.1.bam 

  # Damage

There are methods that can estimate damage over thousands of taxa (LCA mode), references (local mode) or provide a global estimate (global mode).

We used the Local mode to have an overview of the potential references that might show damage in the sample. First, we need to create the config file:
metaDMG config --config-file PRI-TJPGK-CATN-160-162.local.yaml \
  --metaDMG-cpp /usr/local/bin/metaDMG-cpp \
  --parallel-samples 1 \
  --cores-per-sample 5 \
  --output-dir PRI-TJPGK-CATN-160-162.local \
  --max-position 35 \
  --min-similarity-score 0.92 \
  --damage-mode local \   <--- change this to global or lca to run those
  PRI-TJPGK-CATN-160-162.filtered.bam

compute
metaDMG compute PRI-TJPGK-CATN-160-162.local.yaml

convert to csv
metaDMG convert --add-fit-predictions \
  --output PRI-TJPGK-CATN-160-162.csv.gz \
  --results PRI-TJPGK-CATN-160-162.local/results

dashboard

metaDMG dashboard --results PRI-TJPGK-CATN-160-162.local/results/


# Functional Profiling

Workflow
mmseq ----> xfilter ----> anvio


  

  

