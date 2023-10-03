# tmux a
# ctrl + b, d
#!/bin/bash

# Symbolic links
for f in ~/course/data/day2/fastq/PRI-TJPGK*; do ln -s $f .; done

for file in *.fq.gz
do
    vsearch --fastx_uniques $file --fastqout ${file}.vs.fq --minseqlength 30 --strand both
    gzip ${file}.vs.fq
    zcat ${file}.vs.fq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $file.vs.fq.gz.read_length.txt
    bowtie2 --threads 5 -k 100 -x ~/course/data/shared/mapping/db/aegenomics.db -U $file.vs.fq.gz --no-unal | samtools view -bS - > $file.bam
done

# Sort bam files with samtools
conda activate day2
for f in *.fq.gz.bam; do samtools sort -n ${f} -@ 5 > ${f}.sort.bam; done

# Do lca
conda activate metaDMG
for f in *.sort.bam; do metaDMG-cpp lca -bam ${file} -names ~/course/data/shared/mapping/taxonomy/names.dmp -nodes ~/course/data/shared/mapping/taxonomy/nodes.dmp -acc2tax ~/course/data/shared/mapping/taxonomy/acc2taxid.map.gz -weighttype 1 -fix_ncbi 0 -out ${f}; done

# Create config file
metaDMG config *.sort.bam --names ~/course/data/shared/mapping/taxonomy/names.dmp --nodes ~/course/data/shared/mapping/taxonomy/nodes.dmp --acc2tax ~/course/data/shared/mapping/taxonomy/acc2taxid.map.gz -m /usr/local/bin/metaDMG-cpp
vim config.yaml
# Edit config with chose parameters (view in script)

# Compute
metaDMG compute config.yaml 

# View on Dashboard
metaDMG dashboard config.yaml

# LCA Loop Isolated
for f in *.sort.bam;
do metaDMG-cpp lca -bam $f -names ~/course/data/shared/mapping/taxonomy/names.dmp -nodes ~/course/data/shared/mapping/taxonomy/nodes.dmp -acc2tax ~/course/data/shared/mapping/taxonomy/
acc2taxid.map.gz -weighttype 1 -fix_ncbi 0 -out $f.lca;
done
