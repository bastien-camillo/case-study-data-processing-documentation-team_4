
for file in *.fq.gz
do
    vsearch --fastx_uniques $file --fastqout $file.vs.fq --minseqlength 30 --strand both
    gzip $file
    zcat $file.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $file.gz.read_length.txt
    bowtie2 --threads 5 -k 100 -x ~/course/data/shared/mapping/db/aegenomics.db -U $file.gz --no-unal | samtools view -bS - > $file.bam
done

# tmux a
# ctrl + b, d
