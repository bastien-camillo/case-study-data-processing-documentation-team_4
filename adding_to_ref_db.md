Addition of references to reference database

conda activate day1

Nuclear and mito genomes can be downloaded from any source here we use NCBI.

You will need to download Ovis aries mitochondrion, complete genome and save as a file called NC_001941.1.fa in your working directory.

Now lets check the accession number/ID in the header of the faasta

head -1 NC_001941.1.fa
And then check that the accession to taxid matches in the ncbi taxonomy

zgrep NC_001941 ../../../data/shared/taxonomy/acc2taxid.map.gz
Normally these would fit, but for the course here we made a custom taxonomy/database and changed the accession slightly. So let's edit the header to match our taxonomy.

and change it to NC_001941.1

vim NC_001941.1.fa
bowtie2-build NC_001941.1.fa NC_001941.1.fa

for file in *vs.fq.gz
do
bowtie2 -U ../PRI-TJPGK-CATN-224-226.fq.gz.vs.fq.gz -x NC_001941.1.fa --no-unal --threads 5 | samtools view -bS - > PRI-TJPGK-CATN-224-226.fq.gz.vs.fq.gz.ovis.bam 
done &> ovis_map.log
before competitive mapping

samtools depth -a PRI-TJPGK-CATN-224-226.fq.gz.vs.fq.gz.ovis.sort.bam

samtools depth -a PRI-TJPGK-CATN-224-226.fq.gz.vs.fq.gz.ovis.sort.bam | cut -f3

samtools depth -a PRI-TJPGK-CATN-224-226.fq.gz.vs.fq.gz.ovis.sort.bam | cut -f3 | datamash mean 1


if you forget to set the -a option the depth it will be depth per covered base
samtools depth -a PRI-TJPGK-CATN-224-226.fq.gz.vs.fq.gz.ovis.sort.bam | cut -f3 | datamash mean 1
