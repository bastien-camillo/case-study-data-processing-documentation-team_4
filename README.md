[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/-7_RZisP)

# Workflow for Mapping the Professor's Reads
## Team 4

Copied the Day2 Folder into our Working Directory
''' cp -r ~/course/data/day2 ~/course/wdir/mapping'''

ls to list files, you can use -l option for more details. The asterisk is a wildcard that denotes "anything" in this folder.

ls ~/course/data/day1/* .

You can work out the size of files in this directory

du -sk ~/course/data/day1/*

piping the output of cat to the word count command with count lines option

zcat file.fastq.gz | wc -l

This calculates reads as each read consists of four lines

zcat file.fastq.gz | awk '{s++}END{print s/4}'

Fastp can be used to trim the adapter files. We also filter out those less than 30 reads.

fastp -i file.fastq.gz -o file.trimmed.fastq -l 30





