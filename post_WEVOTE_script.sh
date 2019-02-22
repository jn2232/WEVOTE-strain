#!/bin/bash



#First, extract strain reads from the previous FASTA file
#TODO LIST: 
# Change 562 for species of interest.
# Skip the oneline part
# Obviously need to generalize it and change names so they're not ecoli


if [ $# -eq 0 ]
then
    echo "Error: Must provide TAXID as first argument, exiting."
    exit 0
fi

taxid=562

#1. This line turns the fasta into a single line
#echo "Converting fasta file to single-line fasta"
#awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' wvcontigs.fa > oneline_contigs.fa

#Load all the taxids iinto a bash list
input="rumino_species.txt"
i=0
while read F ; do
    taxlist[ $i ]="$F"
    ((i++))
done < "$input"


#taxid=1264 # a rumino sp
outfile="rumino_contigs.fa"
rm $outfile

#2.This awk line will search the WEVOTE details for a given taxid (classified by WEVOTE) and print out the READ_ID. It will then search for that READ_ID in the original fasta file and feed it into a new fasta
for i in "${taxlist[@]}"
do
    echo "Extracting reads from taxid $i"
    awk -v taxid="$i" '{ if ($NF == taxid) {print $1} }' wvdetails.txt | grep -w -A 1 -f - oneline_contigs.fa >> $outfile
echo 'done!'
done

#3 Execute the kmer code

#echo "Executing main script"
#cd $PWD/src/
#./mladen_main ../DATAFILES/ecoli_reads.fa
