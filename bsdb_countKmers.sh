#!/bin/bash -
#===============================================================================
#
#          FILE: bsdb_countKmers.sh
#
#         USAGE: ./bsdb_countKmers.sh
#
#   DESCRIPTION: This script is part 2 in database creation
# GOAL: is to create jellyfish k-mer counted files of each of the genomes downloaded in part 1
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: YOUR NAME (),
#  ORGANIZATION:
#       CREATED: 06/06/2020 19:33
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error



#---  FUNCTION  ----------------------------------------------------------------
#          NAME:  genome_size
#   DESCRIPTION:  optional to get genome_size.dat, Generate Genome Size Counts by TAXID
#    PARAMETERS: used to be a non-function (at the end of the file)
#       RETURNS: N/A but creates a file
#-------------------------------------------------------------------------------
genome_size ()
{
    if [ -e genome_size.tmp ]; then
        rm genome_size.tmp
    fi
    touch genome_size.tmp
    echo -e "taxID\tsize" >> genome_size.tmp

    for fasta in library/added/*.fna;
    do
        seqID=$(head -n1 $fasta | awk '{print $1}' | sed 's/^.\(.*\)..$/\1/')
        taxID=$(grep $seqID seq2tax.txt| awk '{print $2}' )
        genome_size=$(sed 1d $fasta | wc -m)
        echo -e "${taxID}\t${genome_size}" >> genome_size.tmp
    done

    testvar=5

    # Average out genome sizes of multiple seqs w/ same taxid
    if [ -e genome_size.dat ]; then
        rm genome_size.dat
    fi
    touch genome_size.dat
    echo -e "taxID\tsize" >> genome_size.dat
    awk '
        NR>1{
            arr[$1]   += $2
            count[$1] += 1
        }
        END{
            for (a in arr) {
                printf  a "\t" "%.0f\n",arr[a] / count[a]
            }
        }
    ' genome_size.tmp >> genome_size.dat
}	# ----------  end of function genome_size  ----------

#-------------------------------------------------------------------------------
# STEP 0: Arguments, preprocessing, etc
#-------------------------------------------------------------------------------
# Grab first argument as the database directory
dbdir=$1

# k-mer len = 31 by default or argument 2
if [ -z "$2" ]; then
    kmerlen=31
else
    kmerlen=$2
fi
cd $dbdir
# Cleanup
if [ -e seq2tax.txt ]; then
    rm seq2tax.txt
fi

mkdir jelly

#-------------------------------------------------------------------------------
# STEP 1: Grab file ID
#-------------------------------------------------------------------------------
for fasta in library/added/*.fna;
do
    # from >NC_320384028.1 blah blah -------> NC_320384028
    seqID=$(head -n1 $fasta | awk '{print $1}' | sed 's/^.\(.*\)..$/\1/')
    echo -e $seqID'\t'>> seq2tax.tmp
done
# This part is if you dont care about taxid
# cp seq2tax.tmp seq2tax.txt

#-------------------------------------------------------------------------------
# STEP 2: Grab taxonomy from file ID using nucl_gb.acc2taxid
#-------------------------------------------------------------------------------
if [ -f "seq2tax.txt" ]; then
    echo "seq2tax.txt already exists, delete to recreate mapping"
else
    touch seq2tax.txt
    echo -e "seqID\ttaxID" >> seq2tax.txt
    #  Faster to fgrep from a file than to place into loop
    LC_ALL=C fgrep -f seq2tax.tmp ../taxonomy/nucl_gb.accession2taxid | awk '{print $1,$3}'>>seq2tax.txt
fi
rm seq2tax.tmp


#-------------------------------------------------------------------------------
# STEP 3new: Call Jellyfish per file
#-------------------------------------------------------------------------------
# for fasta in library/added/*.fna;
# do
#     seqID=$(head -n1 $fasta | awk '{print $1}' | sed 's/^.\(.*\)..$/\1/')
#     echo $seqID
#     #  Generate Kmers
#     echo  "jellyfish dump for ${fasta} located in jelly/${seqID}.jdb"
#     if [ ! -f jelly/$seqID.jdb ]; then
#         jellyfish count -m $kmerlen -s 100M $fasta #make kmerdb
#         jellyfish dump mer_counts_0 > jelly/$seqID.jdb #  dumb db to a file
#         sed -i '/^>/d' jelly/$seqID.jdb #  remove the counter lines
#     else
#         echo "jdb file for ${seqID} exists, skipping"
#         echo "If you wish to rebuild the jellyfish kmer file, please delete"
#     fi
# done
#-------------------------------------------------------------------------------
# STEP 3old: Call Jellyfish per file
#-------------------------------------------------------------------------------
for fasta in library/added/*.fna;
do
    seqID=$(head -n1 $fasta | awk '{print $1}' | sed 's/^.\(.*\)..$/\1/')
    taxID=$(grep $seqID seq2tax.txt| awk '{print $2}' )
    echo $seqID $taxID
    #  Generate Kmers
    echo  "jellyfish dump for ${fasta} located in jelly/${taxID}.jdb"
    if [ ! -f jelly/$taxID.jdb ]; then
        jellyfish count -m $kmerlen -s 100M $fasta #make kmerdb
        jellyfish dump mer_counts_0 > jelly/$taxID.jdb #  dumb db to a file
        sed -i '/^>/d' jelly/$taxID.jdb #  remove the counter lines
    else
        echo "jdb file for ${taxID} exists, skipping"
        echo "If you wish to rebuild the jellyfish kmer file, please delete"
    fi
done

#-------------------------------------------------------------------------------
# CALL genome_size function (optional)
#-------------------------------------------------------------------------------
# genome_size





