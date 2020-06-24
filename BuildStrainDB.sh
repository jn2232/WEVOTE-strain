#!/bin/bash -
#===============================================================================
#
#          FILE: BuildStrainDB.sh
#
#         USAGE: ./BuildStrainDB.sh [-f taxid.file] dbrootpath taxid(s)
#             EX (taxid file): ./BuildStrainDB.sh -f mytaxidfile databasefolder/
#             EX (single taxa): ./BuildStrainDB.sh databasefolder/ 1305
#             EX (list of taxa): ./BuildStrainDB.sh databasefolder/ 1305 1280 727
#
#   DESCRIPTION:
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: YOUR NAME (),
#  ORGANIZATION:
#       CREATED: 06/09/2020 02:45
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

function read_taxfile {
    taxid=()
    # Read through file
    while read -r line
    do
        # Remember: == is GLOB while =~ is REGEX (different rules!)
        [[ $line == \#* ]] && continue
        taxid+=($line)
    done < $1
}

#-------------------------------------------------------------------------------
### ARGUMENT/PARAMETER HANDLING
#-------------------------------------------------------------------------------
# a. check for any args
if [ $# -eq 0 ]; then
    echo "Missing Options!"
    echo "(run $0 -h for help)"
    echo ""
    exit 1
fi
ARGS=""
# b. process the args (rn just input file)
while getopts 'f:' options
do
    case $options in
        f) infile="$OPTARG"
            ;;
        ?)
            echo "script usage: $(basename $0) [-f taxid_file]" >&2
            exit 1
            ;;
    esac
done
shift $((OPTIND -1))


#-------------------------------------------------------------------------------
# POSITIONAL PARAMETERS
#-------------------------------------------------------------------------------
running_directory=$(pwd)
dbroot=$1
echo $dbroot

# d. Get positional variables
if [ -n "${infile+x}" ]; then # if -f given
    read_taxfile $infile #returns taxid
else #no f given, allow taxid to be list at the end
    shift
    taxid=("$@")
fi
echo 'single'
echo ${taxid[0]}

# e. Check the position variables
# echo "dbroot: " $dbroot
# echo "taxid: " ${taxid[*]}

##TODO: old
kmerlen=31
# f. Optional Kmer parameter
# if [ -z "$4" ]; then
#     kmerlen=31
# else
#     kmerlen=$3
# fi
#########END PARAMS ###############

#-------------------------------------------------------------------------------
# TRAVERSE
#-------------------------------------------------------------------------------
# Go to folder containing this script
# NOTE: Only works in linux, readlink-f is different in mac
scripts_dir="$(dirname "$0")"
scripts_dir="$(readlink -f $scripts_dir)" #make absolute
echo "This script is located in: " $scripts_dir

#-------------------------------------------------------------------------------
#CALL BSDB SCRIPTS FOR EACH TAXA
#-------------------------------------------------------------------------------


# Loop through each taxa possible
for tax in "${taxid[@]}"; do
    echo $tax
    cd $running_directory
    echo "Building Database for: " $tax
    dbname=$tax
    dbfolder="${dbname}_dir"
    dbfolder_jelly="${dbfolder}/jelly"

    # Run
    echo "Step 1: Download Taxonomic and Genomic References"
    echo $tax $dbfolder $dbroot
    echo "in Buildstrain"
    python $scripts_dir"/bsdb_dlGenomes.py"  $tax $dbname $dbroot
    echo "Step 1 Complete: Files stored in ${dbname}_dir"

    pwd

    echo "Step 2: "Generate K-mers for each strain""
    bash $scripts_dir"/bsdb_countKmers.sh" $dbfolder $kmerlen
    echo "Step 2 Complete: K-mer dumps stored in $dbfolder_jelly"

    echo "after step 2 location"
    pwd

    cd $dbfolder_jelly

    echo "Step 3:"K-mer filtering and converting to python dictionary""
    python $scripts_dir"/bsdb_hashKmers.py"
    echo "Step 3 Complete for : $tax"
done
