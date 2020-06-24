import os
import subprocess as sp
import fileinput
import sys
import argparse
import shutil


"""
HOWTO:

This script downloads the required files from the NCBI database, including the taxonomy files (same for all species), followed by the complete genomes in NCBI for the given TAXID
Folder named $dbname_dir and located in cwd

USAGE:

python bsdb_dlGenomes.py $taxid $dbname $dbrootdir

"""


def getArgs():
    """ Function to get input arguments, currently just <taxid> and <dbname>"""
    parser = argparse.ArgumentParser(description='', epilog='')

    # 3 arguments required
    parser.add_argument("taxid", help="This is the taxonomic ID", type=int)
    parser.add_argument("dbname", help="This is the database name")
    parser.add_argument( "dbroot", help="Main folder for all databases, has taxonomy folder")


    args = parser.parse_args()
    taxid = args.taxid
    dbname = args.dbname
    dbroot = args.dbroot
    return taxid, dbname, dbroot


def download_taxonomy_stuff(dbroot, taxpath):
    if not os.path.exists(taxpath):
        os.mkdir(taxpath)
    os.chdir(taxpath)

    # Download accession files
    print("Accession to Taxon")
    if not os.path.isfile("nucl_gb.accession2taxid.gz"):
        os.system(
            "wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz")

    # Download and extract taxonomy dump
    print("Taxonomy Tree")
    if not os.path.isfile("taxdump.tar.gz"):
        print("Extracting... (This takes a while)")
        os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
        os.system("gunzip *.gz")
        os.system("tar -xvf taxdump.tar")

    # Download assembly summary
    print("Download assembly summary")
    if not os.path.isfile("assembly_summary.txt"):
        os.system(
            "wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt")
    # Finished downloading everything for taxonomy folder.
    os.chdir('..')
    return


# STEP 1: Get arguments and path
taxid, dbname, dbroot = getArgs()
print("inside bsdb_dl")
print("taxid",taxid)
print("dbname",dbname)
print("broot",dbroot)
bsdb_folder = os.path.dirname(os.path.abspath(__file__))

""" Locate taxonomy folder or create one """

# STEP 2: Check for path (if no tax, dl it)
print("DBROOT IS: ", dbroot)
dbroot = os.path.abspath(dbroot)
taxpath = os.path.join(dbroot, 'taxonomy')
tax_flag = 0
# if os.path.isdir(taxpath):
for (root, dirs, files) in os.walk(dbroot, topdown=True):
    for dirname in dirs:
        if dirname == 'taxonomy':
            print("Taxonomy already exists.")
            taxpath = os.path.join(root, dirname)
            print(taxpath)
            tax_flag = 1
            break
    if tax_flag == 1:
        break

    break

if tax_flag == 0:
    print("No taxonomy found, need to download taxonomy materials")
    # Execute the taxonomy download
    download_taxonomy_stuff(dbroot, taxpath)


# STEP 3 Clean folder creation

# Generate database folder for species-of-interest
speciesdb = sys.argv[2] + "_dir"
speciesdb = os.path.join(dbroot, speciesdb)
if os.path.exists(speciesdb):
    shutil.rmtree(speciesdb)
    os.mkdir(speciesdb)

# Enter taxpath to grab urls
os.chdir(taxpath)
if os.path.isfile('urls.tmp'):
    os.remove('urls.tmp')
if os.path.isfile('urls.txt'):
    os.remove('urls.txt')


# STEP 4: Extract genome types from assembly
# Grab the genome links prior to download
tmpfile_name = 'urls.tmp'
tmpfile = open('urls.tmp', 'w')

""" Potential Options for genome curation """

# Option 1.
# Complete genome only
# tax_string = 'BEGIN { FS="\t"} { if ($6!=%d && $7==%d && ($12=="Complete Genome" )) print $20;}' % (taxid,taxid)

# Option 2.
# Complete genome and chromosomes
# tax_string = 'BEGIN { FS="\t"} { if ($6!=%d && $7==%d && ($12=="Complete Genome" || $12 =="Chromosome" )) print $20;}' % (taxid,taxid)

# Option 3.
# All types of genomes except when straintax = speciestax pulled (other options at bttm of code)
# tax_string = 'BEGIN { FS="\t"} { if ($6!=%d && $7==%d) print $20;}' % (taxid,taxid) \

# Option 4.
# All genomes pulled

# Option 5.
# Complete genomes but remove the option of strain_id != species_id
# tax_string = 'BEGIN { FS="\t"; OFS="\t'; } { if ($7==%d &&
# $12=="Complete Genome" ) print $6,$7,$8,$9,$20;}' % (taxid)

# Option 6. (this is the taxid-free option)
# Complete genome but strain taxid doesnt have to be different than species
# tax_string = 'BEGIN { FS="\t"; OFS="\t" }{ if ($7==%d && $12=="Complete Genome" ) \
    # print $6,$7,$8,$9,$20;}' % (taxid)


#TODO TODO TODO: SHLEX SPLIT OMG
#TODO: I can adopt 144:150 and do a shlex split on the awk string (it works I checked)
#TODO: Turn the Popen into run below
#TODO: The code below can very easily be looped
#TODO: Alternate option 1 to option 5 can be implemented similarly
#TODO: Improve string formatting
# Grabbing from assembly_summary and putting into urls
# import pdb; pdb.set_trace()  # XXX BREAKPOINT
# tax_string = 'BEGIN { FS="\t"; OFS="\t" }{ if ($7==%d && $12=="Complete Genome" ) \
#     print $6,$7,$8,$9,$20;}' % (taxid)

# STEP 4b: Current selection is to select all COMPLETE GENOMES in REFSEQ
#TODO: Make sure all these tmpfiles or whatever are in their respective folders, not in taxdir
# belonging to the SPECIES (no other criteria necessary)
tax_string = 'BEGIN { FS="\t" ; OFS="\t"} { if ($6!=%d && $7==%d && ($12=="Complete Genome" )) print $6,$7,$8,$9,$20;}' % (taxid,taxid)

# Perform awk and grab info for the partial urls (urls.tmp)
#TODO fix this lol
cmd = ["awk", tax_string, "assembly_summary.txt"]
sp.run(cmd, check=True, stdout=tmpfile)  # only runs in python 3.5+

#### Just a 4am sidenote I could have done this the entire time.
# blah = """
#     awk 'BEGIN { FS="\t"; OFS="\t" } { if ($7==1305 && $12=="Complete Genome" ) print $6,$7,$8,$9,$20; }' assembly_summary.txt
#     """
# xx = shlex.split(blah)
# cc = sp.run(xx)


# STEP 5: Take each column and put it in a file
# Take columns out of tmpfile_name and put them in their correctfile
# TODO: clean this up
with open('strain_taxa.txt', 'w+') as wf:
    col_output = sp.Popen(
        ['cut', '-f1', tmpfile_name],
        stdout=wf, universal_newlines=True).communicate()[0]


with open('species_taxa.txt', 'w+') as wf:
    col_output = sp.Popen(
        ['cut', '-f2', tmpfile_name],
        stdout=wf, universal_newlines=True).communicate()[0]


with open('species_name.txt', 'w+') as wf:
    col_output = sp.Popen(
        ['cut', '-f3', tmpfile_name],
        stdout=wf, universal_newlines=True).communicate()[0]


with open('strain_name.txt', 'w+') as wf:
    col_output = sp.Popen(
        ['cut', '-f4', tmpfile_name],
        stdout=wf, universal_newlines=True).communicate()[0]


with open('url_snippet.txt', 'w+') as wf:
    col_output = sp.Popen(
        ['cut', '-f5', tmpfile_name],
        stdout=wf, universal_newlines=True).communicate()[0]



# STEP 6: Wordcount? idk
# Quick word count of the original
# TODO is this necessary? fix it
count = 0
for line in open('url_snippet.txt'):
    count += 1
print("Grabbed genomes from assembly summary")
print("Found ", count, " results")

#STEP 7: process urls
# Delete header prior to download
findme = "ftp://ftp.ncbi.nlm.nih.gov/genomes/"
with fileinput.input('url_snippet.txt', inplace=True) as f:
    for line in f:
        print(line.replace(findme, '').rstrip())

# Append suffix to url so now its a full URL in urls.txt
f2 = open('urls.txt', 'w')
cmd2 = ["awk", "-F/", '{print $0"/" $6 "_genomic.fna.gz"}', "url_snippet.txt"]
sp.run(cmd2, check=True, stdout=f2)
# Now we have a urls.txt file to rsync from

# STEP 8: Download files
# Call rsync and download all the genomes
lib_added = os.path.join(speciesdb, 'library', 'added')
os.system("mkdir -p " + lib_added)
# rs_command = "rsync --no-motd --files=from=urls.txt rsync://ftp.ncbi.nlm.nih.gov/genomes
sp.run(["rsync", "--no-motd", "--files-from=urls.txt",
        "rsync://ftp.ncbi.nlm.nih.gov/genomes/", "."])

# STEP 9: Extract
""" Fasta extraction and folder rearrangement """
# Extract fasta files
os.system('gunzip -rf all/*')

# STEP 10: Move folders
# Grab all the downloaded+extracted fastas and move to library/added/*.fna
for (root, dirs, files) in os.walk(taxpath, topdown=True):
    for ff in files:
        if ff.endswith(".fna"):
            print(os.path.join(root, ff))
            os.rename(os.path.join(root, ff), lib_added + '/' + ff)

alldir = os.path.join(taxpath, 'all')
os.system('rm -rf ' + os.path.join(taxpath, 'all'))
print("Genomes all downloaded, part 1 complete")


# EXTRA NOTES
##################################

#### 1. Explanation for awk tax_string: ####

# $6 is taxid (if at strain level then shouldn't be equal to species taxid
# $7 is parent taxid (aka should be equal to species)
# $12 is assembly_level which should be a complete genome or chromosome
# $20 is the ftp file location
