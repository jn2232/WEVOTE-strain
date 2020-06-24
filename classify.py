#!/usr/bin/env python

import argparse
import os
import sys
from collections import Counter

try:
    import numpy as np
    import pickle
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    sys.stderr.write("Error: Missing Packages.\n")
    sys.exit(1)



def getArgs():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_fasta", type=str)
    ap.add_argument("--db", help="", type=str, required=True,)
    ap.add_argument("--ftype", help="", \
        choices=["fasta", "fastq"], type=str, required=True)
    ap.add_argument("--athresh", help="", type=float, default=0.001)
    ap.add_argument("--paired", action="store_true", help="")
    ap.add_argument("input_fasta2", type=str, nargs="?")
    ap.add_argument("--kmerLen", help="", type=int, default=31)
    return vars(ap.parse_args())


def paired():
    """function to concatenate two paired
    files into one and return it back into the script"""
    # Important to know that SeqIO parse is a generator function so once you loop through f1 and f2 you change the type of them.

    params = getArgs()
    # Load files and zip
    fasta1 = SeqIO.parse(open(params["input_fasta"]), params["ftype"])
    fasta2 = SeqIO.parse(open(params["input_fasta2"]), params["ftype"])
    pairedFasta = zip(fasta1, fasta2)

    # Regenerate new fasta file with new string
    records = []
    for f1, f2 in pairedFasta:
        # Need to remove phred quality scores in fastq in order to concatenate them
        f1.letter_annotations = {}
        f1.seq = Seq(str(f1.seq) + "NNNNNNNN" + str(f2.seq))
        f1.description = ""
        rec = SeqRecord(f1.seq, id=f1.id)
        records.append(rec)
    SeqIO.write(records, "pairedfasta.tmp", "fasta")
    return 0


annot = {}
params = getArgs()
dbpath = os.path.join(params["db"], "jelly", "strainDB.pickle")
kmerLen = params["kmerLen"]
print("Loading Database from %s" % (dbpath))
with open(dbpath, "rb") as a:
    kdb = pickle.load(a)
print("Database Loaded")
def wvStrainClassify():
    # dbpath,genomesize_dict = wvStrainLoadDB(params['db'])


    if params["paired"] and not params["input_fasta2"]:
        sys.stderr.write("Require two files for paired option\n")
        sys.exit(0)
    if not params["db"]:
        sys.stderr.write("Database .pickle file is required!\n")
        sys.exit(0)
    if not os.path.isfile(params["input_fasta"]):
        sys.stderr.write(
            " Cannot find the input file at: "
            + params["input_fasta"]
            + "\n Please located the correct file location"
        )
        sys.exit(1)

    #  Classify to the Database
    print("Classify")

    numreads = 0

    if params["input_fasta2"]:  # Paired Reads
        print("Obtaining paired read information..")
        print("Reverse reads found at..", params["input_fasta2"])
        paired()
        fastaObject = SeqIO.parse(open("pairedfasta.tmp"), "fasta")
    else:  # Unpaired Reads
        fastaObject = SeqIO.parse(open(params["input_fasta"]), params["ftype"])

    # Read Loop
    for fastaRead in fastaObject:
        numreads += 1
        read_id, read_seq = fastaRead.id, str(fastaRead.seq)
        # debugging
        # if numreads % 1000 != 0:
        #     continue
        readTaxList = []
        # K-mers Loop
        n_kmer = len(read_seq) - kmerLen + 1
        for i in range(0, n_kmer):
            kmer = read_seq[i : i + kmerLen]
            currTaxList = kdb.get(kmer)
            # print('kmer and currtaxlist')
            # print(i,kmer,len(readTaxList))
            if currTaxList is not None:  # Aka k-mer not found
                readTaxList.append(currTaxList)

                # print('broke')
            # else:  # Add the taxid of the k-mer
            # readTaxList.append(currTaxList)

        """Determine taxa designation based on k-mer hits"""

        # mostCommonTuples[Index][0 = taxid, 1 = count ]
        flattenedTaxa = [y for x in readTaxList for y in x]

        # Try to build annot to look at abundance before selecting taxid
        annot[read_id] = flattenedTaxa

        cc = Counter(flattenedTaxa)

        # Case 1: No k-mers found for the entire read
        if not len(flattenedTaxa):
            annot[read_id] = None

        # Case 2: Exactly 1 item for most common
        elif len(cc.most_common()) == 1:
            mostCommonTuples = cc.most_common()
            annot[read_id] = [mostCommonTuples[0][0]]

        # Case 3: Updated Case 2 attempt
        elif len(cc.most_common()) > 1:
            mostCommonTuples = cc.most_common()
            # Add this taxid to the dictionary with the read_id as the key
            annot[read_id] = [mostCommonTuples[0][0]]
            # The first one will have the highest count, set that as the base value
            topCount = mostCommonTuples[0][1]
            # Add taxids tied with topCount
            for i in range(1, len(mostCommonTuples)):
                if mostCommonTuples[i][1] == topCount:
                    annot[read_id].append(mostCommonTuples[i][0])
                else:
                    break

        else:
            print("Error:This shouldn't exist")
            exit(1)

    return 0 #annot, numreads


def wvStrainReadBreakdown(annot):

    unambiguous_reads = 0
    ambig_reads = 0
    unclassif_reads = 0
    for key in annot:
        if isinstance(annot[key], list):
            if len(annot[key]) > 1:
                ambig_reads += 1
            else:
                unambiguous_reads += 1
        else:
            unclassif_reads += 1

    print("Total Number of reads:", len(annot))
    print("Number of reads identified unambiguously: ", unambiguous_reads)
    print("Number of reads with multiple strain hits: ", ambig_reads)
    print("Number of reads not identified with a reference strain: ", unclassif_reads)
    return 0


def wvStrainAbundance(annot, numreads):

    # After classification,
    # each read will have either a max-hit taxid, or an n-way tie, or no hit.
    # This will be read- taxid output eventually
    # print("=========== READ ANNOTATION before filtration==============")
    # for read in annot:
    #     print(read, annot[read])
    # print("=========== end ANNOTATION before filtration==============")

    print("BEFORE WE BREAK TIES")
    # Generate single list of all hits
    priorProb = []
    for x in annot.values():
        if isinstance(x, list) and len(x) == 1:
            for ii in x:
                if not isinstance(ii, type(None)):
                    priorProb.append(ii)

    # Calculation of abundance before breaking ties
    ppCount = Counter(priorProb)
    sumCounts = sum(ppCount.values())
    normCount = ppCount.copy()
    for key in normCount:
        normCount[key] /= sumCounts
    print(ppCount)
    print(normCount)

    # Break ties (Assign taxonomic IDs) according to max probability.
    for read in annot:
        if isinstance(annot[read], list):
            # annot[read] = [t1,t2,t3]
            probList = [ppCount[x] for x in annot[read]]  # [p1, p2, p3]

            #### THIS WILL RETURN 0 IF THERE IS NO MAX
            ind = np.argmax(probList)
            if ind == 0 and probList[ind] == 0:
                annot[read] = "womp"
            else:
                annot[read] = annot[read][ind]  # now its equal to the max prob one

            #### OMG ######

    # Calculation of abundance after breaking ties
    print("AFTER BREAKING TIES")
    ##### COPIED HERE ######
    postProb = []
    for x in annot.values():
        if not isinstance(x, type(None)):
            postProb.append(x)

    postCount = Counter(postProb)
    sumCounts = sum(postCount.values())
    normpostCount = postCount.copy()
    for key in normpostCount:
        normpostCount[key] /= sumCounts
    print(postCount)
    print(normpostCount)

    #### END COPIED HERE #####

    # Calculate abundance after THRESHOLD
    params = getArgs()
    THRESH = params["athresh"]
    print(THRESH)
    print(len(annot))
    thresh_annot = {k: v for k, v in annot.items() if normpostCount[v] > THRESH}
    print(len(thresh_annot))

    ##### COPIED HERE 2######
    postpostProb = []
    for x in thresh_annot.values():
        if not isinstance(x, type(None)):
            postpostProb.append(x)

    postCount = Counter(postpostProb)
    sumCounts = sum(postCount.values())
    normpostCount = postCount.copy()
    for key in normpostCount:
        normpostCount[key] /= sumCounts
    print("after_thresholding")
    print(postCount)
    print(normpostCount)
    return normpostCount

    # Generate frequencies of each taxa identified(MLE)

    # mladen new note:
    # okay but isnt this just gonna be counter
    # yes


#     freq = {}
#     taxList = list(annot.values())
#     for item in taxList:
#         if item in freq:
#             freq[item] += 1
# else:
#     freq[item] = 1

# # Move Unclassified out of dict
# try:
#     nreads_unclassified = freq[None]
#     del freq[None]
# except KeyError:
#     print('None Unclassified')
#     nreads_unclassified=0


# # Filter out low hits
# belowThreshKeys = list()
# params = getArgs()
# THRESH = params['athresh']

# # Generate list of taxa to delete
# for key, value in freq.items():
#     if value < THRESH * numreads:
#         belowThreshKeys.append(key)

# # Go through dict and delete these
# for key in belowThreshKeys:
#     if key in freq:
#         del freq[key]

# # This will be read- taxid output eventually
# print("================ READ ANNOTATION===================")
# # for read in annot:
# #     print(read, annot[read])
# print("================ end ANNOTATION===================")


# return freq, nreads_unclassified


def wvStrainOutput(annot, numreads, relABpost):  # , freq, uReads):
    freq = relABpost
    # Calculate filtered abundances
    print("Abundance Info")
    for key, value in freq.items():
        abund = round(value, 3)
        # abund = round(value / sum(freq.values()), 3)
        print("% s : % s" % (key, abund))

    # print("Classified reads")
    # print(sum(freq.values()))

    # print("Unclassified Reads")
    # print(uReads)

    # print("Total # of Reads")
    # print(numreads)
    return 0


def wvStrainLoadDB(dbfolder):
    dbpath = os.path.join(dbfolder, "jelly", "strainDB.pickle")

    # Genome Sizes
    # gsize = os.path.join(dbfolder,"genome_size.dat")
    # with open(gsize) as infile:
    #     reader = csv.reader(infile, delimiter="\t")
    #     next(reader)
    #     lines = list(reader)
    #     genomesize = dict(lines)
    return dbpath  # , genomesize


def wvStrainAlternateAbundance(annot):
    #### abundance
    clear_annot = {}
    tied_annot = {}
    no_annot = {}

    # Divide the reads
    for k, v in annot.items():
        if not isinstance(v, type(None)) and len(v) == 1:
            clear_annot.update({k: v})
        elif not isinstance(v, type(None)) and len(v) > 1:
            tied_annot.update({k: v})
        else:
            no_annot.update({k: v})

    list_clear_annot = []
    for vals in clear_annot.values():
        for vvals in vals:
            list_clear_annot.append(vvals)

    # Generate hit count of clear reads
    cc_clear_annot = Counter(list_clear_annot)
    total_hits_clear = sum(cc_clear_annot.values())

    # Generate frequency (normalized hits)
    norm_clear_annot = cc_clear_annot.copy()
    for k in norm_clear_annot.keys():
        norm_clear_annot[k] /= total_hits_clear

    print("This is the abundance of only clear hits taxa")
    print(norm_clear_annot)

    print("There are this many ties:", len(tied_annot))
    print("\n")
    PendingDeprecationWarning

    # Loop through ties and add partial counts
    update_abund = cc_clear_annot.copy()
    for k, v in tied_annot.items():
        hitlist = []
        for vvals in v:
            hit = norm_clear_annot.get(str(vvals))
            if hit:
                hitlist.append(vvals)
                # update_abund[vvals] += norm_clear_annot[vvals]
                # print("This is the updated dict after the following k:V pair",k,vvals)
                # print(update_abund)
        if len(hitlist) > 0:
            hitcount = 1.0 / len(hitlist)
            probList = [cc_clear_annot[x] for x in hitlist]
            ind = np.argmax(probList)
            if ind == 0 and probList[ind] == 0:
                print("No max")
            # else:
            #     update_abund[hitlist[ind]] += 1.0

            # Break ties (Assign taxonomic IDs) according to max probability.
            for hh in hitlist:
                update_abund[str(hh)] += hitcount
            # else:
            #     print('taxa not in top clear hits')

    relab_sum = sum(update_abund.values())
    for k in update_abund.keys():
        update_abund[k] /= relab_sum

    print(update_abund)

    # Calculate abundance after THRESHOLD
    params = getArgs()
    THRESH = params["athresh"]
    print(THRESH)
    thresh_annot = {k: v for k, v in update_abund.items() if v > THRESH}
    print(thresh_annot)

    # Sum again
    sum41 = sum(thresh_annot.values())
    for key in thresh_annot:
        thresh_annot[key] /= sum41
    ta = Counter(thresh_annot)
    print("\n")
    print("\n")
    print("after thresh")
    print(ta)
    # print(ta.most_common(5))
    print("\n")
    print("\n")

    ##### COPIED HERE 2######
    # postpostProb = []
    # for x in thresh_annot.values():
    #     if not isinstance(x, type(None)):
    #         postpostProb.append(x)

    # postCount = Counter(postpostProb)
    # sumCounts = sum( postCount.values() )
    # normpostCount = postCount.copy()
    # for key in normpostCount:
    #     normpostCount[key] /= sumCounts
    # print("after_thresholding")
    # print(postCount)
    # print(normpostCount)

    return ta


def wvStrainPrintRawReads(annot, filename):
    with open(filename, "w") as fh:
        for read, taxlist in annot.items():
            if isinstance(taxlist, list):
                for eachtax in taxlist:
                    fh.write(read + "\t" + str(eachtax) + "\n")

            else:
                fh.write(read + "\t" + str(taxlist) + "\n")

    return 0


def main():
    # Call main functions
    params = getArgs()
    print(params)
    dbf = params["db"]
    # annot, numreads = wvStrainClassify()
    wvStrainClassify()
    # wvStrainPrintRawReads(annot)

    # print("NEW ABUNDANCE HERE")
    # newAB_post=wvStrainAlternateAbundance(annot)
    # wvStrainOutput(annot, numreads,newAB_post) # freq, uReads)
    # Look at output of classified vs unclassified

    # print("OG ABUNDANCE HERE")
    # numreads=1

    # relAB_post = wvStrainAbundance(annot, numreads)
    # wvStrainPrintRawReads(annot, dbf)
    # wvStrainReadBreakdown(annot)

    # wvStrainOutput(annot, numreads,relAB_post) # freq, uReads)


    # with open('annot3.pickle','wb') as h3:
    #     pickle.dump(annot,h3)
    # with open('annot3.pickle','rb') as h3:
    #     annot=pickle.load(h3)

if __name__ == "__main__":
    main()
