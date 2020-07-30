#!/usr/bin/env python

import argparse
from datetime import datetime
import csv
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
    # ap.add_argument("--db", help="", type=str ) # DEBUG
    ap.add_argument("--ftype", help="",
                    choices=["fasta", "fastq"], type=str, default = "fasta")
    ap.add_argument("--athresh", help="", type=float, default=0.001)
    ap.add_argument("--paired", action="store_true", help="")
    ap.add_argument("input_fasta2", type=str, nargs="?")
    ap.add_argument("--abund_output", help="", type=str, )
    ap.add_argument("--WV_reads", help="", type=str, ) # comment this line out for debug mode
    ap.add_argument("--species_taxID", help="", type=int, ) # comment this line out for debug mode
    ap.add_argument("--kmerLen", help="", type=int, default=31)

    # debugin= "/Users/lab/Desktop/debugDB1/1strain_R1_trimmed.fastq" #DEBUG
    # debugin_2= "/Users/lab/Desktop/debugDB1/1strain_R2_trimmed.fastq" #DEBUG
    # debugdb= "/Users/lab/Desktop/debugDB1/debug_hflu_dir" #DEBUG
    # args = ap.parse_args(['--db',debugdb, '--paired', debugin_2, debugin]) #DEBUG
    args = ap.parse_args()
    return args

#TODO

# Open a new blank file, {params.species_taxID}_R1.fasta which will be written to with extracted reads
# Open original fwd read file
# Open WV details text file
# Search details text file line by line for species taxID
# If species taxID is found in line, get the name of the read from the file, which shoud be the first element
# Search through old read file for the read name, if found copy that line plus the next line (which should contain the actual nt sequence)
# Write both lines to the new fasta file
# {}_R#_extracted_reads.txt serves as an index containing the read name and taxID information for all extracted reads
def getReads(params):
    pattern = '\t{}\n'.format(params.species_taxID)
    R1_extracted_reads_index = open('{}_R1_extracted_reads.txt'.format(params.species_taxID), 'a+')
    newFasta_R1 = open('{}_R1.fasta'.format(params.species_taxID), 'a+')
    oldFastaR1 = open(params.input_fasta, 'r')
    search_reads = oldFastaR1.readlines()
    WVdetails = open(params.WV_reads, 'r')
    search_details = WVdetails.readlines()
    for i in range(0, len(search_details)):
        line = search_details[i]
        if pattern in line:
            # print("WV_output: ", line)
            readName = line.split('\t')[0]
            # print("WV_output:", readName)
            for j in range(0, len(search_reads)):
                line2 = search_reads[j]
                if readName in line2:
                    newFasta_R1.write(readName + '\n')
                    # print("Input_fasta: ", readName)
                    nextLine = search_reads[j + 1]
                    newFasta_R1.write(nextLine + '\n')
                    R1_extracted_reads_index.write(line)
                    # print("Input_fasta: ", nextLine)
                    j += 2
                    break
        i += 1
    print("Forward reads written to {}_R1.fasta".format(params.species_taxID))
    return 0

def getPairedReads(params):
    if params.paired:
        pattern = '\t{}\n'.format(params.species_taxID)
        R2_extracted_reads_index = open('{}_R2_extracted_reads.txt'.format(params.species_taxID), 'a+')
        newFasta_R2 = open('{}_R2.fasta'.format(params.species_taxID), 'a+')
        oldFastaR2 = open(params.input_fasta2, 'r')
        search_reads = oldFastaR2.readlines()
        WVdetails = open(params.WV_reads, 'r')
        search_details = WVdetails.readlines()
        for i in range(0, len(search_details)):
            line = search_details[i]
            if pattern in line:
                # print("WV_output: ", line)
                readName = line.split('\t')[0]
                # print("WV_output:", readName)
                for j in range(0, len(search_reads)):
                    line2 = search_reads[j]
                    if readName in line2:
                        newFasta_R2.write(readName + '\n')
                        # print("Input_fasta: ", readName)
                        nextLine = search_reads[j + 1]
                        newFasta_R2.write(nextLine + '\n')
                        R2_extracted_reads_index.write(line)
                        # print("Input_fasta: ", nextLine)
                        j += 2
                        break
            i += 1
        print("Reverse reads written to {}_R2.fasta".format(params.species_taxID))
        return 0
    else:
        print("No reverse reads specified")
        return 0

def wvsPaired(params):
    """function to concatenate two paired
    files into one and return it back into the script"""
    # Important to know that SeqIO parse is a generator function so once you loop through f1 and f2 you change the type of them.
    # To enter debug mode (when reads are already extracted, or a strain specific read file was generated) uncomment the commented lines beginning "fasta1" and "fasta2"
    # and comment out the line directly following each

    # Load files and zip
    # fasta1 = SeqIO.parse(open(params.input_fasta), params.ftype)
    fasta1 = SeqIO.parse(open('{params.species_taxID}_R1.fasta'), 'fasta') # New lines for pairing reads that have been extracted for original input fasta files and WEVOTE_details
    # fasta2 = SeqIO.parse(open(params.input_fasta2), params.ftype)
    fasta2 = SeqIO.parse(open('{params.species_taxID}_R2.fasta'), 'fasta')
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

def wvsLoadDB(dbfolder):
    dbpath = os.path.join(dbfolder, "jelly", "strainDB.pickle")

    # Genome Sizes
    # gsize = os.path.join(dbfolder,"genome_size.dat")
    # with open(gsize) as infile:
    #     reader = csv.reader(infile, delimiter="\t")
    #     next(reader)
    #     lines = list(reader)
    #     genomesize = dict(lines)
    return dbpath  # , genomesize

def wvsClassify(params):
    """ hi """
    # Move to load db function
    dbpath = wvsLoadDB(params.db)
    with open(dbpath, "rb") as dbpickle_handle:
        kdb = pickle.load(dbpickle_handle)
    print("Database Loaded")

    # Move to paired function
    if params.paired and not params.input_fasta2:
        sys.stderr.write("Require two files for paired option\n")
        sys.exit(0)
    if not params.db:
        sys.stderr.write("Database .pickle file is required!\n")
        sys.exit(0)
    if not os.path.isfile(params.input_fasta):
        sys.stderr.write(
            " Cannot find the input file at: "
            + params.input_fasta
            + "\n Please located the correct file location"
        )
        sys.exit(1)


    #TODO
    # To turn on debug mode, uncomment the two commented lines of code below, and comment out the line directly following each
    if params.input_fasta2:  # Paired Reads
        print("Obtaining paired read information..")
        # print("Reverse reads found at..", params.input_fasta2)
        print("Reverse reads found at..", 'newFasta_R2.fasta') ########### new line to print the location of extracted reverse reads 
        wvsPaired(params)
        fastaObject = SeqIO.parse(open("pairedfasta.tmp"), "fasta")
    else:  # Unpaired Reads
        # fastaObject = SeqIO.parse(open(params.input_fasta), params.ftype)
        fastaObject = SeqIO.parse(open('newFasta_R1.fasta'), 'fasta') ############## new line to parse sequence from new file containing extracted reads

    #  Classify to the Database
    print("Classify")

    # numreads = 0
    # Read Loop
    for fastaRead in fastaObject:
        params.numreads += 1
        read_id, read_seq = fastaRead.id, str(fastaRead.seq)
        readTaxList = []
        # K-mers Loop
        n_kmer = len(read_seq) - params.kmerLen + 1
        for i in range(0, n_kmer):
            kmer = read_seq[i : i + params.kmerLen]
            currTaxList = kdb.get(kmer)
            if currTaxList is not None:  # Aka k-mer not found
                readTaxList.append(currTaxList)
                # print('broke')
                # else:  # Add the taxid of the k-mer
                # readTaxList.append(currTaxList)

        """Determine taxa designation based on k-mer hits"""

        # mostCommonTuples[Index][0 = taxid, 1 = count ]
        flattenedTaxa = [y for x in readTaxList for y in x]

        # Try to build annot to look at abundance before selecting taxid
        params.annot[read_id] = flattenedTaxa

        cc = Counter(flattenedTaxa)

        # Case 1: No k-mers found for the entire read
        if not len(flattenedTaxa):
            params.annot[read_id] = None

        # Case 2: Exactly 1 item for most common
        elif len(cc.most_common()) == 1:
            mostCommonTuples = cc.most_common()
            params.annot[read_id] = [mostCommonTuples[0][0]]

        # Case 3: Updated Case 2 attempt
        elif len(cc.most_common()) > 1:
            mostCommonTuples = cc.most_common()
            # Add this taxid to the dictionary with the read_id as the key
            params.annot[read_id] = [mostCommonTuples[0][0]]
            # The first one will have the highest count, set that as the base value
            topCount = mostCommonTuples[0][1]
            # Add taxids tied with topCount
            for i in range(1, len(mostCommonTuples)):
                if mostCommonTuples[i][1] == topCount:
                    params.annot[read_id].append(mostCommonTuples[i][0])
                else:
                    break

        else:
            print("Error:This shouldn't exist")
            exit(1)
    print('hi finished classifying')

    return params

def wvsAbundance(params):
    """ hi """
    #hi
    # After classification,
    # each read will have either a max-hit taxid, or an n-way tie, or no hit.
    ##This will be read- taxid output eventually
    # print("=========== READ ANNOTATION before filtration==============")
    #for read in params.annot:
    #    #print(read, params.annot[read])
    #print("=========== end ANNOTATION before filtration==============")
    print("Calculating Abundance\n")

    #### 1. Look at unambiguous reads - use these for prior Probability ####
    # Generate single list of all hits
    priorProb = []
    # Remove list structure for single tax hits and add to priorProb
    for x in params.annot.values():
        if isinstance(x, list) and len(x) == 1:
            for ii in x:
                if not isinstance(ii, type(None)):
                    priorProb.append(ii)

    #### 2. Return this information from 1. into an abundance distribution
    ## Both in hits (ppCount) and percentage (normCount)
    # Calculation of abundance before breaking ties
    ppCount = Counter(priorProb)
    sumCounts = sum(ppCount.values())
    normCount = ppCount.copy()
    for key in normCount:
        normCount[key] /= sumCounts

    # Initial Display
    print("**Abundance calculations**")
    print("+++++++++++++++++++++++++++++++++++++++++++++++")
    print("\n")

    print("Initial distribution after classification (k-mer hit count):")
    print("-------------------------------------------------")
    print(ppCount)
    print("\n")

    print("Initial distribution after classification (percentage):")
    print("-------------------------------------------------")
    print(normCount)

    #### 3. Look at the many reads with ties, and pick the highest probability choice
    # Break ties (Assign taxonomic IDs) according to max probability.
    for read in params.annot:
        if isinstance(params.annot[read], list):
            # params.annot[read] = [t1,t2,t3]
            probList = [ppCount[x] for x in params.annot[read]]  # [p1, p2, p3]

            #### THIS WILL RETURN 0 IF THERE IS NO MAX
            # Aka if argmax does not exist, return womp
            #TODO big error was fixed here but look at it closely again
            ind = np.argmax(probList)
            if ind == 0 and probList[ind] == 0:
                params.annot[read] = "womp"
            else:
                params.annot[read] = params.annot[read][ind]  # now its equal to the max prob one


    #### 4. Calculation of abundance after breaking ties
    # Now we have a final tax for each read, we need to calculate the abundance
    print("AFTER BREAKING TIES")
    ##### COPIED HERE ######
    postProb = []
    for x in params.annot.values():
        if not isinstance(x, type(None)):
            postProb.append(x)

    postCount = Counter(postProb)
    sumCounts = sum(postCount.values())
    normpostCount = postCount.copy()
    for key in normpostCount:
        normpostCount[key] /= sumCounts

    # 4b.  After tiebreaking
    print("\n")
    print("Distribution after resolving ties (hitcount):")
    print("-------------------------------------------------")
    print(postCount)
    print("\n")


    print("Distribution after resolving ties (percentage):")
    print("-------------------------------------------------")
    print(normpostCount)
    print("\n")



    #### 5. Calculate abundance after THRESHOLD
    print("(Thresholding Details)")
    print("----------------------")
    THRESH = params.athresh
    print('threshold value:', THRESH)
    print('length of annot', len(params.annot))
    thresh_annot = {k: v for k, v in params.annot.items() if normpostCount[v] > THRESH}
    print('length of annot after thresholding', len(thresh_annot))

    #### 6. Calculate hit count and percentage after thresholding
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

    # After tiebreaking
    print("\n")
    print("Distribution after resolving ties and applying threshold (hitcount):")
    print("-------------------------------------------------")
    print(postCount)
    print("\n")


    print("Distribution after resolving ties and applying threshold (percentage):")
    print("-------------------------------------------------")
    print(normpostCount)

    def old():
        i = 0
        #### OLD
        #Generate frequencies of each taxa identified(MLE)

        # mladen new note:
        # okay but isnt this just gonna be counter
        # yes

        #     freq = {}
        #     taxList = list(params.annot.values())
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
        # THRESH = params.athresh
        ##

        q=0
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
        # # for read in params.annot:
        # #     print(read, params.annot[read])
        # print("================ end ANNOTATION===================")

        # return freq, nreads_unclassified
        return 0
    return postCount, normpostCount

#TODO
def wvsAlternateAbundance(params):
    """ abundance """
    #abundance
    clear_annot = {}
    tied_annot = {}
    no_annot = {}

    # Divide the reads
    for k, v in params.annot.items():
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
    THRESH = params.athresh
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

def wvsPrintRawReads(params, filename):
    """ Saves per read classification results to a file """
    with open(filename, "w") as fh:
        for read, taxlist in params.annot.items():
            if isinstance(taxlist, list):
                for eachtax in taxlist:
                    fh.write(read + "\t" + str(eachtax) + "\n")

            else:
                fh.write(read + "\t" + str(taxlist) + "\n")
    return 0

def wvsReadBreakdown(params):
    """ Identifies unambigious (classified), ambiguous (multi-class), unclassified """
    unambiguous_reads = 0
    ambig_reads = 0
    unclassif_reads = 0
    for key in params.annot:
        if isinstance(params.annot[key], list):
            if len(params.annot[key]) > 1:
                ambig_reads += 1
            else:
                unambiguous_reads += 1
        else:
            unclassif_reads += 1

    print("Total Number of reads:", len(params.annot))
    print("Number of reads identified unambiguously: ", unambiguous_reads)
    print("Number of reads with multiple strain hits: ", ambig_reads)
    print("Number of reads not identified with a reference strain: ", unclassif_reads)
    return 0

def wvsOutput(params, numreads, relABpost):  # , freq, uReads):
    freq = relABpost
    # Calculate filtered abundances
    print("Abundance Info")
    cwd = os.getcwd()
    cwd_str = "{}\n".format(cwd)
    header2 = "______________________________________________\n______________________________________________\n"
    header = "TaxID\tRA\n"
    time = "{:%a, %d %b %y, %I:%M:%S\n}".format(datetime.now())
    if params.abund_output:
        abund_outfile = True
    else:
        abund_outfile = False
    if os.path.isfile(params.abund_output) and abund_outfile:
        abundance_table = open(params.abund_output, 'a+')
        abundance_table.write(header2)
        abundance_table.write(time)
        abundance_table.write(cwd_str)
        abundance_table.write(header)
        for key, value in freq.items():
            abund = round(value, 3)
        # abund = round(value / sum(freq.values()), 3)
            print("% s : % s" % (key, abund))
            taxid_abund = "{}\t{}\n".format(key, abund)
            abundance_table.write(taxid_abund)        
    else:
        if abund_outfile:
            abundance_table = open(params.abund_output, 'a+')
            abundance_table.write(header2)
            abundance_table.write(time)
            abundance_table.write(cwd_str)
            abundance_table.write(header)
            for key, value in freq.items():
                abund = round(value, 3)
            # abund = round(value / sum(freq.values()), 3)
                print("% s : % s" % (key, abund))
                taxid_abund = "{}\t{}\n".format(key, abund)
                abundance_table.write(taxid_abund)
        else:
            abundance_table = open('strain_abundance_output', 'a+')
            abundance_table.write(header2)
            abundance_table.write(time)
            abundance_table.write(cwd_str)
            abundance_table.write(header)
            for key, value in freq.items():
                abund = round(value, 3)
            # abund = round(value / sum(freq.values()), 3)
                print("% s : % s" % (key, abund))
                taxid_abund = "{}\t{}\n".format(key, abund)
                abundance_table.write(taxid_abund)

    print("Classified reads")
    print(sum(freq.values()))

    # print("Unclassified Reads")
    # print(uReads)

    print("Total # of Reads")
    print(params.numreads)
    return


def main():
    """ hi """
    #call main functions
    params = getArgs()

    # add params
    params.annot = {}
    params.numreads = 0
    print(params)
    print()
    print("Running Classify function from main rn")

    # get reads from WEVOTE
    getReads(params)
    getPairedReads(params)

    # classifying
    params = wvsClassify(params)
    numreads=1

    # abundance calculation post-classification
    hitcount_abundance, pct_abundance = wvsAbundance(params)
    newAB_post=wvsAlternateAbundance(params)

    # per-read classification output to a file
    ########################################
    ############## I didnt change the name of the output file for strain classified reads, but I commented out what I wrote it as to be consistent with what I wrote for getReads, but I'm guessing those need to be edited to be consistent with the rest of the code lol!
    ########################################
    fasta_name = os.path.basename(params.input_fasta) #get filename from filepath
    ############## fasta_name = 'species{}_reads.fasta'.format(params.species_taxID)
    outfile_base = os.path.splitext(fasta_name)[0] #remove extension
    wvsPrintRawReads(params,outfile_base)

    # summary of read annotation (classified, unclassified, etc)
    wvsReadBreakdown(params)

    # print out results
    wvsOutput(params, numreads, pct_abundance) # freq, uReads)

    return


#TODO
# look at 2classify and 3classify and build std model

if __name__ == "__main__":
    print()
    main()
