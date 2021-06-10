# Import libraries
import pysam
import numpy as np
import math
import sys
import os
from math import log
from pathlib import Path

bamName=sys.argv[1]
bamBase = os.path.basename(bamName)
bamBase = bamBase.replace(".bam","_")
print(bamBase)
path = Path(bamName)
svPath = path.parent.absolute()
path = svPath.parent.absolute()

sys.stdout = open(str(path) + '/quickAlStats/' + bamBase + 'quickAlStats.txt', 'w')

# Define variables
readCount=0
unmappedCount=0
suppleCount=0
secondCount=0
suppleBps=0
primaryReads=0
primaryBps=0
bps=0
unmappedReads=0
unmappedbps=0
alignedbps=0
supplealignedbps=0
readLenList=[]
qualityList=[]
alignQuals=[]
percentIdent=0.0

def errs_tab(n):
    """Generate list of error rates for qualities less than equal than n."""
    return [10**(q / -10) for q in range(n+1)]

tab=errs_tab(128)

# Import BAM file
bamfile = pysam.AlignmentFile(bamName, "rb")

# Look for unmapped reads (ONLY WORKS IN THIS LOOP W/ PYSAM UNTIL_EOF
for r in bamfile.fetch(until_eof=True):
    if r.is_unmapped:
        unmappedReads += 1
        unmappedbps += r.query_length                        

# Second read in bamfile to return all relevant variables
for read in bamfile.fetch():
    bps += read.query_length
    readLenList.append(read.query_length)

    if not read.is_unmapped:
        readCount += 1

    if read.is_unmapped:
        pass
    else:
        # Calculate accuracy
        el_count = [0] * 10
        for el in read.cigartuples:
            el_count[el[0]] += el[1]
        nn = 0
        if read.has_tag("nn"):
            nn = read.get_tag("nn")
        nm = read.get_tag("NM")
        dels = el_count[2]
        ins = el_count[1]
        mismatches = nm - dels - ins - nn
        matches = el_count[0] - mismatches
        accuracy = float(matches) / (matches + nm) * 100
        alignQuals.append(accuracy)
        # Calculate quality
        quality = read.query_qualities
        sum_prob = 0.0
        if quality:
            mq = -10 * log(sum([tab[q] for q in quality]) / len(quality), 10)
            qualityList.append(mq)

        
    if read.is_unmapped:
        unmappedCount += 1
    elif read.is_supplementary:
        suppleCount += 1
        suppleBps += read.query_length
        supplealignedbps += read.query_alignment_length
    elif read.is_secondary:
        secondCount += 1
    else:
        primaryReads += 1
        primaryBps += read.query_length
        alignedbps += read.query_alignment_length

# Perform calculations to return relevant statistics
medPercentIdent = np.median(alignQuals)        
meanPercentIdent = np.mean(alignQuals)
totalAligned = alignedbps + supplealignedbps
meanReadLen = np.mean(readLenList)
medianReadLen = np.median(readLenList)
meanReadQual = np.mean(qualityList)
medianReadQual = np.median(qualityList)
readLenList.sort()
n50 = readLenList[np.where(np.cumsum(readLenList) >= 0.5 * np.sum(readLenList))[0][0]]

# Print variables
print("Total/Mapped Reads: " + str(readCount))
print("Unmapped Reads: " + str(unmappedReads))
print("Primary Reads: " + str(primaryReads))
print("Supplementary Reads: " + str(suppleCount))
print("Secondary Reads: " + str(secondCount))
print("Total Base pairs: " + str(bps))
print("Primary Base Pairs: " + str(primaryBps))
print("Supplementary Base Pairs: " + str(suppleBps))
print("Unmapped Base Pairs: " + str(unmappedbps))
print("Aligned Base Pairs: " + str(alignedbps))
print("Aligned Supplementary Base Pairs: " + str(supplealignedbps))
print("Total Aligned Base Pairs: " + str(totalAligned))
print("Mean Read Length: " + str(round(meanReadLen,1)))
print("Median Read Length: " + str(medianReadLen))
print("Mean Read Quality: " + str(round(meanReadQual,1)))
print("Median Read Quality: " + str(round(medianReadQual,1)))
print("N50: " + str(n50))
print("Median Percent Identity: " + str(round(medPercentIdent,1)))
print("Mean Percent Identity: " + str(round(meanPercentIdent,1)))

sys.stdout.close()
