#!/resources/tools/apps/software/lang/Python/3.7.4-GCCcore-8.3.0/bin/python
# Import libraries
import pysam
import numpy as np
import math
import sys
import os
from math import log
from pathlib import Path

bamName=sys.argv[1]

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
alignQualsSNV=[]
alignQualsDel=[]
alignQualsIns=[]
alignQualsNn=[]
percentIdent=0.0
qualityListPrim=[]

def errs_tab(n):
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
        cig_stats_counts = read.get_cigar_stats()[0]
        try:
            nn = read.get_tag('nn')
        except KeyError:
            nn = 0
        nm = cig_stats_counts[10]
        ins = cig_stats_counts[1]
        dels = cig_stats_counts[2]
        mismatches = nm - dels - ins - nn
        matches = cig_stats_counts[0] - mismatches
        accuracy = matches / (matches + nm) * 100
        alignQuals.append(accuracy)
        # Accuracy only SNVs/mismatches and no deletions and insertions
        accuracy_SNV = float(mismatches) / (matches + nm) * 100
        accuracy_Del = float(dels) / (matches + nm) * 100
        accuracy_Ins = float(ins) / (matches + nm) * 100
        accuracy_Nn = float(nn) / (matches + nm) * 100
        alignQualsSNV.append(accuracy_SNV)
        alignQualsDel.append(accuracy_Del)
        alignQualsNn.append(accuracy_Nn)
        alignQualsIns.append(accuracy_Ins)
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
        qualityPrim = read.query_qualities
        sum_prob = 0.0
        if qualityPrim:
            mqPrim = -10 * log(sum([tab[q] for q in qualityPrim]) / len(qualityPrim), 10)
            qualityListPrim.append(mqPrim)


# Perform calculations to return relevant statistics

medPercentIdentIns = np.median(alignQualsIns)
meanPercentIdentIns = np.mean(alignQualsIns)
medPercentIdentNn = np.median(alignQualsNn)
meanPercentIdentNn = np.mean(alignQualsNn)
medPercentIdentDel = np.median(alignQualsDel)
meanPercentIdentDel = np.mean(alignQualsDel)
medPercentIdentSNV = np.median(alignQualsSNV)
meanPercentIdentSNV = np.mean(alignQualsSNV)
medPercentIdent = np.median(alignQuals)        
meanPercentIdent = np.mean(alignQuals)
totalAligned = alignedbps + supplealignedbps
meanReadLen = np.mean(readLenList)
medianReadLen = np.median(readLenList)
meanReadQual = np.mean(qualityList)
medianReadQual = np.median(qualityList)
minReadQual = np.amin(qualityListPrim)
maxReadQual = np.amax(qualityListPrim)
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
print("Lowest Average Read Quality: " + str(round(minReadQual,1)))
print("Highest Average Read Quality: " + str(round(maxReadQual,1)))
print("N50: " + str(n50))
print("Median Percent Identity: " + str(round(medPercentIdent,4)))
print("Mean Percent Identity: " + str(round(meanPercentIdent,4)))
print("Median SNV Error Rate: " + str(round(medPercentIdentSNV,4)))
print("Mean SNV Error Rate: " + str(round(meanPercentIdentSNV,4)))
print("Median INS Error Rate: " + str(round(medPercentIdentIns,4)))
print("Mean INS Error Rate: " + str(round(meanPercentIdentIns,4)))
print("Median DEL Error Rate: " + str(round(medPercentIdentDel,4)))
print("Mean DEL Error Rate: " + str(round(meanPercentIdentDel,4)))
print("Median NN Error Rate: " + str(round(medPercentIdentNn,4)))
print("Mean NN Error Rate: " + str(round(meanPercentIdentNn,4)))
