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
bamBase = os.path.basename(bamName)
bamBase = bamBase.replace(".bam","_")
print(bamBase)
path = Path(bamName)
align_path = path.parent.absolute()
sv_workflow_path = align_path.parent.absolute()

sys.stdout = open(str(sv_workflow_path) + '/quickAlStats/' + bamBase + 'quickAlStats_per_alignment.txt', 'w')

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

# Print column names
print("Read_Name, Alignment_Length, Supplementary_Alignment, Secondary_Alignment, Primary_Alignment, Alignment_Mismatches, Alignment_Ins, Alignment_Dels, Alignment_SNVs, Alignment_Matches, Alignment_Identity, Alignment_SNV_Error_Rate, Alignment_Del_Error_Rate, Alignment_Ins_Error_Rate, Alignment_Nn_Error_Rate, Alignment_Quality")

# Second read in bamfile to return all relevant variables
for read in bamfile.fetch():
    bps += read.query_length
    readLenList.append(read.query_length)

    if not read.is_unmapped:
        readCount += 1

    if read.is_unmapped:
        pass
    else:
        key = read.query_name
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
        mqPrim=0
        print(key, read.query_length, True, False, False, nm, ins, dels, mismatches, matches, accuracy, accuracy_SNV,accuracy_Del,accuracy_Ins,accuracy_Nn, mq,sep=",")
    elif read.is_secondary:
        secondCount += 1
       	mqPrim=0
        print(key, read.query_length, False, True, False, nm, ins, dels, mismatches, matches, accuracy, accuracy_SNV,accuracy_Del,accuracy_Ins,accuracy_Nn, mq,sep=",")
    else:
        primaryReads += 1
        primaryBps += read.query_length
        alignedbps += read.query_alignment_length
        qualityPrim = read.query_qualities
        sum_prob = 0.0
        if qualityPrim:
            mqPrim = -10 * log(sum([tab[q] for q in qualityPrim]) / len(qualityPrim), 10)
            qualityListPrim.append(mqPrim)
        print(key, read.query_length, False, False, True, nm, ins, dels, mismatches, matches, accuracy, accuracy_SNV,accuracy_Del,accuracy_Ins,accuracy_Nn, mq)
