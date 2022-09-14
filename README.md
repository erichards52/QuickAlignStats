# QuickAlignStats

## In accordance with NanoPlot, read length, accuracy and quality is calculated using primary and supplementary reads

An up-to-date version of this script can be found at https://gitlab.com/genomicsengland/bio-pipelines/long-reads-qc/-/blob/master/bin/quick_aln_stats.py.

Quick and dirty Python script which returns alignment metrics from a BAM file

## How to use

MAKE SURE YOU ARE CURRENTLY IN AN INTERACTIVE JOB OR SUBMIT THIS AS AN LSF JOB TO A QUEUE WITH >20000 MB MEMORY AVAILABLE (i.e.: bsub -R rusage[mem>20000] select[mem=20000] -M 20000 -P bio -q pipeline)

Simply run the script by calling `./quickAlignStats.py ${path/to/bamfile}` or `./quickAlignStats_Per_Alignment.py ${path/to/bamfile}`.

If you want to perform this for multiple bamfiles, use the `stats_call.sh` script (make sure the correct quickAlignStats python script is being used as to whether you want per alignment metrics or not):

`./stats_call.sh ${path/to/project/or/experiment/id}`

## Returns:
#### Total/Mapped Reads, Unmapped Reads, Primary Reads, Supplementary Reads, Secondary Reads, Total Base pairs, Primary Base Pairs, Supplementary Base Pairs, Unmapped Base Pairs, Aligned Base Pairs, Aligned Supplementary Base Pairs, Total Aligned Base Pairs, Mean Read Length, Median Read Length, Mean Read Quality, Median Read Quality, N50, Median Percent Identity, Mean Percent Identity

## Calculation for quality metrics
```
qualityList=[]

quality = read.query_qualities
sum_prob = 0.0
if quality:
    mq = -10 * log(sum([tab[q] for q in quality]) / len(quality), 10)
    qualityList.append(mq)

```

## Calculation for mean and median alignment accuracy
```
alignQuals=[]

for read in bamfile.fetch():

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
```

## Calculation for mean and median read lengths
```
readLenList=[]

for read in bamfile.fetch():
    readLenList.append(read.query_length)
    
meanReadLen = np.mean(readLenList)
medianReadLen = np.median(readLenList)
```

## Calculation for mean and median quality scores
```
meanReadQual = np.mean(qualityList)
medianReadQual = np.median(qualityList)
```

## Calculation for N50:
```
readLenList.sort()
n50 = readLenList[np.where(np.cumsum(readLenList) >= 0.5 * np.sum(readLenList))[0][0]]
```
