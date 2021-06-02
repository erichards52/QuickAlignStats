projDir=$1

for d in $(find ${projDir} -type d -name 'aligned_minimap')
do
    svDir="$(dirname $d)"
    mkdir -p ${svDir}/quickAlStats
    echo "${svDir}/quickAlStats"
    bsub -R "select[mem>15000] rusage[mem=15000]" -M 15000 -P bio -q pipeline -o %J.out -e %J.err python quickAlStats.py ${d}/*.bam 
done
