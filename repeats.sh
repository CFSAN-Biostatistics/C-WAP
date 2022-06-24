
echo Initialising the node
srun -N 2 -w compute-dy-r54xlarge-1 compute-dy-r54xlarge-2  date

for i in {1..3}; do
    echo Attempt $i
    time CFSANonly/analyseFolder.sh -i ~/benchmark_data/ covidRefSequences/none.bed >> ./repeats.log
done

