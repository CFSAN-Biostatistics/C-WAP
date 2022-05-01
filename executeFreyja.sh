#! /bin/bash -e


originalPath=$PYTHONPATH
unset PYTHONPATH
freyjaBin=/projects/covidtrakr/software/miniconda3/envs/freyja-env/bin/freyja


echo Submitting Freyja jobs...
for suffix in $(ls freyjaOutput/freyja.depths.* | awk -F '.' '{print $3}'); do
	echo $suffix
	srun -c 2 --mem 8G $freyjaBin demix freyjaOutput/freyja.variants.tsv.$suffix freyjaOutput/freyja.depths.$suffix \
			--output freyjaOutput/freyja.demix.$suffix &
done

echo Waiting for all Freyja calls to finish... 
wait

echo All Freyja computations now over.
export PYTHONPATH=$originalPath

