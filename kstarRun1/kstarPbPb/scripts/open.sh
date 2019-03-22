folder=$1

for left in -1 0 1
do
    for right in -1 0 1
    do
	open ${folder}_l${left}_r${right}/tpc*/*all.png
    done
done
