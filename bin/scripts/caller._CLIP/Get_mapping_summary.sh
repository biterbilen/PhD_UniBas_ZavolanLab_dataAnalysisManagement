
#TODO set
ids=(1341 1342 238 239 240 288 289 308 309)
outdir=Analysis/Project_EWSR1/Mapping_summary
ids=(445 446 308 309 1639 1640 1628 1637 1629 1638)
outdir=Analysis/Project_DIS3L2/Mapping_summary

mkdir -p $outdir; pushd $outdir

for id in ${ids[@]}; do
	rsync -t mirz@web08:~/BITER/Summary/$id.general_summary .
done

