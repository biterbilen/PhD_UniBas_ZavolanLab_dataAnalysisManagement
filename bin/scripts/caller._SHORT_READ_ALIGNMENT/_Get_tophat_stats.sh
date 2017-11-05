for i in TopHat_*/logs/bowtie.left*; do echo $i; cat $i; echo; done | less > TOPHAT.stats
