#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatch_EWSR1_pipeline.sh 
# 
#         USAGE: sbatch_EWSR1_pipeline.sh
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 05/26/2015 17:15
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

config() {
	nc=1
	t=1
	infile=~/Projects/biter_biter_EWSR1/data/clipz/sampleData289.tar.gz
	infile=~/Projects/biter_biter_EWSR1/data/clipz/sampleData5645.tar.gz
	infile=~/Projects/biter_biter_EWSR1/data/clipz/sampleData314.tar.gz

	sampleId=`basename $infile .tar.gz | sed 's/sampleData//g'`
	outdir=clipz_UniqueGenomicAlignments/$sampleId
}

cmd() {
	# TODO run xlinkEnrichment scripts
	cmd_xlinkEnrichment=$(cat <<END

	echo cmd_xlinkEnrichment 
END
)

	cmd_CLIPZ2UniqueGenomicAlignments=$(cat <<END
ln -sf $infile .;
tar -zxvf $infile; 
echo "name chr_start-1 chr_end id copies strand alignment sequence annotation genome_error anno_error five_utr_fract cds_fract three_utr_fract exon_fract intron_fract unknown_fract seq_len" | sed 's/ /\t/g' > $sampleId.bedplus;
less mapped_sequences | awk '\$12==1 { print}' | sort -k 1b,1 > unique_mapped;
less genome_mappings | sort -k 2b,2 > genome_mappings.sorted;
join -t $'\t' -1 1 -2 2 unique_mapped genome_mappings.sorted | \
	awk 'BEGIN{OFS="\t"; }{ print \$18,\$22,\$23,\$1,\$3,\$21,\$24,\$2,\$15,\$13,\$14,\$4,\$5,\$6,\$7,\$8,\$9}' >> $sampleId.bedplus;
	gzip $sampleId.bedplus;
	rm mapped_sequences unmapped_sequences genome_mappings annotation_mappings transcript_expression unique_mapped genome_mappings.sorted;
	mv * ..;
	rm -rf $sampleId
END
)

sbatch -J cuga$sampleId -p normal -t $t:0:0 -c $nc <<SBATCH
#!/bin/bash -l	
echo "$cmd_CLIPZ2UniqueGenomicAlignments"
$cmd_CLIPZ2UniqueGenomicAlignments


echo DONE
SBATCH

}

dummy() {

	echo 
#rm mapped_sequences unmapped_sequences genome_mappings annotation_mappings transcript_expression $sampleId.bedplus;
#less mapped_sequences | awk '\$12==1 || NR==1{ print}' | \
#	~/PI_HOME/Applications/bilebi00/Project_EWSR1/scripts/_PAPD5/innerJoin_stdin.pl genome_mappings 1 2 '19,23,24,1,3,22,25,2,15,13,14,4,5,6,7,8,9' | \
#	awk 'BEGIN{OFS="\t";}{ \$2=\$2-1; print \$0,\$3-\$2 }' >> $sampleId.bedplus;

}

config
mkdir -p $outdir; pushd $outdir
cmd
popd

