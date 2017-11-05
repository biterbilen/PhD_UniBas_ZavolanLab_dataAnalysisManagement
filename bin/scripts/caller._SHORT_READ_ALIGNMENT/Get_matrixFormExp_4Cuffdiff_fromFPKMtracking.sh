#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o LOG.Get_matrixFormExp_4Cuffdiff_fromFPKMtracking.sh_$TASK_ID
#$ -t 1-3

#parallelize
handler() {
	echo "Job $SGE_TASK_ID receives signal : $1"
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

PATH=$HOME/bin:$PATH

indir=Cuffdiff_trusted_RNAseq; tag=CuffdiffTrusted
moutdir=Scheiffele; #new_Scheiffele; #Yoana; 

types=(gene isoform tss_group);
outdir=$moutdir/Expression_matrix
mkdir -p $outdir; pushd $outdir

#	type=gene
type=${types[ $((SGE_TASK_ID-1)) ]}
file=`ls ../$indir/${type}s.fpkm_tracking`;
outf=$tag$type.raw

echo $outf

rand=$RANDOM;

#less $file | perl -e 'while(<>){chomp; @t=split; print $t[0]; $nsam=($#t-9)/3; for $i (1..$nsam) { print "\t", $t[$i*3+7]; } print "\t$t[4]\n"; } print STDERR $nsam; ' > $outf.w_geneSymbol 2> $rand
#cuffdiff v1.1.0 and on
less $file | perl -e 'while(<>){chomp; @t=split; print $t[0]; $nsam=($#t-8)/4; for $i (1..$nsam) { print "\t", $t[5+$i*4]; } print "\t$t[4]\n"; } print STDERR $nsam; ' > $outf.w_geneSymbol 2> $rand
fcount=`cat $rand`;
sed -i 's/_FPKM//g' $outf.w_geneSymbol;
sed -i 's/^tracking_id/id/' $outf.w_geneSymbol;
sed -i 's/gene_short_name$/gid/' $outf.w_geneSymbol;

less $outf.w_geneSymbol | cut -f 1-$((fcount +1)) > $outf
~bilebi00/_DIS3/scripts/matrix2dataframe.pl $outf 1-$((fcount )) exp > $outf.df
~bilebi00/bin/R --no-save --args $outf $outf.df < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/splom_density_exp.R > /dev/null

if [ "`echo $indir | grep clipz`" != "" ]; then
	~bilebi00/bin/R --no-save --args $outf < ~bilebi00/_KNOCKDOWNS/scripts/qn.R > /dev/null
	~bilebi00/_PAPD5/scripts/innerJoin.pl $outf.qn $outf.w_geneSymbol 1 1 1-$((fcount +1)),$((fcount + fcount + 3)) > $outf.qn.w_geneSymbol
	~bilebi00/_DIS3/scripts/matrix2dataframe.pl $outf.qn 1-$((fcount )) exp > $outf.qn.df
	~bilebi00/bin/R --no-save --args $outf.qn $outf.qn.df < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/splom_density_exp.R > /dev/null
fi

rm $rand;

if [ $moutdir == Stepanka ]; then #TODO do u want to compare it to clipz
	for i in ~bilebi00/_DIS3/data/clipz_RNAseq/*.raw*; do
		tag=clipzQNgene;
		if [ "`echo $i | grep ".raw.qn"`" == "" ]; then
			tag=clipzRAWgene;
		fi
		outf=`basename $i`;
		outf=`perl -e '$tag=shift; $outf=shift; $outf =~ s/RNAseq/$tag/; print $outf' $tag $outf`;
		echo $outf;
		ln -sf $i $outf; 
	done
fi

echo DONE
