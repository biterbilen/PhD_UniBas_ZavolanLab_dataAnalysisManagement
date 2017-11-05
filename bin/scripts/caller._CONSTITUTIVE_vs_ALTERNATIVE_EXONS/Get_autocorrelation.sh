#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -l sjpn=1
#$ -q fs_long@@qc_nehalem
# $ -q fs_long@@high_mem_node
#$ -l h_vmem=30G
# $ -l s_vmem=5G
# $ -l h_vmem=1M
# $ -notify
#$ -j y
#$ -cwd
#$ -o log_file
#$ -t 1-98

handler() {
echo "Job $SGE_TASK_ID receives signal : $1" >> log_file
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

chrs=(chr14 chr20 chr6_random chr17_random chr17 chr21 chr19 chr21_random chr15 chr7_random chr7 chr13 chr16 chr2 chr6 chr5 chr19_random chr10 chr8 chr1 chr12 chrX chr11 chr22 chr4_random chr9 chrY chr3 chrM chr18 chr1_random chr4 chr3_random chr15_random chr2_random chr22_random chr5_random chrX_random chr13_random chr10_random chr8_random chr9_random chr16_random chr11_random chr22_h2_hap1 chr6_cox_hap1 chr18_random chr5_h2_hap1 chr6_qbl_hap2);

export PERL5LIB=~bilebi00/lib:\$PERL5LIB;

outdir=output_directory
nhood=nhood
bins=bins

binsize=$(($nhood / $bins));	

mkdir -p $outdir; cd $outdir;

for tag in STRICT_ALTERNATIVE_EXON STRICT_CONSTITUTIVE_EXON; do
	echo $tag
	rm -rf $tag*.corr;
	for chr in ${chrs[@]}; do
		for str in + -; do
			echo $chr$str
			inp=$chr$str.$tag.nhood$nhood.inp
			cls=$inp.cls
			#select regions inside exon nhood
			less profile_distribution | grep $tag | grep -P "$chr\s\\$str" | perl -e 'my $nhood=shift; while(<>){chomp; @t=split; next if ($t[1]<-$nhood or $t[1]>$nhood); $d=$t[1]; $s=1; $s=-1 if($t[4] eq "-"); splice(@t,0,9); print join("\t","SS_0",@t,$d),"\n"; for (1..$nhood) { $t[2]-=$s; $t[3]-=$s; print join("\t","SS_$_",@t,$d),"\n"; }}' $nhood > $inp

			#cluster them
			~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $inp > $cls

			#extract T2C based on nucleotide shift of distance
			~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/parse_cls_for_autocorr_shifts.pl $cls $bins $binsize $tag
			rm $inp $cls
		done
	done
done
