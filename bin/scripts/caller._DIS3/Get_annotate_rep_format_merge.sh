#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o log_file

export PERL5LIB=$HOME/lib:$PERL5LIB

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

outdir=output_directory
tn=repf
gif=gif
filesc=file_count

pushd $outdir

null=NA
annot=$tn.annot

~bilebi00/_DIS3/scripts/select_annotation.pl $filesc $tn.annot > final.$tn.annot
~bilebi00/_DIS3/scripts/assign_gene_id.pl final.$tn.annot $gif 2 NA | perl -e '$tn=shift; while(<>){ ($id, @t)=split; (@a)= ($id =~ /${tn}_(\S+)([+\-])(\d+)\-(\d+)/); print join("\t", $1,$2,$3,$4,@t),"\n"; }' > final2.$tn.annot
~bilebi00/_PAPD5/scripts/outerJoin.pl final2.$tn.annot $tn '1,2,3,4' '1,2,4,5' '' $null | cut -f 5- > $annot
less -S $annot | sort -k 9,9gr > $annot.sorted

rm *\.$tn.annot final*.$tn.annot

source=REP
less $annot.sorted | perl -e '$src=shift; $name=shift; $_=<>; while(<>){@t=split; next if ($t[11] eq "NULL"); my (@s)=split/,/,$t[12]; print join("\t", $t[3],$src,$name,$t[6],$t[7],$t[8],$t[9],$t[4],".",".",".",".",".",$t[0],$t[1],$t[2]),"\n"; }' $source $tn > $tn.data

echo DONE
