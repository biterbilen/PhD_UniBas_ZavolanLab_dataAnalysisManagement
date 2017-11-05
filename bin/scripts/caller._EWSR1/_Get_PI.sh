#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG._Get_PI.sh$TASK_ID
#$ -t 1

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

export PATH=$HOME/bin:$PATH

outdir=Analysis/PI
mkdir -p $outdir; pushd $outdir

p=EWSR1

#TODO set
task=2

#1.
if [ $task == 1 ]; then
	#extract human 
	PIf="PI.human"
	less ~bilebi00/DATA/InteractingProteinDBs/BioGRID/BIOGRID-ALL-3.1.83.tab.txt | grep -w "9606$" > $PIf

	#get direct interactors exclude multiexperiments
	refs="`less $PIf | awk -F $'\t' '$3==p||$4==p {print $0;}' p=$p | cut -f 8 | sort | uniq -c  | perl -e 'while(<>){ $_=~/(\d+)\s+(.+)$/; print $2,"\n" if ($1>10); }'`";
	refs="`echo $refs | sed 's/) /)|/g' | sed 's/)/\\\)/g' | sed 's/(/\\\(/g'`";
	  
	pps=`less $PIf | awk -F $'\t' '$3==p || $4==p { print $0; }' p=$p | grep -v -P "$refs" | awk '{ print $3"\n"$4;}' | grep -v "$p$" | sort | uniq`;	

	for i in $pps; do echo -e "$p\t$i"; done > $p.allPI_singleExp
	#get indirect interactors and exclude multiexperiments
	for pp in $pps; do
		refs="`less $PIf | awk -F $'\t' '$3==p||$4==p {print $0;}' p=$pp | cut -f 8 | sort | uniq -c | perl -e 'while(<>){ $_=~/(\d+)\s+(.+)$/; print $2,"\n" if ($1>10); }'`";
		refs="`echo $refs | sed 's/) /)|/g' | sed 's/)/\\\)/g' | sed 's/(/\\\(/g'`";
		less $PIf | awk -F $'\t' '$3==p || $4==p { print $0; }' p=$pp | grep -v -P "$refs" | awk '{ print $3"\n"$4;}' | grep -v "$pp$" | sort | uniq | \
			while read i; do echo -e "$pp\t$i"; done >> $p.allPI_singleExp;	
	done;
	
fi

#GO w Ontologizer
if [ $task == 2 ]; then

	#for dot
	export PATH=$HOME/bin:/import/bc2/soft/bin/graphviz:$PATH

	dir='GO'
	mkdir -p $dir; pushd $dir

	#GO files
	ln -sf ~bilebi00/_SHORT_READ_ALIGNMENT/Stepanka/SUBSET/GO/gene_association.goa_human .
	ln -sf ~bilebi00/_SHORT_READ_ALIGNMENT/Stepanka/SUBSET/GO/gene_ontology_edit.obo .
	ln -sf ~bilebi00/_SHORT_READ_ALIGNMENT/Stepanka/SUBSET/GO/Ontologizer.jar .

	#population
	ln -sf ~bilebi00/_EWSR1/data/clipz_RNAseq/RNAseq.raw.qn.w_geneSymbol .
	less RNAseq.raw.qn.w_geneSymbol | awk 'NR>1{ print $6; }' > population

	#XXX if more than one GO analysis is required 
	#studyset
	sdir=wo_UBC
	mkdir -p $sdir
	less ../$p.allPI_singleExp | grep -v "^UBC" | cut -f 2 | sort | uniq > $sdir/studyset.all

	#select expressed 
	~bilebi00/_PAPD5/scripts/innerJoin.pl $sdir/studyset.all population 1 1 1 > $sdir/studyset
	
	java -jar Ontologizer.jar -n -a gene_association.goa_human -g gene_ontology_edit.obo -s $sdir/studyset -p population -c Parent-Child-Union -m Westfall-Young-Single-Step -d 0.1 -r 1000 -o $sdir
	dot -Tpng $sdir/*dot -o $sdir/GO_graph.png
fi

#FC
if [ $task == 3 ]; then
	inp=~bilebi00/_EWSR1/Analysis/GO_4clipzRNAseqFCup_down/RNAseq.raw.qn.w_geneSymbol.FC
	ln -sf $inp .
	inp=`basename $inp`;

	#Add PI information to $inp
	ln -sf GO/wo_UBC/studyset . 
	less studyset | perl -e 'print "PIstat\tgid\n"; while(<>){ print "PI\t$_"; }' > PI 
	~bilebi00/_PAPD5/scripts/outerJoin.pl $inp PI 2 2 '1-5' NOTPI > $inp.PI


	#add TF information to INP
	ln -sf ~bilebi00/_EWSR1/data/hg19/gene_info .
	ln -sf /import/bc2/home/nimwegen/GROUP/WMs/Mammals/mat_TF_associations.hg .
	ln -sf ../MARA_combined_from_aj/data_analysis/mergedZvals .

	#change inconsistent names of TFs; make it consistent w clipz gene_info
	less mat_TF_associations.hg | perl -e 'while(<>){ @t=split/\t+/; for(1..$#t) { ($x,$y,$name)=(split /:/,$t[$_]); print "$t[0]\t$name\n" if ($name ne "");  } }' > splitTFs
	missingnames=`~bilebi00/_PAPD5/scripts/leftJoin.pl splitTFs gene_info 2 3 '2-3' NULL | grep NULL | cut -f 1`
	for i in $missingnames; do newname=`grep -w $i gene_info -m 1 | cut -f 3`; echo Changing $i 2 $newname; sed -i "s/\t$i/\t$newname/g" splitTFs ; done

	#add z-scores
	#ask piotr that NR4A2.p2 is missing in the mergedZvals
	echo "TF gid zvalue" | sed 's/ /\t/g' > splitTFs.w_zvalue
	~bilebi00/_PAPD5/scripts/leftJoin.pl splitTFs mergedZvals 1 1 '1,2,4' NULL | grep -v -w NULL >> splitTFs.w_zvalue
	
	~bilebi00/_PAPD5/scripts/leftJoin.pl $inp.PI splitTFs.w_zvalue 2 2 '' NA > $inp.PI.TF

	~bilebi00/bin/R --no-save --args $inp.PI.TF < ~bilebi00/_EWSR1/scripts/plot_PI.R | grep "^Two sample KS test" > KS.stat

fi
