#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG._Get_GO_4clipzRNAseq.sh$TASK_ID
#$ -t 1-4

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

export PATH=$HOME/bin:/import/bc2/soft/bin/graphviz:$PATH

outdir=Analysis/GO_4clipzRNAseqFCup_down/
mkdir -p $outdir; pushd $outdir;
ps=1

if [ ! -e population ]; then
	if [ $SGE_TASK_ID == 1 ]; then
		ln -sf ~bilebi00/_EWSR1/data/clipz_RNAseq/RNAseq.raw.qn.w_geneSymbol .
		less RNAseq.raw.qn.w_geneSymbol | awk 'NR>1{ print $6; }' > population
		less RNAseq.raw.qn.w_geneSymbol  | awk 'BEGIN{OFS="\t";}{ if (NR==1) { print $1,$6,$4."_over_"$2,$5."_over_"$3; } else{ print $1, $6, log($4+ps)/log(2)-log($2+ps)/log(2), log($5+ps)/log(2)-log($3+ps)/log(2); } }' ps=$ps > RNAseq.raw.qn.w_geneSymbol.FC

		#GO files
		ln -sf ~bilebi00/_SHORT_READ_ALIGNMENT/Stepanka/SUBSET/GO/gene_association.goa_human .
		ln -sf ~bilebi00/_SHORT_READ_ALIGNMENT/Stepanka/SUBSET/GO/gene_ontology_edit.obo .
		ln -sf ~bilebi00/_SHORT_READ_ALIGNMENT/Stepanka/SUBSET/GO/Ontologizer.jar .

	else 
		while [ ! -e Ontologizer.jar ]; do 
			sleep 50;
		done
	fi
fi
#R --no-save < a.R > /dev/null

dirs=(batch1_UP batch1_DOWN batch2_UP batch2_DOWN);
dir=${dirs[$((SGE_TASK_ID-1))]}
mkdir -p $dir

echo Doing folder $dir

if [ $dir == batch1_UP ]; then
	less RNAseq.raw.qn.w_geneSymbol.FC | awk 'BEGIN{OFS="\t";} NR>1{ print$2,  $3; }' | sort -k 2,2gr | head -n 200 | cut -f 1 > $dir/studyset
elif [ $dir == batch1_DOWN ]; then
	less RNAseq.raw.qn.w_geneSymbol.FC | awk 'BEGIN{OFS="\t";} NR>1{ print$2,  $3; }' | sort -k 2,2g | head -n 200 | cut -f 1 > $dir/studyset
elif [ $dir == batch2_UP ]; then
	less RNAseq.raw.qn.w_geneSymbol.FC | awk 'BEGIN{OFS="\t";} NR>1{ print $2,  $4; }' | sort -k 2,2gr | head -n 200 | cut -f 1 > $dir/studyset
elif [ $dir == batch2_DOWN ]; then
	less RNAseq.raw.qn.w_geneSymbol.FC | awk 'BEGIN{OFS="\t";} NR>1{ print $2,  $4; }' | sort -k 2,2g | head -n 200 | cut -f 1 > $dir/studyset
fi

java -jar Ontologizer.jar -n -a gene_association.goa_human -g gene_ontology_edit.obo -s $dir/studyset -p population -c Parent-Child-Union -m Westfall-Young-Single-Step -d 0.1 -r 1000 -o $dir
dot -Tpng $dir/*dot -o $dir/GO_graph.png
