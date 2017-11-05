#!/bin/bash
#$ -S /bin/bash
# $ -q fs_long@@qc_nehalem
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG._Get_exonIntronCoverage_fromRNAseq.sh$TASK_ID
# this is for encode RNAseq data (try)
# $ -t 1-70
# this is for coverage computation
#$ -t 1-12
# this is for coverage parsing
# $ -t 1-96
# this is for coverage parsing for encode rnaseq of GIS (x3 * 16)
# $ -t 1-48
# this is for UCSC tracks
# $ -t 1-148

#TODO
#Check all scripts for bedtools intersect for -wa option if the resultant file is unique
exit;

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

outdir=Analysis/RNAseq_exonIntronCoverage
mkdir -p $outdir; pushd $outdir

#TODO set
task=3
indir=~bilebi00/_EWSR1/data/clipz_RNAseq_tags/
ids=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -P "mRNAseq" | cut -f 2`);
libs=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -P "mRNAseq" | cut -f 3`);
cs=(exons introns)

#1.
if [ $task == 1 ]; then
	####mkdir -p $indir; rsync mirz@web07:~/BITER/ChIP/scratch/bed_files/*bed $indir/.
	mkdir -p $indir; 
	rsync -t mirz@web07:~/BITER/RNAseq/scratch/bed_files_hg19/*bed.gz $indir/.
	exit;
fi


#2.
#if [ $task == 2 ] && [ $SGE_TASK_ID == 1 ]; then ##don't submit
if [ $task == 2 ]; then
	sout=regions;
	mkdir -p $sout; pushd $sout;
	#get refseq mRNAs and exons
	less ~bilebi00/_EWSR1/data/hg19_GMAP_mRNAs.gtf | grep "NM_" > hg19_GMAP_mRNAs.refseq.gtf
	less ~bilebi00/_EWSR1/data/hg19_GMAP_EXONS.gtf | grep "NM_" > hg19_GMAP_EXONS.refseq.gtf
	#less ~bilebi00/_EWSR1/data/hg19_GMAP_INTRONS.gtf | grep "NM_" > hg19_GMAP_INTRONS.refseq.gtf

	#this script finds non-overlapping gene loci; if any gene intersect with the loci of another the number of such genes are indicated in the score field of the gtf file
	#
	#set feature as the gene name in mRNA.gtf file
	#use -nms of mergeBed and convert regions to gtf using the number of genes in the cluster as the score of each entry
	less hg19_GMAP_mRNAs.refseq.gtf  | perl -e 'while(<>){$_=~/gene_id \"([^"]*)/; @t=split/\t/; $t[2]=$1; print join ("\t", @t); }' | bedtools merge -s -nms -i stdin | perl -e 'while(<>){ @t=split; @nms=split/;/,$t[3]; %h=(); map { $h{$_}=1; } @nms; $gene_cluster=join(";", sort keys %h); print join ("\t", $t[0], "GMAP", "gene", $t[1]+1, $t[2],scalar keys %h, $t[4], ".", "gene_id \"$gene_cluster\"; transcript_id \"$gene_cluster\";"  ), "\n"; }' > hg19_GMAP_GENE.refseq.gtf

#XXX	
#	# number of genes where sites such that $6>4: 1722
#	#zless ~bilebi00/_ARE/Analysis/Pooled_Reproducibility/EWSR1.gtf.gz | awk '$6>4{ print; }' | annotateBed -s -counts -i _hg19_GMAP_GENE.refseq.gtf -files stdin | perl -e 'while(<>){@t=split/\t/; $t[8] .= " xlink_sites \"$t[9]\";"; splice(@t,9,1); print join("\t", @t); }' > hg19_GMAP_GENE.refseq.gtf
#
#	annotateBed -s -counts -i _hg19_GMAP_GENE.refseq.gtf -files EWSR1.scr_ge$minxlinkscr.gtf | perl -e 'while(<>){@t=split/\t/; $t[8] .= " xlink_sites \"$t[9]\";"; splice(@t,9,1); print join("\t", @t); }' > hg19_GMAP_GENE.refseq.gtf
#
#	rm _hg19_GMAP_GENE.refseq.gtf

	#get non-exon regions as INTRONic regions from GENE.gtf file
	bedtools subtract -s -a hg19_GMAP_GENE.refseq.gtf -b hg19_GMAP_EXONS.refseq.gtf | sed -e 's/\tgene\t/\tintron\t/g' > hg19_GMAP_GENEintrons.refseq.gtf 
#	subtractBed -s -a hg19_GMAP_GENE.refseq.gtf -b hg19_GMAP_EXONS.refseq.gtf | sed -e 's/\tgene\t/\tintron\t/g' | annotateBed -s -counts -i stdin -files EWSR1.scr_ge$minxlinkscr.gtf | perl -e 'while(<>){@t=split/\t/; $t[8] .= " region_xlink_sites \"$t[9]\";"; splice(@t,9,1); print join("\t", @t); }'  > hg19_GMAP_GENEintrons.refseq.gtf
#	#subtractBed -s -a hg19_GMAP_GENE.refseq.gtf -b hg19_GMAP_EXONS.refseq.gtf | sed -e 's/\tgene\t/\tintron\t/g'  > hg19_GMAP_GENEintrons.refseq.gtf
	#get non-intron regions as EXONic regions
	bedtools subtract -s -a hg19_GMAP_GENE.refseq.gtf -b hg19_GMAP_GENEintrons.refseq.gtf | sed -e 's/\tgene\t/\texon\t/g' > hg19_GMAP_GENEexons.refseq.gtf
#	#FIXME These two aproaches are not identical; check
#	subtractBed -s -a hg19_GMAP_GENE.refseq.gtf -b hg19_GMAP_GENEintrons.refseq.gtf | sed -e 's/\tgene\t/\texon\t/g' | annotateBed -s -counts -i stdin -files EWSR1.scr_ge$minxlinkscr.gtf | perl -e 'while(<>){@t=split/\t/; $t[8] .= " region_xlink_sites \"$t[9]\";"; splice(@t,9,1); print join("\t", @t); }' > hg19_GMAP_GENEexons.refseq.gtf
#	##subtractBed -s -a hg19_GMAP_GENE.refseq.gtf -b hg19_GMAP_INTRONS.refseq.gtf | sed -e 's/\tgene\t/\texon\t/g' | sort -k1,1 -k7,7 -k4,5g   > exons

	#linksBed -base http://genome.ucsc.edu -org human -db hg19 -i hg19_GMAP_GENEexons.refseq.gtf > hg19_GMAP_GENEexons.refseq.gtf.htm

#XXX	
#	#get conservation
#	nohup ~bilebi00/_EWSR1/scripts/extractPhastConsFromCoords.pl hg19_GMAP_GENEintrons.refseq.gtf &> intron.log &
#	nohup ~bilebi00/_EWSR1/scripts/extractPhastConsFromCoords.pl hg19_GMAP_GENEexons.refseq.gtf &> exon.log &

#	#rewrite exon and intron gtf file
#	less hg19_GMAP_GENEintrons.refseq.gtf.pc | perl -e 'while(<>){chomp; @t=split/\t/; $t[8] .= " conservation_placentalMammals \"$t[9]\";"; splice(@t,9,1); print join("\t", @t),"\n"; }' > hg19_GMAP_GENEintrons.refseq.gtf
#	less hg19_GMAP_GENEexons.refseq.gtf.pc | perl -e 'while(<>){chomp; @t=split/\t/; $t[8] .= " conservation_placentalMammals \"$t[9]\";"; splice(@t,9,1); print join("\t", @t),"\n"; }' > hg19_GMAP_GENEexons.refseq.gtf
#
#	rm -rf *pc;

#	#TODO other tracks like CTCF polII could be adde here
#	#XXX
#	#encode file doesn't have strand information
#	intersectBed -wao -a hg19_GMAP_GENEintrons.refseq.gtf -b ~bilebi00/DATA/hg19_ucsc_tracks/ENCODE/wgEncodeRegTfbsClustered.bed.gz > intron.w_chipseq

#	#merge to get coverage over
#	cat hg19_GMAP_GENEintrons.refseq.gtf hg19_GMAP_GENEexons.refseq.gtf  | sort -k1,1 -k7,7 -k4,5g  > hg19_GMAP_GENEeis.refseq.gtf
#	~bilebi00/_EWSR1/scripts/sum_length_per_attribute.pl hg19_GMAP_GENEeis.refseq.gtf > hg19_GMAP_GENEeis.id_score_intronLen_exonLen 
	exit;
fi

#3.Get number of read mapping to a region, count the reads that are full length in the region
if [ $task == 3 ]; then
	sout=Erfc;
	mkdir -p $sout; pushd $sout;

	i=`echo $((SGE_TASK_ID-1)) / 2 | bc`;
	id=${ids[$i]}
	file=$indir/$id.bed.gz

	i2=`echo $((SGE_TASK_ID-1)) % 2 | bc`;
	c=${cs[$i2]}
	cfile=../regions/hg19_GMAP_GENE$c.refseq.gtf

	lib=`grep -w $id ~bilebi00/_EWSR1/data/protein_sampleid_list | cut -f 3`;
	echo "Doing $id($lib) $c"
	if [ ! -e $id.${c}_coverage.gz ]; then
		echo `date` coverageBed started
		#zless $file | ~bilebi00/_EWSR1/scripts/get_starts_fromBed.pl | coverageBed -counts -s -a stdin -b $cfile | gzip -c > $id.${c}_coverage.gz
		bedtools intersect -s -a $file -b $cfile -f 1 | bedtools coverage -a - -b $cfile  -counts -s | gzip -c > $id.${c}_coverage.gz
		echo `date` coverageBed finished 
	else
		echo $id.${c}_coverage.gz exists
	fi

	touch $id.$c.coveragestat

	while [ $c == "exons" ] && [ ! -e $id.introns.coveragestat ]; do
		echo waiting for $id.introns.coveragestat file
		sleep 2m;
	done

	if [ $c == "exons" ]; then
		echo Plotting...
		~bilebi00/bin/R --no-save --args $id.introns_coverage.gz $id.exons_coverage.gz count < ~bilebi00/_EWSR1/scripts/plot_intronExpScore_wrtBoth_fromGtfPlus.R > /dev/null
		~bilebi00/bin/R --no-save --args $id.introns_coverage.gz $id.exons_coverage.gz count < ~bilebi00/_EWSR1/scripts/plot_intronExpScore_wrtIntron_fromGtfPlus.R > /dev/null
		rm $id.*.coveragestat
		touch $id.plotstat
	fi

	if [ $SGE_TASK_ID == 1 ]; then
		while [ "`ls *.plotstat | wc -l`" -ne ${#ids[@]} ]; do
			echo waiting for plotstat
			sleep 2m;
		done
		less *.freqstat | sort -r | uniq > all.stat
		~bilebi00/bin/R --no-save --args all.stat < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null
		#clean marks
		rm *freqstat *plotstat

		#library size bias for the Erfc scores
		for i in `seq 0 5`; do 
			id=${ids[$i]}; 
			lib=${libs[$i]}; 
			less $id.introns.txt | awk 'BEGIN{c=0;OFS="\t";}{ if (NR>1) { c=$7+c; print lib, "N", c; } }' lib=$lib | tail -n 1 ; 
		done > N.stat

	fi

	echo DONE
	exit;
fi

#TODO ORGANIZE here after; split the merging scripts
if [ $task == 0 ]; then
	onames="${libs[0]}"
	cp ${ids[0]}.introns.txt all.erfc.introns; 
	for i in `seq 1 $((${#ids[@]}-1))`; do 
		echo $i
		~bilebi00/_PAPD5/scripts/outerJoin.pl all.erfc.introns ${ids[$i]}.introns.txt 5,1,2,3,4 1,2,3,4,5 "" NA >  tmp1; 
		mv tmp1 all.erfc.introns; 
		onames="$onames ${libs[$i]}"
	done
	echo $onames > all.erfc.introns.onames
	~bilebi00/_PAPD5/scripts/outerJoin.pl all.erfc.introns ../xlinkEnrichment/trusted_288_289.introns.xlinkT_mutrate0.05.txt 1,2,3,4,5 5,1,2,3,4 "" NA > all.erfc.introns.wclip; 
	#TODO t-test comparisons are hard coded; generalize
	~bilebi00/bin/R --no-save --args all.erfc.introns erfc.introns ${libs[@]} < ~bilebi00/_EWSR1/scripts/scatter_density_ecdf.R | grep "^Two" > tt.stat

	#clip intersection
	echo "length.cutoff lib method alternative intersection clipped topexpressed total p.value odds.ratio" | sed -e 's/ /\t/g' > intronExpression_clip_intersection.stat
	for id in ${ids[@]}; do
		lib=`grep -w $id ~bilebi00/_EWSR1/data/protein_sampleid_list | cut -f 3`;
		~bilebi00/_PAPD5/scripts/outerJoin.pl $id.introns.txt ../xlinkEnrichment/trusted_288_289.introns.xlinkT_mutrate0.05.txt 1,2,3,4,5 5,1,2,3,4 '' NA > tmp;
		for minlen in 0 100 500 1000 2000 4000 8000; do
			total=`less tmp | awk 'NR>1 && $6 > minlen{print}' minlen=$minlen  | wc -l`
			clipped=`less ../xlinkEnrichment/trusted_288_289.introns.xlinkT_mutrate0.05.txt | awk 'NR > 1 && $3-$2+1 > minlen { print; }' minlen=$minlen | wc -l`;
			topexpr=$clipped
			both=`less tmp | awk 'NR>1 && $6 > minlen {print}' minlen=$minlen | sort -k 10,10gr | head -n $topexpr | grep -v NA$ | wc -l`
			echo -en "$minlen\t$lib\t" 
			~bilebi00/bin/R --no-save --args $both $clipped $topexpr $total < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher" 
		done
		rm tmp
	done >> intronExpression_clip_intersection.stat

	#nuc content
	db=hg19
	genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
	nucBed -pattern GGG -s -fi $genomefa -bed hg19_GMAP_GENEintrons.refseq.gtf > hg19_GMAP_GENEintrons.refseq.gtf.nucbed

	~bilebi00/_PAPD5/scripts/outerJoin.pl hg19_GMAP_GENEintrons.refseq.gtf.nucbed all.erfc.introns.wclip 1,4,5,7 2,3,4,5 '1-19,29,39,49,59,69,79,85' NA | grep "^chr" > hg19_GMAP_GENEintrons.refseq.gtf.nucbed.werfc.wclip

	libs=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -P "mRNAseq" | cut -f 3 | sed 's/$/.erfc/g'`);
	onames="`less *nucbed | head -n 1 | perl -e '$_=<>; $_=~s/#*\d+_usercol|\d+_//g; print; '` ${libs[@]} maxpbeta"
	~bilebi00/bin/R --no-save --args hg19_GMAP_GENEintrons.refseq.gtf.nucbed.werfc.wclip ${onames[@]} < ~bilebi00/_EWSR1/scripts/plot_nucContent_fromGtfPlus.R
	rm hg19_GMAP_GENEintrons.refseq.gtf.nucbed

	#hexamer content
	#XXX

	#ids=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -P "osteo" | cut -f 2`);
	#libs=(`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -P "osteo" | cut -f 3`);
	#cp ${ids[0]}.introns.txt all.erfc.introns.osteo; 
	#for i in `seq 1 $((${#ids[@]}-1))`; do 
	#	echo $i
	#	~bilebi00/_PAPD5/scripts/outerJoin.pl all.erfc.introns.osteo ${ids[$i]}.introns.txt 1,2,3,4,5,6 1,2,3,4,5,6 "" NA >  tmp1; 
	#	mv tmp1 all.erfc.introns.osteo; 
	#done
	#~bilebi00/bin/R --no-save --args all.erfc.introns.osteo erfc.introns ${libs[@]} < ~bilebi00/_EWSR1/scripts/scatter_density_ecdf.R | grep "^Two" > tt.osteo.stat


fi

##TODO ask above
if [ $task == 33 ]; then
	ddir=~bilebi00/DATA/hg19_ucsc_tracks/ENCODE/
	samplesfile=$ddir/RnaSeq_tracks.4download
	ids=(`less $samplesfile | cut -f 2 | sed 's/.bam//g'`);

	i=`echo $((SGE_TASK_ID-1)) / 2 | bc`;
	id=${ids[$i]}
	iddir=`less $samplesfile | grep -w $id | cut -f 1`;
	file=$ddir/$iddir/$id.bam

	i2=`echo $((SGE_TASK_ID-1)) % 2 | bc`;
	c=${cs[$i2]}
	cfile=hg19_GMAP_GENE$c.refseq.gtf

	if [ -e $id.${c}_coverage.gz ]; then
		echo $id.${c}_coverage.gz exists.
		exit;
	fi

	echo "Doing $id(in $iddir) $c"
	echo `date` coverageBed started
	coverageBed -abam -split -d -s -a $file -b $cfile | gzip -c > $id.${c}_coverage.gz
	echo `date` coverageBed finished 

fi

#4.plot region level
#XXX
if [ $task == 44 ] && [ $SGE_TASK_ID == 1 ]; then
	pwd=`pwd`;
	ids=`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep mRNAseq | cut -f 2`; 
	rm -rf topExpressedIntrons_clipIntersection all.intronsExp;
	for id in $ids; do 
		lib=`less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -w $id | cut -f 3`; 
		echo Doing $lib
		~bilebi00/bin/R --no-save --args $id.introns.gtf $id.exons.gtf < ~bilebi00/_EWSR1/scripts/plot_region_fromGtf.R > /dev/null

		#merge files
		if [ ! -e all.intronsExp ]; then
			cat $id.introns.gtf.all.txt > all.intronsExp	
		else
			awk 'NR>1{ print;}' $id.introns.gtf.all.txt >> all.intronsExp	
		fi

		#extract stats
		echo lib=$lib id=$id \
			clipped=`less $id.introns.gtf.all.txt  | head -n 1000 | awk '$7>0{print;}' | wc -l`"/1000" >> topExpressedIntrons_clipIntersection 

		#put them in the web
		pushd ~bilebi00/www/EWSR1;
		ln -sf $pwd/$id.introns.gtf.all.txt $lib.intronsExp 
		popd
	done

	#XXX
	~bilebi00/bin/R --no-save --args all.intronsExp < plot_regionECDF_for_mult_samples.R | grep -B 1 "^Two" > all.intronsExp.stats
	exit;
fi


#5. TF intersection with whole bed file
#XXX automate the number of jobs now it is 148
#Checks if ABx > A'Bx where Bx is all the TFs
if [ $task == 4 ]; then
	jobscount=148
	odir=../wgEncodeRegTfbsClustered/
	mkdir -p $odir; pushd $odir;
	tffile=~bilebi00/DATA/hg19_ucsc_tracks/ENCODE/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClustered.bed.gz
	rfile=hg19_GMAP_GENEintrons.refseq.clip.gtf

	if [ $SGE_TASK_ID == 1 ] && [ ! -e count.stat ]; then
		zless $tffile |	cut -f 4 | sort | uniq -c | awk 'BEGIN{OFS="\t"}{ print $1, $2;}' | sort -k 1,1gr > TFcluster.count
		less ../xlinkEnrichment/trusted_288_289.introns.xlinkT_mutrate0.05.txt | awk 'BEGIN{OFS="\t"} { if (NR>1) { print $1,$2,$3,$5,0,$4;} }' | \
			intersectBed -a ../RNAseq_exonIntronCoverage/hg19_GMAP_GENEintrons.refseq.gtf -b stdin -wao -s | \
			awk -F "\t" 'BEGIN{OFS="\t";} {if ($10 ~ /^chr/) { status=1; } else  { status=0; } $9 = $9" CLIP \""status"\";"; print; }' | cut -f 1-9 > $rfile
		touch count.stat
	else 
		while [ ! -e count.stat ]; do
			echo waiting for count.stat file
			sleep 1m;
		done
	fi

	lib=`less TFcluster.count | head -n $SGE_TASK_ID | tail -n 1 | awk '{ print $2;}'`
	of=$SGE_TASK_ID.TfbsClustered
	zless $tffile | awk '$4==lib{ print; }' lib="$lib" > $SGE_TASK_ID.bed
	#no strand information available in ChIPseq bed tracks
	echo `date` intersectBed started
	intersectBed -a $rfile -b $SGE_TASK_ID.bed  -wao | \
		awk -F "\t" 'BEGIN{OFS="\t";} { print $0,lib,$14*$25; }' lib=$lib | \
		cut -f 1-9,26-27 | ~bilebi00/_EWSR1/scripts/parseIntensity_perlGene_wclip.pl > $of
	echo `date` intersectBed finished

	plot=0
	if [ "`echo $lib | grep -P 'CTCF|Pol'`" != "" ]; then
		plot=1;
	fi

	#clean
	rm $SGE_TASK_ID.bed

	~bilebi00/bin/R --no-save --args $of $lib $plot < ~bilebi00/_EWSR1/scripts/plot_intensity_4ChIPclusters.R | grep "^Two" > $lib.tteststat
	touch $SGE_TASK_ID.plotstat

	if [ $SGE_TASK_ID == 1 ]; then
		while [ "`ls *plotstat | wc -l`" != $jobscount ]; do 
			echo waiting for $jobscount plotstat files
			sleep 1m;
		done
		cat *.tteststat | sort -k 5,5g -t $'\t' > TFintensity.stat
		cat *.TfbsClustered > TFintensity.txt
		gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=TFintensity.pdf *.pdf 
		rm *.plotstat *.tteststat *.TfbsClustered 
	fi
	echo DONE
fi

#TODO set
#RNAseq tracks; miRNA type expression profiling
if [ $task == 5 ]; then
	libindex=3
	ddir=~bilebi00/DATA/hg19_ucsc_tracks/ENCODE/
	samplesfile=$ddir/RnaSeq_tracks.4download; libindex=2
	ids=(`less $samplesfile | grep "lab=GIS;" | cut -f 2 | sed 's/.bam//g'`);

	woTermRegs=(0 1)
	wlog2s=(0 1)
	woZeroCovs=(0 1)

	#TODO generalize the divisors it is the size of the product of following arrays
	i=`echo $((SGE_TASK_ID-1)) / 16 % ${#ids[@]} | bc`;
	id=${ids[$i]}

	i=`echo $((SGE_TASK_ID-1)) / 8 % ${#cs[@]} | bc`;
	c=${cs[$i]}

	i=`echo $((SGE_TASK_ID-1)) / 4 % ${#wlog2s[@]} | bc`;
	wlog2=${wlog2s[$i]}

	i=`echo $((SGE_TASK_ID-1)) / 2 % ${#woZeroCovs[@]} | bc`;
	woZeroCov=${woZeroCovs[$i]}

	i=`echo $((SGE_TASK_ID-1)) / 1 % ${#woTermRegs[@]} | bc`;
	woTermReg=${woTermRegs[$i]}

	tag=${id}.${c}_wlog2_${wlog2}.woZeroCov_${woZeroCov}.woTermReg_${woTermReg}
	if [ $woTermReg == 1 ]; then
		echo woTermRegs=1 is not implemented yet: $tag
		exit;
	fi

	lib=`grep -w $id $samplesfile | cut -f $libindex | sed 's/.bam//g'`;
	infile=$id.${c}_coverage.gz

	echo `date` `uname -a` Copying $tag.gtf
	if [ -e $tag.gtf.gz ]; then
		echo $tag.gtf.gz exists.
		exit;
	fi

	pwd=`pwd`;
	scratchdir=/scratch/bilebi00/$tag
	mkdir -p $scratchdir; 
	rsync $infile $scratchdir/.
	pushd $scratchdir;
	echo `date` Preparing $tag.gtf
	#TODO ~bilebi00/_EWSR1/scripts/get_coverageStats_4geneRegion_from_gtfPlus.pl woTermReg=1 option implementation
	~bilebi00/_EWSR1/scripts/get_coverageStats_4geneRegion_from_gtfPlus.pl $infile $lib $wlog2 $woZeroCov $woTermReg $c | gzip -c > $tag.gtf.gz
	popd
	mv $scratchdir/$tag.gtf.gz .
	rm -rf $scratchdir;
	echo `date` DONE
	exit;
fi


#6.plot gene level
if [ $task == 6 ] && [ $SGE_TASK_ID == 1 ]; then
	rand=$RANDOM
	echo "geneId xlinkSites locusGeneCount intronLength exonLength name intronRNAseqCoverage exonRNAseqCoverage" | sed 's/ /\t/g' > final.plot
	#for id in `less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -P "GROseq|RNAseq" | cut -f 2`; do
	for id in `less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -P "mRNAseq" | cut -f 2`; do
		echo $id
		lib=`grep -w $id ~bilebi00/_EWSR1/data/protein_sampleid_list | cut -f 3`;
		~bilebi00/_PAPD5/scripts/outerJoin.pl hg19_GMAP_GENEeis.id_score_intronLen_exonLen $id.introns_coverage.sum_per_attr 1 1 '' 0 | awk 'BEGIN{OFS="\t";}{print $2,$6,$7,$8,$9,lib,$16,$17; }' lib=$lib > $id.introns
		~bilebi00/_PAPD5/scripts/outerJoin.pl hg19_GMAP_GENEeis.id_score_intronLen_exonLen $id.exons_coverage.sum_per_attr 1 1 '' 0 | awk 'BEGIN{OFS="\t";}{print $16,$2; }' > $id.exons
		~bilebi00/_PAPD5/scripts/outerJoin.pl $id.introns $id.exons 1 2 '' 0 | cut -f 1-8 | sed -e 's/[";]//g' >> final.plot
		rm $id.introns $id.exons
	done
	~bilebi00/bin/R --no-save --args final.plot < ~bilebi00/_EWSR1/scripts/plot_coverage.R > /dev/null

	exit;
fi

#7.plot region level -absolute now
if [ $task == 7 ] && [ $SGE_TASK_ID == 1 ]; then
	h="geneId xlinkSites locusGeneCount name regionXlinkSites regionLocation regionCoverage regionLength regionPMphastCons pexonLocation pexonCoverage pexonLength pexonPMphastCons fexonLocation fexonCoverage fexonLength fexonPMphastCons"
	rand=$RANDOM
	echo $h | sed 's/ /\t/g' > regionFinal.plot
	#for id in `less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -P "GROseq|RNAseq" | cut -f 2`; do
	for id in `less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -P "mRNAseq" | cut -f 2`; do
		echo $id
		lib=`grep -w $id ~bilebi00/_EWSR1/data/protein_sampleid_list | cut -f 3`;
		#*sum_per_attrRegion contain the expressed exons and introns;
		#each intron is 
		#FIXME XXX some introns are skipped because they are consequitive 
		#example SAMD11; check this XXX
		cat $id.*ons_coverage.sum_per_attrRegion | sort -k 1,1 -k4,5g | \
			~bilebi00/_EWSR1/scripts/add_neighbouringExonCovLen2Introns.pl $lib >> regionFinal.plot
	done
	#	~bilebi00/bin/R --no-save --args regionFinal.plot < ~bilebi00/_EWSR1/scripts/plot_region_coverage.R | \
		#		grep "Two sample KS test for CLIPPEDREGION and NOTCLIPPEDREGION RNAseq expression" > KS.stats
	#XXX
	~bilebi00/bin/R --no-save --args regionFinal.plot < plot_region_coverage.R | \
		grep "Two sample KS test for CLIPPEDREGION and NOTCLIPPEDREGION RNAseq expression" > KS.PM.stats

	#	#sRNA track	
	#	#less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -P "GROseq|RNAseq" | cut -f 3 | \
		#		while read i; do  \
		#			echo $i; less regionFinal.plot | awk '$4==lib{ print;}' lib=$i > $i.plot;  \
		#		done

	#	echo "$h $h" | sed 's/ /\t/g' > sRNAseq.joined

	#	#less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -P "GROseq|mRNAseq" | cut -f 3 | \
		#	less ~bilebi00/_EWSR1/data/protein_sampleid_list | grep -P "mRNAseq" | cut -f 3 | \
		#		while read i; do  \
		#			~bilebi00/_PAPD5/scripts/leftJoin.pl sRNAseq.plot $i.plot 1,6 1,6 '' NA; \
		#  done | awk '$NF != "NA"{ print; }' >> sRNAseq.joined
	#
	#	~bilebi00/bin/R --no-save --args sRNAseq.joined < ~bilebi00/_EWSR1/scripts/plot_region_coverage_4sRNAseq.R | \
		#		grep "\"Two sample KS test for CLIPPEDREGION and NOTCLIPPEDREGION RNAseq expression\"" > sRNAseq_KS.stats
	#
	exit;
fi


