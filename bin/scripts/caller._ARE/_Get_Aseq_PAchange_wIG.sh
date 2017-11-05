#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-8
#$ -o LOG._Get_Aseq_PAchange_wIG.sh.$TASK_ID

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

task=3 
indir=~bilebi00/_ARE/data/clipz_Aseq_sites/

#TODO set
#HeLa281
tag=HeLa281
ids=(`less ~bilebi00/_ARE/data/protein_sampleid_list | grep -P "HeLa281" | cut -f 2`);
libs=(`less ~bilebi00/_ARE/data/protein_sampleid_list | grep -P "HeLa281" | cut -f 3`);
#hnRNPQ
tag=hnRNPQold
ids=(`less ~bilebi00/_ARE/data/protein_sampleid_list | grep -P "HEK293wt_253|siCTRL_253|sihnRNPQ_235" | cut -f 2`);
libs=(`less ~bilebi00/_ARE/data/protein_sampleid_list | grep -P "HEK293wt_253|siCTRL_253|sihnRNPQ_235" | cut -f 3`);
#Andreas
tag=253
ids=(`less ~bilebi00/_ARE/data/protein_sampleid_list | grep -P -w "HEK293wt_253|siCTRL_253|siCstf64_253|siCFIm68_253" | cut -f 2`);
libs=(`less ~bilebi00/_ARE/data/protein_sampleid_list | grep -P -w "HEK293wt_253|siCTRL_253|siCstf64_253|siCFIm68_253" | cut -f 3`);
tag=HEK303
ids=(`less ~bilebi00/_ARE/data/protein_sampleid_list | grep -P HEK303 | cut -f 2`);
libs=(`less ~bilebi00/_ARE/data/protein_sampleid_list | grep -P HEK303 | cut -f 3`);
samplecount=${#ids[@]}

outdir=Analysis/PAchange_wIG/$tag/
mkdir -p $outdir; pushd $outdir

db=hg19
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
genomelen=~bilebi00/aux/human.$db.genome

#1.
if [ $task == 1 ]; then
	mkdir -p $indir;
	for id in ${ids[@]}; do
		echo $id
		rsync -t mirz@web08:~/clipzServer/data/samples/$id/aseq.sites $indir/$id.sites
	done
	exit;
fi


if [ $task == 2 ] && [ $SGE_TASK_ID == 1 ]; then
	#adapted from /import/bc2/home/zavolan/grubera/HUMAN-paper/HOWTO.sh
#	i=`echo $((SGE_TASK_ID-1)) / 2 | bc`;
#	id=${ids[$i]}

	# get rid of chrM stuff
	joinparam=''
	for id in ${ids[@]}; do
		grep -v chrM $indir/$id.sites > $id.chrMfiltered
		joinparam="$joinparam --sample=$id.chrMfiltered"
	done

	# auxiliary data generation
	# GFF3 file from the clipZ genome browser
	##perl /import/bc2/home/zavolan/grubera/MOUSE-paper/ag-process-gff-annotation-file.pl --gff=hg19.webserver.gff
	###perl /import/bc2/home/zavolan/grubera/HUMAN-paper/ag-gene2accession-map.pl > hg19.webserver.gff.genes.fromNCBI

	# pool samples into one file
	p='/import/bc2/home/zavolan/grubera/CentralProgs/ag-pool-3p-end-sites-from-multiple-samples.pl';
	$p $joinparam --skip_ip > hg19.sites

	# generating clusters
	# first we call with a very low output value and determine the cutoff with by a FDR of 10%
	p='/import/bc2/home/zavolan/grubera/CentralProgs/ag-single-linkage-clustering-of-sites.pl';
	$p --file=hg19.sites --samples=$samplecount --min_norm_tags_for_output=0.01 > hg19.clusters
	perl /import/bc2/home/zavolan/grubera/MOUSE-paper/ag-filter-clusters.pl --samples=$samplecount --file=hg19.clusters &> fdr

	echo DONE
	exit; 
fi

if [ $task == 3 ] && [ $SGE_TASK_ID == 1 ]; then
	#TODO set
	min_norm_tags=5.50 #poster

	#61.45 0.149607922395062 4526 #HeLa281_CTRLkd253
	min_norm_tags=62.45

	#52.05 0.219919927593496 3861 #HeLa281
	min_norm_tags=52.05

	#3.10 0.0983220506744855 24197 #hnRNPQold
	min_norm_tags=3.10

	#1.90 0.100955814709363 31906 #Andreas
	min_norm_tags=1.90

	#13.60 0.100238917584954 15330 #HEK303
	min_norm_tags=13.60

	echo min_norm_tags: $min_norm_tags

	p='/import/bc2/home/zavolan/grubera/CentralProgs/ag-single-linkage-clustering-of-sites.pl';
	$p --file=hg19.sites --samples=$samplecount --min_norm_tags_for_output=$min_norm_tags > hg19.clusters 2> singleLinkage

	# annotate transcripts
	p='/import/bc2/home/zavolan/grubera/CentralProgs/ag-annotate-transcripts.pl';
	exonsf=/import/bc2/home/zavolan/grubera/HUMAN-paper/hg19.webserver.gff.coding.transcripts.exons
	$p --file=hg19.clusters --exons=$exonsf > hg19.clusters.transcripts

	# annotate genes
	p='/import/bc2/home/zavolan/grubera/CentralProgs/ag-annotate-genes.pl';
	genesf=/import/bc2/home/zavolan/grubera/HUMAN-paper/hg19.webserver.gff.genes
	$p --file=hg19.clusters.transcripts --genes=$genesf > hg19.clusters.transcripts.genes

	# split into subsets
	/import/bc2/home/zavolan/grubera/MOUSE-paper/ag-split-genes-by-number-of-sites.pl --file=hg19.clusters.transcripts.genes --prefix=hg19.clusters.transcripts.genes --splicingindependent 2> PA.stats

	p='/import/bc2/home/zavolan/grubera/HUMAN-paper/ag-calculate-hexamer-scores.pl'
	# calculate hexamer scores
#	perl $p --estimate=hg19.clusters.transcripts.genes.1.sites --file=hg19.clusters.transcripts.genes.1.sites --file=hg19.clusters.transcripts.genes.2.sites --file=hg19.clusters.transcripts.genes.3.sites --file=hg19.clusters.transcripts.genes.4.sites > hg19.clusters.hexamerscore
#XXX bug
#Use of uninitialized value in multiplication (*) at /import/bc2/home/zavolan/grubera/HUMAN-paper/ag-calculate-hexamer-scores.pl line 209.

	rm -rf nucleotide_distance
	libss=`echo ${ids[@]} | sed 's/ /,/g'`;
	rm -rf clean_PA.stats;
	less PA.stats | grep "^Genes" | awk '{ print $3; }' | while read i ; do
		N=`less PA.stats | grep "Genes with $i sites:" | awk '{ print $5; }'`;
		if [ $N -lt 200 ]; then
			break;
		fi
		echo $i

		pa=hg19.clusters.transcripts.genes.$i.sites

		#---------------------
		ti=`echo  "10 + 1 + $((2*samplecount)) + 2" | bc`
		#form gtf file from the sites
		less $pa | cut -f 1-9,$ti > $pa.tmp
		~bilebi00/_PAPD5/scripts/leftJoin.pl $pa.tmp ~bilebi00/_EWSR1/data/hg19/gene_info 1 2 '' NULL | grep -w -v NULL | cut -f 1-10,13 | \
			awk 'BEGIN{OFS="\t";}{ print $5, src, fea, $9, $9, scr, $6, ".", "gene_id \""$1"\"; transcript_id \""$10"\"; gene_name \""$11"\";"; }' src=Aseq fea=cluster scr=0 \
			| ~bilebi00/_ARE/scripts/add_PAcode_2gtf.pl $i | gzip -c > $pa.gtf.gz
		~bilebi00/_PAPD5/scripts/leftJoin.pl $pa.tmp ~bilebi00/_EWSR1/data/hg19/gene_info 1 2 '' NULL | \
			grep -w -v NULL | cut -f 1-$((10 + samplecount)),$((9 + samplecount)),$((9 + 1 + 3 + samplecount)) | \
			awk 'BEGIN{OFS="\t";}{ print $5,src,fea,$9,$9,scr,$6,".","gene_id \""$1"\"; transcript_id \""$(NF-1)"\"; gene_name \""$NF"\";",$0; }' src=Aseq fea=cluster scr=0 | \
			cut -f 1-9,$((9 + 9 + 1 ))-$((9 + 9 + samplecount)) | ~bilebi00/_ARE/scripts/add_PAcode_2gtf.pl $i | gzip -c > $pa.gtfplus.gz
		ids="`less hg19.sites | head -n 100 | grep "^#" | grep "TPM" | perl -e 'while(<>){ $_=~/(\d+)\.chrMfiltered/; $id.=($1.","); } chop $id; print $id,"\n";'`";
		echo "lib.$ids" | sed 's/,/\tlib./g' > $pa.onames
		R --no-save --args $pa.gtfplus.gz $pa.onames <  ~bilebi00/_ARE/scripts/plot_PA_TPM_CI.R > /dev/null
		rm $pa.tmp

		#---------------------
		#get sites that are $d base pairs apart to reduce downstream xlink sites
		d=500;
		bedtools closest -t first -io -D "a" -d -s -a $pa.gtf.gz -b $pa.gtf.gz | \
			awk -F "\t" 'BEGIN{OFS="\t";}{ if ($19<0) { $19=-$19; } if ($19>d) {  print $0;}  }' d=$d | \
			awk '{ print $14; }' | sort | uniq -c | awk '$1==i { print $2; }' i=$i | \
		 	while read i; do zless $pa.gtf.gz | grep -w $i ; done > \
			clean_$pa.gtf;
		N="`zless clean_$pa.gtf | wc -l`";
		perl -e '$i=shift; $N=shift; $n=$N/$i; print "Genes with $i sites: $n\n"; ' $i $N >> clean_PA.stats
		if [ $N -lt 200 ]; then
			rm clean_$pa.gtf
		else
			gzip clean_$pa.gtf;
		fi

		#---------------------
		zless $pa.gtf.gz | bedtools slop -s -l $d -r $d -i stdin -g $genomelen | \
		bedtools nuc "-fi" $genomefa -bed stdin -s -seq | cut -f 1-9,19 | awk 'NR>1{ print; }' | \
			  gzip -c > $pa.gtf.nucbed.gz

		#---------------------
		otag=_$pa.gtf
		nucs="A C G T"
		~bilebi00/_EWSR1/scripts/split_gtf_nucbed.pl $pa.gtf.nucbed.gz "$nucs" $otag 0
		for nuc in $nucs; do
			gzip -f $nuc$otag;
		done

		#get nucleotide profile around PA
		for nuc in $nucs; do 
			zless $pa.gtf.gz | bedtools slop -s -l $d -r $d -i stdin -g $genomelen | \
			bedtools intersect -wao -s -a stdin -b ${nuc}_$pa.gtf.gz | \
			perl -e '$slopl=shift; $slopr=shift; $tag=shift; while(<>){@t=split/\t/; $x=($t[3]+$t[4]+$slopl-$slopr)/2; $d=int($t[12]-$x); if ($t[6] eq "-") { $x=($t[3]+$t[4]+$slopr-$slopl)/2; $d=int($x-$t[12]); } $t[8] =~ /PAcode "([^"]*)"/; print join("\t",$d, $tag.".".$1),"\n"; }' $d $d $nuc.site | \
			sort -k 1,2 | bedtools groupby -g 1,2 -c 1 -o count -i stdin -full; 
		done >> nucleotide_distance
		#---------------------
		#generic column paramaters depending on samplecount
		col1start=11
		col1end=$((col1start + samplecount - 1))
		#col1end=$((col1start + $(($(($((samplecount*samplecount))-samplecount))/2)) ))
		col0indices=""
		for j in `seq 1 $samplecount`; do
			ii=$((j+9))
			col0indices="$col0indices,$ii"
		done
		col0indices=`echo $col0indices | sed 's/^,//'`
	
		#calculate IG score for sites
		~bilebi00/_ARE/scripts/calc_informationGain.pl $pa $col0indices $libss 1 > IG_density.$i.sites
		~bilebi00/bin/R --no-save --args $pa $col1start $col1end IG_density.$i.sites ${ids[@]} < ~bilebi00/_ARE/scripts/plot_IG.R > /dev/null

	done
	R --no-save --args nucleotide_distance < ~bilebi00/_ARE/scripts/plot_profile_4Groups_4Nucleotide.R > /dev/null

	echo DONE
	exit;
fi

if [ $task == 4 ]; then

	#XXX set
	#type="";
	type=clean_

	#get crosslink filename
	sid_prot_f=~bilebi00/_ARE/data/protein_sampleid_list
	tags=(`less $sid_prot_f | grep -w PARCLIP | cut -f 1 | sort | uniq`)
	tag=${tags[$((SGE_TASK_ID-1))]}
	ptag=`less $sid_prot_f | grep -w PARCLIP | grep -w $tag | cut -f 5 | sort | uniq`
	cids=(`less $sid_prot_f | grep -w PARCLIP | grep -w $tag | cut -f 2 | sort`)
	idsXL=(`less $sid_prot_f | grep -w mRNAseq | grep -w background | grep XL | cut -f 2 | sort`)
	idswoXL=(`less $sid_prot_f | grep -w mRNAseq | grep -w background | grep -v XL | cut -f 2 | sort`)

	#-----------------
	#XXX
	c=Nucleotide
	nuc=T
	id1=${cids[0]}
	id2=${cids[1]}

	libpat=`echo ${idsXL[@]} | sed 's/ /|/g'`
	mu_ref_XL=`grep -w -P "$libpat" ../../xlinkEnrichment/basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
	libpat=`echo ${idswoXL[@]} | sed 's/ /|/g'`
	mu_ref_woXL=`grep -w -P "$libpat" ../../xlinkEnrichment/basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
	merged_XL=../../xlinkEnrichment/basedOn$c/$tag/full_xlinked_$tag.${id1}_$id2.*.$nuc.xlink_mutrate$mu_ref_XL.bed

	mkdir -p $tag; pushd $tag;

	#get all crosslink foreground scores for 1000 nucleotide vicinity of the PA region
	for f in `ls ../T_hg19.clusters.transcripts.genes.*.sites.gtf.gz`; do
		i=`echo $f | perl -e '$_ = <>; $_ =~ /(\d+).sites.gtf.gz/; print $1;'`;	
		if [ ! -e merged_XL.$i.sites.gtf.gz ]; then
			echo Doing merged_XL.$i.sites.gtf.gz 
			bedtools intersect -wao -s -a $f -b ../$merged_XL | awk -F "\t" 'BEGIN{OFS="\t";} { if ($10 ~ /^chr/) { $6=1-$14; } print $0; }' | \
				cut -f 1-9 | sort -k 6,6gr | gzip -c > merged_XL.$i.sites.gtf.gz
		else
			echo merged_XL.$i.sites.gtf.gz exists
		fi
	done

	if [ "`echo ${libs[@]} | grep $ptag`" == "" ]; then 
		echo "NO PA infered for $tag -skipping steps about shifts (IG) for $f";
		for f in `ls ../${type}hg19.clusters.transcripts.genes.*.sites.gtf.gz`; do
			i=`echo $f | perl -e '$_ = <>; $_ =~ /(\d+).sites.gtf.gz/; print $1;'`;	
			ln -sf $f $type$i.gtf.gz
		done
	else
		pipelibss=`echo ${libs[@]} | sed 's/ /|/g'`;
		controlAseq_id=`less ~bilebi00/_ARE/data/protein_sampleid_list | grep -P -w $pipelibss | grep siCTRL | cut -f 2`;
		knockdownAseq_id=`less ~bilebi00/_ARE/data/protein_sampleid_list | grep -P -w $pipelibss | grep $ptag | cut -f 2`;
		colname=$controlAseq_id-$knockdownAseq_id
		if [ "$knockdownAseq_id" != "" ] && [ $controlAseq_id -gt $knockdownAseq_id ]; then
			colname=$knockdownAseq_id-$controlAseq_id
		fi
		for f in `ls ../${type}hg19.clusters.transcripts.genes.*.sites.gtf.gz`; do
			i=`echo $f | perl -e '$_ = <>; $_ =~ /(\d+).sites.gtf.gz/; print $1;'`;	
			#form gtf file with IG score
			if [ ! -e $i.gtf.gz ]; then
				echo Doing $i $colname
				~bilebi00/_ARE/scripts/assignIG_score.pl $f ../IG_density.$i.sites $colname | \
					sort -k 6,6gr | gzip -c > $type$i.gtf.gz;
			else
				echo $type$i.gtf.gz exists
			fi
		done
	fi
	
	slopls=(500 500 100 100)
	sloprs=(500 0   100 0)
	for slopi in `seq 0 $((${#slopls[@]}-1))`; do
		slopl=${slopls[$slopi]}
		slopr=${sloprs[$slopi]}
		odir=${type}slopl$slopl-slopr$slopr;
		mkdir -p $odir; pushd $odir;
		rm -rf crosslink_distance
		for f in `ls ../../${type}hg19.clusters.transcripts.genes.*.sites.gtf.gz`; do
			i=`echo $f | perl -e '$_ = <>; $_ =~ /(\d+).sites.gtf.gz/; print $1;'`;	
			#------
			#PA clip intersection
			#--------------------
			if [ ! -e $i.clip_intersection.gz ]; then
				zless ../$type$i.gtf.gz | bedtools slop -s -l $slopl -r $slopr -i stdin -g $genomelen | \
					bedtools intersect -wao -s -a stdin -b ../merged_XL.$i.sites.gtf.gz | \
					gzip -c > $i.clip_intersection.gz
			else
				echo $i.clip_intersection.gz exists
			fi
			#xlink profile around PA site
			zless $i.clip_intersection.gz | \
				perl -e '$slopl=shift; $slopr=shift; $tag=shift; while(<>){@t=split/\t/; $x=($t[3]+$t[4]+$slopl-$slopr)/2; $d=int($t[12]-$x); if ($t[6] eq "-") { $x=($t[3]+$t[4]+$slopr-$slopl)/2; $d=int($x-$t[12]); } $t[8] =~ /PAcode "([^"]*)"/; print join("\t",$d, $t[14], $tag.".".$1),"\n"; }' $slopl $slopr $tag.site | \
				sort -k 1,1g | bedtools groupby -g 1 -c 2,2 -o sum,count -i stdin -full | cut -f 1,3-5 >> crosslink_distance
			#Xlink score around PA site
			if [ ! -e $i.clip_posterior_perGene ]; then
				zless $i.clip_intersection.gz | \
					awk -F "\t" 'BEGIN{OFS="\t";}{ print $9, $15, tag; }' tag=$tag.$i.site | \
					sort -k 1,1 | bedtools groupby -g 1 -c 2,2 -o sum,count -i stdin -full | cut -f 1,3-5 > $i.clip_posterior_perGene
			else
				echo $i.clip_posterior_perGene exists
			fi
			R --no-save --args $i.clip_posterior_perGene Xlink < ~bilebi00/_ARE/scripts/plot_posteriorDensity.R > /dev/null
			if [ "`echo ${libs[@]} | grep $ptag`" == "" ]; then 
				echo No $i.IG_perGene generated
			else
				#IG score of PA shift
				if [ ! -e $i.IG_perGene ]; then
					zless $i.clip_intersection.gz | \
						awk -F "\t" 'BEGIN{OFS="\t";}{ print $9, $6, tag; }' tag=$tag.$i.site | \
						sort -k 1,1 | bedtools groupby -g 1 -c 2,2 -o sum,count -i stdin -full | cut -f 1,3-5 | sort -k 3,3g > $i.IG_perGene
				else
					echo $i.IG_perGene exists
				fi
				R --no-save --args $i.IG_perGene IG < ~bilebi00/_ARE/scripts/plot_posteriorDensity.R > /dev/null
			fi
			#IG and xlink overlap
			if [ $i -gt 1 ] && [ "`echo ${libs[@]} | grep $ptag`" != "" ] && [ ! -e $i.fisher ]; then
				#set cutoffs mean + 1 standard deviation for top
				clip_cut=`less $i.clip_posterior_perGene | awk -F "\t" '{print tag"\t"($3/$4);}' tag=$tag | bedtools groupby -g 1 -c 2,2 -o mean,stdev -i stdin | awk '{print $2+$3}'`
				IG_cut=`less $i.IG_perGene | awk -F "\t" '{print tag"\t"($3/$4);}' tag=$tag | sed -e 's/";*//g' | sort | uniq | bedtools groupby -g 1 -c 2,2 -o mean,stdev -i stdin | awk '{print $2+$3}'`
				#all
				less $i.clip_posterior_perGene | \
					awk -F "\t" 'BEGIN{OFS="\t";}{ cc=$3/$4; if (cc>clip_cut) { print $1, cc; }}' clip_cut=-1 | sort | uniq | sort -k 2,2gr \
					> $i.average_clip_posterior_perGeneSymbol.all
				less $i.IG_perGene | \
					awk -F "\t" 'BEGIN{OFS="\t";}{ igc=$3/$4; if (igc>IG_cut) { print $1, igc }}' IG_cut=-1 | sort | uniq | sort -k 2,2gr \
					> $i.average_IG_perGeneSymbol.all
				~bilebi00/_PAPD5/scripts/innerJoin.pl $i.average_clip_posterior_perGeneSymbol.all $i.average_IG_perGeneSymbol.all 1 1 '1,2,4' | sort | uniq | sort -k 2,2gr |
					sed 's/";*/\t/g' | sed 's/\t+/\t/g' > \
					$i.average_clip_posterior_IG_perGeneSymbol.tab
				R --no-save --args $i.average_clip_posterior_IG_perGeneSymbol.tab $i $IG_cut < ~bilebi00/_ARE/scripts/plot_levelplot_for_PAclip_4IG.R > /dev/null 
				rm -rf $i.average_clip_posterior_perGeneSymbol.all $i.average_IG_perGeneSymbol.all
				#XXX not perGeneSymbol but per Site
				#top
				less $i.clip_posterior_perGene | \
					awk -F "\t" 'BEGIN{OFS="\t";}{ cc=$3/$4; if (cc>clip_cut) { print $1, cc; }}' clip_cut=$clip_cut | sort | uniq | sort -k 2,2gr \
					> $i.average_clip_posterior_perGeneSymbol.top
				less $i.IG_perGene | \
					awk -F "\t" 'BEGIN{OFS="\t";}{ igc=$3/$4; if (igc>IG_cut) { print $1, igc }}' IG_cut=$IG_cut | sort | uniq | sort -k 2,2gr \
					> $i.average_IG_perGeneSymbol.top
				~bilebi00/_PAPD5/scripts/innerJoin.pl $i.average_clip_posterior_perGeneSymbol.top $i.average_IG_perGeneSymbol.top 1 1 '1,2,4' | sort | uniq | sort -k 2,2gr \
					> $i.average_clip_posterior_IG_perGeneSymbol.top
				both=`less $i.average_clip_posterior_IG_perGeneSymbol.top | awk '{ print $6; }' | sort | uniq | wc -l`
				clipped=`less $i.average_clip_posterior_perGeneSymbol.top | awk '{ print $6; }' | sort | uniq | wc -l`
				shifted=`less $i.average_IG_perGeneSymbol.top | awk '{ print $6; }' | sort | uniq | wc -l`
				all=`zless ../$type$i.gtf.gz | cut -f 9 | awk '{ print $6; }' | sort | uniq | wc -l`
				R --no-save --args $both $clipped $shifted $all < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R  | grep "^Fisher" -B 1 > $i.fisher
			else
				echo $i.fisher either exists or no shift data exists for $ptag
			fi
		done
		R --no-save --args crosslink_distance < ~bilebi00/_ARE/scripts/plot_profile.R > /dev/null
		gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=clip_posterior_perGene.PosteriorAverage.pdf ?.clip_posterior_perGene.PosteriorAverage.pdf
		if [ "`echo ${libs[@]} | grep $ptag`" == "" ]; then 
			echo No $i.IG_perGene generated
		else
			gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=IG_perGene.PosteriorAverage.pdf ?.IG_perGene.PosteriorAverage.pdf
		fi
		rm -rf ?.IG_perGene.PosteriorAverage.pdf ?.clip_posterior_perGene.PosteriorAverage.pdf
		popd
	done

	echo DONE
	exit;

	#poster
	~bilebi00/_ARE/scripts/calc_informationGain.pl hg19.clusters.transcripts.genes 10,11,12,13,14,15 $libss 1 > IG_density
	~bilebi00/bin/R --no-save --args hg19.clusters.transcripts.genes 11 16 IG_density ${libs[@]} < ~bilebi00/_ARE/scripts/plot_IG.R > /dev/null

	~bilebi00/_ARE/scripts/calc_informationGain.pl hg19.clusters.transcripts.genes 16,17,18,19,20,21 $libss 1 > IG_rawcount
	~bilebi00/bin/R --no-save --args hg19.clusters.transcripts.genes 17 22 IG_rawcount ${libs[@]} < ~bilebi00/_ARE/scripts/plot_IG.R > /dev/null

	exit;

	#TODO IG scores for densities and raw counts are not totally same
	#check their correlation
	#intersect the data with clip and make profiles for other ARE-BPs
fi

if [ $task == 5 ]; then

	type=clean_
	#XXX
	#	for f in slopl1000-slopr1000 slopl500-slopr500 slopl100-slopr100 slopl100-slopr0; do
	for f in ${type}slopl500-slopr500; do
		ll | grep "^d" | awk '{ print $9; }' | grep -v -P "\.|CFIm68|EWSR1" | while read i; do 
		cat $i/$f/crosslink_distance ; 
		done > $f.crosslink_distance
		R --no-save --args $f.crosslink_distance  < ~bilebi00/_ARE/scripts/plot_profile_4Groups_4Xlink.R > /dev/null
	done

#XXX distance of polyA sites to each other
bedtools closest -t first -io -D "a" -d -s -a 4.gtf.gz -b 4.gtf.gz | \
	awk -F "\t" '{ if ($19<0) { $19=-$19; } print tag"\t"$19; }' tag=gaga | sort -k 2,2 | bedtools groupby -g 1 -c 2,2,2 -o mean,stdev,min

#CFI not shifting anymore!!! XXX
cut -f 1 Analysis/PAchange_wIG/253/_CFIm68/slopl60-slopr0/2.average_clip_posterior_IG_perGeneSymbol.top | while read i; do echo $i; grep -w $i Analysis/PAchange_wIG/253/CFIm68/slopl60-slopr0/2.average_clip_posterior_IG_perGeneSymbol.top; echo; done | less

#levelplot!!!

exit;
for i in 1 2; do
	for f in $i.clip_posterior_perGene $i.extended_clip_posterior_perGene $i.IG_perGene ; do
		gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$tag.$f.pdf */$f.*.pdf
	done
done

exit;
#TODO
innerJoin.pl hnRNPQ/2.average_IG_perGeneSymbol.top ../hnRNPQold/hnRNPQ/_slop_b100/2.average_IG_perGeneSymbol.top  1 1 ''
head */*fisher
head ../hnRNPQold/hnRNPQ/_slop*/*fisher


fi
