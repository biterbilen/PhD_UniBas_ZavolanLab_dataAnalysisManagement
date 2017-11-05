#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
# $ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o log_file

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

cd output_directory

tag=tag
window=polyA_window;
top=topN_polyA_sites;
itop=$((top / 3));

for f in *copies.refined; do grep "^##" $f | sort -k 6,6gr | head -n $itop; done | sort -k 6,6gr | head -n $top | perl -e 'my $i=1; while(<>){ my($a,@t)=split; print join("\t","##$i", @t), "\n"; $i++; }' > maxcopy_Top$top.copies.refined;

rm -rf chr*copies.refined

echo Preparing $window separated peak sites from $tag lib 
~bilebi00/GEORGES_POLYA/scripts/getPolyASites.pl maxcopy_Top$top.copies.refined $window $tag > polyAsites

for i in 100 1000 10000; do 
	genes=`less polyAsites| head -n $i | sort -k1,1 -k4,4g | cut -f 1 | sort -u | wc -l`;
	multgenes=`less polyAsites| head -n $i | sort -k1,1 -k4,4g | cut -f 1 | uniq -D | uniq -d | wc -l`;
	ratio=`perl -e '$g=shift; $mg=shift; print 100*$mg/$g, "\n"' $genes $multgenes`;
	echo among top $i PAS and $genes genes, $multgenes multi-polyA site having genes exist: $ratio percent
done

#annotation event count
~bilebi00/GEORGES_POLYA/scripts/annotate_sites_with_eventCount.pl polyAsites > polyAsites.tmp.ec

##annotations
#TODO make use of old code wrap here
#closest PolyADB distance
#zless ~bilebi00/DATA/PolyA_DB_mammal_human_hg18.gz | awk '{OFS="\t";} { if (NR > 1 ) { print $5, $2, $7, $3+1, $4, "POLYA"; } }'  | ~bilebi00/GEORGES_POLYA/scripts/annotate_sites_with_closestRegionDistance.pl polyAsites.tmp.ec POLYA_$tag ~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel

#closest POLYA_$tag distance (self distance to the gene alternative PAS sites)

#For null annotations check closest distance genes; these could be the novel polyA sites maybe; if not internal
#get their distance


#Check by aligning to AATAAA motif noncanonical polyA sites


