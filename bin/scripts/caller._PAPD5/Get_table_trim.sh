#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o log_file.table_trim

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

r=$RANDOM;

echo region > $r
awk '{ printf("%s_%s_%s_%s\n",$2,$3,$4,$5); }' DB_copies_364/copy* >> $r

tag="table_DB_copies";

for i in DB_copies_*/copy*; do
	echo Processing $i
	f=`perl -e '$f=shift; $f=~/(copies_\d+)/; print $1;' $i`;
	t=`perl -e '$f=shift; $f=~/copies_(\d+)/; print $1;' $i`;
	tag="${tag}_${t}";
	echo $f > $r.$t;
	cut -f 6 $i >> $r.$t;
done  
echo $tag

echo "Pasting $r*";
paste $r* > $tag; 

rm -rf $r*;

echo "~/scripts/trim_min_table_values.pl $tag 0.1 > $tag.trim01";
~/_PAPD5/scripts/trim_min_table_values.pl $tag 0.1 > $tag.trim01;

echo "~/scripts/trim_min_table_values.pl $tag 0.05 > $tag.trim005";
~/_PAPD5/scripts/trim_min_table_values.pl $tag 0.05 > $tag.trim005;

echo Done;
