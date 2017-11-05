#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o log_file

export PERL5LIB=$HOME/lib:$PERL5LIB

#TODO
exit;
#There is one extra field in the input file (check 4WEB_DIS3L2)
#handle that

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

indir=indir
dbname=Differential_Expression
dboutdir=dboutdir

wd=`pwd`;
mkdir -p $dboutdir; cd $dboutdir;

#delete the tables if exists and create new
#sqlite3 $dbname < ~bilebi00/_SAM68/data/db_structure.sql


#load data
for i in $wd/$indir/*metadata; do   
	bn=`basename $i .metadata`; 
	echo -e '.separator "\t"\ndelete from merged_der_metadata where name="'$bn'";\n.import '$wd/$indir'/'$bn'.metadata merged_der_metadata' | sqlite3  $dbname
	echo -e '.separator "\t"\ndelete from merged_der_data where name="'$bn'";\n.import '$wd/$indir'/'$bn'.data merged_der_data' | sqlite3  $dbname
	#added for gbrowse library links
	echo -e 'drop view if exists merged_der; create view merged_der as select * from merged_der_data, merged_der_metadata where merged_der_data.name==merged_der_metadata.name;' | sqlite3 $dbname
	echo -e '.output merged_der\n.separator "\t"\n.header on\nselect * from merged_der;' | sqlite3 $dbname
done

