#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
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
	echo -e '.separator "\t"\ndelete from der_metadata where name="'$bn'";\n.import '$wd/$indir'/'$bn'.metadata der_metadata' | sqlite3  $dbname
	echo -e '.separator "\t"\ndelete from der_data where name="'$bn'";\n.import '$wd/$indir'/'$bn'.data der_data' | sqlite3  $dbname
done

