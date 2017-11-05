#!/bin/bash

ofile=clipz_tables_general.content
rm -rf $ofile

less clipz_tables_general | while read i; do
	echo $i >> $ofile
	/import/bc2/soft/bin/mysql --quick -u mirzadmin --password=wUP4gZMW -h wnzmsql01 -D clipz -A >> $ofile <<END
		select * from $i limit 3; 
		quit
END
	echo >> $ofile
done
