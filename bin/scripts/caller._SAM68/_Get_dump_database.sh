#!/bin/bash

#sqlite3 der
less ~bilebi00/_SAM68/4WEB/a.sql | sqlite3 a 

.show
.separator "\t"
.import '/import/bc2/home/zavolan/bilebi00/_SAM68/4WEB/a.gtf' a
.tables
.schema

#The inverse of the .dump command is the .read command. The syntax would be
#sqlite3> .read export.sqlite3.sql
#or
#echo '.read export.sqlite3.sql' | sqlite3 my_database.sqlite3
#.dump coffees
