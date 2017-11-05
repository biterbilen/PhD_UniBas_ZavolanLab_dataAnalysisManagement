id=1597
date
mapped_sequences=$id.mapped
    /import/bc2/soft/bin/mysql --quick -u mirzadmin --password=wUP4gZMW -h wnzmsql01 -D clipz -A > $mapped_sequences <<END
      select * from t_mapped_sequence_$id where genome_count_total=1;
      quit
END
date
