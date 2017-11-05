for i in DB_*; do echo $i; p=`grep $i ./data/protein_sampleid_list | cut -f 1`; less $i/poly*ec | grep -v NULL | cut -f 1 | sort | uniq -c | head > $p.PAScount;  done
./scripts/outerJoin_two_col_data.pl *PAScount > PAS_count
R --no-save --args PAS_count < ./scripts/plot_corr_scatter.R 
