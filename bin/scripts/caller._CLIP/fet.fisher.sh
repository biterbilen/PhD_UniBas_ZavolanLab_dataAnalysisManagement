less FET.feature.stats  | \
 	awk -F "\t" 'BEGIN{OFS="\t";} { if ($4>0 && NR>1) { $4=1;} print } ' | sort -k 1,4 | bedtools groupby -g 1,2,3,4 -c 4 -o count -i stdin | grep joint
echo EWSR1
R --no-save --args 33 481 94 6286 1 < ~/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher"
echo FUS
R --no-save --args 60 454 215 6165 1 < ~/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher"
echo FUS_mut
R --no-save --args 31 483 57 6323 1 < ~/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher"
echo TAF15
R --no-save --args 24 490 54 6326 1 < ~/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher"
