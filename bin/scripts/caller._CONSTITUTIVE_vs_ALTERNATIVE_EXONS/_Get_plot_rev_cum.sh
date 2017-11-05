#TODO put me at the end of run_getT2C.pl 

./scripts/prep_reverse_cum_forplotting.pl > OUT/reverse_cum_Reproducibility_v05.rep
R --no-save --args OUT/reverse_cum_Reproducibility_v05.rep < ./scripts/plot_reverse_cum.R 

