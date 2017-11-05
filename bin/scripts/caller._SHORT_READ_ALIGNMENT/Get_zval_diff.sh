#/bin/bash

outdir=Shivendra/MARA
pushd $outdir;

for i in mara_report_just_uniq/profiles/*/zvals.dat; do echo $i `less $i | perl -e 'while(<>){chomp;@t=split;$h{$t[2]}=$t[1];} print $h{"mRNAseq_siCTRL-A_HeLa.uniq"}-$h{"mRNAseq_siEWSR1_HeLa.uniq"}, "\n"'`; done | sort -k 2,2g | less | head
echo
for i in mara_report_just_uniq/profiles/*/zvals.dat; do echo $i `less $i | perl -e 'while(<>){chomp;@t=split;$h{$t[2]}=$t[1];} print $h{"mRNAseq_siCTRL-A_HeLa.uniq"}-$h{"mRNAseq_siEWSR1_HeLa.uniq"}, "\n"'`; done | sort -k 2,2gr | less | head
echo
for i in mara_report_uniq/profiles/*/zvals.dat; do echo $i `less $i | perl -e 'while(<>){chomp;@t=split;$h{$t[2]}=$t[1];} print $h{"mRNAseq_siCTRL-A_HeLa.uniq"}-$h{"mRNAseq_siEWSR1_HeLa.uniq"}, "\n"'`; done | sort -k 2,2g | less | head
echo
for i in mara_report_uniq/profiles/*/zvals.dat; do echo $i `less $i | perl -e 'while(<>){chomp;@t=split;$h{$t[2]}=$t[1];} print $h{"mRNAseq_siCTRL-A_HeLa.uniq"}-$h{"mRNAseq_siEWSR1_HeLa.uniq"}, "\n"'`; done | sort -k 2,2gr | less | head
echo
for i in mara_report_all/profiles/*/zvals.dat; do echo $i `less $i | perl -e 'while(<>){chomp;@t=split;$h{$t[2]}=$t[1];} print $h{"mRNAseq_siCTRL-A_HeLa"}-$h{"mRNAseq_siEWSR1_HeLa"}, "\n"'`; done | sort -k 2,2g | less | head
echo
for i in mara_report_all/profiles/*/zvals.dat; do echo $i `less $i | perl -e 'while(<>){chomp;@t=split;$h{$t[2]}=$t[1];} print $h{"mRNAseq_siCTRL-A_HeLa"}-$h{"mRNAseq_siEWSR1_HeLa"}, "\n"'`; done | sort -k 2,2gr | less | head
echo
