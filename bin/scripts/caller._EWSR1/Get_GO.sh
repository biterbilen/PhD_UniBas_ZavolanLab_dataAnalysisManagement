#!/bin/bash

outdir=Analysis
pushd $outdir

soutdir=GO
pushd $soutdir

#TODO write me
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz

less idmapping_selected.tab.gz | cut -f 1,3 > tmp.tab

perl -e 'while (<>) { chomp; @l = split /\t/; next if @l == 1; @k = split /; /, $l[1]; for my $gene (@k) { print join("\t", $l[0], $gene)."\n" } }' < tmp.tab > UniProtKB_GeneID.tab

wget http://www.geneontology.org/ontology/gene_ontology_edit.obo

wget 'http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD'

wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/gene_association.goa_mouse.gz
mv 'gene_association.goa_human.gz?rev=HEAD' gene_association.goa_human.gz

gunzip gene_association.goa_human.gz
gunzip gene_association.goa_mouse.gz

#ID MAPPING
#http://www.uniprot.org/?tab=mapping

perl ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/GO/GeneId2UniProt.pl studyset.geneids > studyset.uniprot

# Background
/import/bc2/home/zavolan/GROUP/RNA-BP_CLIP_Project/Methods/MRNAExpression/expressed_transcripts_No_XL

#grep -v -P '^From' background.mapped | cut -f2 > population.uniprot
less 201108253114T*tab | grep -v -P '^From' | cut -f2 > population.uniprotgrep
grep -v -P '^From' analysis1/studyset.uniprot.pre | cut -f2 > analysis1/studyset.uniprot


java -jar Ontologizer.jar -a gene_association.goa_human -g gene_ontology_edit.obo -s analysis1/studyset.uniprot -p analysis1/population.uniprot -c Parent-Child-Union -m Westfall-Young-Single-Step -d 0.1 -r 1000 -o analysis1

java -jar Ontologizer.jar -a gene_association.goa_human -g gene_ontology_edit.obo -s analysis2/studyset.uniprot -p analysis2population.uniprot -c Parent-Child-Union -m Westfall-Young-Single-Step -d 0.1 -r 1000 -o analysis2


# MOUSE

grep -v -P '^From' MOUSE_bg.mapped | cut -f2 > mouse1/population.uniprot
grep -v -P '^From' MOUSE_fg.mapped | cut -f2 > mouse1/studyset.uniprot

java -jar Ontologizer.jar -a gene_association.goa_mouse -g gene_ontology_edit.obo -s mouse1/studyset.uniprot -p mouse1/population.uniprot -c Parent-Child-Union -m Westfall-Young-Single-Step -d 0.1 -r 1000 -o mouse1
