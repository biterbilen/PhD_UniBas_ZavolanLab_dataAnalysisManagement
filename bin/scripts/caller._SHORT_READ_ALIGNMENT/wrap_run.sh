#!/bin/bash

nohup ./run_isoform_DE.pl hg18 Zavolan 8 /import/bc2/home/zavolan/bilebi00/_SHORT_READ_ALIGNMENT/data/Zavolan_samples /import/bc2/home/zavolan/bilebi00/DATA/Bowtie_DB_INDEX /import/bc2/home/zavolan/bilebi00/DATA/hg18_ucsc_tracks/rmsk.gtf.gz /import/bc2/home/zavolan/bilebi00/_SHORT_READ_ALIGNMENT/data/GMAP_EXONS.gtf /import/bc2/home/zavolan/bilebi00/DATA/hg18_ucsc_tracks/rnaGene.gtf.gz &> zav.log &
