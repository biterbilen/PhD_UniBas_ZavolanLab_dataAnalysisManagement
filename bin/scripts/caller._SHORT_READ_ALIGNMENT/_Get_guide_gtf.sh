#!/bin/bash

~bilebi00/_SHORT_READ_ALIGNMENT/scripts/clsformat2gtf.pl "/import/bc2/home/zavolan/bilebi00/_SAM68/data/coordinates/EXON/chr*" ~bilebi00/_SAM68/data/mm_TR.info.all GMAP exon "type" > ~bilebi00/_SAM68/data/GMAP_EXONS.gtf 2> ~bilebi00/_SAM68/data/GMAP_EXONS.gtf.log
