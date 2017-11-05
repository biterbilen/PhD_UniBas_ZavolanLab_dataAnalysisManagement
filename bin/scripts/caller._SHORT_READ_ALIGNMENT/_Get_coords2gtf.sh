#!/bin/bash 

~bilebi00/_SHORT_READ_ALIGNMENT/scripts/clsformat2gtf.pl "~bilebi00/_SHORT_READ_ALIGNMENT/data/coordinates/EXON/chr*" ~bilebi00/_SHORT_READ_ALIGNMENT/data/hg_TR.info.all GMAP exon > ~bilebi00/_SHORT_READ_ALIGNMENT/data/GMAP_EXONS.gtf 2> ~bilebi00/_SHORT_READ_ALIGNMENT/data/GMAP_EXONS.gtf.log
