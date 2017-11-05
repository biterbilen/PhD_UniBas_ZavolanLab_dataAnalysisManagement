less 1598_1635.annot.bedplus.gz | cut -f 1-6 | bedtools genomecov -trackopts 'name="My Track" visibility=2 color=255,30,30' -bg -strand - -i stdin -g ~/aux/human.hg19.genome > hh.bedgraph
bedGraphToBigWig hh.bedgraph ~/aux/human.hg19.genome hh.bigwig

#TODO do annotation track as bed12
