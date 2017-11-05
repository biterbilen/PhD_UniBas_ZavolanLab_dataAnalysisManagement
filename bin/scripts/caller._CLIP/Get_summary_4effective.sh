#!/bin/bash

sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
typenames=(general annotation mutation)
topNs=(100 100 100)
stacks=(0 0 0)
tag=EWSR1 #project
tag=ARE #project
cell=HEK293
outdir=Analysis/Summary4Effective${cell}_project$tag
mkdir -p $outdir; pushd $outdir;

ids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | grep -w $cell | awk -F "\t" '$8 == tag { print $2; }' tag=$tag | sort`)
echo "${ids[@]}"

for j in `seq 0 $((${#typenames[@]}-1))`; do
	topN=${topNs[$j]}
	typename=${typenames[$j]}
	stack=${stacks[$j]}

	#clean files generated before
	rm -rf $tag.$typename*

	echo "lib type value N" | awk 'BEGIN{OFS="\t";}{ print $1, $2, $3, $4; }' > $tag.$typename

	for i in `seq 0 $((${#ids[@]}-1))`; do
		id=${ids[$i]}
		idn="`less $sid_prot_f | awk -F "\t" '$2==id { print $6}' id=$id`";
		#link data
		ln -sf ~bilebi00/_CLIP/Analysis/UniqueRawData_ALL_woContamination/$id.${typename}_summary .
		#prep for plotting
		t=`less $id.${typename}_summary  | awk -F "\t" 'BEGIN{c=0;}{ if (NR>1) { c= c+$2; print c; }  }' | tail -n 1`
		if [ $typename == general ]; then
			t=`less $id.${typename}_summary | awk -F "\t" 'NR==2{ print $4; }'`
		fi
		#type proportion
		less $id.${typename}_summary | awk -F "\t" 'BEGIN{OFS="\t";}{ if (NR>1) { print lib, $1, $2/t*100,t; }  }' lib="$idn" t=$t >> $tag.$typename;
		#add mutation subtype proportion here
		if [ $typename == mutation ]; then
			pat=MT
			#subtype header
			if [ ! -e $tag.$typename$pat ]; then
				head -n 1 $tag.$typename > $tag.$typename$pat
			fi
			#subtype total
			t=`less $id.${typename}_summary | grep -P "mutation|^$pat" | awk -F "\t" 'BEGIN{c=0;}{ if (NR>1) { c= c+$2; print c; }  }' | tail -n 1`
			#subtype proportion
			less $id.${typename}_summary | grep -P "mutation|^$pat" | \
				awk -F "\t" 'BEGIN{OFS="\t";}{ if (NR>1) { print lib, $1, $2/t*100,t; }  }' lib="$idn" t=$t >> $tag.$typename$pat;
		fi
	done
	#plot samples together
	for f in $tag.$typename*; do
		R --no-save --args $f $topN $stack < ~bilebi00/_EWSR1/scripts/barchart_summary.R > /dev/null
	done
done

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.summary.pdf $tag*pdf
