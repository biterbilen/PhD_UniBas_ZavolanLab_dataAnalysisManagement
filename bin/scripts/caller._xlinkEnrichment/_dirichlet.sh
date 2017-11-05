bedtools intersect -a ../HuR.005.clusters.T.xlink_count.gtf.gz -b ../HuR.006.clusters.T.xlink_count.gtf.gz -s -wb | awk -F "\t" 'BEGIN{OFS="\t";}{ $9=$9" "$18; print; }' > 56.gtf
R --no-save --args 56.gtf count < scatter_density.R
bedtools intersect -a ../HuR.007.clusters.T.xlink_count.gtf.gz -b ../HuR.008.clusters.T.xlink_count.gtf.gz -s -wb | awk -F "\t" 'BEGIN{OFS="\t";}{ $9=$9" "$18; print; }' > 78.gtf
R --no-save --args 78.gtf count < scatter_density.R
less a.gtf | grep -v -P "copies_count \"[0-9]\"" | grep -v -P "copies_count \"1[0-9]\"" > 34.gtf
nohup R --no-save --args 34.gtf count < scatter_density.R & 


less replicate1_XL.bed | awk '$5<0.05{print}' > a
less a | bedtools closest -s -d -t all -D a -a stdin -b a -io > a.dist

less a.dist | awk '{ print type"\t"$13;}' type=distance | sort | uniq -c | sort -k 3,3gr > a.dist.hist
R --no-save --args a.dist.hist < ~bilebi00/_xlinkEnrichment/scripts/self_distance_profile.R

#-------------

#self intersection profile
for N in 1 2; do
	less sample${N}_XL.postcut.sorted.bed | head -n 1000 > $N.1000
	bedtools closest -s -d -t all -D a -a $N.1000 -b $N.1000 -io > $N.selfdist
	less $N.selfdist | awk '{ print type"\t"$13;}' type=distance | sort | uniq -c | sort -k 3,3gr > $N.selfdist.hist
	R --no-save --args $N.selfdist.hist < ~bilebi00/_xlinkEnrichment/scripts/self_distance_profile.R > /dev/null
done

#get intersecting ones within a distance

db=hg18
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
genomelen=~bilebi00/aux/human.$db.genome
distance=1
bedtools closest -s -d -t first -D a -a 1.1000 -b 2.1000 | awk '$13 > -distance && $13 < distance { print ;}' distance=$distance | wc -l

less 1.1000 | bedtools slop -b $distance -i stdin -g $genomelen > tmp.slop.1
less 2.1000 | bedtools slop -b $distance -i stdin -g $genomelen > tmp.slop.2
bedtools closest -s -d -t "first" -a tmp.slop.1 -b tmp.slop.2 | awk '$13==0{print}' | wc -l

