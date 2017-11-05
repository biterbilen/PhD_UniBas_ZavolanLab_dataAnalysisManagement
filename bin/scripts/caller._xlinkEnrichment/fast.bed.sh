#!/bin/bash


s=Get_error_bedgraph.sh; 
for i in `seq 2 2`; do
#for i in `seq 0 7`; do
	echo Doing $i
	sed -e "s/_tagindex/$i/" $s > _$s$i; 
	for k in `seq 1 40`; do 
		sh _$s$i $k &> log.$i.$k & 
	done
	sleep 15m
done
