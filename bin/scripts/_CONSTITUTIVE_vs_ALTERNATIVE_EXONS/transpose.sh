#!/bin/bash

# ./transpose.sh input 8
# ./transpose.sh input

f=$1
for i in $(seq `less $f | head -n 1 | awk '{ print NF}'`); do
	less $f | awk '{print $i; }' i=$i | tr '\n' '\t'
  echo
done

