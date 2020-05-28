#!/bin/bash

r=$1
echo "bin name refdepth genfrac"
for f in `ls mq-*/summary/TXT/Genome_fraction.txt`; do
    x=`grep -v "-" $f`
    bin=`echo $x|cut -d' ' -f2|cut -d '_' -f2`
    echo -n "$bin "
    for i in `echo $x|cut -d' ' -f3,4`; do
        y=`grep $i $r|lastcol.py`
        echo -n "$i $y"
    done
    echo ""
done
