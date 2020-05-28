#!/bin/bash

r=$1
for f in `ls mq-*/summary/TXT/Genome_fraction.txt`; do
    x=`grep -v "-" $f`
    bin=`echo $x|cut -d' ' -f2|cut -d '_' -f2`
    echo -n "$bin "
    for i in `echo $x|cut -d' ' -f3,4`; do
        y=`grep $i $r|lastcol.py`
        echo -n "   $i $y"
    done
    echo ""
    #for i in $x; do echo $i; done
    #echo $x
#    echo -n "$x "
#    y=`echo $x|cut -d' ' -f3`
    #    grep $y /scratch2/shofmeyr/spiked-genomes/synth64d-4/d*.txt|cut -d':' -f1|cut -d'/' -f6
#    grep $y $2
done
