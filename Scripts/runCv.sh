#!/bin/bash

Tr=$1
Bd=$2

for net in `ls Networks/CLASH*.dat`
do
label=`echo $net| cut -f 2 -d / | awk -F '_links' '{print $1}' `;
time python Functions/Cv_analyser.py  $net parameterFile.txt $label . $Tr $Bd;
done
