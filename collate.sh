#! /bin/bash

#Script to collect all eigenM files output by Split_Measure.py and generate a summary text file (by calling some more python)
#The idea behind this design is that it makes it more robust so that if Split_Measure.py crashes then you do not have to repeat all measurements to get a complete text file!.

STAT=$1

echo $STAT
ls ./Data/${STAT}*.sac| grep BHN  > ${STAT}_read_stream.txt
