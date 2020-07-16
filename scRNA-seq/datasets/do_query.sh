#!/bin/bash

sns=$1

echo $sns | tr ',' '\n' | while read sn
do
	cat stat.txt | awk -v ks=$sn '$2~"^"ks"_"'
done | sort | cut -f 1 | sed -e 's/^/"/' -e 's/$/"/' | paste -s -d ','

echo $sns | tr ',' '\n' | while read sn
do
	cat stat.txt | awk -v ks=$sn '$2~"^"ks"_"'
done | sort | cut -f 2 | sed -e 's/^/"/' -e 's/$/"/' | paste -s -d ','

