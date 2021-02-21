#!/bin/bash

cd 00-rawdata

ls | sed 's/.*HCA_\([^-]\+\).*/\1/' | sort -u | while read ida
do
idb=$(echo $ida | sed 's/11/_E/')
rename _${ida}-1- _${idb}-A *
rename _${ida}-1- _${idb}-A */*.gz
rename _${ida}-2- _${idb}-B *
rename _${ida}-2- _${idb}-B */*.gz
rename _${ida}-3- _${idb}-C *
rename _${ida}-3- _${idb}-C */*.gz
rename _${ida}-4- _${idb}-D *
rename _${ida}-4- _${idb}-D */*.gz
done

