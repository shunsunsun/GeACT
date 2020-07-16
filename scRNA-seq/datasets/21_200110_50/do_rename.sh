#!/bin/bash

cd 00-rawdata

function do_xxx {
rename _XXL_ _HCA_ *
rename _XXL_ _HCA_ */*
rename -- _CGM_ -CGM- *
rename -- _CGM_ -CGM- */*
rename -- -CGM- _CGM_ *
rename -- -CGM- _CGM_ */*
mkdir -p tmp
mv *CGM* tmp/
ls -d PR10_HCA_* | while read id; do pid=$(echo $id | sed 's/\([^-]\+\)-\(.*\)/\1_D_\2/'); cd $id; rename $id $pid *; cd -; done
ls -d PR10_HCA_* | while read id; do pid=$(echo $id | sed 's/\([^-]\+\)-\(.*\)/\1_D_\2/'); mv $id $pid; done
mv tmp/* .
rm -rf tmp
}

# 1
#echo CGM C | while read ida idb
# 2
echo D D | while read ida idb
do 
rename _${ida}_1- _${idb}-A *
rename _${ida}_1- _${idb}-A */*.gz
rename _${ida}_2- _${idb}-B *
rename _${ida}_2- _${idb}-B */*.gz
rename _${ida}_3- _${idb}-C *
rename _${ida}_3- _${idb}-C */*.gz
rename _${ida}_4- _${idb}-D *
rename _${ida}_4- _${idb}-D */*.gz
done

