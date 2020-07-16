#!/bin/bash

if [[ $# -ne 1 ]]
then
	echo "Copy rawdata from source"
	echo "Usage: bash $0 source_dir"
	exit
fi

IN=$1
OUT=00-rawdata_copy

rsync -avzP --delete $IN/* $OUT/
