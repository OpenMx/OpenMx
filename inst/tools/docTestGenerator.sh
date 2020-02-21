#!/bin/bash

mkdir -p staging/doctest
for file in docs/source/*.rst
do
	new1=$(basename $file .rst)
	newname=${new1}.R
	inst/tools/sphinxToR.py < $file > staging/doctest/$newname
done
