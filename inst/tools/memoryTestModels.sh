#!/bin/bash

DIR1=demo
DIR2=models/passing

pushd $DIR1

for f in *.R
do
  echo "Running $f..."
  R -d "valgrind --tool=memcheck --leak-check=full --track-origins=yes --suppressions=../inst/tools/OpenMx.supp" --vanilla --slave < $f
done

popd

pushd $DIR2

for f in *.R
do
  echo "Running $f..."
  R -d "valgrind --tool=memcheck --leak-check=full --track-origins=yes --suppressions=../../inst/tools/OpenMx.supp" --vanilla --slave < $f
done

popd
