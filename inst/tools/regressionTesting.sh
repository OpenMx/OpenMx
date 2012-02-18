#!/bin/bash
if [ "$#" -ne "4" ];
then
    echo ""
    echo "expected arguments: [SVNPATH] [OUTDIR] [START-REV] [END-REV]"
    echo ""
    exit
fi

SVN_ROOT=$1
OUTPUT_DIR=$2
START_REVISION=$3
END_REVISION=$4

cd $SVN_ROOT
for ((rev=$START_REVISION; rev <= $END_REVISION; rev++))
do
   svn update -r $rev --non-interactive --force
   svn revert --recursive .
   make install
   make nightly > performance.results
   awk '/Results:/ {flag=1}flag' performance.results | tail -n +2 > performance.data
   cat performance.data | while read line;
   do
      set $line
      echo -e "$rev\t$2" >> $OUTPUT_DIR/$1.results
   done
done

