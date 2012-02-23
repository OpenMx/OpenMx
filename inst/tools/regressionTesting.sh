#!/bin/bash
if [ "$#" -ne "4" ];
then
    echo ""
    echo "expected arguments: [NEWDIR] [OUTDIR] [START-REV] [END-REV]"
    echo ""
    exit
fi

SVN_ROOT=$1
OUTPUT_DIR=$2
START_REVISION=$3
END_REVISION=$4

mkdir -p $SVN_ROOT
svn checkout -r $START_REVISION http://openmx.psyc.virginia.edu/svn/trunk $SVN_ROOT
cd $SVN_ROOT
for ((rev=$START_REVISION; rev <= $END_REVISION; rev++))
do
   svn update -r $rev --non-interactive
   make install
   make nightly > performance.results
   awk '/Runtimes:/ {flag=1}flag' performance.results | tail -n +2 | sed '$d' > performance.data
   cat performance.data | while read line;
   do
      set $line
      location=`basename $1`
      echo -e "$rev\t$2" >> $OUTPUT_DIR/$location.results
   done
done
rm -rf $SVN_ROOT
