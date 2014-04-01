#!/bin/sh

v=$(svnversion -c)
if [ "$v" = 'Unversioned directory' -o "$v" = exported ]; then
  echo 0  # date -u +'%s' changes too often
else
  echo $v | sed -e 's/[MS]//g' -e 's/^[[:digit:]]*://'
fi
