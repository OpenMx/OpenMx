#!/bin/sh

set -o errexit
set -o nounset

v=$(svnversion -c)
if [ "$v" = 'Unversioned directory' -o "$v" = exported ]; then
  git log | grep git-svn-id | head -n 1 | cut -d '@' -f 2 | cut -d ' ' -f 1
else
  echo $v | sed -e 's/[MS]//g' -e 's/^[[:digit:]]*://'
fi
