#!/bin/sh

set -o errexit
set -o nounset

v=$(svnversion -c)
if [ "$v" = 'Unversioned directory' -o "$v" = exported ]; then
  git log -1 | grep git-svn-id | cut -d '@' -f 2 | cut -d ' ' -f 1
else
  echo $v | sed -e 's/[MS]//g' -e 's/^[[:digit:]]*://'
fi
