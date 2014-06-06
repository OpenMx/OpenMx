#!/bin/sh

set -o errexit
set -o nounset

v=$(svnversion -c)
if [ "$v" = 'Unversioned directory' -o "$v" = exported ]; then
  if ! which git >/dev/null; then echo 0; exit; fi
  if ! git log -1 | grep -q git-svn-id; then
      echo "WARNING: Subversion does not have all GIT commits in your current branch" >/dev/fd/2
  fi
  git log | grep git-svn-id | head -n 1 | cut -d '@' -f 2 | cut -d ' ' -f 1
else
  echo $v | sed -e 's/[MS]//g' -e 's/^[[:digit:]]*://'
fi
