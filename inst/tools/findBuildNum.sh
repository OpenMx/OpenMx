#!/bin/sh

set -o errexit
set -o nounset

v=$(svnversion -c)
if [ "$v" = 'Unversioned directory' -o "$v" = exported ]; then
  if ! which git >/dev/null; then echo 0; exit; fi
  git describe | sed -e 's/^v//' -e 's/-[^-]*$//'
else
  # The 1.0 is just a fake placeholder
  echo 1.0-$v | sed -e 's/[MS]//g' -e 's/[[:digit:]]*://'
fi
