#!/bin/sh

set -o errexit
set -o nounset

if ! which git >/dev/null; then
  echo 2.0
else
  git describe | sed -e 's/^v//' -e 's/-[^-]*$//'
fi
