#!/bin/sh

set -o errexit
set -o nounset

if ! which git >/dev/null; then
    echo 2.0
elif ver=$(git describe) 2>/dev/null; then
    echo "$ver" | sed -e 's/^v//' -e 's/-[^-]*$//'
else
    echo 2.0
fi
