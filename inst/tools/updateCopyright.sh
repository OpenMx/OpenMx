#!/bin/sh

fileList=$(comm -12 <(git grep -l 'Copyright' | sort) <(git diff --name-only @{1.year.ago} | sort))

perl -pi -e 's/Copyright (\d+) by the individuals mentioned in the source code history/Copyright $1-2018 by the individuals mentioned in the source code history/' $fileList

perl -pi -e 's/Copyright (\d+)-(\d+) by the individuals mentioned in the source code history/Copyright $1-2018 by the individuals mentioned in the source code history/' $fileList
