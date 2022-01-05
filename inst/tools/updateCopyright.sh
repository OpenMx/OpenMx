#!/bin/bash

anchor=v2.19.8-87-g6916439dd  # Immediately after the previous year's copyright update
year=2022

fileList=$(comm -12 <(git grep -l 'Copyright' | sort) <(git diff --name-only $anchor | sort))

perl -pi -e 's/Copyright (\d+) by the individuals mentioned in the source code history/Copyright $1-'$year' by the individuals mentioned in the source code history/' $fileList

perl -pi -e 's/Copyright (\d+)-(\d+) by the individuals mentioned in the source code history/Copyright $1-'$year' by the individuals mentioned in the source code history/' $fileList
