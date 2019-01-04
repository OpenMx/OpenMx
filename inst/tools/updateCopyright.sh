#!/bin/bash

anchor=v2.9.8-5-gbcdf9f416  # Immediately after the previous year's copyright update
year=2019

fileList=$(comm -12 <(git grep -l 'Copyright' | sort) <(git diff --name-only $anchor | sort))

perl -pi -e 's/Copyright (\d+) by the individuals mentioned in the source code history/Copyright $1-'$year' by the individuals mentioned in the source code history/' $fileList

perl -pi -e 's/Copyright (\d+)-(\d+) by the individuals mentioned in the source code history/Copyright $1-'$year' by the individuals mentioned in the source code history/' $fileList
