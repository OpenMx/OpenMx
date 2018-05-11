#!/bin/sh

perl -pi -e 's/Copyright (\d+) by the individuals mentioned in the source code history/Copyright $1-2018 by the individuals mentioned in the source code history' $(git grep -l 'Copyright')

perl -pi -e 's/Copyright (\d+)-(\d+) by the individuals mentioned in the source code history/Copyright $1-2018 by the individuals mentioned in the source code history' $(git grep -l 'Copyright')
