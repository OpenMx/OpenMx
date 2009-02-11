#!/bin/bash

export R_SCRIPT=$1
export R_DYNLIB="/scratch/projects/tg/SIDGrid/usr/lib/"
shift
export R_SWIFT_ARGS="$*"

R CMD BATCH --vanilla $R_SCRIPT


