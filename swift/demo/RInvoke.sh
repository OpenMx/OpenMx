#!/bin/bash

export R_SCRIPT=$1
shift
export R_SWIFT_ARGS="$*"

R CMD BATCH --vanilla $R_SCRIPT


