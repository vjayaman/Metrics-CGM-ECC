#!/usr/bin/env bash

max=10
for i in `seq 1 $max`
do
	echo "Iteration $i --------------------------------------------------------------------------"
	Rscript tests/testthat/check_eccs.R -t FALSE
done
