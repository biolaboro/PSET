#!/usr/bin/env bash

./app/pool.py snakemake -lag 3 &
shiny run --host 0.0.0.0 ./app/app.py

# kill all background processes
trap 'kill $(jobs -p)' EXIT 

# wait for all the background processes to quit before exit
wait -n 

echo "done"
