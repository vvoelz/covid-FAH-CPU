#!/bin/bash

j=0
for i in LAM*; do
    echo "$i is being moved to RUN$j"
    mv $i RUN$j
    j=$((j+1)); done
