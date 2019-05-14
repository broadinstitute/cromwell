#!/usr/bin/env bash

for i in *.dot
do 
  dot -Tpng -o "$i.png" "$i"
done
