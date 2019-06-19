#!/usr/bin/env bash

while IFS= read -r -d '' file
do
  echo "Rendering graph ${file} into ${file}.png"
  dot -Tpng -o "$file.png" "$file"
done < <(find . -name "*.dot" -print0)
