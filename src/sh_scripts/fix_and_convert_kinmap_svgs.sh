#!/bin/bash

# Loop through all .svg files in the current directory
for file in *.svg; do
  echo "Processing $file"

  # Add xmlns:xlink if it's missing from the <svg> tag
  sed -i '/<svg /{
    /xmlns:xlink=/! s/<svg /<svg xmlns:xlink="http:\/\/www.w3.org\/1999\/xlink" /
  }' "$file"

  # Get the base filename without the extension
  base="${file%.svg}"

  # Convert to PNG with Inkscape
  inkscape -w 1024 -h 1024 "$file" -o "${base}.png"
done