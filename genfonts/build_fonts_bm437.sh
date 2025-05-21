#!/usr/bin/env bash

#
# You will need gcc-14+, bash 4.0+, libfreetype-dev and libharfbuzz-dev
#
# I strongly recomment to check out the options fontbm supports to adapt this script
# for your purposes.
#
# Of note is that this custom version of fontbm has extra features like the ability
# to specify --secondary-font-file which can be used to pull symbols/icons from if
# they don't exist in the primary font.
#

set -euo pipefail
FONT_DIR="${1:-./bm437}"
CHARS=(32-126) # comma separated ranges, accepts 0x1234 type hex values.

git submodule update --init --recursive

mkdir -p fontbm/build
cd fontbm/build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cd ../..

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cd ..

for font in "$FONT_DIR"/*.{ttf,ttf}; do
  [ -e "$font" ] || continue
  base="${font##*/}"
  base="${base%.*}"
  if [[ "$base" =~ _[0-9]+x([0-9]+)$ ]]; then
      height="${BASH_REMATCH[1]}"
  else
      height=32
  fi
  if [[ "$height" -eq 14 ]]; then
    height=16
    echo "Height 14->16"
  fi
  if [[ "$height" -eq 11 ]]; then
    height=12
    echo "Height 11->12"
  fi
  if [[ "$(echo "$base" | tr '[:upper:]' '[:lower:]')" == pxplus* ]]; then
    ./fontbm/build/fontbm --font-file "$font" --output "build/${base}" --font-size ${height} --all-chars --force-auto-hinter --monochrome --extra-info --texture-crop-width --texture-crop-height
    ./build/genfont --font-file "build/${base}.fnt" --out-dir ../fonts
  else
    ./fontbm/build/fontbm --font-file "$font" --output "build/${base}_full" --font-size ${height} --all-chars --force-auto-hinter --monochrome --extra-info --texture-crop-width --texture-crop-height
    ./build/genfont --font-file "build/${base}_full.fnt" --out-dir ../fonts
    ./fontbm/build/fontbm --font-file "$font" --output "build/${base}" --font-size ${height} --chars "${CHARS}" --force-auto-hinter --monochrome --extra-info --texture-crop-width --texture-crop-height
    ./build/genfont --font-file "build/${base}.fnt" --out-dir ../fonts
  fi
done
