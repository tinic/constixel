#!/usr/bin/env bash

#
# You will need libfreetype-dev and libharfbuzz-dev
#
# genfonts.cpp wants gcc-14 or newer due to use of std::print/std::format
# 

set -euo pipefail
FONT_DIR="${1:-./bm437}"
CHARS=(32-126)

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
