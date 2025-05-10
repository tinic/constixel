#!/usr/bin/env bash

#
# You will need libfreetype-dev and libharfbuzz-dev
#
# genfonts.cpp wants gcc-14 or newer due to use of std::print/std::format
# 

set -euo pipefail

FONT_DIR="${1:-./otf}"
SIZES=(12 18 24 32 48)
CHARS=(32-127,65535)

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

for font in "$FONT_DIR"/*.{ttf,otf}; do
  [ -e "$font" ] || continue
  base="${font##*/}"
  base="${base%.*}"
  for size in "${SIZES[@]}"; do
    ./fontbm/build/fontbm --font-file "$font" --font-size "$size" --chars "${CHARS}" --output "build/${base}_${size}" --kerning-pairs extended --force-auto-hinter --extra-info --texture-crop-width --texture-crop-height --tabular-numbers
    ./build/genfont --font-file "build/${base}_${size}.fnt" --out-dir ../fonts
    ./fontbm/build/fontbm --font-file "$font" --font-size "$size" --chars "${CHARS}" --output "build/${base}_${size}" --kerning-pairs extended --monochrome --force-auto-hinter --extra-info --texture-crop-width --texture-crop-height --tabular-numbers
    ./build/genfont --font-file "build/${base}_${size}.fnt" --out-dir ../fonts
  done
done
