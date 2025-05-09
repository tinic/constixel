#!/usr/bin/env bash

#
# You will need libfreetype-dev and libharfbuzz-dev
#
# genfonts.cpp wants gcc-14 or newer due to use of std::print/std::format
# 

set -euo pipefail

ROMAN="32-128"
HIRAGANA="12352-12447"
KATAKANA="12448-12543"
PUNCTUATION="12288-12351,65288,65289"
CHARS="${ROMAN},${HIRAGANA},${KATAKANA},${PUNCTUATION}"
YOYOFILE="yoyo_utf8.txt"
FONT_DIR="${1:-./jp}"
SIZES=(24 32 48)

echo ${CHARS}

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
    ./fontbm/build/fontbm --font-file "$font" --font-size "$size" --output "build/${base}_jp_${size}" --chars 32-128,12352-12447,12448-12543,12288-12351 --chars-file ${YOYOFILE} --texture-crop-width --texture-crop-height --tabular-numbers
    ./build/genfont --font-file "build/${base}_jp_${size}.fnt" --out-dir ../fonts
    ./fontbm/build/fontbm --font-file "$font" --font-size "$size" --output "build/${base}_jp_${size}" --chars 32-128,12352-12447,12448-12543,12288-12351 --chars-file ${YOYOFILE} --monochrome --texture-crop-width --texture-crop-height --tabular-numbers
    ./build/genfont --font-file "build/${base}_jp_${size}.fnt" --out-dir ../fonts
  done
done
