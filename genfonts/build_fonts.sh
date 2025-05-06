#!/usr/bin/env bash
set -euo pipefail

FONT_DIR="${1:-./otf}"
SIZES=(12 18 24 32 48)

mkdir -p fontbm/build
cmake -DCMAKE_BUILD_TYPE=Release fontbm/build
cmake --build fontbm/build

mkdir -p build
cmake -DCMAKE_BUILD_TYPE=Release build
cmake --build build

for font in "$FONT_DIR"/*.{ttf,otf}; do
  [ -e "$font" ] || continue
  base="${font##*/}"
  base="${base%.*}"
  for size in "${SIZES[@]}"; do
    ./fontbm/build/fontbm --font-file "$font" --font-size "$size" --output "build/${base}_${size}" --texture-crop-width --texture-crop-height --tabular-numbers
    ./build/genfont --font-file "build/${base}_${size}.fnt" --out-dir ../fonts
    ./fontbm/build/fontbm --font-file "$font" --font-size "$size" --output "build/${base}_${size}" --monochrome --texture-crop-width --texture-crop-height --tabular-numbers
    ./build/genfont --font-file "build/${base}_${size}.fnt" --out-dir ../fonts
  done
done
