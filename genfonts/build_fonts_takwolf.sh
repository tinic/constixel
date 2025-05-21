#!/usr/bin/env bash

#
# You will need gcc-14+, bash 4.0+, libfreetype-dev and libharfbuzz-dev
#

set -euo pipefail
FONT_DIR="${1:-./takwolf}"
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

for font in "$FONT_DIR"/*.{ttf,otf}; do
  [ -e "$font" ] || continue
  base="${font##*/}"
  base="${base%.*}"
  if [[ "$base" =~ ([0-9]+)px ]]; then
      height="${BASH_REMATCH[1]}"
  fi
  ./fontbm/build/fontbm --font-file "$font" --output "build/${base}" --font-size ${height} --all-chars --force-auto-hinter --monochrome --extra-info --texture-crop-width --texture-crop-height
  ./build/genfont --font-file "build/${base}.fnt" --out-dir ../fonts
done
