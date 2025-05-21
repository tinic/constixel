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
FONT_DIR="${1:-./takwolf}"

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
