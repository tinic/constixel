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
shopt -s nullglob

declare -A meta=(
  [file01:name]=CD-IconsPC.ttf                                 [file01:size]=20       [file01:flags]="--all-chars --monochrome"
  [file02:name]=PixelOperator.ttf                              [file02:size]=16       [file02:flags]="--all-chars --monochrome"
  [file03:name]=PixelOperator-Bold.ttf                         [file03:size]=16       [file03:flags]="--all-chars --monochrome"
  [file04:name]=PixelOperator8-Bold.ttf                        [file04:size]=8        [file04:flags]="--all-chars --monochrome"
  [file05:name]=PixelOperator8.ttf                             [file05:size]=8        [file05:flags]="--all-chars --monochrome"
  [file06:name]=PixelOperatorHB.ttf                            [file06:size]=16       [file06:flags]="--all-chars --monochrome"
  [file07:name]=PixelOperatorHB8.ttf                           [file07:size]=16       [file07:flags]="--all-chars --monochrome"
  [file08:name]=PixelOperatorHBSC.ttf                          [file08:size]=16       [file08:flags]="--all-chars --monochrome"
  [file09:name]=PixelOperatorMono-Bold.ttf                     [file09:size]=16       [file09:flags]="--all-chars --monochrome"
  [file10:name]=PixelOperatorMono8-Bold.ttf                    [file10:size]=8        [file10:flags]="--all-chars --monochrome"
  [file11:name]=PixelOperatorMono8.ttf                         [file11:size]=8        [file11:flags]="--all-chars --monochrome"
  [file12:name]=PixelOperatorMonoHB.ttf                        [file12:size]=16       [file12:flags]="--all-chars --monochrome"
  [file13:name]=PixelOperatorMonoHB8.ttf                       [file13:size]=16       [file13:flags]="--all-chars --monochrome"
  [file14:name]=PixelOperatorSC-Bold.ttf                       [file14:size]=16       [file14:flags]="--all-chars --monochrome"
  [file15:name]=PixelOperatorSC.ttf                            [file15:size]=16       [file15:flags]="--all-chars --monochrome"

  [file15:name]=Bubbly-04B_30.TTF                              [file15:size]=17       [file15:flags]="--all-chars --monochrome"
  [file16:name]=ABSTRACT.TTF                                   [file16:size]=8        [file16:flags]="--all-chars --monochrome"
  [file18:name]=Fipps-Regular.otf                              [file18:size]=8        [file18:flags]="--all-chars --monochrome"
  [file19:name]=hachicro.TTF                                   [file19:size]=8        [file19:flags]="--all-chars --monochrome"
  [file20:name]=m12.TTF                                        [file20:size]=32       [file20:flags]="--all-chars --monochrome"
  [file21:name]=m42.TTF                                        [file21:size]=8        [file21:flags]="--all-chars --monochrome"
  [file22:name]=Mario-Kart-DS.ttf                              [file22:size]=16       [file22:flags]="--all-chars --monochrome"
  [file26:name]=upheavtt.ttf                                   [file26:size]=20       [file26:flags]="--all-chars --monochrome"
  [file27:name]=vhs-gothic.ttf                                 [file27:size]=16       [file27:flags]="--all-chars --monochrome"

  [file29:name]=PxPlus_Cordata_PPC-21.ttf                      [file29:size]=32       [file29:flags]="--all-chars --monochrome"
  [file30:name]=PxPlus_Cordata_PPC-400.ttf                     [file30:size]=32       [file30:flags]="--all-chars --monochrome"

  [file31:name]=retro-pixel-arcade.otf                         [file31:size]=8        [file31:flags]="--all-chars --monochrome"
  [file32:name]=retro-pixel-cute-mono.otf                      [file32:size]=11       [file32:flags]="--all-chars --monochrome"
  [file33:name]=retro-pixel-cute-prop.otf                      [file33:size]=11       [file33:flags]="--all-chars --monochrome"
  [file34:name]=retro-pixel-petty-5h.otf                       [file34:size]=5        [file34:flags]="--all-chars --monochrome"
  [file35:name]=retro-pixel-petty-5x5.otf                      [file35:size]=5        [file35:flags]="--all-chars --monochrome --texture-size 48x32"
  [file36:name]=retro-pixel-thick.otf                          [file36:size]=16       [file36:flags]="--all-chars --monochrome"

  [file37:name]=Bescii-Mono.otf                                [file37:size]=8        [file37:flags]="--all-chars --monochrome --texture-size 256x192"
  [file38:name]=SMW_Whole-Pixel_Spacing.otf                    [file38:size]=8        [file38:flags]="--all-chars --monochrome"
  [file39:name]=SMW_Monospace.otf                              [file39:size]=8        [file39:flags]="--all-chars --monochrome"
  [file40:name]=NotoSansSymbols2-Regular.ttf                   [file40:size]=32       [file40:flags]="--light-hinting --chars 32-126,0x1FBF0-0x1FBF9,0x2300-0x23ff,0x2190-0x21ff,0x25A0-0x25FF,0x1F300-0x1F5FF"
)

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

while IFS= read -r f; do
    font="${meta["$f:name"]}"
    base="${font##*/}"
    base="${base%.*}"
    flags=${meta[$f:flags]//\"}
    ./fontbm/build/fontbm --font-file "special/$font" --output "build/${base}" --font-size "${meta["$f:size"]}" ${flags} --extra-info --texture-crop-width --texture-crop-height
    ./build/genfont --font-file "build/${base}.fnt" --out-dir ../fonts
done < <(printf '%s\n' "${!meta[@]}" | cut -d: -f1 | sort -u)
