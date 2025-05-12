#!/usr/bin/env bash

#
# You will need libfreetype-dev and libharfbuzz-dev
#
# genfonts.cpp wants gcc-14 or newer due to use of std::print/std::format
# 

set -euo pipefail
shopt -s nullglob

declare -A meta=(
  [file1:name]=CD-IconsPC.ttf               [file1:size]=20        [file1:flags]=true
  [file2:name]=PixelOperator.ttf            [file2:size]=16        [file2:flags]=true
  [file3:name]=PixelOperator-Bold.ttf       [file3:size]=16        [file3:flags]=true
  [file4:name]=PixelOperator8-Bold.ttf      [file4:size]=8         [file4:flags]=true
  [file5:name]=PixelOperator8.ttf           [file5:size]=8         [file5:flags]=true
  [file6:name]=PixelOperatorHB.ttf          [file6:size]=16        [file6:flags]=true
  [file7:name]=PixelOperatorHB8.ttf         [file7:size]=16        [file7:flags]=true
  [file8:name]=PixelOperatorHBSC.ttf        [file8:size]=16        [file8:flags]=true
  [file9:name]=PixelOperatorMono-Bold.ttf   [file9:size]=16        [file9:flags]=true
  [file10:name]=PixelOperatorMono8-Bold.ttf [file10:size]=8        [file10:flags]=true
  [file11:name]=PixelOperatorMono8.ttf      [file11:size]=8        [file11:flags]=true
  [file12:name]=PixelOperatorMonoHB.ttf     [file12:size]=16       [file12:flags]=true
  [file13:name]=PixelOperatorMonoHB8.ttf    [file13:size]=16       [file13:flags]=true
  [file14:name]=PixelOperatorSC-Bold.ttf    [file14:size]=16       [file14:flags]=true
  [file15:name]=PixelOperatorSC.ttf         [file15:size]=16       [file15:flags]=true

  [file15:name]=Bubbly-04B_30.TTF           [file15:size]=17       [file15:flags]=true
  [file16:name]=ABSTRACT.TTF                [file16:size]=8        [file16:flags]=true
  [file17:name]=advanced_pixel_lcd-7.ttf    [file17:size]=22       [file17:flags]=true
  [file18:name]=Fipps-Regular.otf           [file18:size]=8        [file18:flags]=true
  [file19:name]=hachicro.TTF                [file19:size]=8        [file19:flags]=true
  [file20:name]=m12.TTF                     [file20:size]=32       [file20:flags]=true
  [file21:name]=m42.TTF                     [file21:size]=8        [file21:flags]=true
  [file22:name]=Mario-Kart-DS.ttf           [file22:size]=16       [file22:flags]=true
  [file23:name]=Paskowy.ttf                 [file23:size]=27       [file23:flags]=true
  [file24:name]=PixArrows.ttf               [file24:size]=10       [file24:flags]=true
  [file26:name]=upheavtt.ttf                [file26:size]=20       [file26:flags]=true
  [file27:name]=VCR_OSD_MONO.ttf            [file27:size]=21       [file27:flags]=true
  [file28:name]=vhs-gothic.ttf              [file28:size]=16       [file28:flags]=true

  [file29:name]=PxPlus_Cordata_PPC-21.ttf   [file29:size]=16       [file29:flags]=true
  [file30:name]=PxPlus_Cordata_PPC-400.ttf  [file30:size]=16       [file30:flags]=true
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
    ./fontbm/build/fontbm --font-file "special/$font" --output "build/${base}" --font-size "${meta["$f:size"]}" --all-chars --force-auto-hinter --monochrome --extra-info --texture-crop-width --texture-crop-height
    ./build/genfont --font-file "build/${base}.fnt" --out-dir ../fonts
done < <(printf '%s\n' "${!meta[@]}" | cut -d: -f1 | sort -u)
