#!/usr/bin/env bash

#
# Fixed font building script 
#

set -euo pipefail

# Default values
FONT_DIR="./otf"
SIZES=(12 18 24 32 48)
CHARSET="ascii"
CUSTOM_CHARS=""
VERBOSE=false

# Predefined character sets
declare -A CHARSETS=(
    ["ascii"]="32-126"
    ["latin1"]="32-255"
    ["symbols"]="32-126,8192-8303,8364,8482,8226,8230,8212,8211,8216-8217,8220-8221"
    ["extended"]="32-126,160-255,8192-8303,8364,8482,8226,8230,8212,8211,8216-8217,8220-8221"
)

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --font-dir)
            FONT_DIR="$2"
            shift 2
            ;;
        --sizes)
            IFS=',' read -ra SIZES <<< "$2"
            shift 2
            ;;
        --charset)
            CHARSET="$2"
            shift 2
            ;;
        --chars)
            CUSTOM_CHARS="$2"
            shift 2
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        --help)
            echo "Usage: $0 [--font-dir DIR] [--sizes SIZES] [--charset PRESET] [--chars RANGES] [--verbose]"
            echo "Character sets: ascii, latin1, symbols, extended"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

# Set character ranges
if [[ "$CHARSET" == "custom" ]]; then
    if [[ -z "$CUSTOM_CHARS" ]]; then
        echo "Error: --charset custom requires --chars option" >&2
        exit 1
    fi
    CHARS="$CUSTOM_CHARS"
elif [[ -n "${CHARSETS[$CHARSET]:-}" ]]; then
    CHARS="${CHARSETS[$CHARSET]}"
else
    echo "Error: Unknown charset '$CHARSET'" >&2
    exit 1
fi

# Show configuration
echo "Font Building Configuration:"
echo "  Font directory: $FONT_DIR"
echo "  Font sizes: ${SIZES[*]}"
echo "  Character set: $CHARSET"
echo "  Character ranges: $CHARS"
echo

# Validate font directory
if [[ ! -d "$FONT_DIR" ]]; then
    echo "Error: Font directory '$FONT_DIR' does not exist" >&2
    exit 1
fi

# Build tools if needed
if [[ ! -f "fontbm/build/fontbm" ]]; then
    echo "Building fontbm..."
    git submodule update --init --recursive
    mkdir -p fontbm/build
    (cd fontbm/build && cmake -DCMAKE_BUILD_TYPE=Release .. && cmake --build .)
    echo "fontbm built successfully"
fi

if [[ ! -f "build/genfont" ]]; then
    echo "Building genfont..."
    mkdir -p build
    (cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && cmake --build .)
    echo "genfont built successfully"
fi

# Create array of font files
font_files=()
for ext in otf ttf OTF TTF; do
    for font in "$FONT_DIR"/*."$ext"; do
        [[ -f "$font" ]] && font_files+=("$font")
    done
done

if [[ ${#font_files[@]} -eq 0 ]]; then
    echo "Error: No font files found in $FONT_DIR"
    exit 1
fi

echo "Found ${#font_files[@]} font files to process"
echo

# Process each font
for i in "${!font_files[@]}"; do
    font="${font_files[$i]}"
    base="$(basename "$font")"
    base="${base%.*}"
    
    echo "Processing font $((i+1))/${#font_files[@]}: $base"
    
    for size in "${SIZES[@]}"; do
        if [[ "$VERBOSE" == "true" ]]; then
            echo "  Size $size (antialiased)..."
        fi
        
        # Antialiased version
        ./fontbm/build/fontbm \
            --font-file "$font" \
            --font-size "$size" \
            --chars "$CHARS" \
            --output "build/${base}_${size}" \
            --kerning-pairs extended \
            --force-auto-hinter \
            --extra-info \
            --texture-crop-width \
            --texture-crop-height \
            --tabular-numbers
        
        ./build/genfont \
            --font-file "build/${base}_${size}.fnt" \
            --out-dir ../fonts
        
        if [[ "$VERBOSE" == "true" ]]; then
            echo "  Size $size (monochrome)..."
        fi
        
        # Monochrome version
        ./fontbm/build/fontbm \
            --font-file "$font" \
            --font-size "$size" \
            --chars "$CHARS" \
            --output "build/${base}_${size}" \
            --kerning-pairs extended \
            --monochrome \
            --force-auto-hinter \
            --extra-info \
            --texture-crop-width \
            --texture-crop-height \
            --tabular-numbers
        
        ./build/genfont \
            --font-file "build/${base}_${size}.fnt" \
            --out-dir ../fonts
    done
    
    if [[ "$VERBOSE" == "false" ]]; then
        echo "  Generated ${#SIZES[@]} sizes (antialiased + monochrome)"
    fi
done

echo
echo "Font generation complete!"
echo "Generated ${#font_files[@]} font families in ${#SIZES[@]} sizes each"
echo "Output directory: ../fonts"