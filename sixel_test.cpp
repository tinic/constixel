#include <cstdio>

#include "sixel_tools.h"

int main() {

    sixel::image<sixel::format_1bit, 256, 256> image;
    image.clear();
    image.fillrect(64,64,128,128,1);
    image.sixel([](uint8_t ch){
        putc(ch,stdout);
    });
} 
