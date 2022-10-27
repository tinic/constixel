#include <cstdio>

#include "sixel_tools.h"

int main() {

    sixel::image<sixel::format_1bit, 256, 256> image;
    image.clear();
    image.fillrect(0,0,16,16,1);
    image.sixel([](uint8_t ch){
        putc(ch,stdout);
    });
} 
