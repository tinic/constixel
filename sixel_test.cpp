#include <cstdio>

#include "sixel_tools.h"

int main() {
    static sixel::image<sixel::format_4bit, 256, 256> image;
    image.clear();
    for (int32_t c=0; c<16; c++) {
        image.fillrect(16+c*16,c*16,64,64,c);
    }
    image.sixel([](uint8_t ch){
        putc(ch,stdout);
    });
}

