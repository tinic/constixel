#include <cstdio>

#include "sixel_tools.h"

static size_t bytesCount = 0;

int main() {
    static sixel::image<sixel::format_4bit, 768, 768> image;
    printf("RAM required: %d bytes\n", int32_t(image.size()));
    image.clear();
    for (int32_t c=0; c<16; c++) {
        image.fillrect(16+c*37,c*32,128,128,c);
    }
    for (int32_t c=0; c<16; c++) {
        image.line(16,16,64+c*42,700,c,c);
    }
    for (int32_t c=0; c<16; c++) {
        image.fillcircle(600,384,256-c*16,c);
    }
    image.sixel([](uint8_t ch){
        putc(ch,stdout);
        bytesCount++;
    });
    printf("Transfer bytes: %d bytes\n", int32_t(bytesCount));
}

