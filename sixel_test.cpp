#include <cstdio>
#include <vector>
#include <chrono>
#include <iostream>
#include <unistd.h>

#include "sixel.h"
#include "sixel_tools.h"

static size_t bytesCount = 0;

int main() {

#if 0
    static sixel::image<sixel::format_4bit, 1024, 256> image;
    sixel::progressbar<sixel::format_4bit, 1024, 256> bar(image);
    bar.start(16,1);
    for (float value = 0.0f; value < 1.0f; value += 0.001f) {
        bar.update(value);
        usleep(1000);
    }
    bar.end();
#endif  // #if 0

#if 1
    static sixel::image<sixel::format_4bit, 768, 768> image;
    printf("RAM required: %d bytes\n", int32_t(image.size()));
    image.clear();

    printf("\033[H\0337");
    for (int32_t c=0; c<16; c++) {
        image.fillrect(16+c*37,c*32,128,128,c);
        printf("\0338");
        image.sixel([](uint8_t ch){
            putc(ch,stdout);
            bytesCount++;
        });
    }
    for (int32_t c=0; c<16; c++) {
        image.line(16,16,64+c*42,700,c,c);
        printf("\0338");
        image.sixel([](uint8_t ch){
            putc(ch,stdout);
            bytesCount++;
        });
    }
    for (int32_t c=0; c<16; c++) {
        image.fillcircle(600,384,256-c*16,c);
        printf("\0338");
        image.sixel([](uint8_t ch){
            putc(ch,stdout);
            bytesCount++;
        });
    }

#if 0
    std::array<uint8_t, 1UL<<21> data;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t c = 0 ; c < 256; c++) {
        size_t d = 0;
        uint8_t *ptr = data.data();
        image.sixel([ptr](uint8_t ch) mutable {
            *ptr++ = ch;
        });
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d0 = t1 - t0;
    std::cout << "time = " << d0.count() << "\n";
#endif // #if 0

    printf("Transfer bytes: %d bytes\n", int32_t(bytesCount));
#endif // #if 0
}

