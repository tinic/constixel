#include <cstdio>

#include "sixel_tools.h"

int main() {

    sixel::image<sixel::format_1bit, 256, 256> image;
    image.plot(128,128,0);

    int a = 0;
    printf("%p\n", &a);
} 
