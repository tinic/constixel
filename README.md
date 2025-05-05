# constixel

constixel is a minimalistic constexpr C++20+ graphics rendering library with the ability to output to a terminal using the sixel protocol.

## Primary features and goals

- Completely constexpr. All operations, including the sixel output stream can be generated during compilation.
- No dynamic allocations. The backbuffer and interal data structures can live as global static variables.
- Minimalistic interface and single header implementation so it can be used without fuzz in any modern C++ project.
- 1, 2, 4 and 8bit palette based back buffers for minimal memory usage. Reasonable standard palettes are provided.
- Blit raw 32-bit RGBA image buffers into the palette backed back buffer (with or without dithering). 

## Applications

- Interface rendering on embedded devices.
- Send graphical data to your terminal through a serial connection from any embedded device.
- Add graphical output to unit tests.
- Programmatically render static graphical assets.
- Help debug dynamic memory handling issues in any complex C++ projects.
- ...
  
## Requirements

- C++20
- gcc 13 or newer, clang 16 or newer, MSVC 17 or newer
- For viewing the sixel image you will need a capabable terminal. Windows Terminal, iTerm2 on MacOS and any Linux terminal will work.

> [!NOTE]
> The Terminal app on MacOS does not support sixel, please use iTerm2.

## Minimal example

```c++
#include "constixel.h"

#include <iostream>

int main() {
    static constixel::image<constixel::format_8bit, 256, 256> image;

    for (int32_t y = 0; y < 16; y++) {
        for (int32_t x = 0; x < 16; x++) {
            image.fillrect(x * 16, y * 16, 16, 16, static_cast<uint8_t>(y * 16 + x));
        }
    }

    std::string out;
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    std::cout << out << std::endl;
}
```

![constixel](./media/constixel.jpg "Example in iTerm")

## Consteval sixel example

As std::vector can not escape consteval (yet), we use std::array. Output of this example should be "Actual byte size: 18537" and the sixel image. The binary will contain the evaluated sixel string.

Compile as such:

```bash
> g++ -fconstexpr-ops-limit=268435456 -std=c++23 constixel.cpp -o constixel -Os
```

```c++
#include "constixel.h"

#include <iostream>
#include <cstring>

consteval auto gen_sixel() {
    constixel::image<constixel::format_8bit, 256, 256> image;
    for (int32_t y = 0; y < 16; y++) {
        for (int32_t x = 0; x < 16; x++) {
            image.fillrect(x * 16, y * 16, 16, 16, static_cast<uint8_t>(y * 16 + x));
        }
    }

    std::array<char, 32767> sixel;
    char *ptr = sixel.data();
    image.sixel([&ptr](char ch) mutable {
        *ptr++ = ch;
    });
    *ptr++ = '\n';
    *ptr++ = 0;
    return sixel;
}

int main() {
    static const auto sixel = gen_sixel();
    std::cout << "Actual byte size: " << strlen(sixel.data()) << "\n";
    std::cout << sixel.data() << std::endl;
}
```

## Consteval embedded image data example

This example will consteval gen_image_1bit(), while dynamically generating the sixel string.

```bash
> g++ -fconstexpr-ops-limit=268435456 -std=c++23 constixel_1bit.cpp -o constixel_1bit -Os
```

```c++
#include "constixel.h"

#include <cstring>

consteval auto gen_image_1bit() {
    constixel::image<constixel::format_1bit, 256, 256, 1> image;
    for (int32_t y = 0; y < 16; y++) {
        for (int32_t x = 0; x < 16; x++) {
            image.fillrect(x * 16, y * 16, 16, 16, static_cast<uint8_t>(y + x) & 1);
        }
    }
    return image;
}

int main() {
    static const auto image = gen_image_1bit();
    printf("image width x height: %d %d x 1bit depth\n", int(image.width()), int(image.height()));
    printf("image instance byte size: %d\n", int(image.size()));
    size_t count = 0;
    image.sixel([&count](char ch) mutable {
        putc(ch, stdout);
        count++;
    });
    putc('\n', stdout);
    printf("sixel byte size: %d\n", int(count));
    return 0;
}
```

![constixel](./media/constixel_1bit.jpg "Example in iTerm")
