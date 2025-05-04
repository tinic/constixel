# constixel

constixel is a constexpr minimal C++20+ graphics library with the ability to output to a terminal using the sixel protocol.

## Primary features and goals

- Completely constexpr. All operations, including the sixel output stream can be generated during compilation.
- No dynamic allocations. The backbuffer and interal data structures can live as global static variables.
- Minimalistic interface and single header implementation so it can be used without fuzz in any modern C++ project.
- Ability to render text from fontbm generated font files. All constexpr safe.
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
- For sixel display you will need a capabable terminal. Windows Terminal, iTerm2 on MacOS and any Linux terminal will work.

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

