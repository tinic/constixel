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

> [!NOTE]
> If you want to consteval either the sixel or image data you likely need to increase the constexpr ops limit. With g++ use '-fconstexpr-ops-limit=268435456' with clang use '-fconstexpr-steps=33554432'. The default limit in MSVC usually seems adequate.

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

## API

The following image formats are available. [Width] is the width in pixels, [Height] is the height in pixels. [Scale] is an optional paramter to specify a scale factor for the sixel output. For instance setting this to 4 would scale the output by a factor of 4.

```c++
    constixel::image<constixel::format_1bit, [Width], [Height], [Scale]>
    constixel::image<constixel::format_2bit, [Width], [Height], [Scale]>
    constixel::image<constixel::format_4bit, [Width], [Height], [Scale]>
    constixel::image<constixel::format_8bit, [Width], [Height], [Scale]>
```

The most important member function of image:

```c++
class image {
    // Size in bytes of the image buffer
    int32_t size();

    // Width in pixel of the image buffer
    int32_t width();

    // Height in pixels of the image buffer
    int32_t height();

    // Return a reference to the internal pixel buffer
    std::array<uint8_t, T<W, H, S>::image_size> &dataRef();

    // Return a clone of this instance
    std::array<uint8_t, T<W, H, S>::image_size> clone();

    // Clear the bitmap, i.e. set everything to color 0
    void clear();

    // Copy another image into this instance. Overwrites the contents, no compositing occurs.
    void copy(const std::array<uint8_t, T<W, H, S>::image_size> &src);

    // Draw a line
    void line(int32_t x0, int32_t y0, int32_t x1, int32_t y1, uint32_t col, uint32_t width = 1, bool clip = true);

    // Draw a filled rectangle
    void fillrect(int32_t x, int32_t y, int32_t w, int32_t h, uint32_t col, bool clip = true);

    // Draw a filled circle
    void fillcircle(int32_t x, int32_t y, int32_t r, uint32_t col, bool clip = true);

    // Get a populated RGBA buffer with the contents of this instance.
    // Color 0 is special will be set to 0 in the returned buffer, while all the other colors which be converted with the alpha channel set to 0xFF.
    std::array<uint32_t, W * H> RGBA();

    // Blit an RGBA buffer into this instance. This is a slow operation.
    void blitRGBA(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride);

    // Blit an RGBA buffer into this instance using line diffusion. This is a slow operation.
    void blitRGBADiffused(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride);

    // Blit an RGBA buffer into this instance using line diffusion in linear color space. This is a very slow operation.
    void blitRGBADiffusedLinear(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride);

    // Convert the current instance in a sixel stream. Provide a lambda function in the form of:
    //
    // sixel([](char ch) mutable {
    //    [do something with the character]
    // });
    //
    void sixel(F &&charOut, const rect<int32_t> &r, bool preserveBackground = true);

}
```
