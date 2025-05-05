# constixel

constixel is a minimalistic constexpr C++20+ graphics rendering library with the ability to output to a terminal using the sixel protocol.

## Primary features and goals

- Completely constexpr. All operations, including the sixel output stream can be generated during compilation.
- No dynamic allocations. The backbuffer and interal data structures can live as global static variables.
- Minimalistic interface and single header implementation so it can be used without fuzz in any modern C++ project.
- 1, 2, 4 and 8bit palette based back buffers for minimal memory usage. Reasonable standard palettes are provided.
- Blit raw 32-bit RGBA image buffers into the palette backed back buffer (with or without dithering).
- Conversion into a RGBA buffer when needed.
- Simple fillrect, line and fillcircle drawing functions.
- Various other simple image manipulation operations.

> [!NOTE]
> This library is not designed for high fidelity graphics generation and should be more thought of a utility library.

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
    return 0;
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
    return 0;
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

The following image formats are available. [Width] is the width in pixels, [Height] is the height in pixels. [Scale] is an optional parameter to specify a scale factor for the sixel output only. For instance setting this to 4 would scale the output of the sixel stream by a factor of 4.

```c++
    constixel::image<constixel::format_1bit, [Width], [Height], [Scale=1]>
    constixel::image<constixel::format_2bit, [Width], [Height], [Scale=1]>
    constixel::image<constixel::format_4bit, [Width], [Height], [Scale=1]>
    constixel::image<constixel::format_8bit, [Width], [Height], [Scale=1]>
```

The most important member function of image:

```c++
class image {

    // Colors in the fixed internal palette
    enum color : uint8_t {
        BLACK_TRANSPARENT = 0,
        BLACK_OPAQUE = 16,
        WHITE = 1,
        RED = 2,
        GREEN = 3,
        BLUE = 4,
        YELLOW = 5,
        CYAN = 6,
        MAGENTA = 7,
        GREY_80 = 8,
        GREY_60 = 9,
        GREY_40 = 10,
        GREY_20 = 11,
        DARK_RED = 12,
        DARK_GREEN = 13,
        DARK_BLUE = 14,
        DARK_YELLOW = 15

        GREY_RAMP_START = 16,
        GREY_RAMP_STOP = GREY_RAMP_START + 15,

        RED_LUMA_RAMP_START = 32,
        RED_LUMA_RAMP_STOP = RED_LUMA_RAMP_START + 15,

        GREEN_LUMA_RAMP_START = 48,
        GREEN_LUMA_RAMP_STOP = GREEN_LUMA_RAMP_START + 15,

        BLUE_LUMA_RAMP_START = 64,
        BLUE_LUMA_RAMP_STOP = BLUE_LUMA_RAMP_START + 15,

        YELLOW_LUMA_RAMP_START = 80,
        YELLOW_LUMA_RAMP_STOP = YELLOW_LUMA_RAMP_START + 15,

        CYAN_LUMA_RAMP_START = 96,
        CYAN_LUMA_RAMP_STOP = CYAN_LUMA_RAMP_START + 15,

        MAGENTA_LUMA_RAMP_START = 112,
        MAGENTA_LUMA_RAMP_STOP = MAGENTA_LUMA_RAMP_START + 15,
    };

    // Size in bytes of the image buffer
    int32_t size();

    // Width in pixel of the image buffer
    int32_t width();

    // Height in pixels of the image buffer
    int32_t height();

    // Return a reference to the internal pixel buffer
    std::array<uint8_t, T<W, H, S>::image_size> &dataRef();

    // Return a clone of the internal pixel buffer
    std::array<uint8_t, T<W, H, S>::image_size> clone();

    // Clear the image, i.e. set everything to color 0
    void clear();

    // Copy another image into this instance. Overwrites the contents, no compositing occurs.
    void copy(const std::array<uint8_t, T<W, H, S>::image_size> &src);

    // Draw a line
    void line(int32_t x0, int32_t y0, int32_t x1, int32_t y1, uint32_t col, uint32_t width = 1, bool clip = true);
    void line(constixel::rect<int32_t> &l, uint32_t col, bool clip = true);

    // Draw a filled rectangle
    void fillrect(int32_t x, int32_t y, int32_t w, int32_t h, uint32_t col, bool clip = true);
    void fillrect(onst constixel::rect<int32_t> &r, uint32_t col, bool clip = true);

    // Draw a filled circle
    void fillcircle(int32_t x, int32_t y, int32_t r, uint32_t col, bool clip = true);

    // Get a populated RGBA buffer with the contents of this instance.
    // Color 0 is special and will be set to 0x0000000 in the returned buffer, 
    // while all the other colors will be converted to a 0xffBBGGRR format.
    std::array<uint32_t, W * H> RGBA();

    // Blit an RGBA (little endian) buffer into this instance while color are quantizied to the internal palette. 
    // NOTE: This is a slow operation due to the brute force color quantization.
    void blitRGBA(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride);
    void blitRGBA(const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride);

    // Blit an RGBA (little endian) buffer into this instance using line diffusion for better quality.
    // NOTE: This is a slow operation due to the brute force color quantization.
    void blitRGBADiffused(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride);
    void blitRGBADiffused(const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride);

    // Blit an RGBA buffer into this instance using line diffusion in linear color space for best quality.
    // NOTE: This is a very slow operation due to the brute force color quantization.
    void blitRGBADiffusedLinear(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride);
    void blitRGBADiffusedLinear(const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride);

    // Convert the current instance into a sixel stream. Provide a lambda function in the form of:
    //
    // image.sixel([](char ch) {
    //    [do something with the character]
    // });
    //
    // Optionally you can provide a rectangle to get a portion of the image only.
    //
    void sixel(F &&charOut, bool preserveBackground = true);
    void sixel(F &&charOut, const constixel::rect<int32_t> &r, bool preserveBackground = true);

}
```
