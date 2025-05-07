/*
MIT License

Copyright (c) 2025 Tinic Uro

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../constixel.h"
#include "../fonts/sf_compact_display_bold_32_mono.h"
#include "../fonts/sf_compact_display_bold_48_mono.h"
#include "../fonts/sf_compact_display_medium_48_mono.h"
#include "../fonts/sf_mono_bold_48_mono.h"
#include "../fonts/sf_mono_regular_18_mono.h"

#if __GNUC__
#pragma GCC diagnostic push

#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma GCC diagnostic ignored "-Wpadded"

#if __clang__
#pragma GCC diagnostic ignored "-Wimplicit-int-conversion"
#pragma GCC diagnostic ignored "-Wshorten-64-to-32"
#pragma GCC diagnostic ignored "-Wsuggest-destructor-override"
#pragma GCC diagnostic ignored "-Wweak-vtables"
#endif  // #if __clang__

#include "../genfonts/fontbm/src/external/lodepng/lodepng.h"

#pragma GCC diagnostic pop

#pragma GCC diagnostic ignored "-Wunused-function"
#endif  // #if __GNUC__

#if 0
static constexpr void hexdump(const uint8_t* data, std::size_t len) {
    for (std::size_t offset = 0; offset < len; offset += 16) {
        std::print("{:08x}  ", offset);

        for (std::size_t i = 0; i < 16; ++i) {
            if (offset + i < len)
                std::print("{:02x} ", data[offset + i]);
            else
                std::print("   ");
            if (i == 7) std::print(" ");
        }

        std::print("|");
        for (std::size_t i = 0; i < 16; ++i) {
            if (offset + i < len) {
                auto c = static_cast<unsigned char>(data[offset + i]);
                std::print("{}", std::isprint(c) ? static_cast<char>(c) : '.');
            } else {
                std::print(" ");
            }
        }
        std::print("|\n");
    }
}
#endif  // #if 0

#if 1
static constexpr std::string test0() {
    std::string out{};
    out.reserve(1UL<<20);
    constixel::image<constixel::format_1bit, 32, 31> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fill_rect(c * 2, c * 2, 8, 8, c & 1);
    }
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test1() {
    std::string out{};
    out.reserve(1UL<<20);
    constixel::image<constixel::format_2bit, 112, 64> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fill_rect(c * 4 - 16, c * 4 - 16, 44, 31, static_cast<uint8_t>(c));
    }
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::string test2() {
    std::string out{};
    out.reserve(1UL<<20);
    constixel::image<constixel::format_4bit, 128, 127> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fill_rect(c * 4 - 16, c * 4 - 16, 44, 31, static_cast<uint8_t>(c));
    }
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test3() {
    std::string out{};
    out.reserve(1UL<<20);
    constixel::image<constixel::format_1bit, 127, 101> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fill_circle(c, c, 256 - c * 16, static_cast<uint8_t>(c));
    }
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test4() {
    std::string out{};
    out.reserve(1UL<<20);
    constixel::image<constixel::format_2bit, 111, 95> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fill_circle(c, c, 256 - c * 16, static_cast<uint8_t>(c));
    }
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test5() {
    std::string out{};
    out.reserve(1UL<<20);
    constixel::image<constixel::format_4bit, 7, 128, 8> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fill_circle(c, c, 256 - c * 16, static_cast<uint8_t>(c));
    }
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test6() {
    std::string out{};
    out.reserve(1UL<<20);
    constixel::image<constixel::format_1bit, 127, 101> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.line(16, 16, c * 8, 64, static_cast<uint8_t>(c), static_cast<uint8_t>(c));
    }
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test7() {
    std::string out{};
    out.reserve(1UL<<20);
    constixel::image<constixel::format_2bit, 111, 95> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.line(-8, -8, c * 8, 64, static_cast<uint8_t>(c), static_cast<uint8_t>(c));
    }
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test8() {
    std::string out{};
    out.reserve(1UL<<20);
    constixel::image<constixel::format_4bit, 7, 7, 8> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.line(-8, -8, c * 8, 64, static_cast<uint8_t>(c), static_cast<uint8_t>(c));
    }
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test9() {
    std::string out{};
    out.reserve(1UL<<20);
    constixel::image<constixel::format_8bit, 64, 64> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.line(-8, -8, c * 8, 48, static_cast<uint8_t>(c), static_cast<uint8_t>(c));
    }
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test10() {
    constixel::image<constixel::format_1bit, 256, 256, 1> image;
    image.fill_rect(0, 0, 256, 256, 2);
    std::string out{};
    out.reserve(1UL<<20);
    image.png([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::array<char, 8192> gen_separator() {
    std::array<char, 8192> sixel{};
    size_t count = 0;
    constixel::image<constixel::format_8bit, 1024, 1> image;
    for (int32_t x = 0; x < 16; x++) {
        image.line(x * 64, 0, x * 64 + 64, 0, static_cast<uint8_t>(constixel::color::GREY_RAMP_STOP - uint32_t(x)));
    }
    image.sixel([&sixel, &count](char ch) mutable {
        sixel[count++] = ch;
    });
    sixel[count++] = 0;
    return sixel;
}

static constexpr void separator() {
    static auto sixel = gen_separator();
    puts(sixel.data());
}

template <typename T>
void draw_image_cut(const std::vector<uint8_t> &rgbaimage, int32_t w, int32_t h) {
    static T image;
    image.blit_RGBA(0, 0, w, h, rgbaimage.data(), w, h, w * 4);
    image.sixel_to_cout();
    printf("%d-bit %dpx %dpx cut\n", int(image.bit_depth()), int(image.width()), int(image.height()));
    separator();
}

template <typename T>
void draw_image_diffused(const std::vector<uint8_t> &rgbaimage, int32_t w, int32_t h) {
    static T image;
    image.blit_RGBA_diffused(0, 0, w, h, rgbaimage.data(), w, h, w * 4);
    image.sixel_to_cout();
    printf("%d-bit %dpx %dpx diffused\n", int(image.bit_depth()), int(image.width()), int(image.height()));
    separator();
}

template <typename T>
void draw_image_linear(const std::vector<uint8_t> &rgbaimage, int32_t w, int32_t h) {
    static T image;
    image.blit_RGBA_diffused_linear(0, 0, w, h, rgbaimage.data(), w, h, w * 4);
    image.sixel_to_cout();
    printf("%d-bit %dpx %dpx diffused linear\n", int(image.bit_depth()), int(image.width()), int(image.height()));
    separator();
}

template <typename T>
void draw_palette() {
    static T image;
    for (int32_t y = 0; y < 16; y++) {
        for (int32_t x = 0; x < 16; x++) {
            image.fill_rect(x, y, 1, 1, static_cast<uint8_t>(y * 16 + x));
        }
    }
    std::string out;
    out.reserve(1UL<<20);
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    puts(out.c_str());
    printf("%d-bit %dpx %dpx\n", int(image.bit_depth()), int(image.width()), int(image.height()));
    separator();
}

template <typename T>
void draw_rgb() {
    static T image;
    std::array<uint32_t, 65536> rgb{};
    for (uint32_t y = 0; y < 256; y++) {
        for (uint32_t x = 0; x < 256; x++) {
            rgb.data()[y * 256 + x] = ((255 - y) << 16) | (x << 0) | (y << 8);
        }
    }
    {
        image.blit_RGBA(0, 0, 256, 256, reinterpret_cast<const uint8_t *>(rgb.data()), 256, 256, 256 * 4);
        std::string out;
        out.reserve(1UL<<20);
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx cut\n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
    {
        image.blit_RGBA_diffused(0, 0, 256, 256, reinterpret_cast<const uint8_t *>(rgb.data()), 256, 256, 256 * 4);
        std::string out;
        out.reserve(1UL<<20);
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx diffused\n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
    {
        image.blit_RGBA_diffused_linear(0, 0, 256, 256, reinterpret_cast<const uint8_t *>(rgb.data()), 256, 256, 256 * 4);
        std::string out;
        out.reserve(1UL<<20);
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx diffused linear\n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
}

template <typename T>
void round_trip() {
    static T image;
    std::array<uint32_t, 65536> rgb{};
    for (uint32_t y = 0; y < 256; y++) {
        for (uint32_t x = 0; x < 256; x++) {
            rgb.data()[y * 256 + x] = ((255 - y) << 16) | (x << 0) | (y << 8);
        }
    }
    image.blit_RGBA_diffused_linear(0, 0, 256, 256, reinterpret_cast<const uint8_t *>(rgb.data()), 256, 256, 256 * 4);
    auto rgba8 = image.RGBA_uint8();
    image.blit_RGBA(0, 0, 256, 256, rgba8.data(), 256, 256, 256 * 4);
    std::string out;
    out.reserve(1UL<<20);
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    puts(out.c_str());
    printf("%d-bit %dpx %dpx round trip\n", int(image.bit_depth()), int(image.width()), int(image.height()));
    separator();
}

template <typename T, size_t I>
void draw_functions() {
    puts("\033[2J\033[H\0337");
    static T image;
    image.clear();
    for (int32_t c = 0; c < static_cast<int32_t>(I); c++) {
        image.fill_rect(16 + c * 37, c * 32, 128, 128, static_cast<uint8_t>(c));
        std::string out("\0338");
        out.reserve(1UL<<20);
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx fill_rect\n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
    for (int32_t c = 0; c < static_cast<int32_t>(I); c++) {
        image.line(16, 16, 64 + c * 42, 700, static_cast<uint8_t>(c), static_cast<uint8_t>(c));
        std::string out("\0338");
        out.reserve(1UL<<20);
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx line\n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
    for (int32_t c = 0; c < static_cast<int32_t>(I); c++) {
        image.fill_circle(600, 384, 256 - c * 16, static_cast<uint8_t>(c));
        std::string out("\0338");
        out.reserve(1UL<<20);
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx circle\n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
}

#endif  // #if 0

int main() {
#if 1
    static_assert(test0().size() == 276);
    static_assert(test1().size() == 952);
    static_assert(test2().size() == 1591);
    static_assert(test3().size() == 1738);
    static_assert(test4().size() == 1716);
    static_assert(test5().size() == 2200);
    static_assert(test6().size() == 1393);
    static_assert(test7().size() == 2372);
    static_assert(test8().size() == 333);
    static_assert(test9().size() == 6037);
    static_assert(test10().size() == 12893);
#endif  // #if 0

#if 1
    puts("\033[2J\033[H\0337");
    puts(test0().c_str());
    puts("\n");
    puts(test1().c_str());
    puts("\n");
    puts(test2().c_str());
    puts("\n");
    puts(test3().c_str());
    puts("\n");
    puts(test4().c_str());
    puts("\n");
    puts(test5().c_str());
    puts("\n");
    puts(test6().c_str());
    puts("\n");
    puts(test7().c_str());
    puts("\n");
    puts(test8().c_str());
    puts("\n");
    puts(test9().c_str());
    puts("\n");
#endif  // #if 1

#if 1
    draw_functions<constixel::image<constixel::format_1bit, 768, 768>, 32>();
    draw_functions<constixel::image<constixel::format_2bit, 768, 768>, 32>();
    draw_functions<constixel::image<constixel::format_4bit, 768, 768>, 32>();
    draw_functions<constixel::image<constixel::format_8bit, 768, 768>, 32>();
    puts("\n");
#endif  // #if 0

#if 1
    draw_palette<constixel::image<constixel::format_1bit, 16, 16, 32>>();
    draw_palette<constixel::image<constixel::format_2bit, 16, 16, 32>>();
    draw_palette<constixel::image<constixel::format_4bit, 16, 16, 32>>();
    draw_palette<constixel::image<constixel::format_8bit, 16, 16, 32>>();
    puts("\n");
#endif  // #if 1

#if 1
    std::vector<uint8_t> rgbaimage;
    uint32_t w = 0;
    uint32_t h = 0;
    constexpr int32_t ow = 1024;
    constexpr int32_t oh = 902;
    if (lodepng::decode(rgbaimage, w, h, "../../media/larikeet.png") == 0) {
        draw_image_cut<constixel::image<constixel::format_1bit, ow, oh>>(rgbaimage, int32_t(w), int32_t(h));
        draw_image_cut<constixel::image<constixel::format_2bit, ow, oh>>(rgbaimage, int32_t(w), int32_t(h));
        draw_image_cut<constixel::image<constixel::format_4bit, ow, oh>>(rgbaimage, int32_t(w), int32_t(h));
        draw_image_cut<constixel::image<constixel::format_8bit, ow, oh>>(rgbaimage, int32_t(w), int32_t(h));

        draw_image_diffused<constixel::image<constixel::format_1bit, ow, oh>>(rgbaimage, int32_t(w), int32_t(h));
        draw_image_diffused<constixel::image<constixel::format_2bit, ow, oh>>(rgbaimage, int32_t(w), int32_t(h));
        draw_image_diffused<constixel::image<constixel::format_4bit, ow, oh>>(rgbaimage, int32_t(w), int32_t(h));
        draw_image_diffused<constixel::image<constixel::format_8bit, ow, oh>>(rgbaimage, int32_t(w), int32_t(h));

        draw_image_linear<constixel::image<constixel::format_1bit, ow, oh>>(rgbaimage, int32_t(w), int32_t(h));
        draw_image_linear<constixel::image<constixel::format_2bit, ow, oh>>(rgbaimage, int32_t(w), int32_t(h));
        draw_image_linear<constixel::image<constixel::format_4bit, ow, oh>>(rgbaimage, int32_t(w), int32_t(h));
        draw_image_linear<constixel::image<constixel::format_8bit, ow, oh>>(rgbaimage, int32_t(w), int32_t(h));
    }
#endif  // #if 1

#if 1
    draw_rgb<constixel::image<constixel::format_1bit, 256, 256, 1>>();
    draw_rgb<constixel::image<constixel::format_2bit, 256, 256, 1>>();
    draw_rgb<constixel::image<constixel::format_4bit, 256, 256, 1>>();
    draw_rgb<constixel::image<constixel::format_8bit, 256, 256, 1>>();
#endif  // #if 1

#if 1
    round_trip<constixel::image<constixel::format_1bit, 256, 256, 1>>();
    round_trip<constixel::image<constixel::format_2bit, 256, 256, 1>>();
    round_trip<constixel::image<constixel::format_4bit, 256, 256, 1>>();
    round_trip<constixel::image<constixel::format_8bit, 256, 256, 1>>();
#endif  // #if 1

#if 1
    {
        // FIXME!
        constixel::image<constixel::format_2bit, 1024, 1024, 1> image;
        image.draw_string_mono<constixel::sf_compact_display_bold_48_mono>(10, 20, "HELLO WORLD!", 1);
        image.sixel_to_cout();
    }
#endif  // #if 0

#if 1
    {
        constixel::image<constixel::format_8bit, 1024, 1024, 1> image;
        image.draw_string_mono<constixel::sf_compact_display_bold_48_mono>(0, 0, "Pack my box with five dozen liquor jugs", 2);
        image.draw_string_mono<constixel::sf_compact_display_bold_32_mono>(0, 50, "Pack my box with five dozen liquor jugs", 3);
        image.draw_string_mono<constixel::sf_mono_regular_18_mono>(0, 100, "Pack my box with five dozen liquor jugs", 5);
        image.draw_string_mono<constixel::sf_mono_bold_48_mono>(0, 200, "Pack my box with five dozen liquor jugs", 6);
        image.sixel_to_cout();
    }
#endif  // #if 0

#if 1
    {
        constixel::image<constixel::format_8bit, 512, 312, 1> image;
        static constexpr std::array<const char *, 5> strings{"ABCDEFGHIJKLM", "NOPQRTSUVWXYZ", "abcdefghijklm", "nopqrstuvwxyz",
                                                             "1234567890&@.,?!'"
                                                             ""};
        for (size_t i = 0; i < strings.size(); i++) {
            uint8_t col = constixel::color::GREY_RAMP_STOP - static_cast<uint8_t>(i * 3);
            image.draw_string_mono<constixel::sf_compact_display_medium_48_mono>(16, 48 * static_cast<int32_t>(i) + 16, strings.at(i), col);
        }
        std::vector<char> out{};
        out.reserve(1UL<<20);
        image.png([&out](char ch) mutable {
            out.push_back(ch);
        });
        std::ofstream file("constixel.png");
        file.write(reinterpret_cast<const char *>(out.data()), static_cast<std::streamsize>(out.size()));
        image.sixel_to_cout();
    }
#endif  // #if 0
}

