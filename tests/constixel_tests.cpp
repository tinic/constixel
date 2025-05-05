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
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "../constixel.h"
#include "../genfonts/fontbm/src/external/lodepng/lodepng.h"

#if 1
constexpr std::string test0() {
    std::string out{};
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

constexpr std::string test1() {
    std::string out{};
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

constexpr std::string test3() {
    std::string out{};
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

constexpr std::string test4() {
    std::string out{};
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

constexpr std::string test5() {
    std::string out{};
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

constexpr std::string test6() {
    std::string out{};
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

constexpr std::string test7() {
    std::string out{};
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

constexpr std::string test8() {
    std::string out{};
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

constexpr std::string test9() {
    std::string out{};
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

template <typename A, typename B>
constexpr int const_memcmp(const A *lhs, const B *rhs, size_t count) {
    const unsigned char *l = reinterpret_cast<const unsigned char *>(lhs);
    const unsigned char *r = reinterpret_cast<const unsigned char *>(rhs);

    for (size_t i = 0; i < count; ++i) {
        if (l[i] != r[i]) {
            return (l[i] < r[i]) ? -1 : 1;
        }
    }
    return 0;
}

constexpr std::array<char, 8192> gen_separator() {
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

void separator() {
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
    for (int32_t c = 0; c < I; c++) {
        image.fill_rect(16 + c * 37, c * 32, 128, 128, static_cast<uint8_t>(c));
        std::string out("\0338");
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx fill_rect\n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
    for (int32_t c = 0; c < I; c++) {
        image.line(16, 16, 64 + c * 42, 700, static_cast<uint8_t>(c), static_cast<uint8_t>(c));
        std::string out("\0338");
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx line\n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
    for (int32_t c = 0; c < I; c++) {
        image.fill_circle(600, 384, 256 - c * 16, static_cast<uint8_t>(c));
        std::string out("\0338");
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
#if 0
    static_assert(test0().size() == 276);
    static_assert(test1().size() == 952);
    static_assert(test2().size() == 1591);
    static_assert(test3().size() == 1738);
    static_assert(test4().size() == 1716);
    static_assert(test5().size() == 2200);
    static_assert(test6().size() == 1393);
    static_assert(test7().size() == 2419);
    static_assert(test8().size() == 333);
    static_assert(test9().size() == 6121);
#endif  // #if 0

#if 0
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

#if 0
    draw_functions<constixel::image<constixel::format_1bit, 768, 768>, 32>();
    draw_functions<constixel::image<constixel::format_2bit, 768, 768>, 32>();
    draw_functions<constixel::image<constixel::format_4bit, 768, 768>, 32>();
    draw_functions<constixel::image<constixel::format_8bit, 768, 768>, 32>();
    puts("\n");
#endif  // #if 0

#if 0
    draw_palette<constixel::image<constixel::format_1bit, 16, 16, 32>>();
    draw_palette<constixel::image<constixel::format_2bit, 16, 16, 32>>();
    draw_palette<constixel::image<constixel::format_4bit, 16, 16, 32>>();
    draw_palette<constixel::image<constixel::format_8bit, 16, 16, 32>>();
    puts("\n");
#endif  // #if 1

#if 0
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

#if 0
    draw_rgb<constixel::image<constixel::format_1bit, 256, 256, 1>>();
    draw_rgb<constixel::image<constixel::format_2bit, 256, 256, 1>>();
    draw_rgb<constixel::image<constixel::format_4bit, 256, 256, 1>>();
    draw_rgb<constixel::image<constixel::format_8bit, 256, 256, 1>>();
#endif  // #if 1

#if 0
    round_trip<constixel::image<constixel::format_1bit, 256, 256, 1>>();
    round_trip<constixel::image<constixel::format_2bit, 256, 256, 1>>();
    round_trip<constixel::image<constixel::format_4bit, 256, 256, 1>>();
    round_trip<constixel::image<constixel::format_8bit, 256, 256, 1>>();
#endif  // #if 1

    constixel::phash<96> hash(
        {constixel::entry{32, 0},   constixel::entry{65, 1},       constixel::entry{66, 2},    constixel::entry{67, 3},       constixel::entry{68, 4},
         constixel::entry{69, 5},   constixel::entry{70, 6},       constixel::entry{71, 7},    constixel::entry{72, 8},       constixel::entry{73, 9},
         constixel::entry{74, 10},  constixel::entry{75, 11},      constixel::entry{76, 12},   constixel::entry{77, 13},      constixel::entry{78, 14},
         constixel::entry{79, 15},  constixel::entry{80, 16},      constixel::entry{81, 17},   constixel::entry{82, 18},      constixel::entry{83, 19},
         constixel::entry{84, 20},  constixel::entry{85, 21},      constixel::entry{86, 22},   constixel::entry{87, 23},      constixel::entry{88, 24},
         constixel::entry{89, 25},  constixel::entry{90, 26},      constixel::entry{97, 27},   constixel::entry{98, 28},      constixel::entry{99, 29},
         constixel::entry{100, 30}, constixel::entry{101, 31},     constixel::entry{102, 32},  constixel::entry{103, 33},     constixel::entry{104, 34},
         constixel::entry{105, 35}, constixel::entry{106, 36},     constixel::entry{107, 37},  constixel::entry{108, 38},     constixel::entry{109, 39},
         constixel::entry{110, 40}, constixel::entry{111, 41},     constixel::entry{112, 42},  constixel::entry{113, 43},     constixel::entry{114, 44},
         constixel::entry{115, 45}, constixel::entry{116, 46},     constixel::entry{117, 47},  constixel::entry{118, 48},     constixel::entry{119, 49},
         constixel::entry{120, 50}, constixel::entry{121, 51},     constixel::entry{122, 52},  constixel::entry{36, 53},      constixel::entry{48, 54},
         constixel::entry{49, 55},  constixel::entry{50, 56},      constixel::entry{51, 57},   constixel::entry{52, 58},      constixel::entry{53, 59},
         constixel::entry{54, 60},  constixel::entry{55, 61},      constixel::entry{56, 62},   constixel::entry{57, 63},      constixel::entry{38, 64},
         constixel::entry{33, 65},  constixel::entry{63, 66},      constixel::entry{40, 67},   constixel::entry{41, 68},      constixel::entry{91, 69},
         constixel::entry{93, 70},  constixel::entry{123, 71},     constixel::entry{125, 72},  constixel::entry{64, 73},      constixel::entry{35, 74},
         constixel::entry{47, 75},  constixel::entry{124, 76},     constixel::entry{92, 77},   constixel::entry{45, 78},      constixel::entry{39, 79},
         constixel::entry{34, 80},  constixel::entry{44, 81},      constixel::entry{46, 82},   constixel::entry{58, 83},      constixel::entry{59, 84},
         constixel::entry{60, 85},  constixel::entry{62, 86},      constixel::entry{61, 87},   constixel::entry{43, 88},      constixel::entry{126, 89},
         constixel::entry{95, 90},  constixel::entry{94, 91},      constixel::entry{42, 92},   constixel::entry{178, 93},     constixel::entry{37, 94},
         constixel::entry{96, 95}});

    printf("%08x\n", hash.scramble);
}
