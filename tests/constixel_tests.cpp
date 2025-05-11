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
#include <memory>
#include <string>
#include <vector>

#define SLOW_TESTS 0
#define MAINLINE_TESTS 1
#define JP_TESTS 0

#include "constixel.h"

#if MAINLINE_TESTS
#include "fonts/ibmplexmono_bold_48_aa.h"
#include "fonts/ibmplexmono_bold_48_mono.h"
#include "fonts/ibmplexmono_regular_18_mono.h"
#include "fonts/ibmplexsans_bold_12_aa.h"
#include "fonts/ibmplexsans_bold_12_mono.h"
#include "fonts/ibmplexsans_bold_32_mono.h"
#include "fonts/ibmplexsans_bold_48_aa.h"
#include "fonts/ibmplexsans_bold_48_mono.h"
#include "fonts/ibmplexsans_medium_48_aa.h"
#include "fonts/ibmplexsans_medium_48_mono.h"
#endif  // #if MAINLINE_TESTS

#if JP_TESTS
#include "fonts/notosansjp_black_jp_24_aa.h"
#include "fonts/notosansjp_black_jp_24_mono.h"
#include "fonts/notosansjp_regular_jp_24_aa.h"
#include "fonts/notosansjp_regular_jp_24_mono.h"
#include "fonts/notosansjp_regular_jp_48_aa.h"
#include "fonts/notosansjp_regular_jp_48_mono.h"
#include "fonts/notosansjp_thin_jp_24_aa.h"
#include "fonts/notosansjp_thin_jp_24_mono.h"
#endif  // #ifdef JP_TESTS

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
#endif  // #if __GNUC__

#include "fontbm/src/external/lodepng/lodepng.h"

#if __GNUC__
#pragma GCC diagnostic pop
#pragma GCC diagnostic ignored "-Wunused-function"
#endif  // #if __GNUC__

#ifndef CMAKE_PROJECT_PATH
#define CMAKE_PROJECT_PATH "."
#endif  // #ifndef CMAKE_PROJECT_PATH

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

#if MAINLINE_TESTS
static constexpr std::string test0() {
    std::string out{};
    out.reserve(1UL << 18);
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
    out.reserve(1UL << 18);
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
    out.reserve(1UL << 18);
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
    out.reserve(1UL << 18);
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
    out.reserve(1UL << 18);
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
    out.reserve(1UL << 18);
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
    out.reserve(1UL << 18);
    constixel::image<constixel::format_1bit, 127, 101> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.draw_line(16, 16, c * 8, 64, static_cast<uint8_t>(c), static_cast<uint8_t>(c));
    }
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test7() {
    std::string out{};
    out.reserve(1UL << 18);
    constixel::image<constixel::format_2bit, 111, 95> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.draw_line(-8, -8, c * 8, 64, static_cast<uint8_t>(c), static_cast<uint8_t>(c));
    }
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test8() {
    std::string out{};
    out.reserve(1UL << 18);
    constixel::image<constixel::format_4bit, 7, 7, 8> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.draw_line(-8, -8, c * 8, 64, static_cast<uint8_t>(c), static_cast<uint8_t>(c));
    }
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test9() {
    std::string out{};
    out.reserve(1UL << 18);
    constixel::image<constixel::format_8bit, 64, 64> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.draw_line(-8, -8, c * 8, 48, static_cast<uint8_t>(c), static_cast<uint8_t>(c));
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
    out.reserve(1UL << 18);
    image.png([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test11() {
    constixel::image<constixel::format_8bit, 64, 64, 1> image;
    image.draw_string_aa<constixel::ibmplexsans_medium_48_aa>(0, 0, "ABCDEFGHIKLMNO", 1);
    std::string out{};
    out.reserve(1UL << 18);
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test12() {
    constixel::image<constixel::format_8bit, 64, 64, 1> image;
    image.draw_string_mono<constixel::ibmplexsans_medium_48_mono>(0, 0, "ABCDEFGHIKLMNO", 1);
    std::string out{};
    out.reserve(1UL << 18);
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test13() {
    constixel::image<constixel::format_4bit, 64, 64, 1> image;
    image.draw_string_aa<constixel::ibmplexsans_medium_48_aa>(0, 0, "ABCDEFGHIKLMNO", 1);
    std::string out{};
    out.reserve(1UL << 18);
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test14() {
    constixel::image<constixel::format_1bit, 64, 64, 1> image;
    image.draw_string_mono<constixel::ibmplexsans_medium_48_mono>(0, 0, "ABCDEFGHIKLMNO", 1);
    std::string out{};
    out.reserve(1UL << 18);
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static constexpr std::string test15() {
    constixel::image<constixel::format_2bit, 64, 64, 1> image;
    image.draw_string_mono<constixel::ibmplexsans_medium_48_mono>(0, 0, "ABCDEFGHIKLMNO", 1);
    std::string out{};
    out.reserve(1UL << 18);
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    return out;
}

static std::array<char, 8192> gen_separator() {
    std::array<char, 8192> sixel{};
    size_t count = 0;
    constixel::image<constixel::format_8bit, 1024, 1> image;
    for (int32_t x = 0; x < 16; x++) {
        image.draw_line(x * 64, 0, x * 64 + 64, 0, static_cast<uint8_t>(constixel::color::GREY_RAMP_STOP - uint32_t(x)));
    }
    image.sixel([&sixel, &count](char ch) mutable {
        sixel[count++] = ch;
    });
    sixel[count++] = 0;
    return sixel;
}

static void separator() {
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
    out.reserve(1UL << 18);
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
        out.reserve(1UL << 18);
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
        out.reserve(1UL << 18);
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
        out.reserve(1UL << 18);
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
    static std::array<uint32_t, 65536> rgb{};
    for (uint32_t y = 0; y < 256; y++) {
        for (uint32_t x = 0; x < 256; x++) {
            rgb.data()[y * 256 + x] = ((255 - y) << 16) | (x << 0) | (y << 8);
        }
    }
    image.blit_RGBA_diffused_linear(0, 0, 256, 256, reinterpret_cast<const uint8_t *>(rgb.data()), 256, 256, 256 * 4);
    auto rgba8 = image.RGBA_uint8();
    image.blit_RGBA(0, 0, 256, 256, rgba8.data(), 256, 256, 256 * 4);
    std::string out;
    out.reserve(1UL << 18);
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    puts(out.c_str());
    printf("%d-bit %dpx %dpx round trip\n", int(image.bit_depth()), int(image.width()), int(image.height()));
    separator();
}

template <typename T, size_t I>
void draw_functions() {
    static T image;
    image.clear();
    puts("\033[2J\033[H");
    for (int32_t c = 0; c < static_cast<int32_t>(I); c++) {
        image.fill_rect(16 + c * 37, c * 32, 128, 128, static_cast<uint8_t>(c));
        puts("\033[H");
        std::string out;
        out.reserve(1UL << 18);
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx fill_rect  \n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
    for (int32_t c = 0; c < static_cast<int32_t>(I); c++) {
        image.draw_line(16, 16, 64 + c * 42, 700, static_cast<uint8_t>(c), static_cast<uint8_t>(c));
        puts("\033[H");
        std::string out;
        out.reserve(1UL << 18);
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx line       \n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
    for (int32_t c = 0; c < static_cast<int32_t>(I); c++) {
        image.fill_circle(600, 384, 256 - c * 16, static_cast<uint8_t>(c));
        puts("\033[H");
        std::string out;
        out.reserve(1UL << 18);
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx circle     \n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
}

template <typename T, size_t I>
void draw_functions_aa() {
    static T image;
    image.clear();
    puts("\033[2J\033[H");
    for (int32_t c = 0; c < static_cast<int32_t>(I); c++) {
        image.fill_rect(16 + c * 37, c * 32, 128, 128, static_cast<uint8_t>(c));
        puts("\033[H");
        std::string out;
        out.reserve(1UL << 18);
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx fill_rect aa \n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
    for (int32_t c = 0; c < static_cast<int32_t>(I); c++) {
        image.draw_line_aa(16, 16, 64 + c * 42, 700, static_cast<uint8_t>(c));
        puts("\033[H");
        std::string out;
        out.reserve(1UL << 18);
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx line aa      \n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
    for (int32_t c = 0; c < static_cast<int32_t>(I); c++) {
        image.fill_circle_aa(600, 384, 256 - c * 16, static_cast<uint8_t>(c));
        puts("\033[H");
        std::string out;
        out.reserve(1UL << 18);
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
        printf("%d-bit %dpx %dpx circle aa    \n", int(image.bit_depth()), int(image.width()), int(image.height()));
        separator();
    }
}

template <typename FONT>
void print_sizeof_font() {
    printf("%s %s %d (%s): %d bytes\n", FONT::name, FONT::style, int(FONT::size), FONT::mono ? "monochrome" : "antialiased",
           int(sizeof(FONT::glyph_tree) + sizeof(FONT::glyph_bitmap) + sizeof(FONT::char_table)));
}

#endif  // #ifdef MAINLINE_TESTS

int main() {
#if SLOW_TESTS
    static_assert(test0().size() == 290);
    static_assert(test1().size() == 981);
    static_assert(test2().size() == 1714);
    static_assert(test3().size() == 1748);
    static_assert(test4().size() == 1722);
    static_assert(test5().size() == 2264);
    static_assert(test6().size() == 1390);
    static_assert(test7().size() == 2422);
    static_assert(test8().size() == 387);
    static_assert(test9().size() == 6098);
    static_assert(test10().size() == 12893);
    static_assert(test11().size() == 5812);
    static_assert(test12().size() == 4434);
    static_assert(test13().size() == 1430);
    static_assert(test14().size() == 579);
    static_assert(test15().size() == 606);
#endif  // #if 0

#if MAINLINE_TESTS
    puts("\033[3J\033[H");
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
    //    puts(test10().c_str());
    //    puts("\n");
    puts(test11().c_str());
    puts("\n");
    puts(test12().c_str());
    puts("\n");
    puts(test13().c_str());
    puts("\n");
    puts(test14().c_str());
    puts("\n");
    puts(test15().c_str());
    puts("\n");
#endif  // #if 1

#if MAINLINE_TESTS
    draw_functions<constixel::image<constixel::format_1bit, 768, 768>, 32>();
    draw_functions<constixel::image<constixel::format_2bit, 768, 768>, 32>();
    draw_functions<constixel::image<constixel::format_4bit, 768, 768>, 32>();
    draw_functions<constixel::image<constixel::format_8bit, 768, 768>, 32>();

    draw_functions<constixel::image<constixel::format_1bit, 768, 768, 1, true>, 32>();
    draw_functions<constixel::image<constixel::format_2bit, 768, 768, 1, true>, 32>();
    draw_functions<constixel::image<constixel::format_4bit, 768, 768, 1, true>, 32>();
    puts("\n");
#endif  // #if 0

#if MAINLINE_TESTS
    draw_functions_aa<constixel::image<constixel::format_4bit, 768, 768>, 32>();
    draw_functions_aa<constixel::image<constixel::format_8bit, 768, 768>, 32>();

    draw_functions_aa<constixel::image<constixel::format_4bit, 768, 768, 1, true>, 32>();
    draw_functions_aa<constixel::image<constixel::format_8bit, 768, 768, 1, true>, 32>();
    puts("\n");
#endif  // #if 0

#if MAINLINE_TESTS
    draw_palette<constixel::image<constixel::format_1bit, 16, 16, 32>>();
    draw_palette<constixel::image<constixel::format_2bit, 16, 16, 32>>();
    draw_palette<constixel::image<constixel::format_4bit, 16, 16, 32>>();
    draw_palette<constixel::image<constixel::format_8bit, 16, 16, 32>>();

    draw_palette<constixel::image<constixel::format_1bit, 16, 16, 32, true>>();
    draw_palette<constixel::image<constixel::format_2bit, 16, 16, 32, true>>();
    draw_palette<constixel::image<constixel::format_4bit, 16, 16, 32, true>>();
    draw_palette<constixel::image<constixel::format_8bit, 16, 16, 32, true>>();
    puts("\n");
#endif  // #if 1

#if MAINLINE_TESTS
    std::vector<uint8_t> rgbaimage;
    uint32_t w = 0;
    uint32_t h = 0;
    constexpr int32_t ow = 1024;
    constexpr int32_t oh = 902;

    std::filesystem::path p{CMAKE_PROJECT_PATH "/../media/larikeet.png"};
    if (lodepng::decode(rgbaimage, w, h, p.lexically_normal().string().c_str()) == 0) {
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

        draw_image_linear<constixel::image<constixel::format_1bit, ow, oh, 1, true>>(rgbaimage, int32_t(w), int32_t(h));
        draw_image_linear<constixel::image<constixel::format_2bit, ow, oh, 1, true>>(rgbaimage, int32_t(w), int32_t(h));
        draw_image_linear<constixel::image<constixel::format_4bit, ow, oh, 1, true>>(rgbaimage, int32_t(w), int32_t(h));
        draw_image_linear<constixel::image<constixel::format_8bit, ow, oh, 1, true>>(rgbaimage, int32_t(w), int32_t(h));
    }
#endif  // #if 1

#if MAINLINE_TESTS
    draw_rgb<constixel::image<constixel::format_1bit, 256, 256, 1>>();
    draw_rgb<constixel::image<constixel::format_2bit, 256, 256, 1>>();
    draw_rgb<constixel::image<constixel::format_4bit, 256, 256, 1>>();
    draw_rgb<constixel::image<constixel::format_8bit, 256, 256, 1>>();

    draw_rgb<constixel::image<constixel::format_1bit, 256, 256, 1, true>>();
    draw_rgb<constixel::image<constixel::format_2bit, 256, 256, 1, true>>();
    draw_rgb<constixel::image<constixel::format_4bit, 256, 256, 1, true>>();
    draw_rgb<constixel::image<constixel::format_8bit, 256, 256, 1, true>>();
#endif  // #if 1

#if MAINLINE_TESTS
    round_trip<constixel::image<constixel::format_1bit, 256, 256, 1>>();
    round_trip<constixel::image<constixel::format_2bit, 256, 256, 1>>();
    round_trip<constixel::image<constixel::format_4bit, 256, 256, 1>>();
    round_trip<constixel::image<constixel::format_8bit, 256, 256, 1>>();
#endif  // #if 1

#if MAINLINE_TESTS
    {
        auto image = std::make_unique<constixel::image<constixel::format_8bit, 1024, 1024, 1>>();
        image->draw_string_mono<constixel::ibmplexsans_bold_48_mono>(0, 0, "Pack my box with five dozen liquor jugs", 2);
        image->draw_string_mono<constixel::ibmplexsans_bold_32_mono>(0, 50, "Pack my box with five dozen liquor jugs", 3);
        image->draw_string_mono<constixel::ibmplexmono_regular_18_mono>(0, 100, "Pack my box with five dozen liquor jugs", 5);
        image->draw_string_mono<constixel::ibmplexmono_bold_48_mono>(0, 200, "Pack my box with five dozen liquor jugs", 6);
        image->sixel_to_cout();
    }
#endif  // #if 0

#if MAINLINE_TESTS
    {
        static constixel::image<constixel::format_8bit, 512, 312, 1> image;
        static constexpr std::array<const char *, 5> strings{"ABCDEFGHIJKLM", "NOPQRTSUVWXYZ", "abcdefghijklm", "nopqrstuvwxyz",
                                                             "1234567890&@.,?!'"
                                                             ""};
        for (size_t i = 0; i < strings.size(); i++) {
            uint8_t col = constixel::color::GREY_RAMP_STOP - static_cast<uint8_t>(i * 3);
            image.draw_string_mono<constixel::ibmplexsans_bold_48_mono>(16, 48 * static_cast<int32_t>(i) + 16, strings.at(i), col);
        }
        std::vector<char> out{};
        out.reserve(1UL << 18);
        image.png([&out](char ch) mutable {
            out.push_back(ch);
        });
        {
            std::ofstream file("constixel_mono.png", std::ios::binary);
            file.write(reinterpret_cast<const char *>(out.data()), static_cast<std::streamsize>(out.size()));
            image.sixel_to_cout();
        }

        image.clear();
        std::array<uint8_t, 5> cols = {0, 4, 3, 125, 1};
        for (size_t i = 0; i < strings.size(); i++) {
            uint8_t col = constixel::color::GREY_RAMP_STOP - static_cast<uint8_t>(i * 3);
            image.fill_rect(0, static_cast<int32_t>(i) * 52 + 22, 512, 52, cols.at(i));
            image.draw_string_aa<constixel::ibmplexsans_bold_48_aa>(16, 52 * static_cast<int32_t>(i) + 16, strings.at(i), col);
        }
        out.clear();
        image.png([&out](char ch) mutable {
            out.push_back(ch);
        });
        {
            std::ofstream file("constixel_aa.png", std::ios::binary);
            file.write(reinterpret_cast<const char *>(out.data()), static_cast<std::streamsize>(out.size()));
            image.sixel_to_cout();
        }
    }
#endif  // #if 0

#if JP_TESTS
    {
        static constixel::image<constixel::format_8bit, 2048, 128, 1> image;
        image.draw_string_aa<constixel::notosansjp_thin_jp_24_aa>(0, 0, "日本語の文章には、漢字(常用漢字)、カタカナ、ひらがな、そして句読点が含まれています。",
                                                                  1);
        image.draw_string_aa<constixel::notosansjp_black_jp_24_aa>(0, 64,
                                                                   "日本語の文章には、漢字(常用漢字)、カタカナ、ひらがな、そして句読点が含まれています。", 1);
        image.sixel_to_cout();
    }
    {
        static constixel::image<constixel::format_1bit, 2048, 128, 1> image;
        image.draw_string_mono<constixel::notosansjp_thin_jp_24_mono>(
            0, 0, "日本語の文章には、漢字(常用漢字)、カタカナ、ひらがな、そして句読点が含まれています。", 1);
        image.draw_string_mono<constixel::notosansjp_black_jp_24_mono>(
            0, 64, "日本語の文章には、漢字(常用漢字)、カタカナ、ひらがな、そして句読点が含まれています。", 1);
        image.sixel_to_cout();
    }
    {
        static constixel::image<constixel::format_1bit, 2048, 128, 1> image;
        image.draw_string_mono<constixel::notosansjp_regular_jp_48_mono>(
            0, 0, "日本語の文章には、漢字(常用漢字)、カタカナ、ひらがな、そして句読点が含まれています。", 1);
        image.draw_string_mono<constixel::notosansjp_regular_jp_24_mono>(
            0, 64, "日本語の文章には、漢字(常用漢字)、カタカナ、ひらがな、そして句読点が含まれています。", 1);
        image.sixel_to_cout();
    }
    {
        static constixel::image<constixel::format_8bit, 2048, 128, 1> image;
        image.draw_string_aa<constixel::notosansjp_regular_jp_48_aa>(0, 0,
                                                                     "日本語の文章には、漢字(常用漢字)、カタカナ、ひらがな、そして句読点が含まれています。", 1);
        image.draw_string_aa<constixel::notosansjp_regular_jp_24_aa>(0, 64,
                                                                     "日本語の文章には、漢字(常用漢字)、カタカナ、ひらがな、そして句読点が含まれています。", 1);
        image.sixel_to_cout();
    }
#endif  // #if 0

#if MAINLINE_TESTS
    {
        auto image = std::make_unique<constixel::image<constixel::format_8bit, 128, 256, 5>>();
        for (int32_t x = 0; x < 16; x++) {
            image->draw_string_aa<constixel::ibmplexsans_bold_48_aa>(x * 7, x * 7, "Pack my box with five dozen liquor jugs", static_cast<uint8_t>(x));
            image->draw_string_aa<constixel::ibmplexmono_bold_48_aa>(x * 7, x * 7 + 100, "Pack my box with five dozen liquor jugs",
                                                                     static_cast<uint8_t>(15 - x));
        }
        image->sixel_to_cout();
    }
#endif  // #if 0

#if MAINLINE_TESTS
    {
        auto image = std::make_unique<constixel::image<constixel::format_8bit, 768, 384, 1>>();
        image->fill_round_rect_aa(32, 32, 400, 100, 32, 3);
        image->fill_round_rect_aa(32, 200, 400, 64, 32, 1);
        image->sixel_to_cout();
    }
#endif  // #if 0

#if MAINLINE_TESTS
    {
        auto image = std::make_unique<constixel::image<constixel::format_8bit, 768, 384, 1>>();
        image->fill_round_rect(32, 32, 400, 100, 32, 3);
        image->fill_round_rect(32, 200, 400, 64, 32, 1);
        image->sixel_to_cout();
    }
#endif  // #if 0

#if MAINLINE_TESTS
    {
        auto image = std::make_unique<constixel::image<constixel::format_8bit, 48, 48, 16>>();

        image->fill_circle_aa(14, 7, 7, 2);
        image->fill_circle(14, 7, 7, 4);

        image->fill_circle_aa(4 + 24, 7, 4, 2);
        image->fill_circle(4 + 24, 7, 4, 4);

        image->fill_circle_aa(4 + 24 + 8, 7, 2, 2);
        image->fill_circle(4 + 24 + 8, 7, 2, 4);

        image->fill_circle_aa(4 + 24 + 8 + 4, 7, 1, 2);
        image->fill_circle(4 + 24 + 8 + 4, 7, 1, 4);

        image->fill_circle_aa(14, 23, 7, 2);
        image->fill_circle_aa(4 + 24, 23, 4, 2);
        image->fill_circle_aa(4 + 24 + 8, 23, 2, 2);
        image->fill_circle_aa(4 + 24 + 8 + 4, 23, 1, 2);

        image->fill_circle(14, 40, 7, 4);
        image->fill_circle(4 + 24, 40, 4, 4);
        image->fill_circle(4 + 24 + 8, 40, 2, 4);
        image->fill_circle(4 + 24 + 8 + 4, 40, 1, 4);

        image->sixel_to_cout();
    }
#endif  // #if 0

#if MAINLINE_TESTS
    {
        using image_type = constixel::image<constixel::format_8bit, 512, 512, 1, false>;

        auto image = std::make_unique<image_type>();
        image->draw_string_aa<constixel::ibmplexsans_bold_48_aa>(0, 0, "Pack my box with five dozen liquor jugs", 1);
        auto image_clone = std::make_unique<image_type>(image->clone());
        auto image_copy = std::make_unique<image_type>();
        image_copy->copy(*image_clone);
        image_copy->sixel_to_cout();
    }
#endif  // #if 0

#if MAINLINE_TESTS
    {
        auto image = std::make_unique<constixel::image<constixel::format_4bit, 512, 96, 1, true>>();
        image->fill_round_rect_aa(0, 0, 512, 96, 32, 15);
        int32_t sw = image->string_width<constixel::ibmplexsans_medium_48_aa>("ABCDEFGHIKLMNO");
        image->draw_string_aa<constixel::ibmplexsans_medium_48_aa>((512 - sw) / 2, 16, "ABCDEFGHIKLMNO", 0);
        image->sixel_to_cout();
    }
#endif  // #if 0

#if MAINLINE_TESTS
    {
        auto image = std::make_unique<constixel::image<constixel::format_4bit, 512, 96, 1, false>>();
        image->fill_round_rect_aa(0, 0, 512, 96, 32, 2);
        int32_t sw = image->string_width<constixel::ibmplexsans_medium_48_aa>("ABCDEFGHIKLMNO");
        image->draw_string_aa<constixel::ibmplexsans_medium_48_aa>((512 - sw) / 2, 16, "ABCDEFGHIKLMNO", 0);
        image->sixel_to_cout();
    }
#endif  // #if 0

#if MAINLINE_TESTS
    {
        auto image = std::make_unique<constixel::image<constixel::format_4bit, 512, 96, 1, false>>();
        image->fill_round_rect_aa(0, 0, 512, 96, 32, constixel::color::GREY_20);
        int32_t sw = image->string_width<constixel::ibmplexsans_medium_48_aa>("ABCDEFGHIKLMNO");
        image->draw_string_aa<constixel::ibmplexsans_medium_48_aa>((512 - sw) / 2, 16, "ABCDEFGHIKLMNO", 2);
        image->sixel_to_cout();
    }
#endif  // #if 0

#if MAINLINE_TESTS
    {
        auto image = std::make_unique<constixel::image<constixel::format_4bit, 512, 96, 1, true>>();
        image->fill_round_rect_aa(0, 0, 512, 96, 32, 3);
        int32_t sw = image->string_width<constixel::ibmplexsans_medium_48_aa>("ABCDEFGHIKLMNO");
        image->draw_string_aa<constixel::ibmplexsans_medium_48_aa>((512 - sw) / 2, 16, "ABCDEFGHIKLMNO", 15);
        image->sixel_to_cout();
    }
#endif  // #if 0

#if 0
    {
        auto image = std::make_unique<constixel::image<constixel::format_1bit, 512, 96, 1, true>>();
        image->fill_round_rect(0, 0, 512, 96, 32, 1);
        int32_t sw = image->string_width<constixel::ibmplexsans_medium_48_aa>("ABCDEFGHIKLMNO");
        image->draw_string_mono<constixel::ibmplexsans_medium_48_aa>((512-sw)/2, 16, "ABCDEFGHIKLMNO", 0);
        image->png_to_iterm();
    }
#endif  // #if 0

#if 0
    {
        auto image = std::make_unique<constixel::image<constixel::format_1bit, 512, 96, 1, true>>();
        image->fill_round_rect(0, 0, 512, 96, 32, 1);
        int32_t sw = image->string_width<constixel::ibmplexsans_medium_48_aa>("ABCDEFGHIKLMNO");
        image->draw_string_mono<constixel::ibmplexsans_medium_48_aa>((512-sw)/2, 16, "ABCDEFGHIKLMNO", 0);
        image->png_to_kitty();
    }
#endif  // #if 0

#if MAINLINE_TESTS
    {
        auto image = std::make_unique<constixel::image<constixel::format_8bit, 512, 512, 1, false>>();
        for (int32_t c = 0; c < 64; c++) {
            image->draw_line_aa(0, 8 * c, 16 * c, 512, static_cast<uint8_t>(c & 0xFF));
        }
        image->sixel_to_cout();
    }
    {
        auto image = std::make_unique<constixel::image<constixel::format_8bit, 512, 512, 1, true>>();
        for (int32_t c = 0; c < 64; c++) {
            image->draw_line_aa(0, 8 * c, 16 * c, 512, static_cast<uint8_t>((c * 2 + 128) & 0xFF));
        }
        image->sixel_to_cout();
    }
#endif  // #if 0

#if MAINLINE_TESTS
    {
        auto image = std::make_unique<constixel::image<constixel::format_8bit, 1024, 512, 1, false>>();

        for (int32_t c = 0; c < 8; c++) {
            const int32_t delta = 8;
            for (int32_t x = 0; x < 1024; x += delta) {
                float xf0 = float(x) / 1024.0f * 2.0f * 3.141f * float(c + 1) / 2.0f;
                float xf1 = float(x + delta) / 1024.0f * 2.0f * 3.141f * float(c + 1) / 2.0f;
                ;
                image->draw_line(x, int32_t(256.0f + std::sin(xf0) * 192.0f), x + delta, int32_t(256.0f + std::sin(xf1) * 192.0f), static_cast<uint8_t>(c), 1);
            }
        }
        image->sixel_to_cout();
    }
    {
        auto image = std::make_unique<constixel::image<constixel::format_8bit, 1024, 512, 1, false>>();
        for (int32_t c = 0; c < 8; c++) {
            const int32_t delta = 8;
            for (int32_t x = 0; x < 1024; x += delta) {
                float xf0 = float(x) / 1024.0f * 2.0f * 3.141f * float(c + 1) / 2.0f;
                float xf1 = float(x + delta) / 1024.0f * 2.0f * 3.141f * float(c + 1) / 2.0f;
                image->draw_line_aa(x, int32_t(256.0f + std::sin(xf0) * 192.0f), x + delta, int32_t(256.0f + std::sin(xf1) * 192.0f), static_cast<uint8_t>(c));
            }
        }
        image->sixel_to_cout();
    }
#endif  // #if MAINLINE_TESTS

#if MAINLINE_TESTS
    {
        using font_mono = constixel::ibmplexsans_bold_12_mono;
        using font_aa = constixel::ibmplexsans_bold_12_aa;
        using namespace constixel;

        static image<format_8bit, 192, 300, 5> i;

        auto draw_bg_for_rotate_test = [=]<typename FONT>(const char *str, int32_t x_pos, int32_t y_pos) mutable {
            int32_t s_wdh = i.string_width<FONT>(str);

            i.fill_rect(x_pos, y_pos, s_wdh, FONT::line_height, 2);
            i.fill_rect(x_pos - s_wdh, y_pos - FONT::line_height, s_wdh, FONT::line_height, 4);
            i.fill_rect(x_pos - FONT::line_height, y_pos, FONT::line_height, s_wdh, 9);
            i.fill_rect(x_pos, y_pos - s_wdh, FONT::line_height, s_wdh, 7);

            i.fill_rect(x_pos, y_pos + FONT::ascent, s_wdh, 1, 12);
            i.fill_rect(x_pos - s_wdh, y_pos - FONT::ascent - 1, s_wdh, 1, 12);
            i.fill_rect(x_pos - FONT::ascent - 1, y_pos, 1, s_wdh, 12);
            i.fill_rect(x_pos + FONT::ascent, y_pos - s_wdh, 1, s_wdh, 12);
        };

        auto str = std::string("Garfunkel");  // random_ascii_string();

        i.clear();

        int32_t y_pos = 96;
        int32_t x_pos = i.width() / 2;
        i.fill_rect(x_pos - 1, 0, 2, i.height(), 1);
        draw_bg_for_rotate_test.template operator()<font_mono>(str.c_str(), x_pos, y_pos);
        y_pos = 224;
        draw_bg_for_rotate_test.template operator()<font_aa>(str.c_str(), x_pos, y_pos);

        y_pos = 96;
        i.draw_string_mono<font_mono, false, text_rotation::DEGREE_0>(x_pos, y_pos, str.c_str(), color::WHITE);
        i.draw_string_mono<font_mono, false, text_rotation::DEGREE_180>(x_pos, y_pos, str.c_str(), color::WHITE);
        i.draw_string_mono<font_mono, false, text_rotation::DEGREE_90>(x_pos, y_pos, str.c_str(), color::WHITE);
        i.draw_string_mono<font_mono, false, text_rotation::DEGREE_270>(x_pos, y_pos, str.c_str(), color::WHITE);

        y_pos = 224;
        i.draw_string_aa<font_aa, false, text_rotation::DEGREE_0>(x_pos, y_pos, str.c_str(), color::WHITE);
        i.draw_string_aa<font_aa, false, text_rotation::DEGREE_180>(x_pos, y_pos, str.c_str(), color::WHITE);
        i.draw_string_aa<font_aa, false, text_rotation::DEGREE_90>(x_pos, y_pos, str.c_str(), color::WHITE);
        i.draw_string_aa<font_aa, false, text_rotation::DEGREE_270>(x_pos, y_pos, str.c_str(), color::WHITE);

        i.sixel_to_cout();
    }
#endif  // #if MAINLINE_TESTS

#if MAINLINE_TESTS
    print_sizeof_font<constixel::ibmplexsans_bold_48_aa>();
    print_sizeof_font<constixel::ibmplexsans_bold_32_mono>();
    print_sizeof_font<constixel::ibmplexsans_bold_48_mono>();
    print_sizeof_font<constixel::ibmplexsans_medium_48_aa>();
    print_sizeof_font<constixel::ibmplexsans_medium_48_mono>();
    print_sizeof_font<constixel::ibmplexmono_bold_48_aa>();
    print_sizeof_font<constixel::ibmplexmono_bold_48_mono>();
    print_sizeof_font<constixel::ibmplexmono_regular_18_mono>();
#endif  // #if MAINLINE_TESTS
}
