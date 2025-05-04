#include <unistd.h>

#include <chrono>
#include <cstdio>
#include <iostream>
#include <mdspan>
#include <string>
#include <vector>

#include "../constixel.h"
#include "../fontbm/src/external/lodepng/lodepng.h"

#if 1
constexpr std::string test0() {
    std::string out{};
    constixel::image<constixel::format_1bit, 32, 31> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fillrect(c * 2, c * 2, 8, 8, c & 1);
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
        image.fillrect(c * 4 - 16, c * 4 - 16, 44, 31, static_cast<uint8_t>(c));
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
        image.fillrect(c * 4 - 16, c * 4 - 16, 44, 31, static_cast<uint8_t>(c));
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
        image.fillcircle(c, c, 256 - c * 16, static_cast<uint8_t>(c));
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
        image.fillcircle(c, c, 256 - c * 16, static_cast<uint8_t>(c));
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
        image.fillcircle(c, c, 256 - c * 16, static_cast<uint8_t>(c));
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

template <typename T>
void draw_image_cut(const std::vector<uint8_t> &rgbaimage, int32_t w, int32_t h) {
    static T image;
    std::string out{};
    image.blitRGBA(0, 0, w, h, rgbaimage.data(), w, h, w * 4);
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    puts(out.c_str());
}

template <typename T>
void draw_image_diffused(const std::vector<uint8_t> &rgbaimage, int32_t w, int32_t h) {
    static T image;
    std::string out{};
    image.blitRGBADiffused(0, 0, w, h, rgbaimage.data(), w, h, w * 4);
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    puts(out.c_str());
}

template <typename T>
void draw_image_linear(const std::vector<uint8_t> &rgbaimage, int32_t w, int32_t h) {
    static T image;
    std::string out{};
    image.blitRGBADiffusedLinear(0, 0, w, h, rgbaimage.data(), w, h, w * 4);
    image.sixel([&out](char ch) mutable {
        out.push_back(ch);
    });
    puts(out.c_str());
}

template <typename T>
void draw_palette() {
    static T image;
    for (int32_t y = 0; y < 16; y++) {
        for (int32_t x = 0; x < 16; x++) {
            image.fillrect(x, y, 1, 1, static_cast<uint8_t>(y * 16 + x));
        }
        std::string out;
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
    }
}

template <typename T>
void draw_rgb() {
    static T image;
    std::array<uint32_t, 65536> rgb{};
    auto m = std::mdspan(rgb.data(), 256, 256);
    for (uint32_t y = 0; y < 256; y++) {
        for (uint32_t x = 0; x < 256; x++) {
            m[x, y] = (x << 0) | (y << 8);
        }
    }
    {
        image.blitRGBA(0, 0, 256, 256, reinterpret_cast<const uint8_t *>(rgb.data()), 256, 256, 256 * 4);
        std::string out;
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
    }
    {
        image.blitRGBADiffused(0, 0, 256, 256, reinterpret_cast<const uint8_t *>(rgb.data()), 256, 256, 256 * 4);
        std::string out;
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
    }
    {
        image.blitRGBADiffusedLinear(0, 0, 256, 256, reinterpret_cast<const uint8_t *>(rgb.data()), 256, 256, 256 * 4);
        std::string out;
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
    }
}

template <typename T, size_t I>
void draw_functions() {
    puts("\033[2J\033[H\0337");
    static T image;
    image.clear();
    for (int32_t c = 0; c < I; c++) {
        image.fillrect(16 + c * 37, c * 32, 128, 128, static_cast<uint8_t>(c));
        std::string out("\0338");
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
    }
    for (int32_t c = 0; c < I; c++) {
        image.line(16, 16, 64 + c * 42, 700, static_cast<uint8_t>(c), static_cast<uint8_t>(c));
        std::string out("\0338");
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
    }
    for (int32_t c = 0; c < I; c++) {
        image.fillcircle(600, 384, 256 - c * 16, static_cast<uint8_t>(c));
        std::string out("\0338");
        image.sixel([&out](char ch) mutable {
            out.push_back(ch);
        });
        puts(out.c_str());
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
    puts("\033[H\0337");
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
    constexpr int32_t oh = 1024;
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
    draw_rgb<constixel::image<constixel::format_2bit, 256, 256, 1>>();
    draw_rgb<constixel::image<constixel::format_4bit, 256, 256, 1>>();
    draw_rgb<constixel::image<constixel::format_8bit, 256, 256, 1>>();
#endif  // #if 1
}
