#include <unistd.h>

#include <chrono>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "../fontbm/src/external/lodepng/lodepng.h"
#include "./sixel_tools.h"
#include "sixel.h"

static size_t bytesCount = 0;

#if 1
constexpr std::string test0() {
    std::string out{};
    sixel::image<sixel::format_1bit, 32, 31> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fillrect(c * 2, c * 2, 8, 8, c & 1);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::string test1() {
    std::string out{};
    sixel::image<sixel::format_2bit, 112, 64> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fillrect(c * 4 - 16, c * 4 - 16, 44, 31, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::string test2() {
    std::string out{};
    sixel::image<sixel::format_4bit, 128, 127> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fillrect(c * 4 - 16, c * 4 - 16, 44, 31, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::string test3() {
    std::string out{};
    sixel::image<sixel::format_1bit, 127, 101> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fillcircle(c, c, 256 - c * 16, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::string test4() {
    std::string out{};
    sixel::image<sixel::format_2bit, 111, 95> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fillcircle(c, c, 256 - c * 16, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::string test5() {
    std::string out{};
    sixel::image<sixel::format_4bit, 7, 128, 8> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fillcircle(c, c, 256 - c * 16, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::string test6() {
    std::string out{};
    sixel::image<sixel::format_1bit, 127, 101> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.line(16, 16, c * 8, 64, c, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::string test7() {
    std::string out{};
    sixel::image<sixel::format_2bit, 111, 95> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.line(-8, -8, c * 8, 64, c, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::string test8() {
    std::string out{};
    sixel::image<sixel::format_4bit, 7, 7, 8> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.line(-8, -8, c * 8, 64, c, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::string test9() {
    std::string out{};
    sixel::image<sixel::format_8bit, 64, 64> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.line(-8, -8, c * 8, 48, c, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

template <typename T>
void draw_image_cut(const std::vector<uint8_t> &rgbaimage, uint32_t w, uint32_t h) {
    static T image;
    std::string out{};
    image.blitRGBA(0, 0, w, h, rgbaimage.data(), w, h, w * 4);
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
        bytesCount++;
    });
    puts(out.c_str());
}

template <typename T>
void draw_image_diffused(const std::vector<uint8_t> &rgbaimage, uint32_t w, uint32_t h) {
    static T image;
    std::string out{};
    image.blitRGBADiffused(0, 0, w, h, rgbaimage.data(), w, h, w * 4);
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
        bytesCount++;
    });
    puts(out.c_str());
}

template <typename T>
void draw_image_linear(const std::vector<uint8_t> &rgbaimage, uint32_t w, uint32_t h) {
    static T image;
    std::string out{};
    image.blitRGBADiffusedLinear(0, 0, w, h, rgbaimage.data(), w, h, w * 4);
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
        bytesCount++;
    });
    puts(out.c_str());
}

template <typename T>
void draw_palette() {
    static T image;
    for (int32_t y = 0; y < 16; y++) {
        for (int32_t x = 0; x < 16; x++) {
            image.fillrect(x, y, 1, 1, y * 16 + x);
        }
    }
    std::string out;
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
        bytesCount++;
    });
    puts(out.c_str());
}
#endif  // #if 0

int main() {
#if 0
    static_assert(test0().size() == 316);
    static_assert(test1().size() == 1221);
    static_assert(test2().size() == 2090);
    static_assert(test3().size() == 1646);
    static_assert(test4().size() == 1639);
    static_assert(test5().size() == 5828);
    static_assert(test6().size() == 1428);
    static_assert(test7().size() == 2630);
    static_assert(test8().size() == 471);
    static_assert(test9().size() == 5916);
#endif  // #if 0

#if 0
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
    printf("\033[H\0337");
    static sixel::image<sixel::format_4bit, 1024, 256> image0;
    sixel::progressbar<sixel::format_4bit, 1024, 256> bar(image0);
    bar.start(16, 1);
    for (float value = 0.0f; value < 1.0f; value += 0.001f) {
        bar.update(value);
        usleep(1000);
    }
    bar.end();
#endif  // #if 0

#if 0
    printf("\033[H\0337");
    static sixel::image<sixel::format_8bit, 768, 768> image1;
    printf("RAM required: %d bytes\n", int32_t(image1.size()));
    image1.clear();

    for (int32_t c = 0; c < 128; c++) {
        image1.fillrect(16 + c * 37, c * 32, 128, 128, c);
        std::string out("\0338");
        image1.sixel([&out](uint8_t ch) mutable {
            out.push_back(ch);
            bytesCount++;
        });
        puts(out.c_str());
    }
    for (int32_t c = 0; c < 128; c++) {
        image1.line(16, 16, 64 + c * 42, 700, c, c);
        std::string out("\0338");
        image1.sixel([&out](uint8_t ch) mutable {
            out.push_back(ch);
            bytesCount++;
        });
        puts(out.c_str());
    }
    for (int32_t c = 0; c < 128; c++) {
        image1.fillcircle(600, 384, 256 - c * 16, c);
        std::string out("\0338");
        image1.sixel([&out](uint8_t ch) mutable {
            out.push_back(ch);
            bytesCount++;
        });
        puts(out.c_str());
    }

    printf("Transfer bytes: %d bytes\n", int32_t(bytesCount));
#endif  // #if 0

#if 0
    draw_palette<sixel::image<sixel::format_1bit, 16, 16, 32>>();
    draw_palette<sixel::image<sixel::format_2bit, 16, 16, 32>>();
    draw_palette<sixel::image<sixel::format_4bit, 16, 16, 32>>();
    draw_palette<sixel::image<sixel::format_8bit, 16, 16, 32>>();
#endif  // #if 1

#if 0
    static_assert(sixel::image<sixel::format_2bit, 1, 1>().octree_memory_length() == sixel::image<sixel::format_2bit, 1, 1>().octree_used_length());
    static_assert(sixel::image<sixel::format_4bit, 1, 1>().octree_memory_length() == sixel::image<sixel::format_4bit, 1, 1>().octree_used_length());
    static_assert(sixel::image<sixel::format_8bit, 1, 1>().octree_memory_length() == sixel::image<sixel::format_8bit, 1, 1>().octree_used_length());
#endif  // #if 1

#if 1
    std::vector<uint8_t> rgbaimage;
    uint32_t w = 0;
    uint32_t h = 0;
    constexpr size_t ow = 1024;
    constexpr size_t oh = 1024;
    constexpr size_t sc = 1;
    if (lodepng::decode(rgbaimage, w, h, "../media/larikeet.png") == 0) {
        draw_image_cut<sixel::image<sixel::format_1bit, ow, oh>>(rgbaimage, w, h);
        draw_image_cut<sixel::image<sixel::format_2bit, ow, oh>>(rgbaimage, w, h);
        draw_image_cut<sixel::image<sixel::format_4bit, ow, oh>>(rgbaimage, w, h);
        draw_image_cut<sixel::image<sixel::format_8bit, ow, oh>>(rgbaimage, w, h);

        draw_image_diffused<sixel::image<sixel::format_1bit, ow, oh>>(rgbaimage, w, h);
        draw_image_diffused<sixel::image<sixel::format_2bit, ow, oh>>(rgbaimage, w, h);
        draw_image_diffused<sixel::image<sixel::format_4bit, ow, oh>>(rgbaimage, w, h);
        draw_image_diffused<sixel::image<sixel::format_8bit, ow, oh>>(rgbaimage, w, h);

        draw_image_linear<sixel::image<sixel::format_1bit, ow, oh>>(rgbaimage, w, h);
        draw_image_linear<sixel::image<sixel::format_2bit, ow, oh>>(rgbaimage, w, h);
        draw_image_linear<sixel::image<sixel::format_4bit, ow, oh>>(rgbaimage, w, h);
        draw_image_linear<sixel::image<sixel::format_8bit, ow, oh>>(rgbaimage, w, h);
    }
#endif  // #if 1
}
