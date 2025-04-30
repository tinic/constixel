#include "sixel.h"

#include <unistd.h>

#include <chrono>
#include <cstdio>
#include <iostream>
#include <vector>
#include <string>

#include "sixel_tools.h"

static size_t bytesCount = 0;

constexpr std::vector<uint8_t> test0() {
    std::vector<uint8_t> out{};
    sixel::image<sixel::format_1bit, 32, 31> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fillrect(16 + c * 37, c * 32, 128, 128, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::vector<uint8_t> test1() {
    std::vector<uint8_t> out{};
    sixel::image<sixel::format_2bit, 112, 64> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fillrect(16 + c * 37, c * 32, 128, 128, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::vector<uint8_t> test2() {
    std::vector<uint8_t> out{};
    sixel::image<sixel::format_4bit, 128, 127> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fillrect(16 + c * 37, c * 32, 128, 128, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::vector<uint8_t> test3() {
    std::vector<uint8_t> out{};
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

constexpr std::vector<uint8_t> test4() {
    std::vector<uint8_t> out{};
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

constexpr std::vector<uint8_t> test5() {
    std::vector<uint8_t> out{};
    sixel::image<sixel::format_4bit, 7, 7> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.fillcircle(c, c, 256 - c * 16, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::vector<uint8_t> test6() {
    std::vector<uint8_t> out{};
    sixel::image<sixel::format_1bit, 127, 101> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.line(16, 16, 64 + c * 42, 700, c, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::vector<uint8_t> test7() {
    std::vector<uint8_t> out{};
    sixel::image<sixel::format_2bit, 111, 95> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.line(16, 16, 64 + c * 42, 700, c, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::vector<uint8_t> test8() {
    std::vector<uint8_t> out{};
    sixel::image<sixel::format_4bit, 7, 7> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.line(0, 0, c * 42, 700, c, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

constexpr std::vector<uint8_t> test9() {
    std::vector<uint8_t> out{};
    sixel::image<sixel::format_8bit, 32, 32> image;
    image.clear();
    for (int32_t c = 0; c < 16; c++) {
        image.line(0, 0, c, 700, c, c);
    }
    image.sixel([&out](uint8_t ch) mutable {
        out.push_back(ch);
    });
    return out;
}

int main() {
#if 0
    static_assert(test0().size() == 144);
    static_assert(test1().size() == 273);
    static_assert(test2().size() == 663);
    static_assert(test3().size() == 1646);
    static_assert(test4().size() == 1638);
    static_assert(test5().size() == 271);
    static_assert(test6().size() == 1164);
    static_assert(test7().size() == 363);
    static_assert(test8().size() == 273);
    static_assert(test9().size() == 4080);
#endif  // #if 0

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

#if 1
    printf("\033[H\0337");

    static sixel::image<sixel::format_8bit, 768, 768> image1;
    printf("RAM required: %d bytes\n", int32_t(image1.size()));
    image1.clear();

    for (int32_t c = 0; c < 16; c++) {
        image1.fillrect(16 + c * 37, c * 32, 128, 128, c);
        std::string out("\0338");
        image1.sixel([&out](uint8_t ch) mutable {
            out.push_back(ch);
            bytesCount++;
        });
        puts(out.c_str());
    }
    for (int32_t c = 0; c < 16; c++) {
        image1.line(16, 16, 64 + c * 42, 700, c, c);
        std::string out("\0338");
        image1.sixel([&out](uint8_t ch) mutable {
            out.push_back(ch);
            bytesCount++;
        });
        puts(out.c_str());
    }
    for (int32_t c = 0; c < 16; c++) {
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
}
