#include "constixel.hpp"

#include "fonts/fipps_regular_mono.hpp"
using font = constixel::fipps_regular_mono;

#include "fonts/ibmplexmono_bold_24_aa.hpp"
using fontaa = constixel::ibmplexmono_bold_24_aa;

#include <cstdio>
#include <format>
#include <limits>
#include <random>
#include <algorithm>

static int32_t random_int32() {
    static std::mt19937 gen(0xDEADBEEF);
    static std::uniform_int_distribution<int> dist(std::numeric_limits<int32_t>::min(),
                                                   std::numeric_limits<int32_t>::max());
    return dist(gen);
}

static int32_t random_limit() {
    static constexpr const std::array<int32_t, 34> limits {{

        0,
        0,
        200,
        300,
        -200,
        -666,
        262144,
        -262144,
        2147483647,
        -2147483647,

        std::numeric_limits<int32_t>::min(),
        std::numeric_limits<int32_t>::max(),
        std::numeric_limits<int16_t>::min(),
        std::numeric_limits<int16_t>::max(),
        std::numeric_limits<uint8_t>::min(),
        std::numeric_limits<uint8_t>::max(),
        static_cast<int32_t>(std::numeric_limits<uint32_t>::min()),
        static_cast<int32_t>(std::numeric_limits<uint32_t>::max()),
        std::numeric_limits<uint16_t>::min(),
        std::numeric_limits<uint16_t>::max(),
        std::numeric_limits<uint8_t>::min(),
        std::numeric_limits<uint8_t>::max(),

        std::numeric_limits<int32_t>::min()+1,
        std::numeric_limits<int32_t>::max()-1,
        std::numeric_limits<int16_t>::min()+1,
        std::numeric_limits<int16_t>::max()-1,
        std::numeric_limits<uint8_t>::min()+1,
        std::numeric_limits<uint8_t>::max()-1,
        static_cast<int32_t>(std::numeric_limits<uint32_t>::min()+1),
        static_cast<int32_t>(std::numeric_limits<uint32_t>::max()-1),
        std::numeric_limits<uint16_t>::min()+1,
        std::numeric_limits<uint16_t>::max()-1,
        std::numeric_limits<uint8_t>::min()+1,
        std::numeric_limits<uint8_t>::max()-1,
    }};

    return limits[static_cast<size_t>(random_int32()) % limits.size()];
}

template <class T, bool AA>
void fuzz_limits() {
    T image;
    for (size_t c = 0; c < 65536; c++) {
        image.draw_line(random_limit(),random_limit(),random_limit(),random_int32(),static_cast<uint8_t>(random_limit()&0xFF));
        image.draw_line(random_limit(),random_limit(),random_limit(),random_int32(),static_cast<uint8_t>(random_limit()&0xFF),random_limit());
        image.fill_circle(random_limit(),random_limit(),random_limit(),static_cast<uint8_t>(random_limit()&0xFF));
        image.fill_rect(random_limit(), random_limit(), random_limit(), random_limit(), static_cast<uint8_t>(random_limit()&0xFF));
        image.stroke_rect(random_limit(), random_limit(), random_limit(), random_int32(), static_cast<uint8_t>(random_limit()&0xFF));
        image.stroke_rect(random_limit(), random_limit(), random_limit(), random_int32(), static_cast<uint8_t>(random_limit()&0xFF), random_limit());
        image.fill_round_rect(random_limit(), random_limit(), random_limit(), random_limit(), random_limit(), static_cast<uint8_t>(random_limit()&0xFF));
        image.stroke_circle(random_limit(), random_limit(), random_limit(), static_cast<uint8_t>(random_limit()&0xFF),  random_limit());
        image.stroke_round_rect(random_limit(), random_limit(), random_limit(), random_limit(), random_limit(), static_cast<uint8_t>(random_limit()&0xFF), random_limit());
        image.plot(random_limit(), random_limit(), static_cast<uint8_t>(random_limit()&0xFF));
        image.plot({.x = random_limit(), .y = random_limit(), .col = static_cast<uint8_t>(random_limit()&0xFF)});
        image.draw_line({.x0 = random_limit(), .y0 = random_limit(), .x1 = random_limit(), .y1 = random_limit(), .col = static_cast<uint8_t>(random_limit()&0xFF), .sw = static_cast<float>(random_limit())});
        image.fill_rect({.x = random_limit(), .y = random_limit(), .w = random_limit(), .h = random_limit(), .col = static_cast<uint8_t>(random_limit()&0xFF)});
        image.stroke_rect({.x = random_limit(), .y = random_limit(), .w = random_limit(), .h = random_limit(), .col = static_cast<uint8_t>(random_limit()&0xFF), .sw = static_cast<int32_t>(random_limit())});
        image.fill_circle({.cx = random_limit(), .cy = random_limit(), .r = random_limit(), .col = static_cast<uint8_t>(random_limit()&0xFF)});
        image.stroke_circle({.cx = random_limit(), .cy = random_limit(), .r = random_limit(), .col = static_cast<uint8_t>(random_limit()&0xFF), .sw = static_cast<int32_t>(random_limit())});
        image.fill_round_rect({.x = random_limit(), .y = random_limit(), .w = random_limit(), .h = random_limit(), .r = random_limit(), .col = static_cast<uint8_t>(random_limit()&0xFF)});
        image.stroke_round_rect({.x = random_limit(), .y = random_limit(), .w = random_limit(), .h = random_limit(), .r = random_limit(), .col = static_cast<uint8_t>(random_limit()&0xFF), .sw = static_cast<int32_t>(random_limit())});
        image.template draw_string_mono<font>(random_limit(), random_limit(), "ABCabc", static_cast<uint8_t>(random_limit()&0xFF));
        image.template draw_string_mono<font, false, constixel::DEGREE_90>(random_limit(), random_limit(), "ABCabc", static_cast<uint8_t>(random_limit()&0xFF));
        image.template draw_string_mono<font, false, constixel::DEGREE_180>(random_limit(), random_limit(), "ABCabc", static_cast<uint8_t>(random_limit()&0xFF));
        image.template draw_string_mono<font, false, constixel::DEGREE_270>(random_limit(), random_limit(), "ABCabc", static_cast<uint8_t>(random_limit()&0xFF));
        image.template draw_string_centered_mono<font>(random_limit(), random_limit(), "ABCabc", static_cast<uint8_t>(random_limit()&0xFF));
        if constexpr (AA) {
            image.stroke_circle_aa(random_limit(), random_limit(), random_limit(), static_cast<uint8_t>(random_limit()&0xFF),  random_limit());
            image.stroke_round_rect_aa(random_limit(), random_limit(), random_limit(), random_limit(), random_limit(), static_cast<uint8_t>(random_limit()&0xFF), random_limit());
            image.draw_line_aa(random_limit(),random_limit(),random_limit(),random_limit(),static_cast<uint8_t>(random_limit()&0xFF));
            image.draw_line_aa(random_limit(),random_limit(),random_limit(),random_limit(),static_cast<uint8_t>(random_limit()&0xFF),random_limit());
            image.fill_circle_aa(random_limit(), random_limit(), random_limit(), static_cast<uint8_t>(random_limit()&0xFF));
            image.fill_round_rect_aa(random_limit(), random_limit(), random_limit(), random_limit(), random_limit(), static_cast<uint8_t>(random_limit()&0xFF));
            image.fill_rect(random_limit(), random_limit(), random_limit(), random_limit(), [](float u, float v, float au, float av) -> std::array<float, 4> { return {u, v, au, av}; });
            image.fill_circle_aa(random_limit(), random_limit(), random_limit(), [](float u, float v, float au, float av) -> std::array<float, 4> { return {u, v, au, av}; });
            image.fill_round_rect_aa(random_limit(), random_limit(), random_limit(), random_limit(), random_limit(), [](float u, float v, float au, float av) -> std::array<float, 4> { return {u, v, au, av}; });
            image.draw_line_aa({.x0 = random_limit(), .y0 = random_limit(), .x1 = random_limit(), .y1 = random_limit(), .col = static_cast<uint8_t>(random_limit()&0xFF), .sw = static_cast<float>(random_limit())});
            image.fill_circle_aa({.cx = random_limit(), .cy = random_limit(), .r = random_limit(), .col = static_cast<uint8_t>(random_limit()&0xFF)});
            image.stroke_circle_aa({.cx = random_limit(), .cy = random_limit(), .r = random_limit(), .col = static_cast<uint8_t>(random_limit()&0xFF), .sw = static_cast<int32_t>(random_limit())});
            image.fill_round_rect_aa({.x = random_limit(), .y = random_limit(), .w = random_limit(), .h = random_limit(), .r = random_limit(), .col = static_cast<uint8_t>(random_limit()&0xFF)});
            image.stroke_round_rect_aa({.x = random_limit(), .y = random_limit(), .w = random_limit(), .h = random_limit(), .r = random_limit(), .col = static_cast<uint8_t>(random_limit()&0xFF), .sw = static_cast<int32_t>(random_limit())});
            image.template draw_string_aa<fontaa>(random_limit(), random_limit(), "ABCabc", static_cast<uint8_t>(random_limit()&0xFF));
            image.template draw_string_aa<fontaa, false, constixel::DEGREE_90>(random_limit(), random_limit(), "ABCabc", static_cast<uint8_t>(random_limit()&0xFF));
            image.template draw_string_aa<fontaa, false, constixel::DEGREE_180>(random_limit(), random_limit(), "ABCabc", static_cast<uint8_t>(random_limit()&0xFF));
            image.template draw_string_aa<fontaa, false, constixel::DEGREE_270>(random_limit(), random_limit(), "ABCabc", static_cast<uint8_t>(random_limit()&0xFF));
            image.template draw_string_centered_aa<fontaa>(random_limit(), random_limit(), "ABCabc", static_cast<uint8_t>(random_limit()&0xFF));
        }
    }
    for (size_t c = 0; c < 256; c++) {
        image.transpose();
    }
}

template <class T, bool AA>
void fuzz_random() {
    T image;
    for (size_t c = 0; c < 65536; c++) {
        image.draw_line(random_int32(),random_int32(),random_int32(),random_int32(),static_cast<uint8_t>(random_int32()&0xFF));
        image.draw_line(random_int32(),random_int32(),random_int32(),random_int32(),static_cast<uint8_t>(random_int32()&0xFF),random_int32());
        image.fill_circle(random_int32(),random_int32(),random_int32(),static_cast<uint8_t>(random_int32()&0xFF));
        image.fill_rect(random_int32(), random_int32(), random_int32(), random_int32(), static_cast<uint8_t>(random_int32()&0xFF));
        image.stroke_rect(random_int32(), random_int32(), random_int32(), random_int32(), static_cast<uint8_t>(random_int32()&0xFF));
        image.stroke_rect(random_int32(), random_int32(), random_int32(), random_int32(), static_cast<uint8_t>(random_int32()&0xFF), random_int32());
        image.fill_round_rect(random_int32(), random_int32(), random_int32(), random_int32(), random_int32(), static_cast<uint8_t>(random_int32()&0xFF));
        image.stroke_circle(random_int32(), random_int32(), random_int32(), static_cast<uint8_t>(random_int32()&0xFF),  random_int32());
        image.stroke_round_rect(random_int32(), random_int32(), random_int32(), random_int32(), random_int32(), static_cast<uint8_t>(random_int32()&0xFF), random_int32());
        image.plot(random_int32(), random_int32(), static_cast<uint8_t>(random_int32()&0xFF));
        image.plot({.x = random_int32(), .y = random_int32(), .col = static_cast<uint8_t>(random_int32()&0xFF)});
        image.draw_line({.x0 = random_int32(), .y0 = random_int32(), .x1 = random_int32(), .y1 = random_int32(), .col = static_cast<uint8_t>(random_int32()&0xFF), .sw = static_cast<float>(random_int32())});
        image.fill_rect({.x = random_int32(), .y = random_int32(), .w = random_int32(), .h = random_int32(), .col = static_cast<uint8_t>(random_int32()&0xFF)});
        image.stroke_rect({.x = random_int32(), .y = random_int32(), .w = random_int32(), .h = random_int32(), .col = static_cast<uint8_t>(random_int32()&0xFF), .sw = static_cast<int32_t>(random_int32())});
        image.fill_circle({.cx = random_int32(), .cy = random_int32(), .r = random_int32(), .col = static_cast<uint8_t>(random_int32()&0xFF)});
        image.stroke_circle({.cx = random_int32(), .cy = random_int32(), .r = random_int32(), .col = static_cast<uint8_t>(random_int32()&0xFF), .sw = static_cast<int32_t>(random_int32())});
        image.fill_round_rect({.x = random_int32(), .y = random_int32(), .w = random_int32(), .h = random_int32(), .r = random_int32(), .col = static_cast<uint8_t>(random_int32()&0xFF)});
        image.stroke_round_rect({.x = random_int32(), .y = random_int32(), .w = random_int32(), .h = random_int32(), .r = random_int32(), .col = static_cast<uint8_t>(random_int32()&0xFF), .sw = static_cast<int32_t>(random_int32())});
        image.template draw_string_mono<font>(random_int32(), random_int32(), "ABCabc", static_cast<uint8_t>(random_int32()&0xFF));
        image.template draw_string_mono<font, false, constixel::DEGREE_90>(random_int32(), random_int32(), "ABCabc", static_cast<uint8_t>(random_int32()&0xFF));
        image.template draw_string_mono<font, false, constixel::DEGREE_180>(random_int32(), random_int32(), "ABCabc", static_cast<uint8_t>(random_int32()&0xFF));
        image.template draw_string_mono<font, false, constixel::DEGREE_270>(random_int32(), random_int32(), "ABCabc", static_cast<uint8_t>(random_int32()&0xFF));
        image.template draw_string_centered_mono<font>(random_int32(), random_int32(), "ABCabc", static_cast<uint8_t>(random_int32()&0xFF));
        if constexpr (AA) {
            image.stroke_circle_aa(random_int32(), random_int32(), random_int32(), static_cast<uint8_t>(random_int32()&0xFF),  random_int32());
            image.stroke_round_rect_aa(random_int32(), random_int32(), random_int32(), random_int32(), random_int32(), static_cast<uint8_t>(random_int32()&0xFF), random_int32());
            image.draw_line_aa(random_int32(),random_int32(),random_int32(),random_int32(),static_cast<uint8_t>(random_int32()&0xFF));
            image.draw_line_aa(random_int32(),random_int32(),random_int32(),random_int32(),static_cast<uint8_t>(random_int32()&0xFF),random_int32());
            image.fill_circle_aa(random_int32(), random_int32(), random_int32(), static_cast<uint8_t>(random_int32()&0xFF));
            image.fill_rect(random_int32(), random_int32(), random_int32(), random_int32(), [](float u, float v, float au, float av) -> std::array<float, 4> { return {u, v, au, av}; });
            image.fill_circle_aa(random_int32(), random_int32(), random_int32(), [](float u, float v, float au, float av) -> std::array<float, 4> { return {u, v, au, av}; });
            image.fill_round_rect_aa(random_int32(), random_int32(), random_int32(), random_int32(), random_int32(), [](float u, float v, float au, float av) -> std::array<float, 4> { return {u, v, au, av}; });
            image.fill_round_rect_aa(random_int32(), random_int32(), random_int32(), random_int32(), random_int32(), static_cast<uint8_t>(random_int32()&0xFF));
            image.draw_line_aa({.x0 = random_int32(), .y0 = random_int32(), .x1 = random_int32(), .y1 = random_int32(), .col = static_cast<uint8_t>(random_int32()&0xFF), .sw = static_cast<float>(random_int32())});
            image.fill_circle_aa({.cx = random_int32(), .cy = random_int32(), .r = random_int32(), .col = static_cast<uint8_t>(random_int32()&0xFF)});
            image.stroke_circle_aa({.cx = random_int32(), .cy = random_int32(), .r = random_int32(), .col = static_cast<uint8_t>(random_int32()&0xFF), .sw = static_cast<int32_t>(random_int32())});
            image.fill_round_rect_aa({.x = random_int32(), .y = random_int32(), .w = random_int32(), .h = random_int32(), .r = random_int32(), .col = static_cast<uint8_t>(random_int32()&0xFF)});
            image.stroke_round_rect_aa({.x = random_int32(), .y = random_int32(), .w = random_int32(), .h = random_int32(), .r = random_int32(), .col = static_cast<uint8_t>(random_int32()&0xFF), .sw = static_cast<int32_t>(random_int32())});
            image.template draw_string_aa<fontaa>(random_int32(), random_int32(), "ABCabc", static_cast<uint8_t>(random_int32()&0xFF));
            image.template draw_string_aa<fontaa, false, constixel::DEGREE_90>(random_int32(), random_int32(), "ABCabc", static_cast<uint8_t>(random_int32()&0xFF));
            image.template draw_string_aa<fontaa, false, constixel::DEGREE_180>(random_int32(), random_int32(), "ABCabc", static_cast<uint8_t>(random_int32()&0xFF));
            image.template draw_string_aa<fontaa, false, constixel::DEGREE_270>(random_int32(), random_int32(), "ABCabc", static_cast<uint8_t>(random_int32()&0xFF));
            image.template draw_string_centered_aa<fontaa>(random_int32(), random_int32(), "ABCabc", static_cast<uint8_t>(random_int32()&0xFF));
        }
    }
    for (size_t c = 0; c < 256; c++) {
        image.transpose();
    }
}

int main() {
    printf("format_1bit 256x256 fuzz_random\n");
    fuzz_random<constixel::image<constixel::format_1bit, 256, 256>, false>();
    printf("format_1bit 256x256 fuzz_limits\n");
    fuzz_limits<constixel::image<constixel::format_1bit, 256, 256>, false>();
    printf("format_1bit 255x7 fuzz_random\n");
    fuzz_random<constixel::image<constixel::format_1bit, 255, 7>, false>();
    printf("format_1bit 255x7 fuzz_limits\n");
    fuzz_limits<constixel::image<constixel::format_1bit, 255, 7>, false>();
    printf("format_2bit 256x256 fuzz_random\n");
    fuzz_random<constixel::image<constixel::format_2bit, 256, 256>, false>();
    printf("format_2bit 256x256 fuzz_limits\n");
    fuzz_limits<constixel::image<constixel::format_2bit, 256, 256>, false>();
    printf("format_2bit 255x7 fuzz_random\n");
    fuzz_random<constixel::image<constixel::format_2bit, 255, 7>, false>();
    printf("format_2bit 255x7 fuzz_limits\n");
    fuzz_limits<constixel::image<constixel::format_2bit, 255, 7>, false>();
    printf("format_4bit 256x256 fuzz_random\n");
    fuzz_random<constixel::image<constixel::format_4bit, 256, 256>, false>();
    printf("format_4bit 256x256 fuzz_limits\n");
    fuzz_limits<constixel::image<constixel::format_4bit, 256, 256>, false>();
    printf("format_4bit 255x7 fuzz_random\n");
    fuzz_random<constixel::image<constixel::format_4bit, 255, 7>, false>();
    printf("format_4bit 255x7 fuzz_limits\n");
    fuzz_limits<constixel::image<constixel::format_4bit, 255, 7>, false>();
    printf("format_8bit 256x256 fuzz_random\n");
    fuzz_random<constixel::image<constixel::format_8bit, 256, 256>, true>();
    printf("format_8bit 256x256 fuzz_limits\n");
    fuzz_limits<constixel::image<constixel::format_8bit, 256, 256>, true>();
    printf("format_8bit 255x7 fuzz_random\n");
    fuzz_random<constixel::image<constixel::format_8bit, 255, 7>, true>();
    printf("format_8bit 255x7 fuzz_limits\n");
    fuzz_limits<constixel::image<constixel::format_8bit, 255, 7>, true>();
    printf("format_24bit 256x256 fuzz_random\n");
    fuzz_random<constixel::image<constixel::format_24bit, 256, 256>, false>();
    printf("format_24bit 256x256 fuzz_limits\n");
    fuzz_limits<constixel::image<constixel::format_24bit, 256, 256>, false>();
    printf("format_24bit 255x7 fuzz_random\n");
    fuzz_random<constixel::image<constixel::format_24bit, 255, 7>, false>();
    printf("format_24bit 255x7 fuzz_limits\n");
    fuzz_limits<constixel::image<constixel::format_24bit, 255, 7>, false>();
    printf("format_32bit 256x256 fuzz_random\n");
    fuzz_random<constixel::image<constixel::format_32bit, 256, 256>, false>();
    printf("format_32bit 256x256 fuzz_limits\n");
    fuzz_limits<constixel::image<constixel::format_32bit, 256, 256>, false>();
    printf("format_32bit 255x7 fuzz_random\n");
    fuzz_random<constixel::image<constixel::format_32bit, 255, 7>, false>();
    printf("format_32bit 255x7 fuzz_limits\n");
    fuzz_limits<constixel::image<constixel::format_32bit, 255, 7>, false>();
}
