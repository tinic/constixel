#define CONSTIXEL_ENABLE_COUT

#include <array>
#include <cmath>
#include <cstdio>
#include <numbers>

#include "../example_output.hpp"
#include "constixel.hpp"
#include "fonts/ibmplexsans_medium_12_mono.hpp"
#include "fonts/ibmplexsans_medium_14_mono.hpp"
#include "fonts/ibmplexsans_semibold_18_mono.hpp"
#include "fonts/ibmplexsans_semibold_24_mono.hpp"
#include "fonts/ibmplexsans_semibold_48_mono.hpp"

using namespace constixel;
using Screen = image<format_1bit, 800, 480>;
using caption = ibmplexsans_medium_12_mono;
using small = ibmplexsans_medium_14_mono;
using heading = ibmplexsans_semibold_18_mono;
using title = ibmplexsans_semibold_24_mono;
using display = ibmplexsans_semibold_48_mono;

namespace {
constexpr uint8_t ink = 0, paper = 1;

// Deliberately monochrome: the pixel buffer is exactly 800 * 480 / 8 = 48 KB.
// No antialiased fonts or RGB staging buffer are needed for an e-paper display.
template <typename Font>
void text(Screen& img, int x, int baseline, const char* str, uint8_t color = ink) {
    img.draw_string_mono<Font, true>(x, baseline - Font::ascent, str, color);
}

template <typename Font>
void centered(Screen& img, int cx, int baseline, const char* str) {
    text<Font>(img, cx - img.string_width<Font, true>(str) / 2, baseline, str);
}

void sun(Screen& img, int cx, int cy, int radius) {
    img.stroke_circle(cx, cy, radius, ink, 2);
    for (int i = 0; i < 8; ++i) {
        const float angle = i * std::numbers::pi_v<float> / 4;
        const int outer = radius + 10;
        img.draw_line(cx + static_cast<int>((radius + 5) * std::cos(angle)),
                      cy + static_cast<int>((radius + 5) * std::sin(angle)),
                      cx + static_cast<int>(outer * std::cos(angle)), cy + static_cast<int>(outer * std::sin(angle)),
                      ink, 2);
    }
}

void cloud(Screen& img, int cx, int cy, bool rainy) {
    img.fill_circle(cx - 14, cy, 12, ink);
    img.fill_circle(cx + 1, cy - 7, 17, ink);
    img.fill_circle(cx + 18, cy + 1, 11, ink);
    img.fill_rect(cx - 14, cy, 33, 13, ink);
    if (rainy) {
        for (int i = 0; i < 3; ++i) {
            img.draw_line(cx - 16 + i * 15, cy + 21, cx - 21 + i * 15, cy + 30, ink, 2);
        }
    }
}

void rain_chart(Screen& img) {
    text<heading>(img, 424, 111, "Rain chance", ink);
    text<caption>(img, 424, 134, "NEXT 12 HOURS", ink);
    constexpr std::array<int, 12> chance = {5, 5, 10, 15, 30, 45, 60, 50, 35, 20, 10, 5};
    constexpr int left = 444, bottom = 237, height = 84;
    for (int i = 0; i < 3; ++i) {
        const int y = bottom - i * height / 2;
        for (int x = left; x < 756; x += 4) {
            img.plot(x, y, ink);
        }
        const std::array<const char*, 3> labels = {"0", "30", "60"};
        text<caption>(img, 416, y + 4, labels[i]);
    }
    for (size_t i = 0; i < chance.size(); ++i) {
        const int x = left + static_cast<int>(i) * 26 + 4;
        const int h = chance[i] * height / 60;
        // A one-bit stipple provides a second visual weight without grayscale.
        for (int y = bottom - h; y < bottom; ++y) {
            for (int px = x; px < x + 18; ++px) {
                if ((px + y) % 2 == 0) {
                    img.plot(px, y, ink);
                }
            }
        }
        img.fill_rect(x, bottom - h, 18, 2, ink);
    }
    for (int i = 0; i < 4; ++i) {
        char label[8];
        std::snprintf(label, sizeof(label), "%02d:00", 10 + i * 3);
        centered<caption>(img, left + 13 + i * 78, 260, label);
    }
    text<small>(img, 424, 283, "Peak: 60% at 16:00");
}
}  // namespace

int main(int argc, char** argv) {
    static Screen img;
    img.fill_rect(0, 0, 800, 480, paper);
    text<title>(img, 32, 49, "Harbor weather");
    text<small>(img, 600, 48, "TUESDAY / 09:41");
    img.draw_line(32, 72, 768, 72, ink);

    text<caption>(img, 32, 111, "OUTSIDE NOW");
    text<display>(img, 32, 170, "22");
    img.stroke_circle(100, 132, 4, ink, 2);
    text<heading>(img, 112, 147, "C");
    sun(img, 282, 153, 28);
    text<title>(img, 32, 208, "Mostly sunny");
    text<small>(img, 32, 235, "High 24 C / Low 17 C");
    text<small>(img, 32, 265, "Wind 12 km/h SW");
    text<small>(img, 202, 265, "Humidity 68%");
    rain_chart(img);
    img.draw_line(390, 96, 390, 282, ink);
    img.draw_line(32, 302, 768, 302, ink);

    struct Forecast {
        const char* day;
        int high, low, condition;
    };
    constexpr std::array<Forecast, 5> days = {
        {{"TODAY", 24, 17, 0}, {"WED", 21, 16, 1}, {"THU", 19, 15, 2}, {"FRI", 20, 15, 1}, {"SAT", 23, 16, 0}}};
    for (size_t i = 0; i < days.size(); ++i) {
        const int cx = 104 + static_cast<int>(i) * 148;
        centered<heading>(img, cx, 332, days[i].day);
        if (days[i].condition == 0) {
            sun(img, cx, 373, 15);
        } else {
            cloud(img, cx, 370, days[i].condition == 2);
        }
        char label[24];
        std::snprintf(label, sizeof(label), "%d / %d C", days[i].high, days[i].low);
        centered<small>(img, cx, 426, label);
        if (i < days.size() - 1) {
            img.draw_line(cx + 74, 322, cx + 74, 430, ink);
        }
    }
    img.draw_line(32, 448, 768, 448, ink);
    text<caption>(img, 32, 469, "DEMO FORECAST / OFFLINE SNAPSHOT");
    text<caption>(img, 602, 469, "NEXT UPDATE 10:00");
    return example::output(img, argc, argv);
}
