#define CONSTIXEL_ENABLE_COUT

#include <array>
#include <cstdio>

#include "../example_output.hpp"
#include "../example_style.hpp"
#include "fonts/ibmplexsans_semibold_18_mono.hpp"
#include "fonts/ibmplexsans_semibold_24_mono.hpp"

using namespace constixel;
using namespace example;
using UI = image<format_8bit, 1024, 864>;
using mono = ibmplexsans_semibold_24_mono;
using mono_small = ibmplexsans_semibold_18_mono;

int main(int argc, char** argv) {
    static UI img;
    const Theme theme(img);
    header(img, "Drawing essentials", "A small vocabulary for building complete interfaces", theme);
    constexpr std::array<const char*, 16> labels = {
        "fill_rect",        "fill_round_rect",  "fill_round_rect_aa", "fill_circle",
        "fill_circle_aa",   "draw_line",        "draw_line_aa",       "stroke_rect",
        "stroke_circle",    "stroke_circle_aa", "stroke_round_rect",  "stroke_round_rect_aa",
        "draw_string_mono", "draw_string_aa",   "rotate / mono",      "rotate / antialiased"};
    for (int i = 0; i < 16; ++i) {
        const int x = 32 + (i % 4) * 244;
        const int y = 108 + (i / 4) * 180;
        const int cx = x + 114, cy = y + 81;
        const uint8_t col = i < 5 ? theme.accent : (i < 12 ? theme.blue : theme.ink);
        card(img, x, y, 228, 164, theme);
        char number[8];
        std::snprintf(number, sizeof(number), "%02d", i + 1);
        text<caption>(img, x + 18, y + 25, number, theme.muted);
        centered_text<small>(img, cx, y + 145, labels[i], theme.ink);
        switch (i) {
            case 0:
                img.fill_rect(cx - 60, cy - 26, 120, 52, col);
                break;
            case 1:
                img.fill_round_rect(cx - 60, cy - 26, 120, 52, 16, col);
                break;
            case 2:
                img.fill_round_rect_aa(cx - 60, cy - 26, 120, 52, 16, col);
                break;
            case 3:
                img.fill_circle(cx, cy, 30, col);
                break;
            case 4:
                img.fill_circle_aa(cx, cy, 30, col);
                break;
            case 5:
                img.draw_line(cx - 60, cy + 24, cx + 60, cy - 24, col, 3);
                break;
            case 6:
                img.draw_line_aa(cx - 60, cy + 24, cx + 60, cy - 24, col, 3);
                break;
            case 7:
                img.stroke_rect(cx - 60, cy - 26, 120, 52, col, 3);
                break;
            case 8:
                img.stroke_circle(cx, cy, 30, col, 4);
                break;
            case 9:
                img.stroke_circle_aa(cx, cy, 30, col, 4);
                break;
            case 10:
                img.stroke_round_rect(cx - 60, cy - 26, 120, 52, 16, col, 4);
                break;
            case 11:
                img.stroke_round_rect_aa(cx - 60, cy - 26, 120, 52, 16, col, 4);
                break;
            case 12:
                img.draw_string_centered_mono<mono, true>(cx, cy + 8 - mono::ascent, "Aa 0123", col);
                break;
            case 13:
                centered_text<title>(img, cx, cy + 8, "Aa 0123", col);
                break;
            case 14:
                img.draw_string_mono<mono_small, true>(cx - 50, cy - 31, "Aa", col);
                img.draw_string_mono<mono_small, true, DEGREE_90>(cx + 48, cy - 31, "Aa", col);
                img.draw_string_mono<mono_small, true, DEGREE_180>(cx + 50, cy + 31, "Aa", col);
                img.draw_string_mono<mono_small, true, DEGREE_270>(cx - 48, cy + 31, "Aa", col);
                break;
            case 15:
                img.draw_string_aa<heading, true>(cx - 50, cy - 31, "Aa", col);
                img.draw_string_aa<heading, true, DEGREE_90>(cx + 48, cy - 31, "Aa", col);
                img.draw_string_aa<heading, true, DEGREE_180>(cx + 50, cy + 31, "Aa", col);
                img.draw_string_aa<heading, true, DEGREE_270>(cx - 48, cy + 31, "Aa", col);
                break;
        }
    }
    text<small>(img, 32, 850, "Fills, strokes and typography", theme.muted);
    right_text<small>(img, 992, 850, "8-bit palette / native resolution", theme.muted);
    return output(img, argc, argv);
}
