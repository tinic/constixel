#pragma once

#include "constixel.hpp"
#include "fonts/ibmplexsans_medium_12_aa.hpp"
#include "fonts/ibmplexsans_medium_14_aa.hpp"
#include "fonts/ibmplexsans_regular_18_aa.hpp"
#include "fonts/ibmplexsans_semibold_18_aa.hpp"
#include "fonts/ibmplexsans_semibold_24_aa.hpp"
#include "fonts/ibmplexsans_semibold_32_aa.hpp"

// Presentation helpers shared by the larger examples. All drawing still uses
// constixel's built-in palette, fonts and primitives; no external assets needed.
namespace example {
using caption = constixel::ibmplexsans_medium_12_aa;
using small = constixel::ibmplexsans_medium_14_aa;
using body = constixel::ibmplexsans_regular_18_aa;
using heading = constixel::ibmplexsans_semibold_18_aa;
using title = constixel::ibmplexsans_semibold_24_aa;
using display = constixel::ibmplexsans_semibold_32_aa;

struct Theme {
    uint8_t background, surface, ink, muted, border, accent, blue, warm, danger;

    template <typename Image>
    explicit constexpr Theme(const Image& img)
        : background(img.get_nearest_color(238, 238, 238)),
          surface(constixel::color::WHITE),
          ink(img.get_nearest_color(34, 34, 34)),
          muted(img.get_nearest_color(102, 102, 102)),
          border(img.get_nearest_color(221, 221, 221)),
          accent(img.get_nearest_color(0, 125, 125)),
          blue(img.get_nearest_color(66, 128, 200)),
          warm(img.get_nearest_color(192, 120, 45)),
          danger(img.get_nearest_color(192, 60, 60)) {
    }
};

// Text coordinates in these examples are baselines, making rows easy to align.
template <typename Font = body, typename Image>
constexpr void text(Image& img, int x, int baseline, const char* str, uint8_t color) {
    img.template draw_string_aa<Font, true>(x, baseline - Font::ascent, str, color);
}

template <typename Font = small, typename Image>
constexpr void right_text(Image& img, int right, int baseline, const char* str, uint8_t color) {
    text<Font>(img, right - img.template string_width<Font, true>(str), baseline, str, color);
}

template <typename Font = body, typename Image>
constexpr void centered_text(Image& img, int center, int baseline, const char* str, uint8_t color) {
    text<Font>(img, center - img.template string_width<Font, true>(str) / 2, baseline, str, color);
}

template <typename Image>
constexpr void card(Image& img, int x, int y, int w, int h, const Theme& theme) {
    img.fill_round_rect_aa(x, y, w, h, 16, theme.surface);
}

template <typename Image>
constexpr void header(Image& img, const char* name, const char* description, const Theme& theme) {
    img.fill_rect(0, 0, img.width(), img.height(), theme.background);
    img.fill_round_rect_aa(32, 30, 6, 48, 3, theme.accent);
    text<display>(img, 52, 54, name, theme.ink);
    text<small>(img, 54, 80, description, theme.muted);
    right_text<small>(img, img.width() - 32, 48, "constixel / examples", theme.muted);
}
}  // namespace example
