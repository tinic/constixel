#include "constixel.hpp"

#include "fonts/ibmplexsans_bold_24_mono.hpp"
#include "fonts/ibmplexsans_bold_32_mono.hpp"
#include "fonts/ibmplexsans_bold_48_mono.hpp"
#include "fonts/ibmplexsans_medium_14_mono.hpp"
#include "fonts/ibmplexsans_regular_14_aa.hpp"
#include "fonts/ibmplexsans_semibold_24_aa.hpp"
#include "fonts/pixeloperator_mono.hpp"
#include "fonts/px437_hp_100lx_6x8_mono.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <memory>
#include <span>
#include <string_view>
#include <vector>

using namespace constixel;

//
// One user interface, three displays.
//
// draw_ui() takes a surface, so the same code draws a 1-bit OLED, a 16 color TFT and a grayscale e-paper screen,
// each with its own size. Every screen then leaves the library the way it would on a device:
//
//   OLED     convert_chunk<MONO_1BIT_VERTICAL_PAGES>() in 32 byte chunks into the display RAM of a SSD1306
//   TFT      convert<RGB565_8BIT_SERIAL>() into a RGB565 frame buffer
//   e-paper  RGBA_uint8()
//
// The picture is rendered from that display data, not from the images.
//

struct theme {
    uint8_t fg;      // text and graph
    uint8_t dim;     // secondary text
    uint8_t fill;    // below the graph
    uint8_t accent;  // title bar
    uint8_t level;   // level bar
};

template <typename PF, typename BIG, typename SMALL>
void draw_ui(surface<PF> s, const theme &t, std::string_view title, std::string_view value, std::string_view unit,
             std::span<const float> history, float level) {
    const int32_t w = s.width();
    const int32_t h = s.height();
    const int32_t pad = std::max(w / 50, int32_t{2});
    const int32_t bar = SMALL::line_height + 1;

    s.fill_rect(0, 0, w, h, 0);

    // Title bar, text knocked out in the background color.
    s.fill_rect(0, 0, w, bar, t.accent);
    s.template draw_string<SMALL>(pad, 0, title, 0);
    s.template draw_string<SMALL>(w - pad, 0, "12:47", 0, ALIGN_RIGHT);

    // Wide screens put the graph next to the value, tall ones below it.
    const bool tall = h * 10 > w * 7;
    const int32_t value_y = bar + pad;
    const int32_t value_x = s.template draw_string<BIG>(pad, value_y, value, t.fg);
    s.template draw_string<SMALL>(value_x + pad, value_y + BIG::ascent - SMALL::ascent, unit, t.dim);

    // Level bar.
    const int32_t level_h = std::max(h / 9, int32_t{7});
    const int32_t level_w = tall ? w - pad * 2 : w / 2 - pad * 2;
    const int32_t level_y = tall ? value_y + BIG::line_height + pad : h - pad - level_h;
    s.stroke_rect(pad, level_y, level_w, level_h, t.fg);
    s.fill_rect(pad + 2, level_y + 2, static_cast<int32_t>(static_cast<float>(level_w - 4) * level), level_h - 4,
                t.level);

    // Graph, stippled below the curve.
    const int32_t gx = tall ? pad : w / 2 + pad;
    const int32_t gy = tall ? level_y + level_h + pad * 2 : bar + pad;
    const int32_t gw = w - pad - gx;
    const int32_t gh = h - pad - gy;
    auto sample = [&](int32_t x) {
        const float pos = static_cast<float>(x) / static_cast<float>(gw - 1) * static_cast<float>(history.size() - 1);
        const auto i = std::min(static_cast<size_t>(pos), history.size() - 2);
        const float f = pos - static_cast<float>(i);
        const float val = history[i] * (1.0f - f) + history[i + 1] * f;
        return gy + gh - 1 - static_cast<int32_t>(val * static_cast<float>(gh - 1));
    };
    for (int32_t x = 0; x < gw; x++) {
        for (int32_t y = sample(x) + 2; y < gy + gh; y++) {
            if (((x + y) & 1) == 0) {
                s.plot(gx + x, y, t.fill);
            }
        }
    }
    for (int32_t x = 1; x < gw; x++) {
        s.draw_line(gx + x - 1, sample(x - 1), gx + x, sample(x), t.fg);
    }
    s.draw_line(gx, gy + gh, gx + gw - 1, gy + gh, t.fg);
}

// Brightness of a point inside a display pixel: the last row and column of every pixel are the gap to the next one.
static float pixel_mask(int32_t px, int32_t py, int32_t scale, float gap) {
    return (px % scale == scale - 1 || py % scale == scale - 1) ? gap : 1.0f;
}

int main() {
    using canvas_type = image<format_32bit, 1020, 690>;
    static auto canvas = std::make_unique<canvas_type>();

    static constexpr std::array<float, 24> history{0.08f, 0.10f, 0.16f, 0.22f, 0.35f, 0.44f, 0.41f, 0.55f,
                                                   0.68f, 0.74f, 0.71f, 0.86f, 0.93f, 0.90f, 0.95f, 0.82f,
                                                   0.63f, 0.70f, 0.52f, 0.40f, 0.31f, 0.18f, 0.12f, 0.07f};

    // ─── the three screens ────────────────────────────────────────────────────
    static image<format_1bit, 128, 64> oled;
    draw_ui<pixel_1bit, ibmplexsans_bold_24_mono, px437_hp_100lx_6x8_mono>(
        oled.view(), {.fg = 1, .dim = 1, .fill = 1, .accent = 1, .level = 1}, "SOLAR", "4.82", "kW", history, 0.72f);

    static image<format_4bit, 160, 80> tft;
    draw_ui<pixel_4bit<false>, ibmplexsans_bold_32_mono, pixeloperator_mono>(
        tft.view(), {.fg = color::WHITE, .dim = color::GRAY_40, .fill = color::BLUE, .accent = color::CYAN, .level = color::YELLOW},
        "SOLAR", "4.82", "kW", history, 0.72f);

    static image<format_2bit, 200, 200, true> paper;
    draw_ui<pixel_2bit<true>, ibmplexsans_bold_48_mono, ibmplexsans_medium_14_mono>(
        paper.view(), {.fg = 3, .dim = 2, .fill = 1, .accent = 3, .level = 2}, "SOLAR", "4.82", "kW", history, 0.72f);

    // ─── out of the library, the way a device would get it ────────────────────
    static std::array<uint8_t, 128 * 64 / 8> ssd1306_ram{};
    {
        std::array<uint8_t, 32> chunk{};
        convert_cursor cursor;
        size_t pos = 0;
        while (const size_t n = oled.convert_chunk<MONO_1BIT_VERTICAL_PAGES>(chunk, cursor)) {
            std::copy_n(chunk.begin(), n, ssd1306_ram.begin() + static_cast<std::ptrdiff_t>(pos));
            pos += n;
        }
    }

    static std::vector<uint8_t> rgb565;
    tft.convert<RGB565_8BIT_SERIAL>([](uint8_t b) {
        rgb565.push_back(b);
    });

    static std::array<uint8_t, 200 * 200 * 4> paper_rgba{};
    paper.RGBA_uint8(paper_rgba);

    // ─── the desk ─────────────────────────────────────────────────────────────
    canvas->fill_rect(0, 0, canvas->width(), canvas->height(), [](float u, float v, int32_t, int32_t) {
        const float d = std::sqrt((u - 0.5f) * (u - 0.5f) + (v - 0.45f) * (v - 0.45f));
        const float l = std::max(0.0f, 1.0f - d * 1.3f);
        return std::array<float, 4>{0.05f + l * 0.09f, 0.06f + l * 0.10f, 0.08f + l * 0.13f, 1.0f};
    });

    auto board = [&](int32_t x, int32_t y, int32_t w, int32_t h, float r, float g, float b) {
        canvas->fill_round_rect_aa(x + 4, y + 8, w, h, 14, [](float, float, int32_t, int32_t) {
            return std::array<float, 4>{0.0f, 0.0f, 0.0f, 0.45f};
        });
        canvas->fill_round_rect_aa(x, y, w, h, 14, [=](float u, float v, int32_t, int32_t) {
            const float l = 1.0f - (u + v) * 0.25f;
            return std::array<float, 4>{r * l, g * l, b * l, 1.0f};
        });
    };

    const uint8_t label_col = canvas->get_nearest_color(190, 200, 215);
    canvas->draw_string<ibmplexsans_semibold_24_aa, true>(40, 22, "One draw_ui(surface<PF>), three displays",
                                                          color::WHITE);

    // SSD1306, 128x64, 3x.
    {
        const int32_t scale = 3, x0 = 62, y0 = 102;
        board(x0 - 22, y0 - 22, 128 * scale + 44, 64 * scale + 64, 0.07f, 0.22f, 0.50f);
        canvas->fill_rect(x0, y0, 128 * scale, 64 * scale, [=](float, float, int32_t px, int32_t py) {
            const int32_t x = (px - x0) / scale;
            const int32_t y = (py - y0) / scale;
            const bool on = ((ssd1306_ram[static_cast<size_t>((y / 8) * 128 + x)] >> (y % 8)) & 1) != 0;
            const float m = pixel_mask(px - x0, py - y0, scale, 0.45f);
            return on ? std::array<float, 4>{0.55f * m, 0.92f * m, 1.0f * m, 1.0f}
                      : std::array<float, 4>{0.01f, 0.02f, 0.03f, 1.0f};
        });
        canvas->draw_string<ibmplexsans_regular_14_aa, true>(x0 - 22, y0 + 64 * scale + 50,
                                                             "SSD1306 OLED, 128x64, 1 bit, 1024 bytes", label_col);
    }

    // RGB565 TFT, 160x80, 3x.
    {
        const int32_t scale = 3, x0 = 54, y0 = 398;
        board(x0 - 14, y0 - 14, 160 * scale + 28, 80 * scale + 28, 0.10f, 0.10f, 0.11f);
        canvas->fill_rect(x0, y0, 160 * scale, 80 * scale, [=](float, float, int32_t px, int32_t py) {
            const auto i = static_cast<size_t>(((py - y0) / scale) * 160 + (px - x0) / scale) * 2;
            const uint32_t a = rgb565[i], b = rgb565[i + 1];
            const float m = pixel_mask(px - x0, py - y0, scale, 0.8f);
            return std::array<float, 4>{static_cast<float>(a >> 3) / 31.0f * m,
                                        static_cast<float>(((a & 7) << 3) | (b >> 5)) / 63.0f * m,
                                        static_cast<float>(b & 31) / 31.0f * m, 1.0f};
        });
        canvas->draw_string<ibmplexsans_regular_14_aa, true>(
            x0 - 14, y0 + 80 * scale + 22, "TFT, 160x80, 16 colors, 6400 bytes, sent as RGB565", label_col);
    }

    // E-paper, 200x200, 2x.
    {
        const int32_t scale = 2, x0 = 582, y0 = 128;
        board(x0 - 26, y0 - 26, 200 * scale + 52, 200 * scale + 70, 0.80f, 0.80f, 0.78f);
        canvas->fill_rect(x0, y0, 200 * scale, 200 * scale, [=](float, float, int32_t px, int32_t py) {
            const auto i = static_cast<size_t>(((py - y0) / scale) * 200 + (px - x0) / scale) * 4;
            // This panel shows index 0 as paper and index 3 as ink.
            const float l = 1.0f - static_cast<float>(paper_rgba[i]) / 255.0f;
            return std::array<float, 4>{0.07f + l * 0.80f, 0.07f + l * 0.79f, 0.08f + l * 0.73f, 1.0f};
        });
        canvas->draw_string<ibmplexsans_regular_14_aa, true>(x0 - 26, y0 + 200 * scale + 52,
                                                             "E-paper, 200x200, 2 bit gray, 10000 bytes", label_col);
    }
    std::ofstream file("example_displays.png", std::ios::binary);
    canvas->png([&](char ch) {
        file.put(ch);
    });
}
