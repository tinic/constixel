#include "constixel.hpp"

#include "fonts/ibmplexsans_bolditalic_48_aa.hpp"
#include "fonts/m42_mono.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <memory>

using namespace constixel;

//
// A retro sunset, built from shader fills, antialiased lines and text.
//

static constexpr int32_t W = 800;
static constexpr int32_t H = 450;
static constexpr int32_t HORIZON = 262;

int main() {
    static auto img = std::make_unique<image<format_32bit, W, H>>();

    // ─── sky ──────────────────────────────────────────────────────────────────
    img->fill_rect(0, 0, W, HORIZON, [](float, float v, int32_t, int32_t) {
        const float glow = v * v * v;
        return std::array<float, 4>{0.05f + 0.25f * v + 0.65f * glow, 0.01f + 0.04f * v + 0.18f * glow,
                                    0.14f + 0.22f * v + 0.20f * glow, 1.0f};
    });

    // ─── stars ────────────────────────────────────────────────────────────────
    uint32_t seed = 0x2545f491;
    auto rnd = [&seed]() {
        seed ^= seed << 13;
        seed ^= seed >> 17;
        seed ^= seed << 5;
        return seed;
    };
    for (int32_t i = 0; i < 220; i++) {
        const auto x = static_cast<int32_t>(rnd() % W);
        const auto y = static_cast<int32_t>(rnd() % 170);
        const float b = 0.35f + static_cast<float>(rnd() % 100) / 150.0f;
        const int32_t size = (rnd() % 12 == 0) ? 2 : 1;
        img->fill_rect(x, y, size, size, [=](float, float, int32_t, int32_t) {
            return std::array<float, 4>{1.0f, 0.95f, 1.0f, b * (1.0f - static_cast<float>(y) / 190.0f)};
        });
    }

    // ─── sun, sliced into bands towards the horizon ───────────────────────────
    img->fill_circle_aa(W / 2, 170, 180, [](float u, float v, int32_t, int32_t) {
        const float d = std::sqrt((u - 0.5f) * (u - 0.5f) + (v - 0.5f) * (v - 0.5f)) * 2.0f;
        const float a = std::max(0.0f, 1.0f - d);
        return std::array<float, 4>{1.0f, 0.35f, 0.55f, a * a * 0.45f};  // halo
    });
    img->fill_circle_aa(W / 2, 170, 112, [](float, float v, int32_t, int32_t py) {
        float alpha = 1.0f;
        if (v > 0.45f) {
            const float band = (v - 0.45f) / 0.55f;  // 0 at the first gap, 1 at the bottom
            const int32_t period = 15;
            const auto gap = static_cast<int32_t>(2.0f + band * 9.0f);
            if ((py + 3) % period < gap) {
                alpha = 0.0f;
            }
        }
        return std::array<float, 4>{1.0f, 0.92f - 0.72f * v, 0.25f + 0.35f * v, alpha};
    });

    // ─── mountains ────────────────────────────────────────────────────────────
    auto ridge = [](float x, float scale, float phase) {
        float h = 0.0f;
        float amp = 1.0f;
        float freq = scale;
        for (int32_t o = 0; o < 5; o++) {
            h += amp * std::abs(std::sin(x * freq + phase * static_cast<float>(o + 1)));
            amp *= 0.5f;
            freq *= 2.13f;
        }
        return h / 1.94f;
    };
    auto mountains = [&](float scale, float phase, float height, float r, float g, float b) {
        for (int32_t x = 0; x < W; x++) {
            // Lower in the middle, so the sun stays visible.
            const float centre = std::abs(static_cast<float>(x) / W - 0.5f) * 2.0f;
            const float shape = 0.25f + 0.75f * centre * centre;
            const auto top = static_cast<int32_t>(ridge(static_cast<float>(x), scale, phase) * height * shape);
            img->fill_rect(x, HORIZON - top, 1, top, [=](float, float v, int32_t, int32_t) {
                const float edge = v < 0.04f ? 1.6f : 1.0f;  // lit ridge line
                return std::array<float, 4>{r * edge * (1.0f - v * 0.6f), g * edge, b * edge * (1.0f - v * 0.4f), 1.0f};
            });
        }
    };
    mountains(0.011f, 1.7f, 120.0f, 0.30f, 0.07f, 0.36f);
    mountains(0.019f, 4.1f, 78.0f, 0.12f, 0.03f, 0.20f);

    // ─── floor ────────────────────────────────────────────────────────────────
    img->fill_rect(0, HORIZON, W, H - HORIZON, [](float, float v, int32_t, int32_t) {
        const float glow = (1.0f - v) * (1.0f - v) * (1.0f - v);
        return std::array<float, 4>{0.04f + 0.50f * glow, 0.00f + 0.08f * glow, 0.10f + 0.30f * glow, 1.0f};
    });

    const uint8_t grid_glow = img->get_nearest_color(150, 20, 160);
    const uint8_t grid_core = img->get_nearest_color(255, 110, 255);
    auto grid_line = [&](int32_t x0, int32_t y0, int32_t x1, int32_t y1, float width) {
        img->draw_line_aa(x0, y0, x1, y1, grid_glow, width * 3.0f);
        img->draw_line_aa(x0, y0, x1, y1, grid_core, width);
    };
    // Lines towards the vanishing point.
    for (int32_t i = -16; i <= 16; i++) {
        grid_line(W / 2 + i * 14, HORIZON, W / 2 + i * 150, H + 40, 1.4f);
    }
    // Lines across, closer together towards the horizon.
    for (int32_t i = 1; i <= 14; i++) {
        const float t = static_cast<float>(i) / 14.0f;
        const auto y = HORIZON + static_cast<int32_t>(std::pow(t, 2.4f) * static_cast<float>(H - HORIZON));
        grid_line(0, y, W, y, 0.8f + t * 1.4f);
    }
    // Horizon line.
    img->fill_rect(0, HORIZON - 1, W, 3, [](float u, float, int32_t, int32_t) {
        const float c = 1.0f - std::abs(u - 0.5f) * 1.4f;
        return std::array<float, 4>{1.0f, 0.55f + 0.4f * c, 0.95f, 0.5f + 0.5f * c};
    });

    // ─── title ────────────────────────────────────────────────────────────────
    using title_font = ibmplexsans_bolditalic_48_aa;
    const char *title = "constixel";
    const int32_t tx = W / 2;
    const int32_t ty = 30;
    const uint8_t shade0 = img->get_nearest_color(90, 0, 130);
    const uint8_t shade1 = img->get_nearest_color(255, 40, 160);
    const uint8_t shade2 = img->get_nearest_color(120, 255, 255);
    for (int32_t d = 6; d >= 3; d--) {
        img->draw_string<title_font, true>(tx + d, ty + d, title, shade0, ALIGN_CENTER);
    }
    img->draw_string<title_font, true>(tx + 2, ty + 2, title, shade1, ALIGN_CENTER);
    img->draw_string<title_font, true>(tx, ty, title, shade2, ALIGN_CENTER);
    img->draw_string<m42_mono>(tx, ty + 66, "CONSTEXPR 2D GRAPHICS", color::WHITE, ALIGN_CENTER);

    // ─── scanlines and vignette ───────────────────────────────────────────────
    img->fill_rect(0, 0, W, H, [](float u, float v, int32_t, int32_t py) {
        const float d = (u - 0.5f) * (u - 0.5f) + (v - 0.5f) * (v - 0.5f);
        return std::array<float, 4>{0.0f, 0.0f, 0.0f, std::min(1.0f, d * 0.9f + ((py & 1) != 0 ? 0.16f : 0.0f))};
    });

    std::ofstream file("example_synthwave.png", std::ios::binary);
    img->png([&](char ch) {
        file.put(ch);
    });
}
