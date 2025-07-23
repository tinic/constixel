#define CONSTIXEL_ENABLE_COUT

#include "constixel.hpp"
#include "fonts/ibmplexsans_regular_18_aa.hpp"
#include "fonts/notosanssymbols2_regular_aa.hpp"
#include <fstream>

using namespace constixel;
using UI = image<format_32bit, 396, 484>;

int main() {
    static UI img;
    img.clear();

    // ─── pick a soft, modern palette via nearest RGB ───────────────────────────
    uint8_t col_text       = img.get_nearest_color(236, 240, 241); // off-white
    uint8_t col_button_fg  = col_text;                             // white

    // ─── full‐screen background with dark gradient ────────────────────────────
    img.fill_rect(0, 0, 396, 484, [&](float /*u*/, float v, float /*au*/, float /*av*/) {
        return std::array<float, 4>{
            (8 + v * 35) / 255.0f,     // very dark to dark blue
            (8 + v * 35) / 255.0f,     // very dark to dark
            (15 + v * 50) / 255.0f,    // very dark to darker blue
            1.0f
        };
    });

    // ─── top header bar with dark gradient ────────────────────────────────────
    img.fill_round_rect_aa(0, 0, 396, 60, 0, [&](float /*u*/, float v, float /*au*/, float /*av*/) {
        return std::array<float, 4>{
            (5 + v * 40) / 255.0f,     // very dark to dark
            (5 + v * 40) / 255.0f,     // very dark to dark
            (12 + v * 55) / 255.0f,    // very dark to darker blue
            1.0f
        };
    });
    img.draw_string_aa<ibmplexsans_regular_18_aa,true>(20, 38, "Now Playing", col_text);

    // ─── album art placeholder ────────────────────────────────────────────────
    const int artSize = 220;
    const int ax = (396 - artSize) / 2;
    const int ay = 80;
    // outer frame
    img.fill_round_rect_aa(ax, ay, artSize, artSize, 20, [&](float u, float v, float /*au*/, float /*av*/) {
        float radial = std::sqrt((u - 0.5f) * (u - 0.5f) + (v - 0.5f) * (v - 0.5f));
        return std::array<float, 4>{
            (15 + radial * 60) / 255.0f,   // dark charcoal gradient
            (20 + radial * 70) / 255.0f,   // dark blue-gray
            (30 + radial * 80) / 255.0f,   // darker blue
            1.0f
        };
    });
    img.stroke_round_rect_aa(ax, ay, artSize, artSize, 20, col_text, 3);
    // inner “art”
    img.fill_round_rect_aa(ax+10, ay+10, artSize-20, artSize-20, 12, [&](float u, float v, float /*au*/, float /*av*/) {
        float center_dist = std::sqrt((u - 0.5f) * (u - 0.5f) + (v - 0.5f) * (v - 0.5f));
        float gradient = 1.0f - center_dist * 1.5f;
        return std::array<float, 4>{
            std::min(1.0f, std::max(0.0f, (10 + gradient * 80) / 255.0f)),  // dark teal
            std::min(1.0f, std::max(0.0f, (60 + gradient * 100) / 255.0f)), // darker teal
            std::min(1.0f, std::max(0.0f, (50 + gradient * 90) / 255.0f)),  // dark blue-green
            1.0f
        };
    });

    // ─── track & artist text ─────────────────────────────────────────────────
    img.draw_string_aa<ibmplexsans_regular_18_aa,true>(20, 320, "Artist Name", col_text);
    img.draw_string_aa<ibmplexsans_regular_18_aa,true>(20, 348, "Track Title", col_text);

    // ─── progress bar ─────────────────────────────────────────────────────────
    {
        const int pbX = 20, pbY = 380, pbW = 356, pbH = 8;
        const float progress = 0.45f;  // 45%
        // background track with subtle gradient
        img.fill_round_rect_aa(pbX, pbY, pbW, pbH, pbH/2, [&](float /*u*/, float v, float /*au*/, float /*av*/) {
            float highlight = (1.0f - v) * 0.7f;
            return std::array<float, 4>{
                (15 + highlight * 50) / 255.0f,   // dark gray gradient
                (20 + highlight * 60) / 255.0f,   // dark blue-gray
                (30 + highlight * 70) / 255.0f,   // darker blue
                1.0f
            };
        });
        // filled portion with warm gradient
        img.fill_round_rect_aa(pbX, pbY, int(pbW * progress), pbH, pbH/2, [&](float u, float v, float /*au*/, float /*av*/) {
            float glow = (1.0f - v) * 0.8f + u * 0.3f;
            return std::array<float, 4>{
                (60 + glow * 100) / 255.0f,    // dark orange
                (30 + glow * 80) / 255.0f,     // dark amber  
                (5 + glow * 25) / 255.0f,      // very dark
                1.0f
            };
        });
    }

    // ─── playback controls ────────────────────────────────────────────────────
    const int btnY    = 440;
    const int btnSize = 60;
    const int spacing = 90;
    const int cx      = 396 / 2;

    // helper to draw a circular button with gradient
    auto drawBtn = [&](int x, const char* label, bool isAccent, uint8_t fg){
        if (isAccent) {
            img.fill_circle_aa(x, btnY, btnSize/2, [&](float u, float v, float /*au*/, float /*av*/) {
                float center_dist = std::sqrt((u - 0.5f) * (u - 0.5f) + (v - 0.5f) * (v - 0.5f));
                float glow = 1.0f - center_dist * 2.0f;
                return std::array<float, 4>{
                    std::min(1.0f, std::max(0.0f, (8 + glow * 70) / 255.0f)),   // dark teal center
                    std::min(1.0f, std::max(0.0f, (40 + glow * 80) / 255.0f)),  // darker teal
                    std::min(1.0f, std::max(0.0f, (30 + glow * 70) / 255.0f)),  // dark blue-green
                    1.0f
                };
            });
        } else {
            img.fill_circle_aa(x, btnY, btnSize/2, [&](float /*u*/, float v, float /*au*/, float /*av*/) {
                float highlight = (1.0f - v) * 0.8f;
                return std::array<float, 4>{
                    (20 + highlight * 60) / 255.0f,   // dark charcoal
                    (25 + highlight * 70) / 255.0f,   // dark blue-gray
                    (35 + highlight * 80) / 255.0f,   // darker blue
                    1.0f
                };
            });
        }
        img.stroke_circle_aa(x, btnY, btnSize/2, fg, 3);
        img.draw_string_centered_aa<notosanssymbols2_regular_aa,true>(x, btnY - 24, label, fg);
    };

    drawBtn(cx - spacing, "⏮", false, col_button_fg);
    drawBtn(cx          , "⏯", true,  col_text);
    drawBtn(cx + spacing, "⏭", false, col_button_fg);

    // ─── save out PNG ─────────────────────────────────────────────────────────
    std::ofstream fout("example_watchui.png", std::ios::binary);
    img.png([&](char c){ fout.put(c); });

    return 0;
}
