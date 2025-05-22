#include "constixel.hpp"
#include "fonts/ibmplexsans_regular_18_aa.hpp"
#include "fonts/notosanssymbols2_regular_aa.hpp"
#include <fstream>

using namespace constixel;
using UI = image<format_32bit, 396, 484>;

int main() {
    UI img;
    img.clear();

    // ─── pick a soft, modern palette via nearest RGB ───────────────────────────
    uint8_t col_bg         = img.get_nearest_color(32, 32, 48);    // dark slate
    uint8_t col_header     = img.get_nearest_color(28, 28, 40);    // slightly lighter
    uint8_t col_text       = img.get_nearest_color(236, 240, 241); // off-white
    uint8_t col_accent     = img.get_nearest_color(26, 188, 156);  // teal
    uint8_t col_track_bg   = img.get_nearest_color(44, 62, 80);    // muted blue
    uint8_t col_progress   = img.get_nearest_color(230, 126, 34);  // orange
    uint8_t col_button_bg  = img.get_nearest_color(60, 76, 100);   // charcoal
    uint8_t col_button_fg  = col_text;                             // white

    // ─── full‐screen background ───────────────────────────────────────────────
    img.fill_rect(0, 0, 396, 484, col_bg);

    // ─── top header bar ───────────────────────────────────────────────────────
    img.fill_round_rect_aa(0, 0, 396, 60, 0, col_header);
    img.draw_string_aa<ibmplexsans_regular_18_aa,true>(20, 38, "Now Playing", col_text);

    // ─── album art placeholder ────────────────────────────────────────────────
    const int artSize = 220;
    const int ax = (396 - artSize) / 2;
    const int ay = 80;
    // outer frame
    img.fill_round_rect_aa(ax, ay, artSize, artSize, 20, col_track_bg);
    img.stroke_round_rect_aa(ax, ay, artSize, artSize, 20, col_text, 3);
    // inner “art”
    img.fill_round_rect_aa(ax+10, ay+10, artSize-20, artSize-20, 12, col_accent);

    // ─── track & artist text ─────────────────────────────────────────────────
    img.draw_string_aa<ibmplexsans_regular_18_aa,true>(20, 320, "Artist Name", col_text);
    img.draw_string_aa<ibmplexsans_regular_18_aa,true>(20, 348, "Track Title", col_text);

    // ─── progress bar ─────────────────────────────────────────────────────────
    {
        const int pbX = 20, pbY = 380, pbW = 356, pbH = 8;
        const float progress = 0.45f;  // 45%
        // background track
        img.fill_round_rect_aa(pbX, pbY, pbW, pbH, pbH/2, col_track_bg);
        // filled portion
        img.fill_round_rect_aa(pbX, pbY, int(pbW * progress), pbH, pbH/2, col_progress);
    }

    // ─── playback controls ────────────────────────────────────────────────────
    const int btnY    = 440;
    const int btnSize = 60;
    const int spacing = 90;
    const int cx      = 396 / 2;

    // helper to draw a circular button with label
    auto drawBtn = [&](int x, const char* label, uint8_t bg, uint8_t fg){
        img.fill_circle_aa(x, btnY, btnSize/2, bg);
        img.stroke_circle_aa(x, btnY, btnSize/2, fg, 3);
        img.draw_string_centered_aa<notosanssymbols2_regular_aa,true>(x, btnY - 24, label, fg);
    };

    drawBtn(cx - spacing, "⏮", col_button_bg, col_button_fg);
    drawBtn(cx          , "⏯", col_accent,    col_text);
    drawBtn(cx + spacing, "⏭", col_button_bg, col_button_fg);

    // ─── save out PNG ─────────────────────────────────────────────────────────
    std::ofstream fout("example_watchui.png", std::ios::binary);
    img.png([&](char c){ fout.put(c); });

    return 0;
}
