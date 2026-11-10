#define CONSTIXEL_ENABLE_COUT

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>

#include "../example_output.hpp"
#include "../example_style.hpp"
#include "fonts/ibmplexsans_semibold_48_aa.hpp"

using namespace constixel;
using namespace example;
using Screen = image<format_8bit, 480, 480>;
using temperature_font = ibmplexsans_semibold_48_aa;

namespace {
struct Climate {
    float measured, target;
    int humidity, co2;
};
constexpr Climate climate{21.8f, 22.5f, 46, 640};

// A circular dial needs only short antialiased line segments. Angles increase
// clockwise on screen, from the lower left through the top to the lower right.
void arc(Screen& img, int cx, int cy, int radius, float start, float sweep, uint8_t color, float width) {
    constexpr float radians = std::numbers::pi_v<float> / 180;
    const int segments = std::max(1, static_cast<int>(std::abs(sweep)));
    auto x = [&](float angle) {
        return cx + static_cast<int>(std::lround(radius * std::cos(angle * radians)));
    };
    auto y = [&](float angle) {
        return cy + static_cast<int>(std::lround(radius * std::sin(angle * radians)));
    };
    for (int i = 1; i <= segments; ++i) {
        const float a = start + sweep * (i - 1) / segments;
        const float b = start + sweep * i / segments;
        img.draw_line_aa(x(a), y(a), x(b), y(b), color, width);
    }
}
}  // namespace

int main(int argc, char** argv) {
    static Screen img;
    const Theme theme(img);
    const auto warm = theme.warm;
    img.fill_rect(0, 0, 480, 480, theme.background);
    text<title>(img, 28, 44, "Living room", theme.ink);
    text<small>(img, 28, 70, "Climate / Following your schedule", theme.muted);
    img.fill_round_rect_aa(362, 25, 90, 28, 14, theme.surface);
    img.fill_circle_aa(378, 39, 4, theme.accent);
    text<caption>(img, 390, 43, "ONLINE", theme.ink);

    card(img, 28, 92, 424, 272, theme);
    constexpr int cx = 240, cy = 230, radius = 102;
    constexpr float minimum = 16, maximum = 28, start = 135, sweep = 270;
    const float fraction = std::clamp((climate.target - minimum) / (maximum - minimum), 0.0f, 1.0f);
    arc(img, cx, cy, radius, start, sweep, theme.background, 10);
    arc(img, cx, cy, radius, start, sweep * fraction, warm, 10);
    for (int i = 0; i <= 24; ++i) {
        const float angle = (start + sweep * i / 24) * std::numbers::pi_v<float> / 180;
        const int inner = i % 4 == 0 ? 117 : 121;
        img.draw_line_aa(cx + static_cast<int>(inner * std::cos(angle)), cy + static_cast<int>(inner * std::sin(angle)),
                         cx + static_cast<int>(126 * std::cos(angle)), cy + static_cast<int>(126 * std::sin(angle)),
                         theme.border, 2);
    }
    const float knob_angle = (start + sweep * fraction) * std::numbers::pi_v<float> / 180;
    const int knob_x = cx + static_cast<int>(std::lround(radius * std::cos(knob_angle)));
    const int knob_y = cy + static_cast<int>(std::lround(radius * std::sin(knob_angle)));
    img.fill_circle_aa(knob_x, knob_y, 9, warm);
    img.fill_circle_aa(knob_x, knob_y, 4, theme.surface);

    centered_text<caption>(img, cx, 187, "TARGET TEMPERATURE", theme.muted);
    char value[32];
    std::snprintf(value, sizeof(value), "%.1f", climate.target);
    centered_text<temperature_font>(img, cx - 8, 243, value, theme.ink);
    // Draw the degree marker explicitly: the bundled Latin font covers ASCII.
    img.stroke_circle_aa(cx + 55, 207, 4, theme.ink, 2);
    text<small>(img, cx + 65, 219, "C", theme.ink);
    const bool heating = climate.measured < climate.target;
    centered_text<small>(img, cx, 274, heating ? "Heating to your target" : "Heating paused",
                         heating ? warm : theme.muted);
    std::snprintf(value, sizeof(value), "Room temperature  %.1f C", climate.measured);
    centered_text<small>(img, cx, 347, value, theme.muted);

    img.fill_circle_aa(75, 314, 20, theme.background);
    img.draw_line_aa(68, 314, 82, 314, theme.ink, 2);
    img.fill_circle_aa(405, 314, 20, theme.background);
    img.draw_line_aa(398, 314, 412, 314, theme.ink, 2);
    img.draw_line_aa(405, 307, 405, 321, theme.ink, 2);
    centered_text<caption>(img, 155, 322, "16", theme.muted);
    centered_text<caption>(img, 325, 322, "28", theme.muted);

    card(img, 28, 380, 204, 68, theme);
    card(img, 248, 380, 204, 68, theme);
    text<caption>(img, 46, 403, "HUMIDITY", theme.muted);
    std::snprintf(value, sizeof(value), "%d%%", climate.humidity);
    text<title>(img, 46, 433, value, theme.ink);
    text<caption>(img, 266, 403, "CO2", theme.muted);
    std::snprintf(value, sizeof(value), "%d ppm", climate.co2);
    text<title>(img, 266, 433, value, theme.ink);
    text<caption>(img, 28, 470, "NEXT: 20.0 C AT 22:00", theme.muted);
    right_text<caption>(img, 452, 470, "DEMO / 8-BIT DISPLAY", theme.muted);
    return output(img, argc, argv);
}
