#define CONSTIXEL_ENABLE_COUT

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <numbers>

#include "../example_output.hpp"
#include "../example_style.hpp"
#include "fonts/ibmplexsans_semibold_48_aa.hpp"

using namespace constixel;
using namespace example;
using Canvas = image<format_32bit, 1024, 720>;
using large_number = ibmplexsans_semibold_48_aa;

namespace {
// Power is in watts. The sample balances: generation = home + battery + export.
struct Power {
    int solar, home, charging, exported, battery_percent;
};
constexpr Power power{4800, 1900, 1600, 1300, 76};
constexpr int battery_capacity_wh = 10000;
static_assert(power.solar == power.home + power.charging + power.exported);

void power_text(Canvas& img, int x, int baseline, int watts, uint8_t color) {
    char label[24];
    std::snprintf(label, sizeof(label), "%.1f kW", watts / 1000.0);
    right_text<heading>(img, x, baseline, label, color);
}

void sun(Canvas& img, int x, int y, uint8_t color) {
    img.fill_circle_aa(x, y, 15, color);
    for (int i = 0; i < 8; ++i) {
        const float angle = i * std::numbers::pi_v<float> / 4;
        img.draw_line_aa(x + static_cast<int>(22 * std::cos(angle)), y + static_cast<int>(22 * std::sin(angle)),
                         x + static_cast<int>(28 * std::cos(angle)), y + static_cast<int>(28 * std::sin(angle)), color,
                         3);
    }
}

void flow_diagram(Canvas& img, const Theme& theme) {
    card(img, 360, 108, 632, 280, theme);
    text<heading>(img, 384, 145, "Where your power goes", theme.ink);
    right_text<small>(img, 968, 145, "LIVE", theme.accent);

    img.draw_line_aa(500, 260, 558, 260, theme.warm, 3);
    img.draw_line_aa(558, 192, 558, 328, theme.border, 3);
    img.fill_round_rect_aa(396, 232, 104, 56, 12, theme.background);
    centered_text<caption>(img, 448, 252, "SOLAR", theme.muted);
    power_text(img, 486, 276, power.solar, theme.ink);
    constexpr std::array<const char*, 3> names = {"Home", "Battery / charging", "Grid / exporting"};
    const std::array<int, 3> watts = {power.home, power.charging, power.exported};
    const std::array<uint8_t, 3> colors = {theme.accent, theme.blue, theme.warm};
    for (int i = 0; i < 3; ++i) {
        const int y = 192 + i * 68;
        img.draw_line_aa(558, y, 612, y, colors[i], 3);
        img.draw_line_aa(606, y - 5, 612, y, colors[i], 2);
        img.draw_line_aa(606, y + 5, 612, y, colors[i], 2);
        img.fill_circle_aa(558, y, 4, colors[i]);
        img.fill_round_rect_aa(624, y - 26, 344, 52, 10, theme.background);
        img.fill_circle_aa(646, y, 5, colors[i]);
        text<small>(img, 664, y + 5, names[i], theme.ink);
        power_text(img, 948, y + 6, watts[i], theme.ink);
    }
}

void history(Canvas& img, const Theme& theme) {
    card(img, 32, 412, 620, 246, theme);
    text<heading>(img, 56, 447, "A day in balance", theme.ink);
    text<small>(img, 56, 472, "Yesterday / hourly average power", theme.muted);
    img.fill_circle_aa(455, 442, 4, theme.warm);
    text<caption>(img, 465, 447, "Solar", theme.muted);
    img.fill_circle_aa(527, 442, 4, theme.blue);
    text<caption>(img, 537, 447, "Home", theme.muted);
    constexpr std::array<float, 24> solar = {0,    0,    0,    0,    0,    0.2f, 0.7f, 1.5f, 2.7f, 3.8f, 4.7f, 5.2f,
                                             5.5f, 5.2f, 4.8f, 4.0f, 2.9f, 1.6f, 0.6f, 0.1f, 0,    0,    0,    0};
    constexpr std::array<float, 24> home = {0.5f, 0.5f, 0.4f, 0.4f, 0.5f, 0.8f, 1.9f, 2.8f, 2.1f, 1.4f, 1.1f, 1.3f,
                                            1.8f, 1.4f, 1.1f, 1.0f, 1.4f, 2.6f, 3.4f, 2.9f, 2.2f, 1.5f, 0.9f, 0.6f};
    constexpr int left = 90, top = 496, bottom = 606, step = 22;
    for (int i = 0; i <= 2; ++i) {
        const int y = top + i * 55;
        img.draw_line(left, y, 618, y, theme.background);
        const std::array<const char*, 3> labels = {"6 kW", "3 kW", "0"};
        right_text<caption>(img, 78, y + 4, labels[i], theme.muted);
    }
    for (size_t i = 0; i < solar.size(); ++i) {
        const int x = left + static_cast<int>(i) * step;
        const int h = static_cast<int>(solar[i] / 6 * (bottom - top));
        if (h > 0) {
            img.fill_round_rect_aa(x + 3, bottom - h, 16, h, 3, theme.warm);
        }
    }
    for (size_t i = 1; i < home.size(); ++i) {
        const int x = left + static_cast<int>(i) * step + 11;
        const int previous = bottom - static_cast<int>(home[i - 1] / 6 * (bottom - top));
        const int current = bottom - static_cast<int>(home[i] / 6 * (bottom - top));
        img.draw_line_aa(x - step, previous, x, current, theme.blue, 3);
    }
    for (int hour : {0, 6, 12, 18, 24}) {
        char label[8];
        std::snprintf(label, sizeof(label), "%02d:00", hour);
        centered_text<caption>(img, left + hour * step, 635, label, theme.muted);
    }
}
}  // namespace

int main(int argc, char** argv) {
    static Canvas img;
    const Theme theme(img);
    header(img, "Home energy", "South roof array / A live view of your power", theme);

    // Linear-light shader color gives the summary panel a soft, warm surface.
    img.fill_round_rect_aa(32, 108, 304, 280, 16, [](float u, float v, int32_t, int32_t) {
        return std::array<float, 4>{0.98f, 0.82f - v * 0.12f, 0.52f + u * 0.12f, 1.0f};
    });
    text<small>(img, 56, 146, "GENERATING NOW", theme.ink);
    sun(img, 279, 178, theme.warm);
    char reading[16];
    std::snprintf(reading, sizeof(reading), "%.1f", power.solar / 1000.0);
    text<large_number>(img, 56, 218, reading, theme.ink);
    text<heading>(img, 146, 216, "kW", theme.ink);
    text<small>(img, 56, 249, "Your home is running on sunshine.", theme.ink);
    img.draw_line(56, 276, 312, 276, theme.warm);
    text<small>(img, 56, 307, "Solar today", theme.ink);
    right_text<title>(img, 312, 309, "24.6 kWh", theme.ink);
    text<caption>(img, 56, 352, "UPDATED 12:40 / DEMO SNAPSHOT", theme.ink);
    flow_diagram(img, theme);
    history(img, theme);

    card(img, 676, 412, 316, 246, theme);
    text<heading>(img, 700, 447, "Stored for later", theme.ink);
    text<small>(img, 700, 472, "Home battery", theme.muted);
    img.stroke_round_rect_aa(700, 504, 60, 100, 8, theme.border, 3);
    img.fill_round_rect_aa(717, 497, 26, 6, 2, theme.border);
    const int charge_height = 86 * std::clamp(power.battery_percent, 0, 100) / 100;
    img.fill_round_rect_aa(707, 597 - charge_height, 46, charge_height, 4, theme.accent);
    std::snprintf(reading, sizeof(reading), "%d%%", power.battery_percent);
    text<large_number>(img, 784, 555, reading, theme.ink);
    char detail[48];
    std::snprintf(detail, sizeof(detail), "Charging at %.1f kW", power.charging / 1000.0);
    text<small>(img, 786, 583, detail, theme.accent);
    std::snprintf(detail, sizeof(detail), "Energy stored: %.1f / %.1f kWh",
                  battery_capacity_wh * power.battery_percent / 100000.0, battery_capacity_wh / 1000.0);
    text<caption>(img, 700, 635, detail, theme.muted);
    text<small>(img, 32, 694, "Generation = home use + battery charging + grid export", theme.muted);
    right_text<small>(img, 992, 694, "100% solar-powered right now", theme.accent);
    return output(img, argc, argv);
}
