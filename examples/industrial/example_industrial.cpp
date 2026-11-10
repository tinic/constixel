#define CONSTIXEL_ENABLE_COUT

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>

#include "../example_output.hpp"
#include "../example_style.hpp"

using namespace constixel;
using namespace example;
using UI = image<format_8bit, 1024, 768>;

namespace {
struct Equipment {
    const char* name;
    const char* status;
    float temperature;
    float pressure;
    float flow;
    bool alarm;
};

void metric(UI& img, const Theme& theme, int x, int y, const char* label, float value, float maximum, const char* unit,
            uint8_t color) {
    char value_text[32];
    std::snprintf(value_text, sizeof(value_text), "%.1f %s", value, unit);
    text<small>(img, x, y, label, theme.muted);
    right_text<small>(img, x + 180, y, value_text, theme.ink);
    img.fill_round_rect_aa(x, y + 9, 180, 5, 2, theme.background);
    // Out-of-range values retain their numeric reading, but never overflow a track.
    const int width = static_cast<int>(180 * std::clamp(value / maximum, 0.0f, 1.0f));
    if (width > 0) {
        img.fill_rect(x, y + 9, width, 5, color);
    }
}

void equipment_card(UI& img, const Theme& theme, int x, const Equipment& equipment) {
    card(img, x, 138, 228, 224, theme);
    text<heading>(img, x + 24, 173, equipment.name, theme.ink);
    const auto status_color = equipment.alarm ? theme.danger : theme.accent;
    img.fill_circle_aa(x + 28, 194, 4, status_color);
    text<caption>(img, x + 40, 198, equipment.status, status_color);
    metric(img, theme, x + 24, 233, "Temperature", equipment.temperature, 100, "C", theme.blue);
    metric(img, theme, x + 24, 278, "Pressure", equipment.pressure, 10, "bar",
           equipment.alarm ? theme.danger : theme.accent);
    metric(img, theme, x + 24, 323, "Flow", equipment.flow, 50, "L/s", theme.blue);
}

void trend(UI& img, const Theme& theme) {
    card(img, 32, 386, 620, 246, theme);
    text<heading>(img, 56, 420, "Flow over time", theme.ink);
    text<small>(img, 56, 444, "Combined output / last 24 hours", theme.muted);
    right_text<heading>(img, 626, 421, "106.9 L/s", theme.accent);

    constexpr int left = 90, top = 470, width = 534, height = 106;
    for (int i = 0; i <= 2; ++i) {
        const int y = top + i * height / 2;
        img.draw_line(left, y, left + width, y, theme.background);
        const std::array<const char*, 3> labels = {"120", "80", "40"};
        right_text<caption>(img, left - 12, y + 4, labels[i], theme.muted);
    }
    auto sample = [](float t) {
        return 106.9f + 16 * (std::sin(t * 7) - std::sin(7.0f)) + 6 * (std::sin(t * 19) - std::sin(19.0f));
    };
    auto screen_y = [](float value) {
        return top + height - static_cast<int>((value - 40) / 80 * height);
    };
    int previous_y = screen_y(sample(0));
    for (int i = 1; i <= width; ++i) {
        const int y = screen_y(sample(static_cast<float>(i) / width));
        img.draw_line_aa(left + i - 1, previous_y, left + i, y, theme.accent, 3);
        previous_y = y;
    }
    img.fill_circle_aa(left + width, previous_y, 4, theme.accent);
    const std::array<const char*, 5> times = {"00:00", "06:00", "12:00", "18:00", "24:00"};
    for (int i = 0; i < 5; ++i) {
        centered_text<caption>(img, left + i * width / 4, 604, times[i], theme.muted);
    }
}
}  // namespace

int main(int argc, char** argv) {
    static UI img;
    const Theme theme(img);
    header(img, "Plant overview", "North facility / Process water", theme);
    text<small>(img, 32, 117, "EQUIPMENT", theme.muted);
    right_text<small>(img, 992, 117, "3 running   /   1 needs attention", theme.muted);

    constexpr std::array<Equipment, 4> equipment = {{
        {"Pump station 01", "Running", 72.5f, 8.2f, 35.7f, false},
        {"Pump station 02", "Pressure fault", 45.3f, 12.8f, 0.0f, true},
        {"Compressor A", "Running", 88.1f, 6.4f, 42.3f, false},
        {"Cooling tower", "Running", 34.7f, 2.1f, 28.9f, false},
    }};
    for (size_t i = 0; i < equipment.size(); ++i) {
        equipment_card(img, theme, 32 + static_cast<int>(i) * 244, equipment[i]);
    }
    trend(img, theme);

    card(img, 676, 386, 316, 246, theme);
    text<heading>(img, 700, 420, "System health", theme.ink);
    const std::array<const char*, 4> systems = {"Primary power", "Backup power", "Network", "Safety interlock"};
    const std::array<const char*, 4> states = {"Online", "Standby", "Connected", "Armed"};
    for (int i = 0; i < 4; ++i) {
        const int y = 461 + i * 42;
        text<small>(img, 700, y, systems[i], theme.muted);
        right_text<small>(img, 968, y, states[i], i == 1 ? theme.muted : theme.accent);
        if (i < 3) {
            img.draw_line(700, y + 14, 968, y + 14, theme.background);
        }
    }

    card(img, 32, 656, 960, 64, theme);
    img.fill_round_rect_aa(32, 656, 6, 64, 3, theme.danger);
    img.fill_circle_aa(62, 688, 10, theme.danger);
    centered_text<heading>(img, 62, 695, "!", theme.surface);
    text<heading>(img, 86, 683, "Pump station 02: pressure limit exceeded", theme.ink);
    text<small>(img, 86, 706, "12.8 bar measured / 10.0 bar limit. Inspect the outlet valve.", theme.muted);
    img.fill_round_rect_aa(824, 672, 144, 32, 8, theme.background);
    centered_text<small>(img, 896, 693, "View equipment", theme.ink);
    text<caption>(img, 32, 748, "PROCESS MONITOR", theme.muted);
    right_text<caption>(img, 992, 748, "Demo snapshot / 14:32:08", theme.muted);
    return output(img, argc, argv);
}
