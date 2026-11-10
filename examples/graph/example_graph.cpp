#define CONSTIXEL_ENABLE_COUT

#include <cmath>
#include <cstdio>

#include "../example_output.hpp"
#include "../example_style.hpp"

using namespace constixel;
using namespace example;
using Plot = image<format_8bit, 1024, 720>;

int main(int argc, char** argv) {
    static Plot img;
    const Theme theme(img);
    header(img, "Signal explorer", "Two waveforms, one shared time base", theme);
    card(img, 32, 108, 960, 548, theme);
    text<heading>(img, 60, 145, "Amplitude", theme.ink);
    img.draw_line_aa(608, 138, 632, 138, theme.accent, 3);
    text<small>(img, 642, 143, "Damped sine", theme.ink);
    img.draw_line_aa(796, 138, 820, 138, theme.blue, 3);
    text<small>(img, 830, 143, "Modulated wave", theme.ink);

    constexpr int left = 100, top = 184, width = 848, height = 380;
    constexpr float xmin = 0.0f, xmax = 6.0f, ymin = -1.2f, ymax = 1.2f;
    auto screen_y = [](float y) {
        return top + static_cast<int>((ymax - y) / (ymax - ymin) * height);
    };

    for (int i = 0; i <= 6; ++i) {
        const int x = left + i * width / 6;
        const int y = top + i * height / 6;
        img.draw_line(left, y, left + width, y, i == 3 ? theme.border : theme.background);
        img.draw_line(x, top, x, top + height, theme.background);
        char label[16];
        std::snprintf(label, sizeof(label), "%.1f", ymax - static_cast<float>(i) * 0.4f);
        right_text<small>(img, left - 16, y + 5, label, theme.muted);
        std::snprintf(label, sizeof(label), "%d", i);
        centered_text<small>(img, x, top + height + 29, label, theme.muted);
    }
    centered_text<small>(img, left + width / 2, 630, "Time (seconds)", theme.muted);

    auto damped_sine = [](float x) {
        return std::exp(-0.4f * x) * std::sin(5.0f * x);
    };
    auto modulated_wave = [](float x) {
        return 0.6f * std::sin(2.0f * x + 1.0f) * std::cos(0.5f * x);
    };
    auto draw_wave = [&](auto sample, uint8_t color) {
        int previous_y = screen_y(sample(xmin));
        for (int i = 1; i <= width; ++i) {
            const float x = xmin + (xmax - xmin) * static_cast<float>(i) / width;
            const int y = screen_y(sample(x));
            img.draw_line_aa(left + i - 1, previous_y, left + i, y, color, 3);
            previous_y = y;
        }
        img.fill_circle_aa(left + width, previous_y, 4, color);
    };
    draw_wave(damped_sine, theme.accent);
    draw_wave(modulated_wave, theme.blue);

    text<small>(img, 32, 691, "Damped sine   exp(-0.4t) sin(5t)", theme.muted);
    right_text<small>(img, 992, 691, "Modulated wave   0.6 sin(2t + 1) cos(0.5t)", theme.muted);
    return output(img, argc, argv);
}
