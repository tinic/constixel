#define CONSTIXEL_ENABLE_COUT

#include <array>
#include <cmath>

#include "../example_output.hpp"
#include "../example_style.hpp"

using namespace constixel;
using namespace example;
using UI = image<format_32bit, 396, 484>;

int main(int argc, char** argv) {
    static UI img;
    const auto background = img.get_nearest_color(17, 17, 17);
    const auto surface = img.get_nearest_color(34, 34, 34);
    const auto muted = img.get_nearest_color(170, 170, 170);
    const auto accent = img.get_nearest_color(130, 220, 190);
    img.fill_rect(0, 0, 396, 484, background);

    // A compact player with generous touch targets and a clear title hierarchy.
    img.draw_line_aa(31, 32, 25, 38, muted, 2);
    img.draw_line_aa(25, 38, 31, 44, muted, 2);
    centered_text<small>(img, 198, 43, "ON REPEAT", muted);
    text<small>(img, 330, 43, "09:41", muted);

    // Procedural cover art: a warm sun over layered, cool-toned hills.
    // Shader channels are linear-light values, so low values keep the art rich.
    img.fill_round_rect_aa(90, 68, 216, 216, 16, [](float u, float v, int32_t, int32_t) {
        const float dx = u - 0.66f, dy = v - 0.30f;
        if (dx * dx + dy * dy < 0.021f) {
            return std::array<float, 4>{0.98f, 0.59f, 0.28f, 1.0f};
        }
        if (v > 0.80f + 0.08f * std::sin(u * 6.0f + 2.0f)) {
            return std::array<float, 4>{0.025f, 0.09f, 0.12f, 1.0f};
        }
        if (v > 0.65f + 0.10f * std::sin(u * 5.0f - 1.0f)) {
            return std::array<float, 4>{0.04f, 0.22f, 0.25f, 1.0f};
        }
        if (v > 0.48f + 0.12f * std::sin(u * 5.0f + 1.0f)) {
            return std::array<float, 4>{0.18f, 0.39f, 0.39f, 1.0f};
        }
        return std::array<float, 4>{0.64f - v * 0.12f, 0.76f - v * 0.10f, 0.73f - v * 0.08f, 1.0f};
    });
    centered_text<title>(img, 198, 318, "Distant shores", color::WHITE);
    centered_text<small>(img, 198, 342, "Low Tide / Coastal Sessions", muted);

    constexpr int progress_x = 32, progress_y = 363, progress_width = 332;
    constexpr int elapsed_seconds = 102, duration_seconds = 238;
    const int filled = progress_width * elapsed_seconds / duration_seconds;
    img.fill_round_rect_aa(progress_x, progress_y, progress_width, 4, 2, surface);
    img.fill_round_rect_aa(progress_x, progress_y, filled, 4, 2, accent);
    img.fill_circle_aa(progress_x + filled, progress_y + 2, 4, accent);
    text<caption>(img, 32, 388, "1:42", muted);
    right_text<caption>(img, 364, 388, "3:58", muted);

    // Draw playback symbols with primitives, independent of Unicode coverage.
    img.fill_circle_aa(198, 433, 29, accent);
    img.fill_round_rect_aa(187, 422, 7, 22, 2, background);
    img.fill_round_rect_aa(202, 422, 7, 22, 2, background);
    for (int direction : {-1, 1}) {
        const int cx = 198 + direction * 88;
        img.fill_rect(cx + direction * 12 - 1, 423, 3, 20, color::WHITE);
        for (int i = 0; i < 17; ++i) {
            const int x = cx + direction * (i - 8);
            const int half_height = (17 - i) * 10 / 17;
            img.draw_line(x, 433 - half_height, x, 433 + half_height, color::WHITE);
        }
    }
    return output(img, argc, argv, "example_watchui.png");
}
