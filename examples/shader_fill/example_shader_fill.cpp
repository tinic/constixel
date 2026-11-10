#define CONSTIXEL_ENABLE_COUT

#include <algorithm>
#include <array>
#include <cmath>

#include "../example_output.hpp"
#include "../example_style.hpp"

using namespace constixel;
using namespace example;
using Canvas = image<format_32bit, 1024, 908>;
using RGBA = std::array<float, 4>;

namespace {
// Interpolate a coordinated set of colors in linear light.
RGBA gradient(float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    constexpr RGBA navy{0.015f, 0.05f, 0.14f, 1.0f};
    constexpr RGBA teal{0.035f, 0.34f, 0.36f, 1.0f};
    constexpr RGBA sand{0.93f, 0.66f, 0.38f, 1.0f};
    const auto& from = t < 0.5f ? navy : teal;
    const auto& to = t < 0.5f ? teal : sand;
    const float blend = t < 0.5f ? t * 2 : (t - 0.5f) * 2;
    return {std::lerp(from[0], to[0], blend), std::lerp(from[1], to[1], blend), std::lerp(from[2], to[2], blend), 1.0f};
}

RGBA shader(int pattern, float u, float v) {
    const float dx = u - 0.5f, dy = v - 0.5f;
    const float radius = std::sqrt(dx * dx + dy * dy);
    const float angle = std::atan2(dy, dx);
    switch (pattern) {
        case 0:
            return gradient(u);
        case 1:
            return gradient(std::max(0.0f, 1.0f - radius * 1.8f));
        case 2:
            return gradient((static_cast<int>(u * 12) + static_cast<int>(v * 6)) % 2 ? 0.8f : 0.2f);
        case 3:
            return gradient(std::max(0.0f, 1.0f - radius * 2));
        case 4:
            return gradient(0.5f + 0.5f * std::sin((u + v) * 16));
        case 5:
            return gradient(0.5f + 0.5f * std::sin(angle * 3 + radius * 18));
        case 6:
            return {0.08f + u * 0.55f, 0.16f + v * 0.45f, 0.30f + (1 - u) * 0.22f, 1.0f};
        case 7:
            return gradient(0.5f + (std::sin(u * 9) + std::sin(v * 9)) * 0.25f);
        default:
            return gradient(0.5f + std::sin(radius * 22) * std::exp(-radius * 3) * 0.5f);
    }
}
}  // namespace

int main(int argc, char** argv) {
    static Canvas img;
    const Theme theme(img);
    header(img, "Color in motion", "Procedural fills / one function, every pixel", theme);
    constexpr std::array<const char*, 9> names = {"Linear gradient", "Radial light",   "Checkerboard",
                                                  "Soft sphere",     "Diagonal bands", "Spiral",
                                                  "Color field",     "Interference",   "Damped ripple"};
    for (int i = 0; i < 9; ++i) {
        const int x = 32 + i % 3 * 328;
        const int y = 108 + i / 3 * 184;
        card(img, x, y, 304, 168, theme);
        auto fill = [i](float u, float v, int32_t, int32_t) {
            return shader(i, u, v);
        };
        if (i < 3) {
            img.fill_rect(x + 20, y + 18, 264, 104, fill);
        } else if (i < 6) {
            img.fill_circle_aa(x + 152, y + 70, 52, fill);
        } else {
            img.fill_round_rect_aa(x + 20, y + 18, 264, 104, 16, fill);
        }
        text<small>(img, x + 20, y + 149, names[i], theme.ink);
        right_text<caption>(img, x + 284, y + 148, i < 3 ? "RECT" : (i < 6 ? "CIRCLE" : "ROUND RECT"), theme.muted);
    }

    card(img, 32, 660, 632, 196, theme);
    img.fill_round_rect_aa(52, 680, 592, 126, 12, [](float u, float v, int32_t, int32_t) {
        const float x = (u - 0.5f) * 8, y = (v - 0.5f) * 4;
        const float radius = std::sqrt(x * x + y * y);
        const float ripple = std::sin(radius * 4 - 8) * std::exp(-radius * 0.25f);
        const float waves = std::sin(u * 16) * std::sin(v * 12);
        return gradient(0.5f + ripple * 0.3f + waves * 0.2f);
    });
    text<heading>(img, 52, 835, "Wave study", theme.ink);
    right_text<small>(img, 644, 835, "Layered ripples + interference", theme.muted);

    card(img, 688, 660, 304, 196, theme);
    img.fill_circle_aa(840, 744, 64, [](float u, float v, int32_t, int32_t) {
        const float x = u - 0.5f, y = v - 0.5f;
        const float radius = std::sqrt(x * x + y * y);
        return gradient(0.5f + 0.5f * std::sin(std::atan2(y, x) * 6 + radius * 18));
    });
    text<small>(img, 708, 835, "Radial symmetry", theme.ink);
    right_text<caption>(img, 972, 835, "POLAR", theme.muted);
    text<small>(img, 32, 890, "Normalized coordinates / linear-light RGBA", theme.muted);
    right_text<small>(img, 992, 890, "32-bit color / antialiased edges", theme.muted);
    // Keep the original iTerm preview when no PNG filename is supplied.
    if (argc == 1) {
#ifdef _MSC_VER
        img.sixel_to_cout<1>();
#else
        img.png_to_iterm();
#endif
        return 0;
    }
    return output(img, argc, argv);
}
