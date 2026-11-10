#define CONSTIXEL_ENABLE_COUT

#include <algorithm>
#include <array>
#include <cmath>

#include "../example_output.hpp"
#include "../example_style.hpp"

using namespace constixel;
using namespace example;
using Canvas = image<format_8bit, 1024, 768>;

namespace {
struct Point3 {
    float x, y, z;
};
struct Point2 {
    int x, y;
};

float surface(float x, float z) {
    const float radius = std::sqrt(x * x + z * z);
    return 0.6f * std::sin(radius * 2) * std::exp(-0.3f * radius) + 0.3f * std::sin(x * 3) * std::cos(z * 2);
}

Point2 project(Point3 point) {
    // A fixed orthographic camera keeps the entire [-3, 3] domain in view.
    constexpr float yaw = 0.65f, tilt = 0.62f, scale = 88.0f;
    const float x = point.x * std::cos(yaw) + point.z * std::sin(yaw);
    const float z = -point.x * std::sin(yaw) + point.z * std::cos(yaw);
    const float y = point.y * std::cos(tilt) - z * std::sin(tilt);
    return {512 + static_cast<int>(x * scale), 386 - static_cast<int>(y * scale)};
}

void segment(Canvas& img, Point3 a, Point3 b, uint8_t color) {
    const auto p = project(a), q = project(b);
    img.draw_line_aa(p.x, p.y, q.x, q.y, color);
}
}  // namespace

int main(int argc, char** argv) {
    static Canvas img;
    const Theme theme(img);
    header(img, "Surface study", "Damped radial waves with a cross-wave perturbation", theme);
    card(img, 32, 108, 960, 568, theme);
    text<small>(img, 56, 140, "ORTHOGRAPHIC VIEW", theme.muted);
    right_text<small>(img, 968, 140, "x, z in [-3, 3]", theme.muted);

    // A subtle floor grid establishes orientation without competing with the mesh.
    for (int i = -3; i <= 3; ++i) {
        const float t = static_cast<float>(i);
        segment(img, {-3, -0.8f, t}, {3, -0.8f, t}, theme.background);
        segment(img, {t, -0.8f, -3}, {t, -0.8f, 3}, theme.background);
    }
    constexpr int resolution = 36;
    constexpr float step = 6.0f / resolution;
    float min_height = 1, max_height = -1;
    for (int i = 0; i <= resolution; ++i) {
        for (int j = 0; j <= resolution; ++j) {
            const float y = surface(-3 + i * step, -3 + j * step);
            min_height = std::min(min_height, y);
            max_height = std::max(max_height, y);
        }
    }
    const auto low = img.get_nearest_color(70, 110, 180);
    const auto middle = theme.accent;
    const auto high = theme.warm;
    auto height_color = [&](float y) {
        const float t = (y - min_height) / (max_height - min_height);
        return t < 0.35f ? low : (t < 0.7f ? middle : high);
    };
    for (int i = 0; i <= resolution; ++i) {
        for (int j = 0; j <= resolution; ++j) {
            const float x = -3 + i * step, z = -3 + j * step;
            const Point3 a{x, surface(x, z), z};
            if (i < resolution) {
                const Point3 b{x + step, surface(x + step, z), z};
                segment(img, a, b, height_color((a.y + b.y) / 2));
            }
            if (j < resolution) {
                const Point3 b{x, surface(x, z + step), z + step};
                segment(img, a, b, height_color((a.y + b.y) / 2));
            }
        }
    }

    // Keep the orientation key separate from the surface and its title.
    img.draw_line_aa(76, 625, 111, 641, theme.muted, 2);
    img.draw_line_aa(76, 625, 102, 607, theme.muted, 2);
    img.draw_line_aa(76, 625, 76, 589, theme.muted, 2);
    text<caption>(img, 117, 647, "z", theme.muted);
    text<caption>(img, 109, 608, "x", theme.muted);
    centered_text<caption>(img, 76, 582, "y", theme.muted);
    text<small>(img, 652, 643, "Height", theme.muted);
    const std::array<const char*, 3> labels = {"Low", "Mid", "High"};
    const std::array<uint8_t, 3> colors = {low, middle, high};
    for (int i = 0; i < 3; ++i) {
        img.fill_circle_aa(728 + i * 84, 638, 5, colors[i]);
        text<small>(img, 740 + i * 84, 643, labels[i], theme.ink);
    }
    text<small>(img, 32, 711, "f(x,z) = 0.6 sin(2r) exp(-0.3r) + 0.3 sin(3x) cos(2z)", theme.ink);
    text<caption>(img, 32, 739, "r = sqrt(x*x + z*z)", theme.muted);
    right_text<small>(img, 992, 711, "36 x 36 cells", theme.muted);
    right_text<caption>(img, 992, 739, "2D lines / 3D projection", theme.muted);
    return output(img, argc, argv);
}
