#define CONSTIXEL_ENABLE_COUT

#include "constixel.hpp"

#include "../example_output.hpp"

#include "fonts/ibmplexsans_regular_14_aa.hpp"
#include "fonts/ibmplexsans_semibold_24_aa.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using namespace constixel;

//
// A small ray tracer renders a scene into a RGBA buffer. The buffer is then blitted into every palette format with
// blit_RGBA_diffused_linear(), which shows what the same picture costs and looks like at 1, 2, 4 and 8 bits per pixel.
//

static constexpr int32_t W = 320;
static constexpr int32_t H = 240;

struct vec {
    float x = 0, y = 0, z = 0;
    constexpr vec operator+(const vec &o) const { return {x + o.x, y + o.y, z + o.z}; }
    constexpr vec operator-(const vec &o) const { return {x - o.x, y - o.y, z - o.z}; }
    constexpr vec operator*(float s) const { return {x * s, y * s, z * s}; }
    constexpr vec operator*(const vec &o) const { return {x * o.x, y * o.y, z * o.z}; }
    [[nodiscard]] constexpr float dot(const vec &o) const { return x * o.x + y * o.y + z * o.z; }
    [[nodiscard]] vec norm() const { return *this * (1.0f / std::sqrt(dot(*this))); }
};

struct sphere {
    vec center;
    float radius;
    vec albedo;
    float mirror;
};

static constexpr std::array<sphere, 4> spheres{{
    {{-1.9f, 0.9f, 5.6f}, 0.9f, {0.95f, 0.25f, 0.20f}, 0.15f},
    {{0.0f, 1.2f, 6.6f}, 1.2f, {0.90f, 0.90f, 0.95f}, 0.75f},
    {{1.9f, 0.8f, 5.2f}, 0.8f, {0.20f, 0.55f, 0.95f}, 0.15f},
    {{0.5f, 0.35f, 4.1f}, 0.35f, {0.98f, 0.80f, 0.15f}, 0.10f},
}};

static const vec light = vec{-0.5f, 0.9f, -0.45f}.norm();

static bool hit(const vec &o, const vec &d, float &t, vec &n, vec &albedo, float &mirror) {
    t = 1e9f;
    if (d.y < -1e-4f) {  // floor, checkered
        t = -o.y / d.y;
        const vec p = o + d * t;
        const bool c = ((static_cast<int32_t>(std::floor(p.x)) + static_cast<int32_t>(std::floor(p.z))) & 1) != 0;
        n = {0, 1, 0};
        albedo = c ? vec{0.85f, 0.85f, 0.85f} : vec{0.18f, 0.20f, 0.25f};
        mirror = 0.25f;
    }
    for (const auto &s : spheres) {
        const vec oc = o - s.center;
        const float b = oc.dot(d);
        const float disc = b * b - oc.dot(oc) + s.radius * s.radius;
        if (disc > 0) {
            const float ts = -b - std::sqrt(disc);
            if (ts > 1e-3f && ts < t) {
                t = ts;
                n = (o + d * t - s.center).norm();
                albedo = s.albedo;
                mirror = s.mirror;
            }
        }
    }
    return t < 1e8f;
}

static vec trace(const vec &o, const vec &d, int depth) {
    float t = 0, mirror = 0;
    vec n, albedo;
    if (!hit(o, d, t, n, albedo, mirror)) {
        const float h = std::max(d.y, 0.0f);
        const float k = std::pow(1.0f - h, 6.0f);
        return vec{1.00f, 0.50f, 0.22f} * k + vec{0.04f, 0.10f, 0.38f} * (1.0f - k);  // sunset sky
    }
    const vec p = o + d * t;
    float ts = 0, ms = 0;
    vec ns, as;
    const bool shadow = hit(p + n * 1e-3f, light, ts, ns, as, ms);
    const float diffuse = shadow ? 0.0f : std::max(n.dot(light), 0.0f);
    const vec r = d - n * (2.0f * d.dot(n));
    const float spec = shadow ? 0.0f : std::pow(std::max(r.dot(light), 0.0f), 40.0f);
    vec col = albedo * (0.10f + diffuse * 0.90f) + vec{1, 1, 1} * (spec * 0.6f);
    if (depth > 0 && mirror > 0) {
        col = col * (1.0f - mirror) + trace(p + n * 1e-3f, r, depth - 1) * mirror;
    }
    const float fog = std::min(t / 70.0f, 1.0f);
    return col * (1.0f - fog) + vec{1.00f, 0.50f, 0.22f} * fog;
}

static std::vector<uint8_t> render() {
    std::vector<uint8_t> rgba(static_cast<size_t>(W * H * 4));
    const vec eye{0.0f, 1.6f, 0.0f};
    for (int32_t y = 0; y < H; y++) {
        for (int32_t x = 0; x < W; x++) {
            vec sum;
            for (int32_t s = 0; s < 4; s++) {  // 2x2 supersampling
                const float u = (static_cast<float>(x) + 0.25f + 0.5f * static_cast<float>(s & 1)) / W - 0.5f;
                const float v = 0.5f - (static_cast<float>(y) + 0.25f + 0.5f * static_cast<float>(s >> 1)) / H;
                sum = sum + trace(eye, vec{u * 1.33f, v - 0.12f, 1.0f}.norm(), 3);
            }
            auto *px = &rgba[static_cast<size_t>((y * W + x) * 4)];
            const auto to8 = [](float l) {
                return static_cast<uint8_t>(std::pow(std::clamp(l * 0.25f, 0.0f, 1.0f), 1.0f / 2.2f) * 255.0f + 0.5f);
            };
            px[0] = to8(sum.x);
            px[1] = to8(sum.y);
            px[2] = to8(sum.z);
            px[3] = 255;
        }
    }
    return rgba;
}

int main(int argc, char **argv) {
    using canvas_type = image<format_32bit, 3 * W + 80, 2 * H + 150>;
    static auto canvas = std::make_unique<canvas_type>();
    canvas->fill_rect(0, 0, canvas->width(), canvas->height(), [](float, float v, int32_t, int32_t) {
        return std::array<float, 4>{0.07f + v * 0.03f, 0.07f + v * 0.03f, 0.09f + v * 0.05f, 1.0f};
    });
    canvas->draw_string<ibmplexsans_semibold_24_aa, true>(20, 14, "One picture, six pixel formats", color::WHITE);

    const std::vector<uint8_t> scene = render();
    const uint8_t label_col = canvas->get_nearest_color(190, 200, 215);

    int32_t slot = 0;
    auto show = [&](const uint8_t *rgba, const char *name, size_t bytes) {
        const int32_t x = 20 + (slot % 3) * (W + 20);
        const int32_t y = 62 + (slot / 3) * (H + 44);
        slot++;
        canvas->blit_RGBA(x, y, W, H, rgba, W, H, W * 4);
        canvas->stroke_rect(x - 1, y - 1, W + 2, H + 2, color::GRAY_40);
        std::array<char, 96> label{};
        std::snprintf(label.data(), label.size(), "%s, %d bytes", name, static_cast<int>(bytes));
        canvas->draw_string<ibmplexsans_regular_14_aa, true>(x, y + H + 6, label.data(), label_col);
    };

    // Dither the scene into a palette image and show what the image holds.
    auto dithered = [&](auto &img, const char *name, bool diffuse = true) {
        if (diffuse) {
            img.blit_RGBA_diffused_linear(0, 0, W, H, scene.data(), W, H, W * 4);
        } else {
            img.blit_RGBA(0, 0, W, H, scene.data(), W, H, W * 4);
        }
        static std::array<uint8_t, W * H * 4> rgba;
        img.RGBA_uint8(rgba);
        for (size_t i = 3; i < rgba.size(); i += 4) {
            rgba[i] = 255;  // palette index 0 is exported as transparent
        }
        show(rgba.data(), name, img.size());
    };

    show(scene.data(), "Ray traced source, 32 bit", scene.size());
    static image<format_8bit, W, H> img8;
    dithered(img8, "8 bit, 256 colors");
    static image<format_4bit, W, H> img4;
    dithered(img4, "4 bit, 16 colors");
    static image<format_2bit, W, H, true> img2;
    dithered(img2, "2 bit, 4 grays");
    static image<format_1bit, W, H> img1;
    dithered(img1, "1 bit, error diffusion");
    static image<format_1bit, W, H> img1t;
    dithered(img1t, "1 bit, threshold", false);

    return example::output(*canvas, argc, argv);
}
