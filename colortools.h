#ifndef _COLORTOOLS_H_
#define _COLORTOOLS_H_

#include <cstdint>
#include <array>

#include "./fastmath.h"

namespace colortools {

struct node {
    int16_t child[8]{-1, -1, -1, -1, -1, -1, -1, -1};
    int16_t palette = -1;
};

template <size_t S, size_t N>
class octree {
    size_t node_idx = 0;
    static constexpr size_t palette_size = S;
    static constexpr size_t node_size = N;
    std::array<node, node_size> nodes;
    const std::array<uint32_t, palette_size> &pal;

    static constexpr int child_index(uint8_t r, uint8_t g, uint8_t b, int level) {
        return ((r >> level) & 1) << 2 | ((g >> level) & 1) << 1 | ((b >> level) & 1);
    }

    constexpr void insert(uint32_t c, size_t idx) {
        uint32_t n = 0;
        for (int32_t lvl = 7; lvl >= 0; --lvl) {
            int32_t ci = child_index((c >> 16) & 0xFF, (c >> 8) & 0xFF, (c >> 0) & 0xFF, lvl);
            if (nodes[n].child[ci] == -1) {
                nodes[n].child[ci] = static_cast<int16_t>(node_idx);
                node_idx++;
            }
            n = nodes[n].child[ci];
        }
        nodes[n].palette = static_cast<int16_t>(idx);
    }

   public:

    constexpr octree(const std::array<uint32_t, palette_size> &palette) : pal(palette) {
        for (size_t i = 0; i < palette_size; ++i) {
            insert(palette[i], i);
        }
    }

    constexpr uint8_t nearest(uint8_t r, uint8_t g, uint8_t b) const {
        int n = 0;
        for (int lvl = 7; lvl >= 0; --lvl) {
            int ci = child_index(r, g, b, lvl);
            int next = nodes[n].child[ci];
            if (next == -1)
                break;
            n = next;
            if (nodes[n].palette != -1)
                return static_cast<uint8_t>(nodes[n].palette);
        }
        // fallback: brute-force
        int best = 0, bestd = INT_MAX;
        for (size_t i = 0; i < pal.size(); ++i) {
            int dr = int(r) - ((pal[i] >> 16) & 0xFF);
            int dg = int(g) - ((pal[i] >> 8) & 0xFF);
            int db = int(b) - ((pal[i] >> 0) & 0xFF);
            int d = dr * dr + dg * dg + db * db;
            if (d < bestd) {
                bestd = d;
                best = int(i);
            }
        }
        return static_cast<uint8_t>(best);
    }

    constexpr size_t nodes_used_length() const { return node_idx; }
    constexpr size_t nodes_memory_length() const { return node_size; }
};

static constexpr double cos(double x, int terms = 10) {
    x = x - 6.283185307179586 * int(x / 6.283185307179586);  // wrap x to [0, 2π)
    double res = 1.0, term = 1.0;
    double x2 = x * x;
    for (int i = 1; i < terms; ++i) {
        term *= -x2 / ((2 * i - 1) * (2 * i));
        res += term;
    }
    return res;
}

static consteval double sin(double x, int terms = 10) {
    x = x - 6.283185307179586 * int(x / 6.283185307179586);  // wrap x to [0, 2π)
    double res = x, term = x;
    double x2 = x * x;
    for (int i = 1; i < terms; ++i) {
        term *= -x2 / ((2 * i) * (2 * i + 1));
        res += term;
    }
    return res;
}

static consteval double pow(double base, double exp, int terms = 10) {
    if (base <= 0.0)
        return (base == 0.0) ? 0.0 : 0.0 / 0.0;  // NaN for negative base
    double ln = 0.0, y = (base - 1) / (base + 1);
    double y2 = y * y, num = y;
    for (int i = 1; i <= terms; ++i) {
        ln += num / (2 * i - 1);
        num *= y2;
    }
    ln *= 2;
    double res = 1.0, term = 1.0, x = exp * ln;
    for (int i = 1; i < terms; ++i) {
        term *= x / i;
        res += term;
    }
    return res;
}

struct oklch {
    double l, c, h;
};

struct oklab {
    double l, a, b;
};

struct srgb {
    double r, g, b;
};

static consteval double linear_to_srgb(double c) {
    if (c <= 0.0031308) {
        return 12.92 * c;
    } else {
        return 1.055 * pow(c, 1.0 / 2.4) - 0.055;
    }
}

static consteval double srgb_to_linear(double s) {
    if (s <= 0.040449936) {
        return s / 12.92;
    } else {
        return pow((s + 0.055) / 1.055, 2.4);
    }
}

static constexpr float linear_to_srgb(float c) {
    if (c <= 0.0031308f) {
        return 12.92f * c;
    } else {
        return 1.055f * fast_pow(c, 1.0f / 2.4f) - 0.055f;
    }
}

static constexpr float srgb_to_linear(float s) {
    if (s <= 0.040449936f) {
        return s / 12.92f;
    } else {
        return fast_pow((s + 0.055f) / 1.055f, 2.4f);
    }
}

static consteval srgb oklab_to_srgb(const oklab &oklab) {
    double l = oklab.l;
    double a = oklab.a;
    double b = oklab.b;

    double l_ = l + 0.3963377774 * a + 0.2158037573 * b;
    double m_ = l - 0.1055613458 * a - 0.0638541728 * b;
    double s_ = l - 0.0894841775 * a - 1.2914855480 * b;

    double r = 4.0767416621 * l_ - 3.3077115913 * m_ + 0.2309699292 * s_;
    double g = -1.2684380046 * l_ + 2.6097574011 * m_ - 0.3413193965 * s_;
    double bl = -0.0041960863 * l_ - 0.7034186168 * m_ + 1.7076147031 * s_;

    return {linear_to_srgb(std::max(0.0, std::min(1.0, r))), linear_to_srgb(std::max(0.0, std::min(1.0, g))), linear_to_srgb(std::max(0.0, std::min(1.0, bl)))};
}

static consteval oklab oklch_to_oklab(const oklch &oklch) {
    return {oklch.l, oklch.c * cos(oklch.h * M_PI / 180.0), oklch.c * sin(oklch.h * M_PI / 180.0)};
}

};  // namespace colortools

#endif  // #ifndef _COLORTOOLS_H_
