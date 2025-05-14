/*
MIT License

Copyright (c) 2025 Tinic Uro

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#ifndef CONSTIXEL_H_
#define CONSTIXEL_H_

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

#if defined(__ARM_NEON)
#include <arm_neon.h>
#endif  // #if defined(__ARM_NEON)

#if defined(__AVX2__)
#include <immintrin.h>
#endif  // #if defined(__AVX2__)

namespace constixel {

/// @cond DOXYGEN_EXCLUDE
[[nodiscard]] static constexpr float fast_exp2(const float p) {
    const float offset = (p < 0) ? 1.0f : 0.0f;
    const float clipp = (p < -126) ? -126.0f : p;
    const float z = clipp - static_cast<float>(static_cast<int32_t>(clipp)) + offset;
    return std::bit_cast<float>(static_cast<uint32_t>((1 << 23) * (clipp + 121.2740575f + 27.7280233f / (4.84252568f - z) - 1.49012907f * z)));
}

[[nodiscard]] static constexpr float fast_log2(const float x) {
    uint32_t xi = std::bit_cast<uint32_t>(x);
    float xf = std::bit_cast<float>((xi & 0x007FFFFF) | 0x3f000000);
    const float y = static_cast<float>(xi) * 1.1920928955078125e-7f;
    return y - 124.22551499f - 1.498030302f * xf - 1.72587999f / (0.3520887068f + xf);
}

static constexpr float fast_sqrtf(const float x) {
    int32_t i = std::bit_cast<int>(x);
    int k = i & 0x00800000;
    float y;
    if (k != 0) {
        i = 0x5ed9d098 - (i >> 1);
        y = std::bit_cast<float>(i);
        y = 2.33139729f * y * ((-x * y * y) + 1.07492042f);
    } else {
        i = 0x5f19d352 - (i >> 1);
        y = std::bit_cast<float>(i);
        y = 0.82420468f * y * ((-x * y * y) + 2.14996147f);
    }
    float c = x * y;
    float r = ((y * -c) + 1.0f);
    y = ((0.5f * c * r) + c);
    return y;
}

[[nodiscard]] static constexpr float fast_pow(const float x, const float p) {
    return fast_exp2(p * fast_log2(x));
}

static constexpr double m_pi_d = 3.14159265358979323846;

[[nodiscard]] static consteval double cos(double x, int32_t terms = 10) {
    x = x - 6.283185307179586 * static_cast<int32_t>(x / 6.283185307179586);  // wrap x to [0, 2π)
    double res = 1.0, term = 1.0;
    double x2 = x * x;
    for (int32_t i = 1; i < terms; ++i) {
        term *= -x2 / ((2 * i - 1) * (2 * i));
        res += term;
    }
    return res;
}

[[nodiscard]] static consteval double sin(double x, int32_t terms = 10) {
    x = x - 6.283185307179586 * static_cast<int32_t>(x / 6.283185307179586);  // wrap x to [0, 2π)
    double res = x, term = x;
    double x2 = x * x;
    for (int32_t i = 1; i < terms; ++i) {
        term *= -x2 / ((2 * i) * (2 * i + 1));
        res += term;
    }
    return res;
}

[[nodiscard]] static consteval double pow(double base, double exp, int32_t terms = 10) {
    if (base <= 0.0)
        return (base == 0.0) ? 0.0 : 0.0 / 0.0;  // NaN for negative base
    double ln = 0.0, y = (base - 1) / (base + 1);
    double y2 = y * y, num = y;
    for (int32_t i = 1; i <= terms; ++i) {
        ln += num / (2 * i - 1);
        num *= y2;
    }
    ln *= 2;
    double res = 1.0, term = 1.0, x = exp * ln;
    for (int32_t i = 1; i < terms; ++i) {
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

[[nodiscard]] static consteval double linear_to_srgb(double c) {
    if (c <= 0.0031308) {
        return 12.92 * c;
    } else {
        return 1.055 * pow(c, 1.0 / 2.4) - 0.055;
    }
}

[[nodiscard]] static consteval double srgb_to_linear(double s) {
    if (s <= 0.040449936) {
        return s / 12.92;
    } else {
        return pow((s + 0.055) / 1.055, 2.4);
    }
}

[[nodiscard]] static constexpr float linear_to_srgb(float c) {
    if (c <= 0.0031308f) {
        return 12.92f * c;
    } else {
        if (std::is_constant_evaluated()) {
            return 1.055f * fast_pow(c, 1.0f / 2.4f) - 0.055f;
        } else {
            return 1.055f * std::pow(c, 1.0f / 2.4f) - 0.055f;
        }
    }
}

[[nodiscard]] static constexpr float srgb_to_linear(float s) {
    if (s <= 0.040449936f) {
        return s / 12.92f;
    } else {
        if (std::is_constant_evaluated()) {
            return fast_pow((s + 0.055f) / 1.055f, 2.4f);
        } else {
            return std::pow((s + 0.055f) / 1.055f, 2.4f);
        }
    }
}

[[nodiscard]] static consteval srgb oklab_to_srgb(const oklab &oklab) {
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

[[nodiscard]] static consteval oklab oklch_to_oklab(const oklch &oklch) {
    return {oklch.l, oklch.c * cos(oklch.h * m_pi_d / 180.0), oklch.c * sin(oklch.h * m_pi_d / 180.0)};
}

static constexpr float epsilon_low = static_cast<float>(srgb_to_linear(0.5 / 255.0));
static constexpr float epsilon_high = static_cast<float>(srgb_to_linear(254.5 / 255.0));
/// @endcond

/// @cond PRIVATE_CLASS
template <size_t S>
class quantize {
    static constexpr size_t palette_size = S;

#if defined(__ARM_NEON)
    std::array<float, palette_size * 3> linearpal_neon{};
#endif  // #if defined(__ARM_NEON)
#if defined(__AVX2__)
    alignas(32) std::array<float, palette_size * 3> linearpal_avx2{};
#endif  // #if defined(__ARM_NEON)

    const std::array<uint32_t, palette_size> &pal;

 public:
    explicit constexpr quantize(const std::array<uint32_t, palette_size> &palette) : pal(palette) {
        for (size_t i = 0; i < pal.size(); i++) {
            linearpal.at(i * 3 + 0) = srgb_to_linear(static_cast<float>((pal[i] >> 16) & 0xFF) * (1.0f / 255.0f));
            linearpal.at(i * 3 + 1) = srgb_to_linear(static_cast<float>((pal[i] >> 8) & 0xFF) * (1.0f / 255.0f));
            linearpal.at(i * 3 + 2) = srgb_to_linear(static_cast<float>((pal[i] >> 0) & 0xFF) * (1.0f / 255.0f));
        }
#if defined(__ARM_NEON)
        if (pal.size() >= 4) {
            for (size_t i = 0; i < pal.size(); i += 4) {
                for (size_t j = 0; j < 4; j++) {
                    linearpal_neon.at(i * 3 + j) = linearpal.at((i + j) * 3 + 0);
                    linearpal_neon.at(i * 3 + j + 4) = linearpal.at((i + j) * 3 + 1);
                    linearpal_neon.at(i * 3 + j + 8) = linearpal.at((i + j) * 3 + 2);
                }
            }
        }
#endif  // #if defined(__ARM_NEON)
#if defined(__AVX2__)
        if (pal.size() >= 8) {
            for (size_t i = 0; i < pal.size(); i += 8) {
                for (size_t j = 0; j < 8; j++) {
                    linearpal_avx2.at(i * 3 + j) = linearpal.at((i + j) * 3 + 0);
                    linearpal_avx2.at(i * 3 + j + 8) = linearpal.at((i + j) * 3 + 1);
                    linearpal_avx2.at(i * 3 + j + 16) = linearpal.at((i + j) * 3 + 2);
                }
            }
        }
#endif  // #if defined(__AVX2__)
    }

    [[nodiscard]] constexpr uint8_t nearest(int32_t r, int32_t g, int32_t b) const {
        int32_t best = 0, bestd = 1UL << 30;
        for (size_t i = 0; i < pal.size(); ++i) {
            int32_t dr = r - static_cast<int32_t>((pal[i] >> 16) & 0xFF);
            int32_t dg = g - static_cast<int32_t>((pal[i] >> 8) & 0xFF);
            int32_t db = b - static_cast<int32_t>((pal[i] >> 0) & 0xFF);
            int32_t d = dr * dr + dg * dg + db * db;
            if (d < bestd) {
                bestd = d;
                best = static_cast<int32_t>(i);
            }
        }
        return static_cast<uint8_t>(best);
    }

    std::array<float, palette_size * 3> linearpal{};

    [[nodiscard]] constexpr uint8_t nearest_linear(float r, float g, float b) const {
#if defined(__ARM_NEON)
        if (!std::is_constant_evaluated() && pal.size() >= 4) {
            const float32x4_t vR = vdupq_n_f32(r);
            const float32x4_t vG = vdupq_n_f32(g);
            const float32x4_t vB = vdupq_n_f32(b);

            float best = std::numeric_limits<float>::infinity();
            std::size_t bestIdx = 0;

            for (size_t i = 0; i < pal.size(); i += 4) {
                float32x4_t dr = vsubq_f32(vld1q_f32(&linearpal_neon[i * 3]), vR);
                float32x4_t dg = vsubq_f32(vld1q_f32(&linearpal_neon[i * 3 + 4]), vG);
                float32x4_t db = vsubq_f32(vld1q_f32(&linearpal_neon[i * 3 + 8]), vB);

#if defined(__aarch64__) && defined(__ARM_FEATURE_FMA)
                float32x4_t dist = vfmaq_f32(vfmaq_f32(vmulq_f32(dr, dr), dg, dg), db, db);
#else
                float32x4_t dist = vaddq_f32(vaddq_f32(vmulq_f32(dr, dr), vmulq_f32(dg, dg)), vmulq_f32(db, db));
#endif

                float32x2_t lo = vget_low_f32(dist), hi = vget_high_f32(dist);
                float d0 = vget_lane_f32(lo, 0);
                if (d0 < best) {
                    best = d0;
                    bestIdx = i;
                }
                float d1 = vget_lane_f32(lo, 1);
                if (d1 < best) {
                    best = d1;
                    bestIdx = i + 1;
                }
                float d2 = vget_lane_f32(hi, 0);
                if (d2 < best) {
                    best = d2;
                    bestIdx = i + 2;
                }
                float d3 = vget_lane_f32(hi, 1);
                if (d3 < best) {
                    best = d3;
                    bestIdx = i + 3;
                }
            }
            return static_cast<uint8_t>(bestIdx);
        }
#endif  // #if defined(__ARM_NEON)
#if defined(__AVX2__)
        if (!std::is_constant_evaluated() && pal.size() >= 8) {
            const __m256 vR = _mm256_set1_ps(r);
            const __m256 vG = _mm256_set1_ps(g);
            const __m256 vB = _mm256_set1_ps(b);

            float best = std::numeric_limits<float>::infinity();
            std::uint8_t bestIdx = 0;

            for (size_t i = 0; i < pal.size(); i += 8) {
                __m256 pr = _mm256_load_ps(&linearpal_avx2[i * 3]);
                __m256 pg = _mm256_load_ps(&linearpal_avx2[i * 3 + 8]);
                __m256 pb = _mm256_load_ps(&linearpal_avx2[i * 3 + 16]);

                __m256 dist =
                    _mm256_fmadd_ps(_mm256_sub_ps(pr, vR), _mm256_sub_ps(pr, vR),
                                    _mm256_fmadd_ps(_mm256_sub_ps(pg, vG), _mm256_sub_ps(pg, vG), _mm256_mul_ps(_mm256_sub_ps(pb, vB), _mm256_sub_ps(pb, vB))));

                alignas(32) float d[8];
                _mm256_store_ps(d, dist);

                for (size_t lane = 0; lane < 8; lane++) {
                    if (d[lane] < best) {
                        best = d[lane];
                        bestIdx = static_cast<uint8_t>(i + lane);
                    }
                }
            }

            return static_cast<uint8_t>(bestIdx);
        }
#endif  // #if defined(__AVX2__)
        size_t best = 0;
        float bestd = 100.0f;
        for (size_t i = 0; i < pal.size(); ++i) {
            float dr = r - linearpal.at(i * 3 + 0);
            float dg = g - linearpal.at(i * 3 + 1);
            float db = b - linearpal.at(i * 3 + 2);
            float d = dr * dr + dg * dg + db * db;
            if (d < bestd) {
                bestd = d;
                best = i;
            }
        }
        return static_cast<uint8_t>(best);
    }
};
/// @endcond

/// @cond PRIVATE_CLASS
template <size_t N, typename T>
class hextree {
    static constexpr T bitslices = ((sizeof(T) * 8) / 4) - 1;

    static constexpr size_t child_nodes_n = 1UL << 4;
    static constexpr size_t child_nodes_n_mask = child_nodes_n - 1;

    struct node {
        T child[child_nodes_n]{};
        constexpr node() {
            for (auto &c : child) {
                c = invalid;
            }
        }
    };

 public:
    static constexpr T invalid = std::numeric_limits<T>::max();
    std::array<node, N> nodes{};

    [[nodiscard]] consteval size_t byte_size() const {
        return sizeof(node) * nodes.size();
    }

    hextree(const hextree &) = delete;
    hextree &operator=(const hextree &) = delete;

    explicit consteval hextree() {
    }

    template <std::size_t NS>
    explicit consteval hextree(const std::array<std::pair<T, T>, NS> &in) {
        nodes[0] = node{};
        T node_cnt = 1;
        for (auto [key, val] : in) {
            T idx = 0;
            for (T d = 0; d < bitslices; d++) {
                T nib = (key >> ((bitslices - d) * 4)) & child_nodes_n_mask;
                T next = nodes.at(idx).child[nib];
                if (next == invalid) {
                    next = node_cnt;
                    nodes.at(idx).child[nib] = next;
                    node_cnt++;
                }
                idx = next;
            }
            nodes.at(idx).child[key & child_nodes_n_mask] = val;
        }
    }

    template <std::size_t NS>
    [[nodiscard]] static consteval T size(const std::array<std::pair<T, T>, NS> &in) {
        std::vector<node> vnodes{1};
        vnodes.assign(1, node{});
        for (auto [key, val] : in) {
            T idx = 0;
            for (T d = 0; d < bitslices; d++) {
                T nib = (key >> ((bitslices - d) * 4)) & child_nodes_n_mask;
                T next = vnodes.at(idx).child[nib];
                if (next == invalid) {
                    next = static_cast<T>(vnodes.size());
                    vnodes.at(idx).child[nib] = next;
                    vnodes.emplace_back();
                }
                idx = next;
            }
            vnodes.at(idx).child[key & child_nodes_n_mask] = val;
        }
        return static_cast<T>(vnodes.size());
    }

    [[nodiscard]] constexpr T lookup(T key) const {
        T idx = 0;
        for (T d = 0; d < bitslices; ++d) {
            T nib = (key >> ((bitslices - d) * 4)) & child_nodes_n_mask;
            T next = nodes.at(idx).child[nib];
            if (next == invalid) {
                return invalid;
            }
            idx = next;
        }
        return nodes.at(idx).child[key & child_nodes_n_mask];
    }
};
/// @endcond

/// @cond PRIVATE_CLASS
template <typename T>
struct char_info {
    T x;
    T y;
    T width;
    T height;
    T xadvance;
    T xoffset;
    T yoffset;
};
/// @endcond

/// @brief basic rectangle structure
/// @tparam T coordinate number type
template <typename T>
struct rect {
    T x = 0;  ///< x coordinate
    T y = 0;  ///< y coordinate
    T w = 0;  ///< width
    T h = 0;  ///< height

    /// @brief intersects one rect with another
    /// @param other the other rectangle
    constexpr rect operator&(const rect &other) const {
        T nax = std::max(x, other.x);
        T nay = std::max(y, other.y);
        T nix = std::min(x + w, other.x + other.w);
        T niy = std::min(x + h, other.y + other.h);
        return {nax, nay, nix - nax, niy - nay};
    }

    /// @brief intersects one rect with another
    /// @param other the other rectangle
    constexpr rect &operator&=(const rect &other) {
        T nax = std::max(x, other.x);
        T nay = std::max(y, other.y);
        T nix = std::min(x + w, other.x + other.w);
        T niy = std::min(x + h, other.y + other.h);
        x = nax;
        y = nay;
        w = nix - nax;
        h = niy - nay;
        return *this;
    }
};

/// @cond PRIVATE_CLASS
class format {
 public:
    [[nodiscard]] static constexpr uint32_t adler32(const uint8_t *data, std::size_t len, uint32_t adler32_sum) {
        uint32_t adler32_s1 = adler32_sum & 0xFFFF;
        uint32_t adler32_s2 = adler32_sum >> 16;
        for (size_t c = 0; c < len; c++) {
            adler32_s1 = (adler32_s1 + data[c]) % 65521;
            adler32_s2 = (adler32_s2 + adler32_s1) % 65521;
        }
        return ((adler32_s2 << 16) | adler32_s1);
    }

    //@private
    template <typename F>
    static constexpr void png_write_be(F &&char_out, uint32_t value) {
        char_out(static_cast<char>((value >> 24) & 0xFF));
        char_out(static_cast<char>((value >> 16) & 0xFF));
        char_out(static_cast<char>((value >> 8) & 0xFF));
        char_out(static_cast<char>((value >> 0) & 0xFF));
    }

    template <typename F, typename A>
    static constexpr void png_write_crc32(F &&char_out, const A &array, size_t bytes) {
        size_t idx = 0;
        uint32_t crc = 0xFFFFFFFF;
        size_t len = bytes;
        while (len--) {
            uint8_t d = static_cast<uint8_t>(array.at(idx++));
            crc ^= d;
            for (int i = 0; i < 8; ++i) {
                crc = (crc & 1) ? (crc >> 1) ^ 0xEDB88320u : crc >> 1;
            }
        }
        png_write_be(char_out, crc ^ 0xFFFFFFFF);
    }

    template <typename F, typename A>
    static constexpr void png_write_array(F &&char_out, const A &array, size_t bytes) {
        for (size_t c = 0; c < bytes; c++) {
            char_out(static_cast<char>(array.at(c)));
        }
    }

    template <typename F>
    static constexpr void png_marker(F &&char_out) {
        char_out(static_cast<char>(0x89));
        char_out(0x50);
        char_out(0x4E);
        char_out(0x47);
        char_out(0x0D);
        char_out(0x0A);
        char_out(0x1A);
        char_out(0x0A);
    }

    template <typename F>
    static constexpr void png_header(F &&char_out, size_t w, size_t h, size_t depth) {
        const size_t chunkLength = 17;
        std::array<char, chunkLength> header;
        uint32_t i = 0;
        header.at(i++) = 'I';
        header.at(i++) = 'H';
        header.at(i++) = 'D';
        header.at(i++) = 'R';
        header.at(i++) = static_cast<char>((w >> 24) & 0xFF);
        header.at(i++) = static_cast<char>((w >> 16) & 0xFF);
        header.at(i++) = static_cast<char>((w >> 8) & 0xFF);
        header.at(i++) = static_cast<char>((w >> 0) & 0xFF);
        header.at(i++) = static_cast<char>((h >> 24) & 0xFF);
        header.at(i++) = static_cast<char>((h >> 16) & 0xFF);
        header.at(i++) = static_cast<char>((h >> 8) & 0xFF);
        header.at(i++) = static_cast<char>((h >> 0) & 0xFF);
        header.at(i++) = static_cast<char>(depth);
        header.at(i++) = 3;
        header.at(i++) = 0;
        header.at(i++) = 0;
        header.at(i++) = 0;
        png_write_be(char_out, i - 4);
        png_write_array(char_out, header, i);
        png_write_crc32(char_out, header, i);
    }

    template <typename F, typename P>
    static constexpr void png_palette(F &&char_out, const P &palette) {
        std::array<char, 256 * 3 + 4> header;
        uint32_t i = 0;
        header.at(i++) = 'P';
        header.at(i++) = 'L';
        header.at(i++) = 'T';
        header.at(i++) = 'E';
        for (size_t c = 0; c < palette.size(); c++) {
            header.at(i++) = static_cast<char>((palette[c] >> 16) & 0xFF);
            header.at(i++) = static_cast<char>((palette[c] >> 8) & 0xFF);
            header.at(i++) = static_cast<char>((palette[c] >> 0) & 0xFF);
        }
        png_write_be(char_out, i - 4);
        png_write_array(char_out, header, i);
        png_write_crc32(char_out, header, i);
    }

    template <typename F>
    static constexpr void png_end(F &&char_out) {
        std::array<char, 4> header;
        uint32_t i = 0;
        header.at(i++) = 'I';
        header.at(i++) = 'E';
        header.at(i++) = 'N';
        header.at(i++) = 'D';
        png_write_be(char_out, i - 4);
        png_write_array(char_out, header, i);
        png_write_crc32(char_out, header, i);
    }

    template <typename F>
    static constexpr void png_idat_zlib_header(F &&char_out) {
        std::array<char, 6> header;
        uint32_t i = 0;
        header.at(i++) = 'I';
        header.at(i++) = 'D';
        header.at(i++) = 'A';
        header.at(i++) = 'T';
        header.at(i++) = 0x78;
        header.at(i++) = 0x01;
        png_write_be(char_out, i - 4);
        png_write_array(char_out, header, i);
        png_write_crc32(char_out, header, i);
    }

    template <typename F>
    [[nodiscard]] static constexpr uint32_t png_idat_zlib_stream(F &&char_out, const uint8_t *line, size_t bytes, uint32_t adler32_sum) {
        const size_t max_data_use = 1024;
        const size_t extra_data = 24;
        const size_t max_stack_use = max_data_use + extra_data;
        std::array<uint8_t, max_stack_use> header;
        while (bytes > 0) {
            uint32_t i = 0;
            header.at(i++) = 'I';
            header.at(i++) = 'D';
            header.at(i++) = 'A';
            header.at(i++) = 'T';
            header.at(i++) = 0x00;

            size_t bytes_to_copy = std::min(static_cast<size_t>(max_data_use), bytes);
            header.at(i++) = (((bytes_to_copy + 1) >> 0) & 0xFF);
            header.at(i++) = (((bytes_to_copy + 1) >> 8) & 0xFF);
            header.at(i++) = ((((bytes_to_copy + 1) ^ 0xffff) >> 0) & 0xFF);
            header.at(i++) = ((((bytes_to_copy + 1) ^ 0xffff) >> 8) & 0xFF);

            uint32_t adlersum32_start_pos = i;
            header.at(i++) = 0;
            for (size_t c = 0; c < bytes_to_copy; c++) {
                header.at(i++) = line[c];
            }
            adler32_sum = adler32(&header[adlersum32_start_pos], i - adlersum32_start_pos, adler32_sum);

            png_write_be(char_out, i - 4);
            png_write_array(char_out, header, i);
            png_write_crc32(char_out, header, i);

            bytes -= bytes_to_copy;
        }
        return adler32_sum;
    }

    template <typename F>
    static constexpr void png_idat_zlib_trailer(F &&char_out, uint32_t adler32_sum) {
        std::array<char, 8> header;
        uint32_t i = 0;
        header.at(i++) = 'I';
        header.at(i++) = 'D';
        header.at(i++) = 'A';
        header.at(i++) = 'T';
        header.at(i++) = static_cast<char>((adler32_sum >> 24) & 0xFF);
        header.at(i++) = static_cast<char>((adler32_sum >> 16) & 0xFF);
        header.at(i++) = static_cast<char>((adler32_sum >> 8) & 0xFF);
        header.at(i++) = static_cast<char>((adler32_sum >> 0) & 0xFF);
        png_write_be(char_out, i - 4);
        png_write_array(char_out, header, i);
        png_write_crc32(char_out, header, i);
    }

    template <typename F>
    static constexpr void sixel_header(F &&char_out) {
        char_out(0x1b);
        char_out('P');
        char_out('0');
        char_out(';');
        char_out('1');
        char_out(';');
        char_out('0');
        char_out('q');
    }

    template <size_t W, size_t H, size_t S, typename F>
    static constexpr void sixel_raster_attributes(F &&char_out) {
        char_out('\"');
        sixel_number(char_out, 2);
        char_out(';');
        sixel_number(char_out, 2);
        char_out(';');
        sixel_number(char_out, W * S);
        char_out(';');
        sixel_number(char_out, H * S);
    }

    template <typename F>
    static constexpr void sixel_number(F &&char_out, uint16_t u) {
        if (u < 10) {
            char_out(static_cast<char>('0' + u));
        } else if (u < 100) {
            char_out(static_cast<char>('0' + (((u / 10) % 10))));
            char_out(static_cast<char>('0' + (u % 10)));
        } else if (u < 1000) {
            char_out(static_cast<char>('0' + ((u / 100) % 10)));
            char_out(static_cast<char>('0' + ((u / 10) % 10)));
            char_out(static_cast<char>('0' + (u % 10)));
        } else if (u < 10000) {
            char_out(static_cast<char>('0' + ((u / 1000) % 10)));
            char_out(static_cast<char>('0' + ((u / 100) % 10)));
            char_out(static_cast<char>('0' + ((u / 10) % 10)));
            char_out(static_cast<char>('0' + (u % 10)));
        } else {
            char_out(static_cast<char>('0' + ((u / 10000) % 10)));
            char_out(static_cast<char>('0' + ((u / 1000) % 10)));
            char_out(static_cast<char>('0' + ((u / 100) % 10)));
            char_out(static_cast<char>('0' + ((u / 10) % 10)));
            char_out(static_cast<char>('0' + (u % 10)));
        }
    }

    template <typename F>
    static constexpr void sixel_color(F &&char_out, uint16_t i, uint32_t col) {
        char_out('#');
        sixel_number(char_out, i);
        char_out(';');
        char_out('2');
        char_out(';');
        for (size_t c = 0; c < 3; c++) {
            sixel_number(char_out, static_cast<uint16_t>((((col >> (8 * (2 - c))) & 0xFF) * 100) / 255));
            char_out(';');
        }
    }

    template <typename F>
    static constexpr void sixel_end(F &&char_out) {
        char_out(0x1b);
        char_out('\\');
    }

    template <typename PBT, size_t PBS>
    struct palette_bitset {
        constexpr void mark(PBT col) {
            PBT idx = (col >> 5) & set_idx_mask;
            set[idx] |= 1UL << (col & 0x1F);
        }

        constexpr void clear() {
            set.fill(0);
        }

        [[nodiscard]] constexpr size_t genstack(std::array<PBT, PBS> &stack) const {
            size_t count = 0;
            for (size_t c = 0; c < PBS / 32; c++) {
                for (size_t d = 0; d < 32; d++) {
                    if ((set[c] >> d) & 1) {
                        stack.at(count++) = static_cast<PBT>(c * 32 + d);
                    }
                }
            }
            return count;
        }

        std::array<uint32_t, PBS / 32> set{};
        static constexpr size_t set_idx_mask = (1UL << sizeof(set) / 4) - 1;
    };

    template <size_t W, size_t H, int32_t S, typename PBT, size_t PBS, typename P, typename F, typename L>
    static constexpr void png_image(const uint8_t *data, const P &palette, F &&char_out, const L &line_ptr) {
        png_marker(char_out);
        png_header(char_out, W, H, PBS);
        png_palette(char_out, palette);
        png_idat_zlib_header(char_out);
        uint32_t adler32_sum = 1;
        for (size_t y = 0; y < H; y++) {
            size_t bpl = 0;
            const uint8_t *ptr = line_ptr(data, y, bpl);
            adler32_sum = png_idat_zlib_stream(char_out, ptr, bpl, adler32_sum);
        }
        png_idat_zlib_trailer(char_out, adler32_sum);
        png_end(char_out);
    }

    template <size_t W, size_t H, int32_t S, typename PBT, size_t PBS, typename P, typename F, typename C, typename D>
    static constexpr void sixel_image(const uint8_t *data, const P &palette, F &&char_out, const rect<int32_t> &_r, const C &collect6, const D &set6) {
        sixel_header(char_out);
        sixel_raster_attributes<W, H, S>(char_out);
        for (size_t c = 0; c < palette.size(); c++) {
            sixel_color(char_out, static_cast<uint16_t>(c), palette[c]);
        }
        const auto r = rect<int32_t>{_r.x * S, _r.y * S, _r.w * S, _r.h * S} & rect<int32_t>{0, 0, W * S, H * S};
        std::array<PBT, std::max(32UL, 1UL << PBS)> stack{};
        palette_bitset<PBT, std::max(32UL, 1UL << PBS)> pset{};
        for (size_t y = static_cast<size_t>(r.y); y < static_cast<size_t>(r.y + r.h); y += 6) {
            pset.clear();
            set6(data, static_cast<size_t>(r.x), static_cast<size_t>(r.w), y, pset);
            size_t stack_count = pset.genstack(stack);
            for (size_t s = 0; s < stack_count; s++) {
                PBT col = stack[s];
                if (col != 0) {
                    char_out('$');
                }
                char_out('#');
                sixel_number(char_out, static_cast<uint16_t>(col));
                for (size_t x = static_cast<size_t>(r.x); x < static_cast<size_t>(r.x + r.w); x++) {
                    PBT bits6 = collect6(data, x, col, y);
                    uint16_t repeat_count = 0;
                    for (size_t xr = (x + 1); xr < (std::min(x + size_t{255}, W * size_t{S})); xr++) {
                        if (bits6 == collect6(data, xr, col, y)) {
                            repeat_count++;
                            continue;
                        }
                        break;
                    }
                    if (repeat_count > 3) {
                        char_out('!');
                        sixel_number(char_out, repeat_count + 1);
                        x += repeat_count;
                    }
                    char_out(static_cast<char>('?' + bits6));
                }
            }
            char_out('-');
        }
        sixel_end(char_out);
    }
};
/// @endcond

/// @brief 1-bit format, just b/w. Use as template parameter for image.
/// @tparam W Width in pixels.
/// @tparam H Height in pixels.
/// @tparam S scale of sixel output.
/// @tparam GR Grayscale palette.
template <size_t W, size_t H, int32_t S, bool GR>
class format_1bit : public format {
 public:
    /// @cond DOXYGEN_EXCLUDE
    static constexpr bool grayscale = GR;
    static constexpr size_t bits_per_pixel = 1;
    static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
    static constexpr size_t image_size = H * bytes_per_line;
    static constexpr std::array<uint32_t, (1UL << bits_per_pixel)> palette = {0x00000000, 0x00ffffff};

    static constexpr uint8_t reverse(uint8_t b) {
        b = static_cast<uint8_t>((b & uint8_t{0xF0}) >> 4 | (b & uint8_t{0x0F}) << 4);
        b = static_cast<uint8_t>((b & uint8_t{0xCC}) >> 2 | (b & uint8_t{0x33}) << 2);
        b = static_cast<uint8_t>((b & uint8_t{0xAA}) >> 1 | (b & uint8_t{0x55}) << 1);
        return b;
    }

    static constexpr void transpose8x8(std::array<uint8_t, 8> &a) {
        // clang-format off
        uint32_t x = (static_cast<uint32_t>(a[0]) << 24) | 
                     (static_cast<uint32_t>(a[1]) << 16) | 
                     (static_cast<uint32_t>(a[2]) <<  8) | 
                     (static_cast<uint32_t>(a[3])      );
        uint32_t y = (static_cast<uint32_t>(a[4]) << 24) | 
                     (static_cast<uint32_t>(a[5]) << 16) | 
                     (static_cast<uint32_t>(a[6]) <<  8) | 
                     (static_cast<uint32_t>(a[7])      );

        uint32_t t = (x ^ (x >>  7)) & 0x00AA00AA; x = x ^ t ^ (t <<  7);
                 t = (y ^ (y >>  7)) & 0x00AA00AA; y = y ^ t ^ (t <<  7);
                 t = (x ^ (x >> 14)) & 0x0000CCCC; x = x ^ t ^ (t << 14);
                 t = (y ^ (y >> 14)) & 0x0000CCCC; y = y ^ t ^ (t << 14);

        t = ((y >> 4) & 0x0F0F0F0F) | (x & 0xF0F0F0F0);
        y = ((x << 4) & 0xF0F0F0F0) | (y & 0x0F0F0F0F);

        a[0] = static_cast<uint8_t>(t >> 24);
        a[1] = static_cast<uint8_t>(t >> 16);
        a[2] = static_cast<uint8_t>(t >>  8);
        a[3] = static_cast<uint8_t>(t      );
        a[4] = static_cast<uint8_t>(y >> 24);
        a[5] = static_cast<uint8_t>(y >> 16);
        a[6] = static_cast<uint8_t>(y >>  8);
        a[7] = static_cast<uint8_t>(y      );
        // clang-format on
    }

    static constexpr void transpose(const uint8_t *src, uint8_t *dst) {
        std::array<uint8_t, 8> tmp;
        size_t src_stride = ((W + 7) / 8);
        size_t dst_stride = ((H + 7) / 8);
        for (size_t y = 0; y < dst_stride; y++) {
            for (size_t x = 0; x < src_stride; x++) {
                size_t xl = std::min(size_t{8}, H - (y * 8));
                for (size_t c = 0; c < xl; c++) {
                    tmp[c] = src[(y * 8 + c) * src_stride + x];
                }
                for (size_t c = xl; c < 8; c++) {
                    tmp[c] = 0;
                }
                transpose8x8(tmp);
                size_t yl = std::min(size_t{8}, W - (x * 8));
                for (size_t c = 0; c < yl; c++) {
                    dst[(x * 8 + c) * dst_stride + y] = tmp[c];
                }
            }
        }
    }

    static constexpr void compose(std::array<uint8_t, image_size> &, size_t, size_t, float, float, float, float) {
#ifndef _MSC_VER
        static_assert(false, "composing not supported on 1-bit format, use a mono font.");
#endif  // #ifndef _MSC_VER
    }

    static constexpr void plot(std::array<uint8_t, image_size> &data, size_t x0, size_t y, uint8_t col) {
        col &= (1UL << bits_per_pixel) - 1;
        size_t x8 = x0 / 8;
        x0 %= 8;
        uint8_t *yptr = &data.data()[y * bytes_per_line];
        yptr[x8] &= ~static_cast<uint8_t>(1UL << (7 - x0));
        yptr[x8] |= static_cast<uint8_t>(col << (7 - x0));
    }

    static constexpr void span(std::array<uint8_t, image_size> &data, size_t xl0, size_t xr0, size_t y, uint8_t col) {
        col &= (1UL << bits_per_pixel) - 1;
        size_t xl8 = xl0 / 8;
        xl0 %= 8;
        size_t xr8 = xr0 / 8;
        xr0 %= 8;
        size_t xs8 = xr8 - xl8;
        uint8_t c8 = static_cast<uint8_t>(col << 7 | col << 6 | col << 5 | col << 4 | col << 3 | col << 2 | col << 1 | col << 0);
        constexpr uint8_t ml[] = {0b11111111, 0b01111111, 0b00111111, 0b00011111, 0b00001111, 0b00000111, 0b00000011, 0b00000001};
        constexpr uint8_t mr[] = {0b00000000, 0b10000000, 0b11000000, 0b11100000, 0b11110000, 0b11111000, 0b11111100, 0b11111110};
        uint8_t *yptr = &data.data()[y * bytes_per_line];
        if (xs8 > 0) {
            yptr[xl8] &= ~ml[xl0];
            yptr[xl8] |= ml[xl0] & c8;
            for (size_t x = xl8 + 1; x < xr8; x++) {
                yptr[x] = c8;
            }
            yptr[xr8] &= ~mr[xr0];
            yptr[xr8] |= mr[xr0] & c8;
        } else {
            yptr[xl8] &= ~(ml[xl0] & mr[xr0]);
            yptr[xl8] |= (ml[xl0] & mr[xr0] & c8);
        }
    }

    [[nodiscard]] static constexpr uint8_t get_col(const uint8_t *line, size_t x) {
        size_t x8 = x / 8;
        size_t xb = x % 8;
        return (line[x8] >> (7 - xb)) & 1;
    }

    [[nodiscard]] static constexpr auto RGBA_uint32(const std::array<uint8_t, image_size> &data) {
        std::array<uint32_t, W * H> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = get_col(ptr, x);
                rgba[y * W + x] = palette.at(col) | (col ? 0xFF000000 : 0x00000000);
            }
            ptr += bytes_per_line;
        }
        return rgba;
    }

    [[nodiscard]] static constexpr auto RGBA_uint8(const std::array<uint8_t, image_size> &data) {
        std::array<uint8_t, W * H * 4> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = get_col(ptr, x);
                rgba[y * W * 4 + x * 4 + 0] = static_cast<uint8_t>((palette.at(col) >> 16) & 0xFF);
                rgba[y * W * 4 + x * 4 + 1] = static_cast<uint8_t>((palette.at(col) >> 8) & 0xFF);
                rgba[y * W * 4 + x * 4 + 2] = static_cast<uint8_t>((palette.at(col) >> 0) & 0xFF);
                rgba[y * W * 4 + x * 4 + 3] = (col ? 0xFF : 0x00);
            }
            ptr += bytes_per_line;
        }
        return rgba;
    }

    static constexpr void blit_RGBA(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                uint32_t R = ptr[y * static_cast<size_t>(stride) + x * 4 + 0];
                uint32_t G = ptr[y * static_cast<size_t>(stride) + x * 4 + 1];
                uint32_t B = ptr[y * static_cast<size_t>(stride) + x * 4 + 2];
                plot(data, x + static_cast<size_t>(r.x), y + static_cast<size_t>(r.y), (R * 2 + G * 3 + B * 1) > 768 ? 1 : 0);
            }
        }
    }

    static constexpr void blit_RGBA_diffused(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        int32_t err = 0;
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                int32_t R = ptr[y * static_cast<size_t>(stride) + x * 4 + 0];
                int32_t G = ptr[y * static_cast<size_t>(stride) + x * 4 + 1];
                int32_t B = ptr[y * static_cast<size_t>(stride) + x * 4 + 2];
                int32_t V = (R * 2 + G * 3 + B * 1) + err;
                uint8_t n = V > 768 ? 1 : 0;
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                err = std::clamp(V - (n ? int32_t{0xFF * 6} : int32_t{0x00}), int32_t{-0xFF * 6}, int32_t{0xFF * 6});
            }
        }
    }

    static constexpr void blit_RGBA_diffused_linear(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            float err_r = 0;
            float err_g = 0;
            float err_b = 0;
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                float R = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 0]) * (1.0f / 255.0f);
                float G = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 1]) * (1.0f / 255.0f);
                float B = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 2]) * (1.0f / 255.0f);
                float Rl = constixel::srgb_to_linear(R);
                float Gl = constixel::srgb_to_linear(G);
                float Bl = constixel::srgb_to_linear(B);
                Rl = Rl + err_r;
                Gl = Gl + err_g;
                Bl = Bl + err_b;
                uint8_t n = (Rl * 2 * +Gl * 3 + Bl * 1) > 3.0f ? 1 : 0;
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                const float c = 0.75;
                err_r = std::clamp(Rl - (n ? 1.0f : 0.0f), -c, c);
                err_g = std::clamp(Gl - (n ? 1.0f : 0.0f), -c, c);
                err_b = std::clamp(Bl - (n ? 1.0f : 0.0f), -c, c);
            }
        }
    }

    template <typename F>
    static constexpr void png(const std::array<uint8_t, image_size> &data, F &&char_out) {
        png_image<W, H, S, uint8_t, bits_per_pixel>(data.data(), palette, char_out, [](const uint8_t *data_raw, size_t y, size_t &bpl) {
            bpl = bytes_per_line;
            return data_raw + y * bytes_per_line;
        });
    }

    template <typename F>
    static constexpr void sixel(const std::array<uint8_t, image_size> &data, F &&char_out, const rect<int32_t> &r) {
        sixel_image<W, H, S, uint8_t, bits_per_pixel>(
            data.data(), palette, char_out, r,
            [](const uint8_t *data_raw, size_t x, size_t col, size_t y) {
                const uint8_t *ptr = &data_raw[(y / S) * bytes_per_line + (x / S) / 8];
                size_t x8 = (x / S) % 8;
                uint8_t out = 0;
                size_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    out >>= 1;
                    if ((y + y6) < H * S) {
                        out |= ((((*ptr) >> (7 - x8)) & 1) == col) ? (1UL << 5) : 0;
                        if (y6 != 5) {
                            if (++inc >= S) {
                                inc = 0;
                                ptr += bytes_per_line;
                            }
                        }
                    }
                }
                return out;
            },
            [](const uint8_t *data_raw, size_t x, size_t w, size_t y, palette_bitset<uint8_t, 32> &set) {
                size_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    if ((y + y6) < H * S) {
                        for (size_t xx = 0; xx < (w + S - 1) / S; xx++) {
                            const uint8_t *ptr = &data_raw[((y + y6) / S) * bytes_per_line + (xx + x) / 8];
                            size_t x8 = (xx + x) % 8;
                            set.mark((((*ptr) >> (7 - x8)) & 1));
                        }
                        if (++inc >= S) {
                            inc = 0;
                        }
                    }
                }
            });
    }
    /// @endcond
};

/// @brief 2-bit color format, 4 colors total. Use as template parameter for image.
/// @tparam W Width in pixels.
/// @tparam H Height in pixels.
/// @tparam S scale of sixel output.
/// @tparam GR Grayscale palette.
template <size_t W, size_t H, int32_t S, bool GR>
class format_2bit : public format {
 public:
    /// @cond DOXYGEN_EXCLUDE
    static constexpr bool grayscale = GR;
    static constexpr size_t bits_per_pixel = 2;
    static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
    static constexpr size_t image_size = H * bytes_per_line;
    static consteval const std::array<uint32_t, (1UL << bits_per_pixel)> gen_palette() {
        if (GR) {
            return {0x000000, 0x444444, 0x888888, 0xffffff};
        } else {
            return {0x000000, 0xffffff, 0xff0000, 0x0077ff};
        }
    }
    static constexpr std::array<uint32_t, (1UL << bits_per_pixel)> palette = gen_palette();

    static constexpr const constixel::quantize<1UL << bits_per_pixel> gen_quant() {
        return constixel::quantize<1UL << bits_per_pixel>(palette);
    }

    static constexpr const constixel::quantize<1UL << bits_per_pixel> quant = gen_quant();

    static constexpr uint8_t reverse(uint8_t b) {
        b = static_cast<uint8_t>((b & uint8_t{0xF0}) >> 4 | (b & uint8_t{0x0F}) << 4);
        b = static_cast<uint8_t>((b & uint8_t{0xCC}) >> 2 | (b & uint8_t{0x33}) << 2);
        return b;
    }

    static constexpr void transpose(const uint8_t *, uint8_t *) {
#ifndef _MSC_VER
        static_assert(false, "Not implemented yet.");
#endif  // #ifndef _MSC_VER
    }

    static constexpr void compose(std::array<uint8_t, image_size> &, size_t, size_t, float, float, float, float) {
#ifndef _MSC_VER
        static_assert(false, "composing not supported on 2-bit format, use a mono font.");
#endif  // #ifndef _MSC_VER
    }

    static constexpr void plot(std::array<uint8_t, image_size> &data, size_t x0, size_t y, uint8_t col) {
        col &= (1UL << bits_per_pixel) - 1;
        size_t x4 = x0 / 4;
        x0 %= 4;
        uint8_t *yptr = &data.data()[y * bytes_per_line];
        yptr[x4] &= ~static_cast<uint8_t>(3UL << (6 - x0 * 2));
        yptr[x4] |= static_cast<uint8_t>(col << (6 - x0 * 2));
    }

    static constexpr void span(std::array<uint8_t, image_size> &data, size_t xl0, size_t xr0, size_t y, uint8_t col) {
        col &= (1UL << bits_per_pixel) - 1;
        size_t xl4 = xl0 / 4;
        xl0 %= 4;
        size_t xr4 = xr0 / 4;
        xr0 %= 4;
        size_t xs4 = xr4 - xl4;
        uint8_t c4 = static_cast<uint8_t>(col << 6 | col << 4 | col << 2 | col << 0);
        constexpr uint8_t ml[] = {0b11111111, 0b00111111, 0b00001111, 0b00000011};
        constexpr uint8_t mr[] = {0b00000000, 0b11000000, 0b11110000, 0b11111100};
        uint8_t *yptr = &data.data()[y * bytes_per_line];
        if (xs4 > 0) {
            yptr[xl4] &= ~ml[xl0];
            yptr[xl4] |= ml[xl0] & c4;
            for (size_t x = xl4 + 1; x < xr4; x++) {
                yptr[x] = c4;
            }
            yptr[xr4] &= ~mr[xr0];
            yptr[xr4] |= mr[xr0] & c4;
        } else {
            yptr[xl4] &= ~(ml[xl0] & mr[xr0]);
            yptr[xl4] |= (ml[xl0] & mr[xr0] & c4);
        }
    }

    static constexpr uint8_t get_col(const uint8_t *line, size_t x) {
        size_t x4 = x / 4;
        size_t xb = x % 4;
        return (line[x4] >> ((3 - xb) * 2)) & 0x3;
    }

    [[nodiscard]] static constexpr auto RGBA_uint32(const std::array<uint8_t, image_size> &data) {
        std::array<uint32_t, W * H> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = get_col(ptr, x);
                rgba[y * W + x] = palette.at(col) | (col ? 0xFF000000 : 0x00000000);
            }
            ptr += bytes_per_line;
        }
        return rgba;
    }

    [[nodiscard]] static constexpr auto RGBA_uint8(const std::array<uint8_t, image_size> &data) {
        std::array<uint8_t, W * H * 4> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = get_col(ptr, x);
                rgba[y * W * 4 + x * 4 + 0] = static_cast<uint8_t>((palette.at(col) >> 16) & 0xFF);
                rgba[y * W * 4 + x * 4 + 1] = static_cast<uint8_t>((palette.at(col) >> 8) & 0xFF);
                rgba[y * W * 4 + x * 4 + 2] = static_cast<uint8_t>((palette.at(col) >> 0) & 0xFF);
                rgba[y * W * 4 + x * 4 + 3] = (col ? 0xFF : 0x00);
            }
            ptr += bytes_per_line;
        }
        return rgba;
    }

    static constexpr void blit_RGBA(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                plot(data, (x + static_cast<uint32_t>(r.x)), (y + static_cast<uint32_t>(r.y)),
                     quant.nearest(ptr[y * static_cast<size_t>(stride) + x * 4 + 0], ptr[y * static_cast<size_t>(stride) + x * 4 + 1],
                                   ptr[y * static_cast<size_t>(stride) + x * 4 + 2]));
            }
        }
    }

    static constexpr void blit_RGBA_diffused(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            int32_t err_r = 0;
            int32_t err_g = 0;
            int32_t err_b = 0;
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                int32_t R = ptr[y * static_cast<size_t>(stride) + x * 4 + 0];
                int32_t G = ptr[y * static_cast<size_t>(stride) + x * 4 + 1];
                int32_t B = ptr[y * static_cast<size_t>(stride) + x * 4 + 2];
                R = R + err_r;
                G = G + err_g;
                B = B + err_b;
                uint8_t n = quant.nearest(R, G, B);
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                err_r = std::clamp(R - static_cast<int32_t>((palette.at(n) >> 16) & 0xFF), int32_t{-255}, int32_t{255});
                err_g = std::clamp(G - static_cast<int32_t>((palette.at(n) >> 8) & 0xFF), int32_t{-255}, int32_t{255});
                err_b = std::clamp(B - static_cast<int32_t>((palette.at(n) >> 0) & 0xFF), int32_t{-255}, int32_t{255});
            }
        }
    }

    static constexpr void blit_RGBA_diffused_linear(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            float err_r = 0;
            float err_g = 0;
            float err_b = 0;
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                float R = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 0]) * (1.0f / 255.0f);
                float G = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 1]) * (1.0f / 255.0f);
                float B = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 2]) * (1.0f / 255.0f);
                float Rl = constixel::srgb_to_linear(R);
                float Gl = constixel::srgb_to_linear(G);
                float Bl = constixel::srgb_to_linear(B);
                Rl = Rl + err_r;
                Gl = Gl + err_g;
                Bl = Bl + err_b;
                uint8_t n = quant.nearest_linear(Rl, Gl, Bl);
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                err_r = std::clamp(Rl - quant.linearpal.at(n * 3 + 0), -1.0f, 1.0f);
                err_g = std::clamp(Gl - quant.linearpal.at(n * 3 + 1), -1.0f, 1.0f);
                err_b = std::clamp(Bl - quant.linearpal.at(n * 3 + 2), -1.0f, 1.0f);
            }
        }
    }

    template <typename F>
    static constexpr void png(const std::array<uint8_t, image_size> &data, F &&char_out) {
        png_image<W, H, S, uint8_t, bits_per_pixel>(data.data(), palette, char_out, [](const uint8_t *data_raw, size_t y, size_t &bpl) {
            bpl = bytes_per_line;
            return data_raw + y * bytes_per_line;
        });
    }

    template <typename F>
    static constexpr void sixel(const std::array<uint8_t, image_size> &data, F &&char_out, const rect<int32_t> &r) {
        sixel_image<W, H, S, uint8_t, bits_per_pixel>(
            data.data(), palette, char_out, r,
            [](const uint8_t *data_raw, size_t x, size_t col, size_t y) {
                const uint8_t *ptr = &data_raw[(y / S) * bytes_per_line + (x / S) / 4];
                size_t x4 = (x / S) % 4;
                uint8_t out = 0;
                size_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    out >>= 1;
                    if ((y + y6) < H * S) {
                        out |= ((((*ptr) >> (6 - x4 * 2)) & 3) == col) ? (1UL << 5) : 0;
                        if (y6 != 5) {
                            if (++inc >= S) {
                                inc = 0;
                                ptr += bytes_per_line;
                            }
                        }
                    }
                }
                return out;
            },
            [](const uint8_t *data_raw, size_t x, size_t w, size_t y, palette_bitset<uint8_t, 32> &set) {
                size_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    if ((y + y6) < H * S) {
                        for (size_t xx = 0; xx < (w + S - 1) / S; xx++) {
                            const uint8_t *ptr = &data_raw[((y + y6) / S) * bytes_per_line + (xx + x) / 4];
                            size_t x4 = (xx + x) % 4;
                            set.mark((((*ptr) >> (6 - x4 * 2)) & 3));
                        }
                        if (++inc >= S) {
                            inc = 0;
                        }
                    }
                }
            });
    }
    /// @endcond
};

/// @brief 4-bit color format, 16 colors total. Use as template parameter for image.
/// @tparam W Width in pixels.
/// @tparam H Height in pixels.
/// @tparam S scale of sixel output.
/// @tparam GR Grayscale palette.
template <size_t W, size_t H, int32_t S, bool GR>
class format_4bit : public format {
 public:
    /// @cond DOXYGEN_EXCLUDE
    static constexpr bool grayscale = GR;
    static constexpr size_t bits_per_pixel = 4;
    static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
    static constexpr size_t image_size = H * bytes_per_line;

    static consteval const std::array<uint32_t, (1UL << bits_per_pixel)> gen_palette() {
        if (GR) {
            return {0x000000, 0x111111, 0x222222, 0x333333, 0x444444, 0x555555, 0x666666, 0x777777,
                    0x888888, 0x999999, 0xaaaaaa, 0xbbbbbb, 0xcccccc, 0xdddddd, 0xeeeeee, 0xffffff};
        } else {
            return {0x000000, 0xffffff, 0xff0000, 0x00ff00, 0x0000ff, 0xffff00, 0x00ffff, 0xff00ff,
                    0x333333, 0x666666, 0x999999, 0xcccccc, 0x7f0000, 0x007f00, 0x00007f, 0x7f7f00};
        }
    }

    static constexpr std::array<uint32_t, (1UL << bits_per_pixel)> palette = gen_palette();

    static constexpr const constixel::quantize<1UL << bits_per_pixel> gen_quant() {
        return constixel::quantize<1UL << bits_per_pixel>(palette);
    }

    static constexpr const constixel::quantize<1UL << bits_per_pixel> quant = gen_quant();

    static constexpr uint8_t reverse(uint8_t b) {
        b = static_cast<uint8_t>((b & uint8_t{0xF0}) >> 4 | (b & uint8_t{0x0F}) << 4);
        return b;
    }

    static constexpr void transpose(const uint8_t *, uint8_t *) {
#ifndef _MSC_VER
        static_assert(false, "Not implemented yet.");
#endif  // #ifndef _MSC_VER
    }

    static constexpr void compose(std::array<uint8_t, image_size> &data, size_t x, size_t y, float cola, float colr, float colg, float colb) {
        size_t bg = static_cast<size_t>(get_col(data, x, y));
        float Rl = colr + quant.linearpal.at(bg * 3 + 0) * (1.0f - cola);
        float Gl = colg + quant.linearpal.at(bg * 3 + 1) * (1.0f - cola);
        float Bl = colb + quant.linearpal.at(bg * 3 + 2) * (1.0f - cola);
        plot(data, x, y, quant.nearest_linear(Rl, Gl, Bl));
    }

    static constexpr void plot(std::array<uint8_t, image_size> &data, size_t x0, size_t y, uint8_t col) {
        col &= (1UL << bits_per_pixel) - 1;
        size_t x2 = x0 / 2;
        x0 %= 2;
        uint8_t *yptr = &data.data()[y * bytes_per_line];
        yptr[x2] &= ~static_cast<uint8_t>(0xFUL << (4 - x0 * 4));
        yptr[x2] |= static_cast<uint8_t>(col << (4 - x0 * 4));
    }

    static constexpr void span(std::array<uint8_t, image_size> &data, size_t xl0, size_t xr0, size_t y, uint8_t col) {
        col &= (1UL << bits_per_pixel) - 1;
        size_t xl2 = xl0 / 2;
        xl0 %= 2;
        size_t xr2 = xr0 / 2;
        xr0 %= 2;
        size_t xs2 = xr2 - xl2;
        uint8_t c2 = static_cast<uint8_t>(col << 4 | col << 0);
        constexpr uint8_t ml[] = {0b11111111, 0b00001111};
        constexpr uint8_t mr[] = {0b00000000, 0b11110000};
        uint8_t *yptr = &data.data()[y * bytes_per_line];
        if (xs2 > 0) {
            yptr[xl2] &= ~ml[xl0];
            yptr[xl2] |= ml[xl0] & c2;
            for (size_t x = xl2 + 1; x < xr2; x++) {
                yptr[x] = c2;
            }
            yptr[xr2] &= ~mr[xr0];
            yptr[xr2] |= mr[xr0] & c2;
        } else {
            yptr[xl2] &= ~(ml[xl0] & mr[xr0]);
            yptr[xl2] |= (ml[xl0] & mr[xr0] & c2);
        }
    }

    [[nodiscard]] static constexpr uint8_t get_col(const uint8_t *line, size_t x) {
        size_t x2 = x / 2;
        size_t xb = x % 2;
        return (line[x2] >> ((1 - xb) * 4)) & 0xF;
    }

    [[nodiscard]] static constexpr uint8_t get_col(const std::array<uint8_t, image_size> &data, size_t x, size_t y) {
        const uint8_t *ptr = data.data() + y * bytes_per_line;
        size_t x2 = x / 2;
        size_t xb = x % 2;
        return (ptr[x2] >> ((1 - xb) * 4)) & 0xF;
    }

    [[nodiscard]] static constexpr auto RGBA_uint32(const std::array<uint8_t, image_size> &data) {
        std::array<uint32_t, W * H> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = get_col(ptr, x);
                rgba[y * W + x] = palette.at(col) | (col ? 0xFF000000 : 0x00000000);
            }
            ptr += bytes_per_line;
        }
        return rgba;
    }

    [[nodiscard]] static constexpr auto RGBA_uint8(const std::array<uint8_t, image_size> &data) {
        std::array<uint8_t, W * H * 4> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = get_col(ptr, x);
                rgba[y * W * 4 + x * 4 + 0] = static_cast<uint8_t>((palette.at(col) >> 16) & 0xFF);
                rgba[y * W * 4 + x * 4 + 1] = static_cast<uint8_t>((palette.at(col) >> 8) & 0xFF);
                rgba[y * W * 4 + x * 4 + 2] = static_cast<uint8_t>((palette.at(col) >> 0) & 0xFF);
                rgba[y * W * 4 + x * 4 + 3] = (col ? 0xFF : 0x00);
            }
            ptr += bytes_per_line;
        }
        return rgba;
    }

    static constexpr void blit_RGBA(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)),
                     quant.nearest(ptr[y * static_cast<size_t>(stride) + x * 4 + 0], ptr[y * static_cast<size_t>(stride) + x * 4 + 1],
                                   ptr[y * static_cast<size_t>(stride) + x * 4 + 2]));
            }
        }
    }

    static constexpr void blit_RGBA_diffused(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            int32_t err_r = 0;
            int32_t err_g = 0;
            int32_t err_b = 0;
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                int32_t R = ptr[y * static_cast<size_t>(stride) + x * 4 + 0];
                int32_t G = ptr[y * static_cast<size_t>(stride) + x * 4 + 1];
                int32_t B = ptr[y * static_cast<size_t>(stride) + x * 4 + 2];
                R = R + err_r;
                G = G + err_g;
                B = B + err_b;
                uint8_t n = quant.nearest(R, G, B);
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                err_r = std::clamp(R - static_cast<int32_t>((palette.at(n) >> 16) & 0xFF), int32_t{-255}, int32_t{255});
                err_g = std::clamp(G - static_cast<int32_t>((palette.at(n) >> 8) & 0xFF), int32_t{-255}, int32_t{255});
                err_b = std::clamp(B - static_cast<int32_t>((palette.at(n) >> 0) & 0xFF), int32_t{-255}, int32_t{255});
            }
        }
    }

    static constexpr void blit_RGBA_diffused_linear(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            float err_r = 0;
            float err_g = 0;
            float err_b = 0;
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                float R = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 0]) * (1.0f / 255.0f);
                float G = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 1]) * (1.0f / 255.0f);
                float B = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 2]) * (1.0f / 255.0f);
                float Rl = constixel::srgb_to_linear(R);
                float Gl = constixel::srgb_to_linear(G);
                float Bl = constixel::srgb_to_linear(B);
                Rl = Rl + err_r;
                Gl = Gl + err_g;
                Bl = Bl + err_b;
                uint8_t n = quant.nearest_linear(Rl, Gl, Bl);
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                err_r = std::clamp(Rl - quant.linearpal.at(n * 3 + 0), -1.0f, 1.0f);
                err_g = std::clamp(Gl - quant.linearpal.at(n * 3 + 1), -1.0f, 1.0f);
                err_b = std::clamp(Bl - quant.linearpal.at(n * 3 + 2), -1.0f, 1.0f);
            }
        }
    }

    template <typename F>
    static constexpr void png(const std::array<uint8_t, image_size> &data, F &&char_out) {
        png_image<W, H, S, uint8_t, bits_per_pixel>(data.data(), palette, char_out, [](const uint8_t *data_raw, size_t y, size_t &bpl) {
            bpl = bytes_per_line;
            return data_raw + y * bytes_per_line;
        });
    }

    template <typename F>
    static constexpr void sixel(const std::array<uint8_t, image_size> &data, F &&char_out, const rect<int32_t> &r) {
        sixel_image<W, H, S, uint8_t, bits_per_pixel>(
            data.data(), palette, char_out, r,
            [](const uint8_t *data_raw, size_t x, size_t col, size_t y) {
                const uint8_t *ptr = &data_raw[(y / S) * bytes_per_line + (x / S) / 2];
                size_t x2 = (x / S) % 2;
                uint8_t out = 0;
                size_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    out >>= 1;
                    if ((y + y6) < H * S) {
                        out |= ((((*ptr) >> (4 - x2 * 4)) & 0xF) == col) ? (1UL << 5) : 0;
                        if (y6 != 5) {
                            if (++inc >= S) {
                                inc = 0;
                                ptr += bytes_per_line;
                            }
                        }
                    }
                }
                return out;
            },
            [](const uint8_t *data_raw, size_t x, size_t w, size_t y, palette_bitset<uint8_t, 32> &set) {
                size_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    if ((y + y6) < H * S) {
                        for (size_t xx = 0; xx < (w + S - 1) / S; xx++) {
                            const uint8_t *ptr = &data_raw[((y + y6) / S) * bytes_per_line + (xx + x) / 2];
                            size_t x2 = (xx + x) % 2;
                            set.mark((((*ptr) >> (4 - x2 * 4)) & 0xF));
                        }
                        if (++inc >= S) {
                            inc = 0;
                        }
                    }
                }
            });
    }
    /// @endcond
};

/// @brief 8-bit format, 256 colors total. Use as template parameter for image.
/// @tparam W Width in pixels.
/// @tparam H Height in pixels.
/// @tparam S scale of sixel output.
/// @tparam GR Grayscale palette.
template <size_t W, size_t H, int32_t S, bool GR>
class format_8bit : public format {
 public:
    /// @cond DOXYGEN_EXCLUDE
    static constexpr bool grayscale = GR;
    static constexpr size_t bits_per_pixel = 8;
    static constexpr size_t bytes_per_line = W;
    static constexpr size_t image_size = H * bytes_per_line;

    static consteval const std::array<uint32_t, (1UL << bits_per_pixel)> gen_palette() {
        std::array<uint32_t, (1UL << bits_per_pixel)> pal{};
        if (GR) {
            for (size_t c = 0; c < 256; c++) {
                pal[c] = static_cast<uint32_t>((c << 16) | (c << 8) | c);
            }
        } else {
            pal[0] = 0x000000;
            pal[1] = 0xffffff;
            pal[2] = 0xff0000;
            pal[3] = 0x00ff00;
            pal[4] = 0x0000ff;
            pal[5] = 0xffff00;
            pal[6] = 0x00ffff;
            pal[7] = 0xff00ff;
            pal[8] = 0x333333;
            pal[9] = 0x666666;
            pal[10] = 0x999999;
            pal[11] = 0xcccccc;
            pal[12] = 0x7f0000;
            pal[13] = 0x007f00;
            pal[14] = 0x00007f;
            pal[15] = 0x7f7f00;
            for (size_t c = 0; c < 16; c++) {
                uint32_t y = (0xff * static_cast<uint32_t>(c)) / 15;
                pal[0x10 + c] = (y << 16) | (y << 8) | (y << 0);
            }
            for (size_t c = 0; c < 8; c++) {
                uint32_t y = (0xff * static_cast<uint32_t>(c)) / 7;
                uint32_t x = (0xff * (static_cast<uint32_t>(c) + 1)) / 8;
                pal[0x20 + c + 0] = (y << 16) | (0 << 8) | (0 << 0);
                pal[0x20 + c + 8] = (255 << 16) | (x << 8) | (x << 0);
                pal[0x30 + c + 0] = (0 << 16) | (y << 8) | (0 << 0);
                pal[0x30 + c + 8] = (x << 16) | (255 << 8) | (x << 0);
                pal[0x40 + c + 0] = (0 << 16) | (0 << 8) | (y << 0);
                pal[0x40 + c + 8] = (x << 16) | (x << 8) | (255 << 0);
                pal[0x50 + c + 0] = (y << 16) | (y << 8) | (0 << 0);
                pal[0x50 + c + 8] = (255 << 16) | (255 << 8) | (x << 0);
                pal[0x60 + c + 0] = (0 << 16) | (y << 8) | (y << 0);
                pal[0x60 + c + 8] = (x << 16) | (255 << 8) | (255 << 0);
                pal[0x70 + c + 0] = (y << 16) | (0 << 8) | (y << 0);
                pal[0x70 + c + 8] = (255 << 16) | (x << 8) | (255 << 0);
            }
            for (size_t c = 0; c < 8; c++) {
                constixel::oklab lft{static_cast<double>(c) / 7 - 0.2, 0.2, 0.0};
                constixel::oklab rgh{static_cast<double>(c) / 7 - 0.2, 0.2, 337.5};
                for (size_t d = 0; d < 16; d++) {
                    auto res = constixel::oklab_to_srgb(constixel::oklch_to_oklab(constixel::oklch{
                        std::lerp(lft.l, rgh.l, static_cast<double>(d) / 15.0),
                        std::lerp(lft.a, rgh.a, static_cast<double>(d) / 15.0),
                        std::lerp(lft.b, rgh.b, static_cast<double>(d) / 15.0),
                    }));
                    pal[0x80 + c * 16 + d] = (static_cast<uint32_t>(std::max(0.0, std::min(1.0, res.r)) * 255.0) << 16) |
                                             (static_cast<uint32_t>(std::max(0.0, std::min(1.0, res.g)) * 255.0) << 8) |
                                             (static_cast<uint32_t>(std::max(0.0, std::min(1.0, res.b)) * 255.0) << 0);
                }
            }
        }
        return pal;
    }

    static constexpr std::array<uint32_t, 1UL << bits_per_pixel> palette = gen_palette();

    static constexpr const constixel::quantize<1UL << bits_per_pixel> gen_quant() {
        return constixel::quantize<1UL << bits_per_pixel>(palette);
    }

    static constexpr const constixel::quantize<1UL << bits_per_pixel> quant = gen_quant();

    static constexpr uint8_t reverse(uint8_t b) {
        return b;
    }

    static constexpr void transpose(const uint8_t *src, uint8_t *dst) {
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                dst[x * H + y] = *src++;
            }
        }
    }

    static constexpr void plot(std::array<uint8_t, image_size> &data, size_t x, size_t y, uint8_t col) {
        data.data()[y * bytes_per_line + x] = static_cast<uint8_t>(col);
    }

    static constexpr void span(std::array<uint8_t, image_size> &data, size_t xl0, size_t xr0, size_t y, uint8_t col) {
        uint8_t *yptr = &data.data()[y * bytes_per_line];
        for (size_t x = xl0; x < xr0; x++) {
            yptr[x] = static_cast<uint8_t>(col);
        }
    }

    static constexpr void compose(std::array<uint8_t, image_size> &data, size_t x, size_t y, float cola, float colr, float colg, float colb) {
        size_t bg = static_cast<size_t>(data.data()[y * bytes_per_line + x]);
        float Rl = colr + quant.linearpal.at(bg * 3 + 0) * (1.0f - cola);
        float Gl = colg + quant.linearpal.at(bg * 3 + 1) * (1.0f - cola);
        float Bl = colb + quant.linearpal.at(bg * 3 + 2) * (1.0f - cola);
        plot(data, x, y, quant.nearest_linear(Rl, Gl, Bl));
    }

    [[nodiscard]] static constexpr auto RGBA_uint32(const std::array<uint8_t, image_size> &data) {
        std::array<uint32_t, W * H> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = ptr[x];
                rgba[y * W + x] = palette.at(col) | (col ? 0xFF000000 : 0x00000000);
            }
            ptr += bytes_per_line;
        }
        return rgba;
    }

    [[nodiscard]] static constexpr auto RGBA_uint8(const std::array<uint8_t, image_size> &data) {
        std::array<uint8_t, W * H * 4> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = ptr[x];
                rgba[y * W * 4 + x * 4 + 0] = static_cast<uint8_t>((palette.at(col) >> 16) & 0xFF);
                rgba[y * W * 4 + x * 4 + 1] = static_cast<uint8_t>((palette.at(col) >> 8) & 0xFF);
                rgba[y * W * 4 + x * 4 + 2] = static_cast<uint8_t>((palette.at(col) >> 0) & 0xFF);
                rgba[y * W * 4 + x * 4 + 3] = (col ? 0xFF : 0x00);
            }
            ptr += bytes_per_line;
        }
        return rgba;
    }

    static constexpr void blit_RGBA(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                data.data()[(y + static_cast<size_t>(r.y)) * bytes_per_line + (x + static_cast<size_t>(r.x))] =
                    quant.nearest(ptr[y * static_cast<size_t>(stride) + x * 4 + 0], ptr[y * static_cast<size_t>(stride) + x * 4 + 1],
                                  ptr[y * static_cast<size_t>(stride) + x * 4 + 2]);
            }
        }
    }

    static constexpr void blit_RGBA_diffused(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            int32_t err_r = 0;
            int32_t err_g = 0;
            int32_t err_b = 0;
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                int32_t R = ptr[y * static_cast<size_t>(stride) + x * 4 + 0];
                int32_t G = ptr[y * static_cast<size_t>(stride) + x * 4 + 1];
                int32_t B = ptr[y * static_cast<size_t>(stride) + x * 4 + 2];
                R = R + err_r;
                G = G + err_g;
                B = B + err_b;
                uint8_t n = quant.nearest(R, G, B);
                data.data()[(y + static_cast<size_t>(r.y)) * bytes_per_line + (x + static_cast<size_t>(r.x))] = n;
                err_r = std::clamp(R - static_cast<int32_t>((palette.at(n) >> 16) & 0xFF), int32_t{-255}, int32_t{255});
                err_g = std::clamp(G - static_cast<int32_t>((palette.at(n) >> 8) & 0xFF), int32_t{-255}, int32_t{255});
                err_b = std::clamp(B - static_cast<int32_t>((palette.at(n) >> 0) & 0xFF), int32_t{-255}, int32_t{255});
            }
        }
    }

    static constexpr void blit_RGBA_diffused_linear(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            float err_r = 0;
            float err_g = 0;
            float err_b = 0;
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                float R = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 0]) * (1.0f / 255.0f);
                float G = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 1]) * (1.0f / 255.0f);
                float B = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 2]) * (1.0f / 255.0f);
                float Rl = constixel::srgb_to_linear(R);
                float Gl = constixel::srgb_to_linear(G);
                float Bl = constixel::srgb_to_linear(B);
                Rl = Rl + err_r;
                Gl = Gl + err_g;
                Bl = Bl + err_b;
                uint8_t n = quant.nearest_linear(Rl, Gl, Bl);
                data.data()[(y + static_cast<size_t>(r.y)) * bytes_per_line + (x + static_cast<size_t>(r.x))] = n;
                err_r = std::clamp(Rl - quant.linearpal.at(n * 3 + 0), -1.0f, 1.0f);
                err_g = std::clamp(Gl - quant.linearpal.at(n * 3 + 1), -1.0f, 1.0f);
                err_b = std::clamp(Bl - quant.linearpal.at(n * 3 + 2), -1.0f, 1.0f);
            }
        }
    }

    template <typename F>
    static constexpr void png(const std::array<uint8_t, image_size> &data, F &&char_out) {
        png_image<W, H, S, uint8_t, bits_per_pixel>(data.data(), palette, char_out, [](const uint8_t *data_raw, size_t y, size_t &bpl) {
            bpl = bytes_per_line;
            return data_raw + y * bytes_per_line;
        });
    }

    template <typename F>
    static constexpr void sixel(const std::array<uint8_t, image_size> &data, F &&char_out, const rect<int32_t> &r) {
        sixel_image<W, H, S, uint8_t, bits_per_pixel>(
            data.data(), palette, char_out, r,
            [](const uint8_t *data_raw, size_t x, size_t col, size_t y) {
                const uint8_t *ptr = &data_raw[(y / S) * bytes_per_line + x / S];
                uint8_t out = 0;
                size_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    out >>= 1;
                    if ((y + y6) < H * S) {
                        out |= (*ptr == col) ? (1UL << 5) : 0;
                        if (y6 != 5) {
                            if (++inc >= S) {
                                inc = 0;
                                ptr += bytes_per_line;
                            }
                        }
                    }
                }
                return out;
            },
            [](const uint8_t *data_raw, size_t x, size_t w, size_t y, palette_bitset<uint8_t, 1UL << bits_per_pixel> &set) {
                const uint8_t *ptr = &data_raw[(y / S) * bytes_per_line + x / S];
                size_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    if ((y + y6) < H * S) {
                        for (size_t xx = 0; xx < (w + S - 1) / S; xx++) {
                            set.mark(ptr[xx + x]);
                        }
                        if (y6 != 5) {
                            if (++inc >= S) {
                                inc = 0;
                                ptr += bytes_per_line;
                            }
                        }
                    }
                }
            });
    }
    /// @endcond
};

/*! Enum class which contains a list of predefined color name and values to use with the drawing API. */
enum color : uint8_t {
    BLACK = 0,         /*!< Black */
    TRANSPARENT = 0,   /*!< Transparent */
    BLACK_OPAQUE = 16, /*!< Black, only valid for format_8bit */
    WHITE = 1,         /*!< White */
    RED = 2,           /*!< Red */
    GREEN = 3,         /*!< Green */
    BLUE = 4,          /*!< Blue */
    YELLOW = 5,        /*!< Yellow */
    CYAN = 6,          /*!< Cyan */
    MAGENTA = 7,       /*!< Magenta */
    GRAY_80 = 8,       /*!< Gray 80% */
    GRAY_60 = 9,       /*!< Gray 60% */
    GRAY_40 = 10,      /*!< Gray 40% */
    GRAY_20 = 11,      /*!< Gray 20% */
    DARK_RED = 12,     /*!< Dark Red (50%) */
    DARK_GREEN = 13,   /*!< Dark Green (50%) */
    DARK_BLUE = 14,    /*!< Dark Blue (50%) */
    DARK_YELLOW = 15,  /*!< Dark Yellow (50%) */

    GRAY_RAMP_START = 16,
    GRAY_RAMP_COUNT = 16,
    GRAY_RAMP_STOP = GRAY_RAMP_START + 15,

    RED_LUMA_RAMP_START = 32,
    RED_LUMA_RAMP_COUNT = 16,
    RED_LUMA_RAMP_STOP = RED_LUMA_RAMP_START + 15,

    GREEN_LUMA_RAMP_START = 48,
    GREEN_LUMA_RAMP_COUNT = 16,
    GREEN_LUMA_RAMP_STOP = GREEN_LUMA_RAMP_START + 15,

    BLUE_LUMA_RAMP_START = 64,
    BLUE_LUMA_RAMP_COUNT = 16,
    BLUE_LUMA_RAMP_STOP = BLUE_LUMA_RAMP_START + 15,

    YELLOW_LUMA_RAMP_START = 80,
    YELLOW_LUMA_RAMP_COUNT = 16,
    YELLOW_LUMA_RAMP_STOP = YELLOW_LUMA_RAMP_START + 15,

    CYAN_LUMA_RAMP_START = 96,
    CYAN_LUMA_RAMP_COUNT = 16,
    CYAN_LUMA_RAMP_STOP = CYAN_LUMA_RAMP_START + 15,

    MAGENTA_LUMA_RAMP_START = 112,
    MAGENTA_LUMA_RAMP_COUNT = 16,
    MAGENTA_LUMA_RAMP_STOP = MAGENTA_LUMA_RAMP_START + 15
};

/// @brief Text rotation
enum text_rotation {
    DEGREE_0,
    DEGREE_90,
    DEGREE_180,
    DEGREE_270
};

/*! Data formats for the convert function */
enum device_format {
    STRAIGHT_THROUGH,      //!< Just copy the data as is.
    RGB565_8BIT_SERIAL,    //!< RGB565 pixel data is stored from left to right, each two bytes containing 1 pixel value in the x direction.
                           //   Byte encoding: 0xRRRRRGGG 0xGGGBBBBB
    RGB666_8BIT_SERIAL_1,  //!< RGB565 pixel data is stored from left to right, each three bytes containing 1 pixel values in the x direction.
                           //   Byte encoding: 0x00RRRRRR 0x00GGGGGG 0x00BBBBBB
    RGB666_8BIT_SERIAL_2   //!< RGB565 pixel data is stored from left to right, each three bytes containing 1 pixel values in the x direction.
                           //   Byte encoding: 0xRRRRRR00 0xGGGGGG00 0xBBBBBB00
};

/// @brief Core class of constixel, holds a buffer of an image of a certain size and format.
/// @tparam T Type of the image buffer. One of format_1bit, format_2bit, format_4bit or format_8bit
/// @tparam W Width in pixels
/// @tparam H Height in pixels
/// @tparam S Scale factor for sixel output
/// @tparam GR boolean to indicate if palette should be grayscale
template <template <size_t, size_t, int32_t, bool> class T, size_t W, size_t H, int32_t S = 1, bool GR = false>
class image {
    static_assert(sizeof(W) >= sizeof(uint32_t));
    static_assert(sizeof(H) >= sizeof(uint32_t));

    static_assert(W > 0 && H > 0);
    static_assert((W * S) <= 16384 && (H * S) <= 16384);
    static_assert(S >= 1 && S <= 256);

 public:
    /**
     * \brief Is this image grayscale?
     * \return Is this image grayscale?
     */
    [[nodiscard]] static constexpr size_t grayscale() {
        return T<W, H, S, GR>::grayscale;
    }

    /**
     * \brief Bit depth of the image.
     * \return Bit depth of the image.
     */
    [[nodiscard]] static constexpr size_t bit_depth() {
        return T<W, H, S, GR>::bits_per_pixel;
    }

    /**
     * \brief Size in bytes of the image data.
     * \return Size in bytes of the image data.
     */
    [[nodiscard]] static constexpr size_t size() {
        return T<W, H, S, GR>::image_size;
    }

    /**
     * \brief Bytes per line in the image data.
     * \return Bytes per line in the image data.
     */
    [[nodiscard]] static constexpr size_t bytes_per_line() {
        return T<W, H, S, GR>::bytes_per_line;
    }

    /**
     * \brief Width in pixels of the image.
     * \return Width in pixels of the image.
     */
    [[nodiscard]] static constexpr int32_t width() {
        return static_cast<int32_t>(W);
    }

    /**
     * \brief Height in pixels of the image.
     * \return Height in pixels of the image.
     */
    [[nodiscard]] static constexpr int32_t height() {
        return static_cast<int32_t>(H);
    }

    /**
     * \brief Clear the image by setting everything to zero.
     */
    constexpr void clear() {
        data.fill(0);
    }

    /**
     * \brief Get a reference to the underlying raw data of the image.
     * \return Reference to the data array which contains the raw image data.
     */
    [[nodiscard]] constexpr std::array<uint8_t, T<W, H, S, GR>::image_size> &data_ref() {
        return data;
    }

    /**
     * \brief Returns a clone of this image. Data is copied.
     * \return cloned image instance.
     */
    [[nodiscard]] constexpr image<T, W, H, S, GR> clone() const {
        return *this;
    }
    /**
     * \brief Copy source image into this instance. No compositing occurs.
     * \param src source image.
     */
    constexpr void copy(const image<T, W, H, S, GR> &src) {
        data = src.data;
    }

    /**
     * \brief Copy raw source data into this instance. No compositing occurs.
     * \tparam BYTE_SIZE Amount data in the source data. Typically a sizeof() of an array. Must match size()
     * \param src source data.
     */
    template <size_t BYTE_SIZE>
    constexpr void copy(const uint8_t *src) {
        static_assert(size() == BYTE_SIZE, "Copied length much match the image size");
        for (size_t c = 0; c < BYTE_SIZE; c++) {
            data.data()[c] = src[c];
        }
    }

    /**
     * \brief Draw a line with the specified color and thickness.
     * \param x0 starting x-coordinate in pixels.
     * \param y0 starting x-coordinate in pixels.
     * \param x1 ending x-coordinate in pixels.
     * \param y1 ending x-coordinate in pixels.
     * \param col Color palette index to use.
     * \param stroke_width width of the stroke in pixels.
     */
    constexpr void draw_line(int32_t x0, int32_t y0, int32_t x1, int32_t y1, uint8_t col, uint32_t stroke_width = 1) {
        int32_t steep = abs(y1 - y0) > abs(x1 - x0);

        if (steep) {
            std::swap(x0, y0);
            std::swap(x1, y1);
        }

        if (x0 > x1) {
            std::swap(x0, x1);
            std::swap(y0, y1);
        }

        int32_t dx, dy;
        dx = x1 - x0;
        dy = abs(y1 - y0);

        int32_t err = dx / 2;
        int32_t ystep;
        if (y0 < y1) {
            ystep = 1;
        } else {
            ystep = -1;
        }

        if (stroke_width == 1) {
            for (; x0 <= x1; x0++) {
                if (steep) {
                    plot(y0, x0, col);
                } else {
                    plot(x0, y0, col);
                }
                err -= dy;
                if (err < 0) {
                    y0 += ystep;
                    err += dx;
                }
            }
        } else if (stroke_width > 1) {
            for (; x0 <= x1; x0++) {
                if (steep) {
                    fill_circle(y0, x0, (stroke_width + 1) / 2, col);
                } else {
                    fill_circle(x0, y0, (stroke_width + 1) / 2, col);
                }
                err -= dy;
                if (err < 0) {
                    y0 += ystep;
                    err += dx;
                }
            }
        }
    }

    /**
     * \brief Draw a line with the specified color and thickness.
     * \param rect rectangle containing the line coordinates in pixels.
     * \param col Color palette index to use.
     * \param stroke_width width of the stroke in pixels.
     */
    constexpr void draw_line(const rect<int32_t> &rect, uint8_t col, uint32_t stroke_width = 1) {
        draw_line(rect.x, rect.y, rect.x + rect.w, rect.y + rect.h, col, stroke_width);
    }

    /**
     * \brief Draw a 1-pixel wide antialiased line with the specified color and thickness.
     * \param x0 starting x-coordinate in pixels.
     * \param y0 starting x-coordinate in pixels.
     * \param x1 ending x-coordinate in pixels.
     * \param y1 ending x-coordinate in pixels.
     * \param col Color palette index to use.
     */
    constexpr void draw_line_aa(int32_t x0, int32_t y0, int32_t x1, int32_t y1, uint8_t col) {
        auto ipart = [](float x) {
            return std::floor(x);
        };
        auto fpart = [](float x) {
            return x - std::floor(x);
        };
        auto rfpart = [&](float x) {
            return 1.0f - fpart(x);
        };

        bool steep = std::fabs(y1 - y0) > std::fabs(x1 - x0);
        if (steep) {
            std::swap(x0, y0);
            std::swap(x1, y1);
        }
        if (x0 > x1) {
            std::swap(x0, x1);
            std::swap(y0, y1);
        }

        float Rl = format.quant.linearpal.at((col & ((1UL << format.bits_per_pixel) - 1)) * 3 + 0);
        float Gl = format.quant.linearpal.at((col & ((1UL << format.bits_per_pixel) - 1)) * 3 + 1);
        float Bl = format.quant.linearpal.at((col & ((1UL << format.bits_per_pixel) - 1)) * 3 + 2);

        float dx = static_cast<float>(x1 - x0);
        float dy = static_cast<float>(y1 - y0);
        float gradient = dx == 0.0f ? 1.0f : dy / dx;

        auto color_compose = [&](int32_t x, int32_t y, float a) {
            if (a < epsilon_low) {
                return;
            } else if (a >= epsilon_high) {
                plot(x, y, col);
            } else {
                compose(x, y, a, Rl * a, Gl * a, Bl * a);
            }
        };

        // first endpoint
        float xend = static_cast<float>(x0);
        float yend = static_cast<float>(y0) + gradient * (xend - static_cast<float>(x0));
        float xgap = rfpart(static_cast<float>(x0) + 0.5f);
        int32_t xpxl1 = static_cast<int32_t>(xend);
        int32_t ypxl1 = static_cast<int32_t>(ipart(yend));
        if (steep) {
            color_compose(ypxl1, xpxl1, rfpart(yend) * xgap);
            color_compose(ypxl1 + 1, xpxl1, fpart(yend) * xgap);
        } else {
            color_compose(xpxl1, ypxl1, rfpart(yend) * xgap);
            color_compose(xpxl1, ypxl1 + 1, fpart(yend) * xgap);
        }
        float intery = yend + gradient;

        // second endpoint
        xend = static_cast<float>(x1);
        yend = static_cast<float>(y1) + gradient * (xend - static_cast<float>(x1));
        xgap = fpart(static_cast<float>(x1) + 0.5f);
        int32_t xpxl2 = static_cast<int32_t>(xend);
        int32_t ypxl2 = static_cast<int32_t>(ipart(yend));
        if (steep) {
            color_compose(ypxl2, xpxl2, rfpart(yend) * xgap);
            color_compose(ypxl2 + 1, xpxl2, fpart(yend) * xgap);
        } else {
            color_compose(xpxl2, ypxl2, rfpart(yend) * xgap);
            color_compose(xpxl2, ypxl2 + 1, fpart(yend) * xgap);
        }

        // main loop
        for (int32_t x = xpxl1 + 1; x < xpxl2; ++x) {
            int32_t y = static_cast<int32_t>(ipart(intery));
            if (steep) {
                color_compose(y, x, rfpart(intery));
                color_compose(y + 1, x, fpart(intery));
            } else {
                color_compose(x, y, rfpart(intery));
                color_compose(x, y + 1, fpart(intery));
            }
            intery += gradient;
        }
    }

    /**
     * \brief Draw an antialiased line with the specified color and thickness.
     * \param rect rectangle containing the line coordinates in pixels.
     * \param col Color palette index to use.
     */
    constexpr void draw_line_aa(const rect<int32_t> &rect, uint8_t col) {
        draw_line_aa(rect.x, rect.y, rect.x + rect.w, rect.y + rect.h, col);
    }

    /**
     * \brief Plot a single pixel at the specified coordinates used the supplied color.
     * \param x x-coordinate in pixels.
     * \param y y-coordinate in pixels.
     * \param col Color palette index to use.
     */
    constexpr void plot(int32_t x, int32_t y, uint8_t col) {
        if (static_cast<uint32_t>(x) >= static_cast<uint32_t>(W) || static_cast<uint32_t>(y) >= static_cast<uint32_t>(H)) {
            return;
        }
        T<W, H, S, GR>::plot(data, static_cast<uint32_t>(x), static_cast<uint32_t>(y), col);
    }

    /**
     * \brief Return closest match in the color palette based on the supplied red, green and blue values.
     * \param r Red value (0-255)
     * \param g Green value (0-255)
     * \param b Blue value (0-255)
     * \return The closest matching color palette index.
     */
    [[nodiscard]] constexpr uint8_t get_nearest_color(uint8_t r, uint8_t g, uint8_t b) const {
        return format.quant.nearest(r, g, b);
    }

    /**
     * \brief Return the width of a string using the specified font in the template parameter.
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available.
     * \param str UTF-8 string.
     * \return width of the string in pixels
     */
    template <typename FONT, bool KERNING = false>
    [[nodiscard]] constexpr int32_t string_width(const char *str) {
        int32_t x = 0;
        while (*str != 0) {
            uint32_t utf32 = 0;
            str = get_next_utf32(str, &utf32);
            uint32_t index = 0;
            if (lookup_glyph<FONT>(utf32, &index)) {
                const char_info<typename FONT::char_info_type> &ch_info = FONT::char_table.at(index);
                if (*str == 0) {
                    x += static_cast<int32_t>(ch_info.width);
                } else {
                    x += static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        x += get_kerning<FONT>(utf32, str);
                    }
                }
            }
        }
        return x;
    }

    /**
     * \brief Draw text at the specified coordinate. The template parameter selects which mono font to use. Only format_8bit targets are supported.
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available.
     * \tparam ROTATION Rotation around the x/y coordinate. Can be text_rotation::DEGREE_0, text_rotation::DEGREE_90, text_rotation::DEGREE_180 or
     * text_rotation::DEGREE_270
     * \param x Starting x-coordinate in pixels.
     * \param y Starting y-coordinate in pixels.
     * \param str UTF-8 string.
     * \param col Color palette index to use.
     * \return Returns the new caret position.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_mono(int32_t x, int32_t y, const char *str, uint8_t col) {
        static_assert(FONT::mono == true, "Can't use an antialiased font to draw mono/pixelized text.");
        while (*str != 0) {
            uint32_t utf32 = 0;
            str = get_next_utf32(str, &utf32);
            uint32_t index = 0;
            if (lookup_glyph<FONT>(utf32, &index)) {
                const char_info<typename FONT::char_info_type> &ch_info = FONT::char_table.at(index);
                draw_char_mono<FONT, ROTATION>(x, y, ch_info, col);
                if constexpr (ROTATION == DEGREE_0) {
                    x += static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        x += get_kerning<FONT>(utf32, str);
                    }
                } else if constexpr (ROTATION == DEGREE_90) {
                    y += static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        y += get_kerning<FONT>(utf32, str);
                    }
                } else if constexpr (ROTATION == DEGREE_180) {
                    x -= static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        x -= get_kerning<FONT>(utf32, str);
                    }
                } else if constexpr (ROTATION == DEGREE_270) {
                    y -= static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        y -= get_kerning<FONT>(utf32, str);
                    }
                }
            }
        }
        if constexpr (ROTATION == DEGREE_90 || ROTATION == DEGREE_270) {
            return y;
        } else {
            return x;
        }
    }

    /**
     * \brief Draw text centered at the specified coordinate. The template parameter selects which mono font to use. Only format_8bit targets are supported.
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available.
     * \tparam ROTATION Rotation around the x/y coordinate. Can be text_rotation::DEGREE_0, text_rotation::DEGREE_90, text_rotation::DEGREE_180 or
     * text_rotation::DEGREE_270
     * \param cx Center x-coordinate in pixels.
     * \param y Starting y-coordinate in pixels.
     * \param str UTF-8 string.
     * \param col Color palette index to use.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr void draw_string_centered_mono(int32_t cx, int32_t y, const char *str, uint8_t col) {
        if constexpr (ROTATION == DEGREE_0) {
            draw_string_mono<FONT, KERNING, ROTATION>(cx - string_width<FONT, KERNING>(str) / 2, y, str, col);
        } else if constexpr (ROTATION == DEGREE_180) {
            draw_string_mono<FONT, KERNING, ROTATION>(cx + string_width<FONT, KERNING>(str) / 2, y, str, col);
        } else if constexpr (ROTATION == DEGREE_90) {
            draw_string_mono<FONT, KERNING, ROTATION>(cx, y + string_width<FONT, KERNING>(str) / 2, str, col);
        } else if constexpr (ROTATION == DEGREE_270) {
            draw_string_mono<FONT, KERNING, ROTATION>(cx, y - string_width<FONT, KERNING>(str) / 2, str, col);
        }
    }

    /**
     * \brief Draw antialiased text at the specified coordinate. The template parameter selects which antialiased font to use. Only format_8bit targets are
     * supported.
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available.
     * \tparam ROTATION Rotation around the x/y coordinate. Can be text_rotation::DEGREE_0, text_rotation::DEGREE_90, text_rotation::DEGREE_180 or
     * text_rotation::DEGREE_270
     * \param x Starting x-coordinate in pixels.
     * \param y Starting y-coordinate in pixels.
     * \param str UTF-8 string.
     * \param col Color palette index to use.
     * \return Returns the new caret position
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_aa(int32_t x, int32_t y, const char *str, uint8_t col) {
        static_assert(FONT::mono == false, "Can't use a mono font to draw antialiased text.");
        while (*str != 0) {
            uint32_t utf32 = 0;
            str = get_next_utf32(str, &utf32);
            uint32_t index = 0;
            if (lookup_glyph<FONT>(utf32, &index)) {
                const char_info<typename FONT::char_info_type> &ch_info = FONT::char_table.at(index);
                draw_char_aa<FONT, ROTATION>(x, y, ch_info, col);
                if constexpr (ROTATION == DEGREE_0) {
                    x += static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        x += get_kerning<FONT>(utf32, str);
                    }
                } else if constexpr (ROTATION == DEGREE_90) {
                    y += static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        y += get_kerning<FONT>(utf32, str);
                    }
                } else if constexpr (ROTATION == DEGREE_180) {
                    x -= static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        x -= get_kerning<FONT>(utf32, str);
                    }
                } else if constexpr (ROTATION == DEGREE_270) {
                    y -= static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        y -= get_kerning<FONT>(utf32, str);
                    }
                }
            }
        }
        if constexpr (ROTATION == DEGREE_90 || ROTATION == DEGREE_270) {
            return y;
        } else {
            return x;
        }
    }

    /**
     * \brief Draw antialiased text centered at the specified coordinate. The template parameter selects which antialiased font to use. Only format_8bit
     * targets are supported.
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available.
     * \tparam ROTATION Rotation around the x/y coordinate. Can be text_rotation::DEGREE_0, text_rotation::DEGREE_90, text_rotation::DEGREE_180 or
     * text_rotation::DEGREE_270
     * \param cx Senter x-coordinate in pixels.
     * \param y Starting y-coordinate in pixels.
     * \param str UTF-8 string.
     * \param col Color palette index to use.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr void draw_string_centered_aa(int32_t cx, int32_t y, const char *str, uint8_t col) {
        if constexpr (ROTATION == DEGREE_0) {
            draw_string_aa<FONT, KERNING, ROTATION>(cx - string_width<FONT, KERNING>(str) / 2, y, str, col);
        } else if constexpr (ROTATION == DEGREE_180) {
            draw_string_aa<FONT, KERNING, ROTATION>(cx + string_width<FONT, KERNING>(str) / 2, y, str, col);
        } else if constexpr (ROTATION == DEGREE_90) {
            draw_string_aa<FONT, KERNING, ROTATION>(cx, y + string_width<FONT, KERNING>(str) / 2, str, col);
        } else if constexpr (ROTATION == DEGREE_270) {
            draw_string_aa<FONT, KERNING, ROTATION>(cx, y - string_width<FONT, KERNING>(str) / 2, str, col);
        }
    }

    /**
     * \brief Fill a rectangle with the specified color.
     * \param x Starting x-coordinate in pixels.
     * \param y Starting y-coordinate in pixels.
     * \param w Sidth of the rectangle.
     * \param h Seight of the rectangle.
     * \param col Color palette index to use.
     */
    constexpr void fill_rect(int32_t x, int32_t y, int32_t w, int32_t h, uint8_t col) {
        if (y < 0) {
            h += y;
            y = 0;
        }
        if (y + h >= static_cast<int32_t>(H)) {
            h = static_cast<int32_t>(H) - y;
        }
        h += y;
        for (; y < h; y++) {
            span(x, w, y, col);
        }
    }

    /**
     * \brief Fill a rectangle with the specified color.
     * \param rect Rectangle containing the line coordinates in pixels.
     * \param col Color palette index to use.
     */
    constexpr void fill_rect(const rect<int32_t> &rect, uint8_t col) {
        fill_rect(rect.x, rect.y, rect.w, rect.h, col);
    }

    /**
     * \brief Draw a stroked rectangle with the specified color and stroke width.
     * \param x Starting x-coordinate in pixels.
     * \param y Starting y-coordinate in pixels.
     * \param w Width of the rectangle.
     * \param h Height of the rectangle.
     * \param col Color palette index to use.
     * \param stroke_width Width of the stroke in pixels.
     */
    constexpr void stroke_rect(int32_t x, int32_t y, int32_t w, int32_t h, uint8_t col, uint32_t stroke_width = 1) {
        draw_line(x, y, x + w, y, col, stroke_width);
        draw_line(x + w, y, x + w, y + h, col, stroke_width);
        draw_line(x + w, y + h, x, y + h, col, stroke_width);
        draw_line(x, y + h, x, y, col, stroke_width);
    }

    /**
     * \brief Draw a stroked rectangle with the specified color and stroke width.
     * \param rect Rectangle containing the line coordinates in pixels.
     * \param col Color palette index to use.
     * \param stroke_width Width of the stroke in pixels.
     */
    constexpr void stroke_rect(const rect<int32_t> &rect, uint8_t col, uint32_t stroke_width = 1) {
        stroke_rect(rect.x, rect.y, rect.w, rect.h, col, stroke_width);
    }

    /**
     * \brief Fill a circle with the specified radius and color.
     * \param cx Center x-coordinate of the circle in pixels.
     * \param cy Center y-coordinate of the circle in pixels.
     * \param radius radius of the circle in pixels.
     * \param col Color palette index to use.
     */
    constexpr void fill_circle(int32_t cx, int32_t cy, int32_t radius, uint8_t col) {
        if (radius == 1) {
            fill_rect(cx - 1, cy - 1, 2, 2, col);
            return;
        }
        fill_arc(cx, cy - 1, radius, 9, 0, col);
        fill_arc(cx, cy, radius, 10, 0, col);
        fill_arc(cx, cy - 1, radius, 1, -1, col);
        fill_arc(cx, cy, radius, 2, -1, col);
    }

    /**
     * \brief Fill a circle using antialiasing with with the specified radius and color. Only format_8bit targets are supported.
     * \param cx Center x-coordinate of the circle in pixels.
     * \param cy Center y-coordinate of the circle in pixels.
     * \param radius Radius of the circle in pixels.
     * \param col Color palette index to use.
     */
    constexpr void fill_circle_aa(int32_t cx, int32_t cy, int32_t radius, uint8_t col) {
        fill_circle_aa_int(cx, cy, radius, 0, 0, col);
    }

    /**
     * \brief Fill a rounded rectangle with the specified color.
     * \param x Starting x-coordinate in pixels.
     * \param y Starting y-coordinate in pixels.
     * \param w Width of the rectangle.
     * \param h Height of the rectangle.
     * \param radius Radius of the rounded corners.
     * \param col Color palette index to use.
     */
    constexpr void fill_round_rect(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius, uint8_t col) {
        int32_t cr = std::min((w) / 2, std::min((w) / 2, radius));
        int32_t dx = w - cr * 2;
        int32_t dy = h - cr * 2;
        fill_arc(x + cr, y + h - cr - 1, cr, 9, 0, col);
        fill_arc(x + cr, y + cr, cr, 10, 0, col);
        fill_arc(x + w - cr, y + h - cr - 1, cr, 1, -1, col);
        fill_arc(x + w - cr, y + cr, cr, 2, -1, col);
        fill_rect(x, y + cr, cr, dy, col);
        fill_rect(x + w - cr, y + cr, cr, dy, col);
        fill_rect(x + cr, y, dx, h, col);
    }

    /**
     * \brief Fill a rounded rectangle with the specified color.
     * \param rect Rectangle containing the line coordinates in pixels
     * \param radius Radius of the rounded corners.
     * \param col Color palette index to use.
     */
    constexpr void fill_round_rect(const rect<int32_t> &rect, int32_t radius, uint8_t col) {
        stroke_rect(rect.x, rect.y, rect.w, rect.h, radius, col);
    }

    /**
     * \brief Fill a rounded rectangle using antialiasing with the specified color. Only format_8bit targets are supported.
     * \param x Starting x-coordinate in pixels.
     * \param y Starting y-coordinate in pixels.
     * \param w Width of the rectangle.
     * \param h Height of the rectangle.
     * \param radius Radius of the rounded corners.
     * \param col Color palette index to use.
     */
    constexpr void fill_round_rect_aa(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius, uint8_t col) {
        int32_t cr = std::min((w) / 2, std::min((w) / 2, radius));
        int32_t dx = w - cr * 2;
        int32_t dy = h - cr * 2;
        fill_circle_aa_int(x + cr, y + cr, cr, dx, dy, col);
        fill_rect(x, y + cr, cr, dy, col);
        fill_rect(x + w - cr, y + cr, cr, dy, col);
        fill_rect(x + cr, y, dx, h, col);
    }

    /**
     * \brief Fill a rounded rectangle using antialiasing with the specified color. Only format_8bit targets are supported.
     * \param rect Rectangle containing the line coordinates in pixels.
     * \param radius Radius of the rounded corners.
     * \param col Color palette index to use.
     */
    constexpr void fill_round_rect_aa(const rect<int32_t> &rect, int32_t radius, uint8_t col) {
        fill_round_rect_aa(rect.x, rect.y, rect.w, rect.h, radius, col);
    }

    /**
     * \brief Convert this instance to an equivalent RGBA8 data array.
     * \return RGBA8 array made of uint32_t values.
     */
    [[nodiscard]] constexpr std::array<uint32_t, W * H> RGBA_uint32() const {
        return T<W, H, S, GR>::RGBA_uint32(data);
    }

    /**
     * \brief Convert this instance to an equivalent RGBA8 data array.
     * \return RGBA8 array made of uint8_t values.
     */
    [[nodiscard]] constexpr std::array<uint8_t, W * H * 4> RGBA_uint8() const {
        return T<W, H, S, GR>::RGBA_uint8(data);
    }

    /**
     * \brief Flip the contents of this image horizontally.
     */
    constexpr void flip_h() {
        uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < bytes_per_line() / 2; x++) {
                uint8_t a = T<W, H, S, GR>::reverse(ptr[x]);
                uint8_t b = T<W, H, S, GR>::reverse(ptr[bytes_per_line() - x - 1]);
                ptr[x] = b;
                ptr[bytes_per_line() - x - 1] = a;
            }
            ptr += bytes_per_line();
        }
    }

    /**
     * \brief Flip the contents of this image vertically.
     */
    constexpr void flip_v() {
        uint8_t *ptr = data.data();
        for (size_t x = 0; x < bytes_per_line(); x++) {
            for (size_t y = 0; y < H / 2; y++) {
                uint8_t *ptr_a = ptr + y * bytes_per_line();
                uint8_t *ptr_b = ptr + (H - y - 1) * bytes_per_line();
                uint8_t a = ptr_a[x];
                uint8_t b = ptr_b[x];
                ptr_a[x] = b;
                ptr_b[x] = a;
            }
        }
    }

    /**
     * \brief Flip the contents of this image horizontally & vertically.
     */
    constexpr void flip_hv() {
        uint8_t *ptr = data.data();
        for (size_t x = 0; x < bytes_per_line(); x++) {
            for (size_t y = 0; y < H / 2; y++) {
                uint8_t *ptr_a = ptr + y * bytes_per_line();
                uint8_t *ptr_b = ptr + (H - y - 1) * bytes_per_line();
                uint8_t a = T<W, H, S, GR>::reverse(ptr_a[x]);
                uint8_t b = T<W, H, S, GR>::reverse(ptr_b[bytes_per_line() - x - 1]);
                ptr_a[x] = b;
                ptr_b[bytes_per_line() - x - 1] = a;
            }
        }
    }

    /**
     * \brief Return a transposed version of this image.
     */
    constexpr image<T, H, W, S, GR> transpose() const {
        image<T, H, W, S, GR> transposed;
        static_assert(T<W, H, S, GR>::bits_per_pixel != 1 || ((H + 7) / 8) == transposed.bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 1 || ((W + 7) / 8) == bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 2 || ((H + 3) / 4) == transposed.bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 2 || ((W + 3) / 4) == bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 4 || ((H + 1) / 2) == transposed.bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 4 || ((W + 1) / 2) == bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 8 || H == transposed.bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 8 || W == bytes_per_line());
        T<W, H, S, GR>::transpose(data.data(), transposed.data_ref().data());
        return transposed;
    }

    /**
     * \brief Transpose this image into another.
     */
    constexpr void transpose(image<T, H, W, S, GR> &dst) const {
        static_assert(T<W, H, S, GR>::bits_per_pixel != 1 || ((H + 7) / 8) == dst.bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 1 || ((W + 7) / 8) == bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 2 || ((H + 3) / 4) == dst.bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 2 || ((W + 3) / 4) == bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 4 || ((H + 1) / 2) == dst.bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 4 || ((W + 1) / 2) == bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 8 || H == dst.bytes_per_line());
        static_assert(T<W, H, S, GR>::bits_per_pixel != 8 || W == bytes_per_line());
        T<W, H, S, GR>::transpose(data, dst.data);
    }

    /**
     * \brief Blit an RGBA8 buffer into this instance using brute force color mapping.
     * \param x Starting x-coordinate position in the target instance in pixels
     * \param y Starting y-coordinate position in the target instance in pixels
     * \param w Width of the rectangle. If the width is smaller than the RGBA8 buffer content will be clipped.
     * \param h Weight of the rectangle. If the height is smaller than the RGBA8 buffer content will be clipped.
     * \param ptr Pointer to the RGBA8 buffer.
     * \param iw Width in pixels of the RGBA8 buffer.
     * \param ih Height in pixels of the RGBA8 buffer.
     * \param stride Bytes per line of pixels of the RGBA8 buffer.
     */
    constexpr void blit_RGBA(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{x, y, w, h};
        blitrect &= {0, 0, W, H};
        blitrect &= {x, y, iw, ih};
        T<W, H, S, GR>::blit_RGBA(data, blitrect, ptr, stride);
    }

    /**
     * \brief Blit an RGBA8 buffer into this instance using brute force color mapping.
     * \param dstrect Rectangular area in the target buffer to blit into. If the rectangle is smaller than the RGBA8 buffer, clipping occurs.
     * \param ptr Pointer to the RGBA8 buffer.
     * \param iw Width in pixels of the RGBA8 buffer.
     * \param ih Weight in pixels of the RGBA8 buffer.
     * \param stride Bytes per line of pixels of the RGBA8 buffer.
     */
    constexpr void blit_RGBA(const rect<int32_t> &dstrect, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{dstrect};
        blitrect &= {0, 0, W, H};
        blitrect &= {dstrect.x, dstrect.y, iw, ih};
        T<W, H, S, GR>::blit_RGBA(data, blitrect, ptr, stride);
    }

    /**
     * \brief Blit an RGBA8 buffer into this instance using brute force color mapping. Simple integer based diffusion is applied.
     * \param x Starting x-coordinate position in the target instance in pixels
     * \param y Starting y-coordinate position in the target instance in pixels
     * \param w Width of the rectangle. If the width is smaller than the RGBA8 buffer content will be clipped.
     * \param h Weight of the rectangle. If the height is smaller than the RGBA8 buffer content will be clipped.
     * \param ptr Pointer to the RGBA8 buffer.
     * \param iw Width in pixels of the RGBA8 buffer.
     * \param ih Height in pixels of the RGBA8 buffer.
     * \param stride Bytes per line of pixels of the RGBA8 buffer.
     */
    constexpr void blit_RGBA_diffused(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{x, y, w, h};
        blitrect &= {0, 0, W, H};
        blitrect &= {x, y, iw, ih};
        T<W, H, S, GR>::blit_RGBA_diffused(data, blitrect, ptr, stride);
    }

    /**
     * \brief Blit an RGBA8 buffer into this instance using brute force color mapping. Simple integer based diffusion is applied.
     * \param dstrect Rectangular area in the target buffer to blit into. If the rectangle is smaller than the RGBA8 buffer, clipping occurs.
     * \param ptr Pointer to the RGBA8 buffer.
     * \param iw Width in pixels of the RGBA8 buffer.
     * \param ih Weight in pixels of the RGBA8 buffer.
     * \param stride Bytes per line of pixels of the RGBA8 buffer.
     */
    constexpr void blit_RGBA_diffused(const rect<int32_t> &dstrect, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{dstrect};
        blitrect &= {0, 0, W, H};
        blitrect &= {dstrect.x, dstrect.y, iw, ih};
        T<W, H, S, GR>::blit_RGBA_diffused(data, blitrect, ptr, stride);
    }

    /**
     * \brief Blit an RGBA8 buffer into this instance using brute force color mapping. Diffusion in linear color space is applied.
     * \param x Starting x-coordinate position in the target instance in pixels
     * \param y Starting y-coordinate position in the target instance in pixels
     * \param w Width of the rectangle. If the width is smaller than the RGBA8 buffer content will be clipped.
     * \param h Weight of the rectangle. If the height is smaller than the RGBA8 buffer content will be clipped.
     * \param ptr Pointer to the RGBA8 buffer.
     * \param iw Width in pixels of the RGBA8 buffer.
     * \param ih Height in pixels of the RGBA8 buffer.
     * \param stride Bytes per line of pixels of the RGBA8 buffer.
     */
    constexpr void blit_RGBA_diffused_linear(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{x, y, w, h};
        blitrect &= {0, 0, W, H};
        blitrect &= {x, y, iw, ih};
        T<W, H, S, GR>::blit_RGBA_diffused_linear(data, blitrect, ptr, stride);
    }

    /**
     * \brief Blit an RGBA8 buffer into this instance using brute force color mapping. Diffusion in linear color space is applied.
     * \param dstrect Rectangular area in the target buffer to blit into. If the rectangle is smaller than the RGBA8 buffer, clipping occurs.
     * \param ptr Pointer to the RGBA8 buffer.
     * \param iw Width in pixels of the RGBA8 buffer.
     * \param ih Weight in pixels of the RGBA8 buffer.
     * \param stride Bytes per line of pixels of the RGBA8 buffer.
     */
    constexpr void blit_RGBA_diffused_linear(const rect<int32_t> &dstrect, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{dstrect};
        blitrect &= {0, 0, W, H};
        blitrect &= {dstrect.x, dstrect.y, iw, ih};
        T<W, H, S, GR>::blit_RGBA_diffused_linear(data, blitrect, ptr, stride);
    }

    /**
     * \brief Convert the current instance into a png image.
     * \param char_out A lambda function which consumes the png image data one byte at a time.
     */
    template <typename F>
    constexpr void png(F &&char_out) const {
        T<W, H, S, GR>::png(data, char_out);
    }

    /**
     * \brief Convert the current instance into a sixel stream.
     * \param char_out A lambda function which consumes the sixel stream data one byte at a time.
     */
    template <typename F>
    constexpr void sixel(F &&char_out) const {
        T<W, H, S, GR>::sixel(data, char_out, {0, 0, W, H});
    }

    /**
     * \brief Convert the current instance into a sixel stream.
     * \param char_out A lambda function which consumes the sixel stream data one byte at a time
     * \param rect Clipping rectangle; to only show a portion of the image.
     */
    template <typename F>
    constexpr void sixel(F &&char_out, const rect<int32_t> &rect) const {
        T<W, H, S, GR>::sixel(data, char_out, rect);
    }

    /**
     * \brief Convert the current instance into a sixel stream and output it to std::cout.
     */
    void sixel_to_cout() const {
        std::string out;
        T<W, H, S, GR>::sixel(data,
                              [&out](char ch) mutable {
                                  out.push_back(ch);
                              },
                              {0, 0, W, H});
        std::cout << out << std::endl;
    }

    /**
     * \brief Send a escape command to std::cout to clear the screen and scroll buffer of a vt100 compatible terminal.
     */
    void vt100_clear() const {
        std::cout << "\033[2J\033[H" << std::endl;
    }

    /**
     * \brief Send a escape command to std::cout to home the cursor of a vt100 compatible terminal.
     */
    void vt100_home() const {
        std::cout << "\033[H" << std::endl;
    }

    /**
     * \brief Convert the current instance into a png and display it in iTerm.
     */
    void png_to_iterm() const {
        std::string output;
        output.append("\033]1337;File=inline=1:");
        append_png_as_base64(output);
        output.append("\07");
        std::cout << output << std::endl;
    }

    /**
     * \brief Convert the current instance into a png and display it in a terminal with kitty graphics support.
     */
    void png_to_kitty() const {
        std::string base64{};
        std::string output{};
        append_png_as_base64(base64);
        bool first = true;
        for (; base64.length();) {
            if (first) {
                first = false;
                output.append("\033_Ga=T,f=100,");
            } else {
                output.append("\033_G");
            }
            output.append(base64.length() <= 4096 ? "m=0;" : "m=1;");
            size_t bytes_to_append = std::min(base64.length(), static_cast<size_t>(4096));
            output.append(base64.substr(0, bytes_to_append));
            base64.erase(0, bytes_to_append);
            output.append("\033\\");
        }
        std::cout << output << std::endl;
    }

    /**
     * \brief Convert the current instance into a byte stream formatted for embedded displays. This function will write chunk_length of data into dst and
     * update chunk_index on each call. You should pass chunk_index = 0 at the beginnig of the sequence. Returns true of there is more data, or false if
     * there is no data left. This function is typically used to drive interrupt driven DMA transfers.
     * \tparam dst_format The desired data format.
     * \param dst The buffer the data will be written to.
     * \param chunk_size The requested chunk size in bytes. This value has to be kept the same during a full conversion sequence.
     * \param chunk_actual The actual amount of bytes which were written into ptr during this call of the function.
     * \param chunk_index This value will increment after each call of the function. Has to be set to 0 on the first call.
     * \return true if there is more data to convert, false if we reached the end.
     */
    template <device_format dst_format>
    bool convert_chunk(char *dst, size_t chunk_size, size_t &chunk_actual, size_t &chunk_index) {
        (void)dst;
        (void)chunk_size;
        (void)chunk_actual;
        (void)chunk_index;
        (void)dst_format;
        return false;
    }

    /**
     * \brief Convert the current instance into a byte stream formatted for embedded displays.
     * \tparam dst_format The desired data format.
     * \param char_out A lambda function which consumes the data stream one byte at a time.
     */
    template <device_format dst_format, typename F>
    constexpr void convert(F &&char_out) {
        if constexpr (dst_format == STRAIGHT_THROUGH) {
            for (auto c : data) {
                char_out(static_cast<char>(c));
            }
        } else if constexpr (dst_format == RGB565_8BIT_SERIAL) {
            const uint8_t *ptr = data.data();
            for (size_t y = 0; y < H; y++) {
                for (size_t x = 0; x < W; x++) {
                    uint32_t col = format.get_col(ptr, x);
                    uint32_t a = ((((col >> 16) & 0xff) >> 3) << 3) | ((((col >> 8) & 0xff) >> 2) >> 3);
                    uint32_t b = (((((col >> 8) & 0xff) >> 2) & 0x7) << 5) | ((col & 0xff) >> 3);
                    char_out(static_cast<char>(a));
                    char_out(static_cast<char>(b));
                }
                ptr += format.bytes_per_line;
            }
        } else if constexpr (dst_format == RGB666_8BIT_SERIAL_1) {
            const uint8_t *ptr = data.data();
            for (size_t y = 0; y < H; y++) {
                for (size_t x = 0; x < W; x++) {
                    uint32_t col = format.get_col(ptr, x);
                    char_out(static_cast<char>(((col >> 16) & 0xff) >> 2));
                    char_out(static_cast<char>(((col >> 8) & 0xff) >> 2));
                    char_out(static_cast<char>(((col >> 0) & 0xff) >> 2));
                }
                ptr += format.bytes_per_line;
            }
        } else if constexpr (dst_format == RGB666_8BIT_SERIAL_2) {
            const uint8_t *ptr = data.data();
            for (size_t y = 0; y < H; y++) {
                for (size_t x = 0; x < W; x++) {
                    uint32_t col = format.get_col(ptr, x);
                    char_out(static_cast<char>(((col >> 16) & 0xff) >> 2) << 2);
                    char_out(static_cast<char>(((col >> 8) & 0xff) >> 2) << 2);
                    char_out(static_cast<char>(((col >> 0) & 0xff) >> 2) << 2);
                }
                ptr += format.bytes_per_line;
            }
        }
    }

 private:
    template <typename abs_T>
    [[nodiscard]] static constexpr abs_T abs(abs_T v) {
        return v < 0 ? -v : v;
    }

    template <typename FONT>
    constexpr bool lookup_glyph(uint32_t utf32, uint32_t *glyph_index) {
        auto index = FONT::glyph_tree.lookup(static_cast<FONT::lookup_type>(utf32));
        if (index == FONT::glyph_tree.invalid) {
            index = FONT::glyph_tree.lookup(static_cast<FONT::lookup_type>(0xFFFD));
            if (index == FONT::glyph_tree.invalid) {
                index = FONT::glyph_tree.lookup(static_cast<FONT::lookup_type>(0x0000));
                if (index == FONT::glyph_tree.invalid) {
                    return false;
                }
            }
        }
        *glyph_index = index;
        return true;
    }

    template <typename FONT>
    constexpr int32_t get_kerning(uint32_t utf32, const char *str) const {
        if (FONT::kerning_tree.byte_size() > 0) {
            uint32_t utf_l = utf32;
            uint32_t utf_r = 0;
            get_next_utf32(str, &utf_r);
            auto amount = FONT::kerning_tree.lookup(static_cast<FONT::kerning_lookup_type>(utf_l << FONT::kerning_code_shift | utf_r));
            if (amount != FONT::kerning_tree.invalid) {
                return static_cast<int32_t>(static_cast<FONT::kerning_amount_type>(amount) -
                                            static_cast<FONT::kerning_amount_type>(FONT::kerning_amount_offset));
            }
        }
        return 0;
    }

    constexpr const char *get_next_utf32(const char *str, uint32_t *utf32) const {
        *utf32 = 0;
        uint32_t lead = static_cast<uint32_t>(*str) & 0xFF;
        if (lead < 0x80) {
            *utf32 = lead;
            str += 1;
        } else if ((lead >> 5) == 0x06 && str[1] != 0) {
            *utf32 = ((lead & 0x1F) << 6) | (static_cast<uint32_t>(str[1]) & 0x3F);
            str += 2;
        } else if ((lead >> 4) == 0x0E && str[1] != 0 && str[2] != 0) {
            *utf32 = ((lead & 0x0F) << 12) | ((static_cast<uint32_t>(str[1]) & 0x3F) << 6) | (static_cast<uint32_t>(str[2]) & 0x3F);
            str += 3;
        } else if ((lead >> 3) == 0x1E && str[1] != 0 && str[2] != 0 && str[3] != 0) {
            *utf32 = ((lead & 0x07) << 18) | ((static_cast<uint32_t>(str[1]) & 0x3F) << 12) | ((static_cast<uint32_t>(str[2]) & 0x3F) << 6) |
                     (static_cast<uint32_t>(str[3]) & 0x3F);
            str += 4;
        } else {
            str += 1;
        }
        return str;
    }

    constexpr void compose(int32_t x, int32_t y, float cola, float colr, float colg, float colb) {
        if (static_cast<uint32_t>(x) >= static_cast<uint32_t>(W) || static_cast<uint32_t>(y) >= static_cast<uint32_t>(H)) {
            return;
        }
        T<W, H, S, GR>::compose(data, static_cast<uint32_t>(x), static_cast<uint32_t>(y), cola, colr, colg, colb);
    }

    void append_png_as_base64(std::string &output) const {
        size_t buffer = 0;
        size_t bits_collected = 0;
        T<W, H, S, GR>::png(data, [&buffer, &bits_collected, &output](char byte) mutable {
            static constexpr char base64_chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

            buffer = (buffer << 8) | static_cast<uint8_t>(byte);
            bits_collected += 8;
            while (bits_collected >= 6) {
                bits_collected -= 6;
                output.push_back(base64_chars[(buffer >> bits_collected) & 0x3F]);
            }
            return [&output]() mutable {
                if (output.capacity() == 0)
                    return;
                if (output.size() % 4) {
                    size_t padding = 4 - output.size() % 4;
                    while (padding--) output.push_back('=');
                }
            };
        });
    }

    constexpr void fill_circle_aa_int(int32_t cx, int32_t cy, int32_t r, int32_t ox, int32_t oy, uint8_t col) {
        int32_t x0 = cx - r - 1;
        int32_t y0 = cy - r - 1;
        float Rl = format.quant.linearpal.at((col & ((1UL << format.bits_per_pixel) - 1)) * 3 + 0);
        float Gl = format.quant.linearpal.at((col & ((1UL << format.bits_per_pixel) - 1)) * 3 + 1);
        float Bl = format.quant.linearpal.at((col & ((1UL << format.bits_per_pixel) - 1)) * 3 + 2);
        for (int32_t y = y0; y <= cy; ++y) {
            for (int32_t x = x0; x <= cx; ++x) {
                float dx = (static_cast<float>(x) + 0.5f) - static_cast<float>(cx);
                float dy = (static_cast<float>(y) + 0.5f) - static_cast<float>(cy);
                float dist_sq = dx * dx + dy * dy;
                if (dist_sq > (r + 0.5f) * (r + 0.5f)) {
                    continue;
                }
                if (dist_sq < (r - 0.5f) * (r - 0.5f)) {
                    int32_t lx = x;
                    int32_t ly = y;
                    int32_t rx = cx + (x0 - x) + r + ox;
                    int32_t ry = cy + (y0 - y) + r + oy;
                    plot(lx, ly, col);
                    plot(rx, ly, col);
                    plot(rx, ry, col);
                    plot(lx, ry, col);
                    continue;
                }
                float a = static_cast<float>(r);
                if (std::is_constant_evaluated()) {
                    a -= fast_sqrtf(dist_sq);
                } else {
                    a -= std::sqrt(dist_sq);
                }
                a = std::clamp(a + 0.5f, 0.0f, 1.0f);
                if (a >= epsilon_low) {
                    int32_t lx = x;
                    int32_t ly = y;
                    int32_t rx = cx + (x0 - x) + r + ox;
                    int32_t ry = cy + (y0 - y) + r + oy;
                    if (a >= epsilon_high) {
                        plot(lx, ly, col);
                        plot(rx, ly, col);
                        plot(rx, ry, col);
                        plot(lx, ry, col);
                    } else {
                        compose(lx, ly, a, Rl * a, Gl * a, Bl * a);
                        compose(rx, ly, a, Rl * a, Gl * a, Bl * a);
                        compose(rx, ry, a, Rl * a, Gl * a, Bl * a);
                        compose(lx, ry, a, Rl * a, Gl * a, Bl * a);
                    }
                }
            }
        }
    }

    template <typename FONT, text_rotation ROTATION>
    constexpr void draw_char_mono(int32_t x, int32_t y, const char_info<typename FONT::char_info_type> &ch, uint8_t col) {
        static_assert(FONT::mono == true, "Can't use an antialiased font to draw mono/pixelized text.");
        int32_t ch_data_off = static_cast<int32_t>(ch.y) * static_cast<int32_t>(FONT::glyph_bitmap_stride) + static_cast<int32_t>(ch.x) / 8;
        if constexpr (ROTATION == DEGREE_0) {
            x += ch.xoffset;
            y += ch.yoffset;
            const int32_t x2 = x + static_cast<int32_t>(ch.width);
            const int32_t y2 = y + static_cast<int32_t>(ch.height);
            for (int32_t yy = y; yy < y2; yy++) {
                for (int32_t xx = x; xx < x2; xx++) {
                    const int32_t x_off = (xx - x) + ch.x % 8;
                    const int32_t bit_index = 7 - (x_off % 8);
                    const size_t byte_index = static_cast<size_t>(ch_data_off + x_off / 8);
                    if (byte_index < (FONT::glyph_bitmap_stride * FONT::glyph_bitmap_height)) {
                        const uint8_t a = (FONT::glyph_bitmap[byte_index] >> bit_index) & 1;
                        if (a) {
                            plot(xx, yy, col);
                        }
                    }
                }
                ch_data_off += FONT::glyph_bitmap_stride;
            }
        } else if constexpr (ROTATION == DEGREE_180) {
            x -= ch.xoffset;
            y -= FONT::ascent;
            const int32_t x2 = x - static_cast<int32_t>(ch.width);
            const int32_t y2 = y + static_cast<int32_t>(ch.height);
            for (int32_t yy = y2 - 1; yy >= y; yy--) {
                for (int32_t xx = x2; xx < x; xx++) {
                    const int32_t x_off = (ch.width - (xx - x2) - 1) + ch.x % 8;
                    const int32_t bit_index = 7 - (x_off % 8);
                    const size_t byte_index = static_cast<size_t>(ch_data_off + x_off / 8);
                    if (byte_index < (FONT::glyph_bitmap_stride * FONT::glyph_bitmap_height)) {
                        const uint8_t a = (FONT::glyph_bitmap[byte_index] >> bit_index) & 1;
                        if (a) {
                            plot(xx, yy, col);
                        }
                    }
                }
                ch_data_off += FONT::glyph_bitmap_stride;
            }
        } else if constexpr (ROTATION == DEGREE_90) {
            x -= FONT::ascent;
            y += ch.xoffset;
            const int32_t x2 = x + static_cast<int32_t>(ch.height);
            const int32_t y2 = y + static_cast<int32_t>(ch.width);
            for (int32_t xx = x2 - 1; xx >= x; xx--) {
                for (int32_t yy = y; yy < y2; yy++) {
                    const int32_t x_off = (yy - y) + ch.x % 8;
                    const int32_t bit_index = 7 - (x_off % 8);
                    const size_t byte_index = static_cast<size_t>(ch_data_off + x_off / 8);
                    if (byte_index < (FONT::glyph_bitmap_stride * FONT::glyph_bitmap_height)) {
                        const uint8_t a = (FONT::glyph_bitmap[byte_index] >> bit_index) & 1;
                        if (a) {
                            plot(xx, yy, col);
                        }
                    }
                }
                ch_data_off += FONT::glyph_bitmap_stride;
            }
        } else if constexpr (ROTATION == DEGREE_270) {
            x += ch.yoffset;
            y -= ch.xoffset;
            const int32_t x2 = x + static_cast<int32_t>(ch.height);
            const int32_t y2 = y - static_cast<int32_t>(ch.width);
            for (int32_t xx = x; xx < x2; xx++) {
                for (int32_t yy = y2; yy < y; yy++) {
                    const int32_t x_off = (ch.width - (yy - y2) - 1) + ch.x % 8;
                    const int32_t bit_index = 7 - (x_off % 8);
                    const size_t byte_index = static_cast<size_t>(ch_data_off + x_off / 8);
                    if (byte_index < (FONT::glyph_bitmap_stride * FONT::glyph_bitmap_height)) {
                        const uint8_t a = (FONT::glyph_bitmap[byte_index] >> bit_index) & 1;
                        if (a) {
                            plot(xx, yy, col);
                        }
                    }
                }
                ch_data_off += FONT::glyph_bitmap_stride;
            }
        }
    }

    static consteval auto gen_a2al() {
        std::array<float, 16> a2alg{};
        for (size_t c = 0; c < 16; c++) {
            a2alg[c] = constixel::srgb_to_linear(static_cast<float>(c) * (1.0f / 15.0f));
        }
        return a2alg;
    }

    static constexpr std::array<float, 16> a2al = gen_a2al();

    template <typename FONT, text_rotation ROTATION>
    constexpr void draw_char_aa(int32_t x, int32_t y, const char_info<typename FONT::char_info_type> &ch, uint8_t col) {
        static_assert(FONT::mono == false, "Can't use a mono font to draw antialiased text.");
        int32_t ch_data_off = static_cast<int32_t>(ch.y) * static_cast<int32_t>(FONT::glyph_bitmap_stride) + static_cast<int32_t>(ch.x) / 2;
        float Rl = format.quant.linearpal.at((col & ((1UL << format.bits_per_pixel) - 1)) * 3 + 0);
        float Gl = format.quant.linearpal.at((col & ((1UL << format.bits_per_pixel) - 1)) * 3 + 1);
        float Bl = format.quant.linearpal.at((col & ((1UL << format.bits_per_pixel) - 1)) * 3 + 2);
        if constexpr (ROTATION == DEGREE_0) {
            x += ch.xoffset;
            y += ch.yoffset;
            const int32_t x2 = x + static_cast<int32_t>(ch.width);
            const int32_t y2 = y + static_cast<int32_t>(ch.height);
            for (int32_t yy = y; yy < y2; yy++) {
                for (int32_t xx = x; xx < x2; xx++) {
                    const int32_t x_off = (xx - x) + ch.x % 2;
                    const int32_t bit_index = (1 - (x_off % 2)) * 4;
                    const size_t byte_index = static_cast<size_t>(ch_data_off + x_off / 2);
                    if (byte_index < (FONT::glyph_bitmap_stride * FONT::glyph_bitmap_height)) {
                        const uint8_t a = (FONT::glyph_bitmap[byte_index] >> bit_index) & 0xF;
                        if (a != 0) {
                            if (a == 0xF) {
                                plot(xx, yy, col);
                            } else {
                                float Al = a2al[a];
                                compose(xx, yy, Al, Rl * Al, Gl * Al, Bl * Al);
                            }
                        }
                    }
                }
                ch_data_off += FONT::glyph_bitmap_stride;
            }
        } else if constexpr (ROTATION == DEGREE_180) {
            x -= ch.xoffset;
            y -= FONT::ascent;
            const int32_t x2 = x - static_cast<int32_t>(ch.width);
            const int32_t y2 = y + static_cast<int32_t>(ch.height);
            for (int32_t yy = y2 - 1; yy >= y; yy--) {
                for (int32_t xx = x2; xx < x; xx++) {
                    const int32_t x_off = (ch.width - (xx - x2) - 1) + ch.x % 2;
                    const int32_t bit_index = (1 - (x_off % 2)) * 4;
                    const size_t byte_index = static_cast<size_t>(ch_data_off + x_off / 2);
                    if (byte_index < (FONT::glyph_bitmap_stride * FONT::glyph_bitmap_height)) {
                        const uint8_t a = (FONT::glyph_bitmap[byte_index] >> bit_index) & 0xF;
                        if (a != 0) {
                            if (a == 0xF) {
                                plot(xx, yy, col);
                            } else {
                                float Al = a2al[a];
                                compose(xx, yy, Al, Rl * Al, Gl * Al, Bl * Al);
                            }
                        }
                    }
                }
                ch_data_off += FONT::glyph_bitmap_stride;
            }
        } else if constexpr (ROTATION == DEGREE_90) {
            x -= FONT::ascent;
            y += ch.xoffset;
            const int32_t x2 = x + static_cast<int32_t>(ch.height);
            const int32_t y2 = y + static_cast<int32_t>(ch.width);
            for (int32_t xx = x2 - 1; xx >= x; xx--) {
                for (int32_t yy = y; yy < y2; yy++) {
                    const int32_t x_off = (yy - y) + ch.x % 2;
                    const int32_t bit_index = (1 - (x_off % 2)) * 4;
                    const size_t byte_index = static_cast<size_t>(ch_data_off + x_off / 2);
                    if (byte_index < (FONT::glyph_bitmap_stride * FONT::glyph_bitmap_height)) {
                        const uint8_t a = (FONT::glyph_bitmap[byte_index] >> bit_index) & 0xF;
                        if (a != 0) {
                            if (a == 0xF) {
                                plot(xx, yy, col);
                            } else {
                                float Al = a2al[a];
                                compose(xx, yy, Al, Rl * Al, Gl * Al, Bl * Al);
                            }
                        }
                    }
                }
                ch_data_off += FONT::glyph_bitmap_stride;
            }
        } else if constexpr (ROTATION == DEGREE_270) {
            x += ch.yoffset;
            y -= ch.xoffset;
            const int32_t x2 = x + static_cast<int32_t>(ch.height);
            const int32_t y2 = y - static_cast<int32_t>(ch.width);
            for (int32_t xx = x; xx < x2; xx++) {
                for (int32_t yy = y2; yy < y; yy++) {
                    const int32_t x_off = (ch.width - (yy - y2) - 1) + ch.x % 2;
                    const int32_t bit_index = (1 - (x_off % 2)) * 4;
                    const size_t byte_index = static_cast<size_t>(ch_data_off + x_off / 2);
                    if (byte_index < (FONT::glyph_bitmap_stride * FONT::glyph_bitmap_height)) {
                        const uint8_t a = (FONT::glyph_bitmap[byte_index] >> bit_index) & 0xF;
                        if (a != 0) {
                            if (a == 0xF) {
                                plot(xx, yy, col);
                            } else {
                                float Al = a2al[a];
                                compose(xx, yy, Al, Rl * Al, Gl * Al, Bl * Al);
                            }
                        }
                    }
                }
                ch_data_off += FONT::glyph_bitmap_stride;
            }
        }
    }

    constexpr void span(int32_t x, int32_t w, int32_t y, uint8_t col) {
        if (x < 0) {
            w += x;
            x = 0;
        }
        if ((x >= static_cast<int32_t>(W)) || (x + w < 0) || (y >= static_cast<int32_t>(H)) || (y < 0)) {
            return;
        }
        if (x + w >= static_cast<int32_t>(W)) {
            w = static_cast<int32_t>(W) - x;
        }
        size_t _xl = static_cast<size_t>(x);
        size_t _xr = static_cast<size_t>(x + w);
        size_t _y = static_cast<size_t>(y);
        T<W, H, S, GR>::span(data, _xl, _xr, _y, col);
    }

    constexpr void fill_arc(int32_t x0, int32_t y0, int32_t r, uint8_t corners, int32_t delta, uint8_t col) {
        int32_t f = 1 - r;
        int32_t ddx = -2 * r;
        int32_t ddy = 1;
        int32_t x = r;
        int32_t y = 0;
        int32_t px = x;
        int32_t py = y;
        delta++;
        int32_t hl = corners & 4 ? 2 : 1;
        int32_t hr = corners & 8 ? 1 : 0;
        while (y < x) {
            if (f >= 0) {
                x--;
                ddx += 2;
                f += ddx;
            }
            ddy += 2;
            f += ddy;
            if (++y < (x + 1)) {
                if (corners & 1)
                    span(x0 - hr * x, hl * x + delta, y0 + y, col);
                if (corners & 2)
                    span(x0 - hr * x, hl * x + delta, y0 - y, col);
            }
            if (x != px) {
                if (corners & 1)
                    span(x0 - hr * py, hl * py + delta, y0 + px, col);
                if (corners & 2)
                    span(x0 - hr * py, hl * py + delta, y0 - px, col);
                px = x;
            }
            py = y;
        }
    }

    std::array<uint8_t, T<W, H, S, GR>::image_size> data{};
    T<W, H, S, GR> format{};
};

}  // namespace constixel

#endif  // CONSTIXEL_H_
