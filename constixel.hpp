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
#ifndef CONSTIXEL_HPP_
#define CONSTIXEL_HPP_

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef CONSTIXEL_ENABLE_COUT
#include <iostream>
#include <string>
#endif  // #ifdef CONSTIXEL_ENABLE_COUT

#if defined(__ARM_NEON)
#include <arm_neon.h>
#endif  // #if defined(__ARM_NEON)

#if defined(__AVX2__)
#include <immintrin.h>
#endif  // #if defined(__AVX2__)

// An empty member takes no storage. MSVC ignores the standard spelling.
#if defined(_MSC_VER) && !defined(__clang__)
#define CONSTIXEL_NO_UNIQUE_ADDRESS [[msvc::no_unique_address]]
#else  // #if defined(_MSC_VER) && !defined(__clang__)
#define CONSTIXEL_NO_UNIQUE_ADDRESS [[no_unique_address]]
#endif  // #if defined(_MSC_VER) && !defined(__clang__)

namespace constixel {

/// @cond DOXYGEN_EXCLUDE
namespace hidden {

[[nodiscard]] inline consteval double consteval_cos(double x, int32_t terms = 10) {
    x = x - 6.283185307179586 * static_cast<int32_t>(x / 6.283185307179586);  // wrap x to [0, 2π)
    double res = 1.0;
    double term = 1.0;
    const double x2 = x * x;
    for (int32_t i = 1; i < terms; ++i) {
        term *= -x2 / ((2 * i - 1) * (2 * i));
        res += term;
    }
    return res;
}

[[nodiscard]] inline consteval double consteval_sin(double x, int32_t terms = 10) {
    x = x - 6.283185307179586 * static_cast<int32_t>(x / 6.283185307179586);  // wrap x to [0, 2π)
    double res = x;
    double term = x;
    const double x2 = x * x;
    for (int32_t i = 1; i < terms; ++i) {
        term *= -x2 / ((2 * i) * (2 * i + 1));
        res += term;
    }
    return res;
}

inline constexpr float fast_sqrtf(const float x) {
    auto i = std::bit_cast<int32_t>(x);
    const int k = i & 0x00800000;
    float y = 0.0f;
    if (k != 0) {
        i = 0x5ed9d098 - (i >> 1);
        y = std::bit_cast<float>(i);
        y = 2.33139729f * y * ((-x * y * y) + 1.07492042f);
    } else {
        i = 0x5f19d352 - (i >> 1);
        y = std::bit_cast<float>(i);
        y = 0.82420468f * y * ((-x * y * y) + 2.14996147f);
    }
    const float c = x * y;
    const float r = ((y * -c) + 1.0f);
    y = ((0.5f * c * r) + c);
    return y;
}

[[nodiscard]] inline constexpr float fast_exp2(const float p) {
    const float offset = (p < 0) ? 1.0f : 0.0f;
    const float clipp = (p < -126) ? -126.0f : p;
    const float z = clipp - static_cast<float>(static_cast<int32_t>(clipp)) + offset;
    return std::bit_cast<float>(
        static_cast<uint32_t>((1 << 23) * (clipp + 121.2740575f + 27.7280233f / (4.84252568f - z) - 1.49012907f * z)));
}

[[nodiscard]] inline constexpr float fast_log2(const float x) {
    const auto xi = std::bit_cast<uint32_t>(x);
    const auto xf = std::bit_cast<float>((xi & 0x007FFFFF) | 0x3f000000);
    const float y = static_cast<float>(xi) * 1.1920928955078125e-7f;
    return y - 124.22551499f - 1.498030302f * xf - 1.72587999f / (0.3520887068f + xf);
}

[[nodiscard]] inline constexpr float fast_pow(const float x, const float p) {
    return fast_exp2(p * fast_log2(x));
}

[[nodiscard]] inline constexpr float linear_to_srgb(float c) {
    if (c <= 0.0031308f) {
        return 12.92f * c;
    }
    if (std::is_constant_evaluated()) {
        return 1.055f * fast_pow(c, 1.0f / 2.4f) - 0.055f;
    }
    return 1.055f * std::pow(c, 1.0f / 2.4f) - 0.055f;
}

[[nodiscard]] inline constexpr float srgb_to_linear(float s) {
    if (s <= 0.040449936f) {
        return s / 12.92f;
    }
    if (std::is_constant_evaluated()) {
        return fast_pow((s + 0.055f) / 1.055f, 2.4f);
    }
    return std::pow((s + 0.055f) / 1.055f, 2.4f);
}

#if defined(__ARM_NEON)
inline float32x4_t fast_log2_f32_neon(float32x4_t x) {
    uint32x4_t xi = vreinterpretq_u32_f32(x);
    float32x4_t y = vmulq_n_f32(vcvtq_f32_u32(xi), 1.1920928955078125e-7f);
    uint32x4_t mantissa = vorrq_u32(vandq_u32(xi, vdupq_n_u32(0x007FFFFF)), vdupq_n_u32(0x3f000000));
    float32x4_t xf = vreinterpretq_f32_u32(mantissa);
    float32x4_t a = vdupq_n_f32(124.22551499f);
    float32x4_t b = vdupq_n_f32(1.498030302f);
    float32x4_t c = vdupq_n_f32(1.72587999f);
    float32x4_t d = vdupq_n_f32(0.3520887068f);
    return vsubq_f32(vsubq_f32(y, vmlaq_f32(a, b, xf)), vdivq_f32(c, vaddq_f32(d, xf)));
}

inline float32x4_t fast_exp2_f32_neon(float32x4_t p) {
    float32x4_t one = vdupq_n_f32(1.0f);
    float32x4_t zero = vdupq_n_f32(0.0f);
    float32x4_t neg126 = vdupq_n_f32(-126.0f);
    float32x4_t offset = vbslq_f32(vcltq_f32(p, zero), one, zero);
    float32x4_t clipp = vmaxq_f32(p, neg126);
    int32x4_t ipart = vcvtq_s32_f32(clipp);
    float32x4_t fpart = vsubq_f32(clipp, vcvtq_f32_s32(ipart));
    float32x4_t z = vaddq_f32(fpart, offset);
    float32x4_t a = vdupq_n_f32(121.2740575f);
    float32x4_t b = vdupq_n_f32(27.7280233f);
    float32x4_t c = vdupq_n_f32(4.84252568f);
    float32x4_t d = vdupq_n_f32(1.49012907f);
    float32x4_t t = vaddq_f32(clipp, vsubq_f32(a, vmulq_f32(d, z)));
    t = vaddq_f32(t, vdivq_f32(b, vsubq_f32(c, z)));
    int32x4_t res = vcvtq_s32_f32(vmulq_n_f32(t, static_cast<float>(1 << 23)));
    return vreinterpretq_f32_s32(res);
}

inline float32x4_t pow_1_over_2_4_f32_neon(float32x4_t x) {
    return fast_exp2_f32_neon(vmulq_n_f32(fast_log2_f32_neon(x), 1.0f / 2.4f));
}

static inline float32x4_t linear_to_srgb_approx_neon(float32x4_t l) {
    const float32x4_t cutoff = vdupq_n_f32(0.0031308f);
    const float32x4_t scale = vdupq_n_f32(12.92f);
    const float32x4_t a = vdupq_n_f32(1.055f);
    const float32x4_t b = vdupq_n_f32(-0.055f);
    uint32x4_t mask = vcltq_f32(l, cutoff);
    float32x4_t srgb = vmulq_f32(l, scale);
    float32x4_t approx = vmlaq_f32(b, a, pow_1_over_2_4_f32_neon(l));
    return vbslq_f32(mask, srgb, approx);
}
#endif  // #if defined(__ARM_NEON)

#if defined(__AVX2__)

inline __m128 fast_exp2_ps(__m128 p) {
    const __m128 one = _mm_set1_ps(1.0f);
    const __m128 zero = _mm_setzero_ps();
    const __m128 neg126 = _mm_set1_ps(-126.0f);
    __m128 offset = _mm_and_ps(_mm_cmplt_ps(p, zero), one);
    __m128 clipp = _mm_max_ps(p, neg126);
    __m128i ipart = _mm_cvttps_epi32(clipp);
    __m128 fpart = _mm_sub_ps(clipp, _mm_cvtepi32_ps(ipart));
    __m128 z = _mm_add_ps(fpart, offset);
    const __m128 c1 = _mm_set1_ps(121.2740575f);
    const __m128 c2 = _mm_set1_ps(27.7280233f);
    const __m128 c3 = _mm_set1_ps(4.84252568f);
    const __m128 c4 = _mm_set1_ps(1.49012907f);
    const __m128 bias = _mm_set1_ps(1 << 23);
    __m128 t = _mm_add_ps(clipp, _mm_sub_ps(c1, _mm_mul_ps(c4, z)));
    t = _mm_add_ps(t, _mm_div_ps(c2, _mm_sub_ps(c3, z)));
    __m128i result = _mm_cvtps_epi32(_mm_mul_ps(t, bias));
    return _mm_castsi128_ps(result);
}

inline __m128 fast_log2_ps(__m128 x) {
    const __m128i xi = _mm_castps_si128(x);
    const __m128i mant_mask = _mm_set1_epi32(0x007FFFFF);
    const __m128i one_bits = _mm_set1_epi32(0x3f000000);
    const __m128 y = _mm_mul_ps(_mm_cvtepi32_ps(xi), _mm_set1_ps(1.1920928955078125e-7f));  // 1/(1<<23)
    __m128i mant_bits = _mm_or_si128(_mm_and_si128(xi, mant_mask), one_bits);
    __m128 xf = _mm_castsi128_ps(mant_bits);
    const __m128 c0 = _mm_set1_ps(124.22551499f);
    const __m128 c1 = _mm_set1_ps(1.498030302f);
    const __m128 c2 = _mm_set1_ps(1.72587999f);
    const __m128 c3 = _mm_set1_ps(0.3520887068f);
    return _mm_sub_ps(_mm_sub_ps(y, _mm_add_ps(c0, _mm_mul_ps(c1, xf))), _mm_div_ps(c2, _mm_add_ps(c3, xf)));
}

inline __m128 fast_pow_ps(__m128 x, __m128 p) {
    return fast_exp2_ps(_mm_mul_ps(p, fast_log2_ps(x)));
}

inline __m128 pow_1_over_2_4_ps(__m128 x) {
    return fast_pow_ps(x, _mm_set1_ps(1.0f / 2.4f));
}

inline __m128 linear_to_srgb_approx_sse(__m128 l) {
    const __m128 cutoff = _mm_set1_ps(0.0031308f);
    const __m128 scale = _mm_set1_ps(12.92f);
    const __m128 a = _mm_set1_ps(1.055f);
    const __m128 b = _mm_set1_ps(-0.055f);
    __m128 below = _mm_mul_ps(l, scale);
    __m128 powed = pow_1_over_2_4_ps(l);
    __m128 above = _mm_add_ps(_mm_mul_ps(a, powed), b);
    __m128 mask = _mm_cmplt_ps(l, cutoff);
    return _mm_or_ps(_mm_and_ps(mask, below), _mm_andnot_ps(mask, above));
}
#endif  // #if defined(__AVX2__)

struct oklch {
    double l = 0.0;
    double c = 0.0;
    double h = 0.0;
};

struct oklab {
    double l = 0.0;
    double a = 0.0;
    double b = 0.0;
};

struct srgb {
    double r = 0.0;
    double g = 0.0;
    double b = 0.0;
};

[[nodiscard]] inline consteval srgb oklab_to_srgb_consteval(const oklab &oklab) {
    const double l = oklab.l;
    const double a = oklab.a;
    const double b = oklab.b;

    const double l_ = l + 0.3963377774 * a + 0.2158037573 * b;
    const double m_ = l - 0.1055613458 * a - 0.0638541728 * b;
    const double s_ = l - 0.0894841775 * a - 1.2914855480 * b;

    const double r = 4.0767416621 * l_ - 3.3077115913 * m_ + 0.2309699292 * s_;
    const double g = -1.2684380046 * l_ + 2.6097574011 * m_ - 0.3413193965 * s_;
    const double bl = -0.0041960863 * l_ - 0.7034186168 * m_ + 1.7076147031 * s_;

    return {.r = static_cast<double>(linear_to_srgb(std::max(0.0f, std::min(1.0f, static_cast<float>(r))))),
            .g = static_cast<double>(linear_to_srgb(std::max(0.0f, std::min(1.0f, static_cast<float>(g))))),
            .b = static_cast<double>(linear_to_srgb(std::max(0.0f, std::min(1.0f, static_cast<float>(bl)))))};
}

[[nodiscard]] inline consteval oklab oklch_to_oklab_consteval(const oklch &oklch) {
    return {.l = oklch.l,
            .a = oklch.c * consteval_cos(oklch.h * 3.14159265358979323846 / 180.0),
            .b = oklch.c * consteval_sin(oklch.h * 3.14159265358979323846 / 180.0)};
}

inline constexpr const float epsilon_low = srgb_to_linear(0.5f / 255.0f);
inline constexpr const float epsilon_high = srgb_to_linear(254.5f / 255.0f);

inline consteval auto gen_a2al_4bit_consteval() {
    std::array<float, 16> a2al{};
    for (size_t c = 0; c < 16; c++) {
        a2al[c] = srgb_to_linear(static_cast<float>(c) * (1.0f / 15.0f));
    }
    return a2al;
}

inline constexpr const std::array<float, 16> a2al_4bit = gen_a2al_4bit_consteval();

inline consteval auto gen_a2al_8bit_consteval() {
    std::array<float, 256> a2al{};
    for (size_t c = 0; c < 256; c++) {
        a2al[c] = srgb_to_linear(static_cast<float>(c) * (1.0f / 255.0f));
    }
    return a2al;
}

inline constexpr const std::array<float, 256> a2al_8bit = gen_a2al_8bit_consteval();

template <size_t S>
class quantize {
    static constexpr size_t palette_size = S;

    std::array<float, palette_size * 3> linearpal{};
#if defined(__ARM_NEON)
    // Note: ordering after linearpal critical for overreads
    std::array<float, palette_size * 3> linearpal_neon{};
#endif  // #if defined(__ARM_NEON)
#if defined(__AVX2__)
    // Note: ordering after linearpal critical for overreads
    alignas(32) std::array<float, palette_size * 3> linearpal_avx2{};
#endif  // #if defined(__ARM_NEON)

    std::array<uint32_t, palette_size> pal{};

 public:
    ~quantize() = default;
    quantize(const quantize &&) = delete;
    quantize(const quantize &) = delete;
    quantize &operator=(const quantize &) = delete;
    quantize &operator=(quantize &&) = delete;

    explicit consteval quantize(const std::array<uint32_t, palette_size> &palette) : pal(palette) {
        for (size_t i = 0; i < pal.size(); i++) {
            pal[i] = (pal[i] & 0xFF000000) | (pal[i] & 0x0000FF00) | ((pal[i] >> 16) & 0x000000FF) |
                     ((pal[i] << 16) & 0x00FF0000);
        }
        for (size_t i = 0; i < pal.size(); i++) {
            linearpal.at(i * 3 + 0) =
                hidden::srgb_to_linear(static_cast<float>((pal[i] >> 0) & 0xFF) * (1.0f / 255.0f));
            linearpal.at(i * 3 + 1) =
                hidden::srgb_to_linear(static_cast<float>((pal[i] >> 8) & 0xFF) * (1.0f / 255.0f));
            linearpal.at(i * 3 + 2) =
                hidden::srgb_to_linear(static_cast<float>((pal[i] >> 16) & 0xFF) * (1.0f / 255.0f));
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
        int32_t best = 0;
        int32_t bestd = 1UL << 30;
        for (size_t i = 0; i < pal.size(); ++i) {
            const int32_t dr = r - static_cast<int32_t>((pal[i] >> 0) & 0xFF);
            const int32_t dg = g - static_cast<int32_t>((pal[i] >> 8) & 0xFF);
            const int32_t db = b - static_cast<int32_t>((pal[i] >> 16) & 0xFF);
            const int32_t d = dr * dr + dg * dg + db * db;
            if (d < bestd) {
                bestd = d;
                best = static_cast<int32_t>(i);
            }
        }
        return static_cast<uint8_t>(best);
    }

    [[nodiscard]] constexpr auto &palette() const {
        return pal;
    }

    [[nodiscard]] constexpr auto &linear_palette() const {
        return linearpal;
    }

    [[nodiscard]] constexpr uint8_t nearest_linear(float r, float g, float b) const {
        const float epsilon = 0.0001f;
#if defined(__ARM_NEON)
        if (!std::is_constant_evaluated() && pal.size() >= 4) {
            const float32x4_t vR = vdupq_n_f32(r);
            const float32x4_t vG = vdupq_n_f32(g);
            const float32x4_t vB = vdupq_n_f32(b);

            float best = std::numeric_limits<float>::infinity();
            std::size_t bestIdx = 0;

            for (size_t i = 0; i < pal.size(); i += 4) {
                const float32x4_t dr = vsubq_f32(vld1q_f32(&linearpal_neon[i * 3]), vR);
                const float32x4_t dg = vsubq_f32(vld1q_f32(&linearpal_neon[i * 3 + 4]), vG);
                const float32x4_t db = vsubq_f32(vld1q_f32(&linearpal_neon[i * 3 + 8]), vB);

#if defined(__aarch64__) && defined(__ARM_FEATURE_FMA)
                const float32x4_t dist = vfmaq_f32(vfmaq_f32(vmulq_f32(dr, dr), dg, dg), db, db);
#else
                const float32x4_t dist = vaddq_f32(vaddq_f32(vmulq_f32(dr, dr), vmulq_f32(dg, dg)), vmulq_f32(db, db));
#endif

                const float32x2_t lo = vget_low_f32(dist);
                const float32x2_t hi = vget_high_f32(dist);
                const float d0 = vget_lane_f32(lo, 0);
                if (d0 < best) {
                    best = d0;
                    bestIdx = i;
                    if (d0 <= epsilon) {
                        break;
                    }
                }
                const float d1 = vget_lane_f32(lo, 1);
                if (d1 < best) {
                    best = d1;
                    bestIdx = i + 1;
                    if (d1 <= epsilon) {
                        break;
                    }
                }
                const float d2 = vget_lane_f32(hi, 0);
                if (d2 < best) {
                    best = d2;
                    bestIdx = i + 2;
                    if (d2 <= epsilon) {
                        break;
                    }
                }
                const float d3 = vget_lane_f32(hi, 1);
                if (d3 < best) {
                    best = d3;
                    bestIdx = i + 3;
                    if (d3 <= epsilon) {
                        break;
                    }
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
                                    _mm256_fmadd_ps(_mm256_sub_ps(pg, vG), _mm256_sub_ps(pg, vG),
                                                    _mm256_mul_ps(_mm256_sub_ps(pb, vB), _mm256_sub_ps(pb, vB))));

                alignas(32) float d[8];
                _mm256_store_ps(d, dist);

                for (size_t lane = 0; lane < 8; lane++) {
                    if (d[lane] < best) {
                        best = d[lane];
                        bestIdx = static_cast<uint8_t>(i + lane);
                        if (d[lane] <= epsilon) {
                            return static_cast<uint8_t>(bestIdx);
                        }
                    }
                }
            }

            return static_cast<uint8_t>(bestIdx);
        }
#endif  // #if defined(__AVX2__)
        size_t best = 0;
        float bestd = 100.0f;
        for (size_t i = 0; i < pal.size(); ++i) {
            const float dr = r - linearpal.at(i * 3 + 0);
            const float dg = g - linearpal.at(i * 3 + 1);
            const float db = b - linearpal.at(i * 3 + 2);
            const float d = dr * dr + dg * dg + db * db;
            if (d < bestd) {
                bestd = d;
                best = i;
                if (d <= epsilon) {
                    break;
                }
            }
        }
        return static_cast<uint8_t>(best);
    }
};

}  // namespace hidden
/// @endcond // DOXYGEN_EXCLUDE

/// @cond PRIVATE_CLASS
template <size_t N, typename T>
class hextree {
    static constexpr T bitslices = ((sizeof(T) * 8) / 4) - 1;

    static constexpr size_t child_nodes_n = 1UL << 4;
    static constexpr size_t child_nodes_n_mask = child_nodes_n - 1;

    struct node {
        std::array<T, child_nodes_n> child{};
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

    ~hextree() = default;
    hextree(const hextree &&) = delete;
    hextree(const hextree &) = delete;
    hextree &operator=(const hextree &) = delete;
    hextree &operator=(hextree &&) = delete;

    explicit consteval hextree() = default;
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
                    vnodes[idx].child[nib] = next;
                    vnodes.emplace_back();
                }
                idx = next;
            }
            vnodes[idx].child[key & child_nodes_n_mask] = val;
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
/// @endcond // PRIVATE_CLASS

/// @cond PRIVATE_CLASS
template <typename T>
struct char_info {
    T x = 0;
    T y = 0;
    T width = 0;
    T height = 0;
    T xadvance = 0;
    T xoffset = 0;
    T yoffset = 0;
};
/// @endcond // PRIVATE_CLASS

/**
 * @brief basic rectangle structure
 * @tparam T coordinate number type
 */
template <typename T>
struct rect {
    T x = 0;  ///< x coordinate
    T y = 0;  ///< y coordinate
    T w = 0;  ///< width
    T h = 0;  ///< height

    /// @brief intersects one rect with another
    /// @param other the other rectangle
    constexpr rect operator&(const rect &other) const {
        T x1 = std::max(x, other.x);
        T y1 = std::max(y, other.y);
        T x2 = std::min(x + w, other.x + other.w);
        T y2 = std::min(y + h, other.y + other.h);

        T nw = x2 - x1;
        T nh = y2 - y1;

        if (nw <= 0 || nh <= 0) {
            return {T{0}, T{0}, T{0}, T{0}};
        }

        return {x1, y1, nw, nh};
    }

    /// @brief intersects one rect with another
    /// @param other the other rectangle
    constexpr rect &operator&=(const rect &other) {
        T x1 = std::max(x, other.x);
        T y1 = std::max(y, other.y);
        T x2 = std::min(x + w, other.x + other.w);
        T y2 = std::min(y + h, other.y + other.h);

        T nw = x2 - x1;
        T nh = y2 - y1;

        if (nw <= 0 || nh <= 0) {
            *this = {T{0}, T{0}, T{0}, T{0}};
            return *this;
        }

        x = x1;
        y = y1;
        w = nw;
        h = nh;

        return *this;
    }
};

/// @cond PRIVATE_CLASS
namespace hidden {

// Pixel buffer as seen by the format classes: stride is the byte distance between lines.
template <typename BYTE>
struct raster {
    BYTE *data = nullptr;
    size_t width = 0;
    size_t height = 0;
    size_t stride = 0;

    constexpr operator raster<const BYTE>() const
        requires(!std::is_const_v<BYTE>)
    {
        return {data, width, height, stride};
    }
};

}  // namespace hidden
/// @endcond

/// @cond PRIVATE_CLASS
class format {
 public:
    using raster = hidden::raster<uint8_t>;
    using const_raster = hidden::raster<const uint8_t>;

    [[nodiscard]] static constexpr uint32_t adler32(const uint8_t *data, std::size_t len, uint32_t adler32_sum) {
        uint32_t adler32_s1 = adler32_sum & 0xFFFF;
        uint32_t adler32_s2 = adler32_sum >> 16;
        for (size_t c = 0; c < len; c++) {
            adler32_s1 = (adler32_s1 + data[c]) % 65521;
            adler32_s2 = (adler32_s2 + adler32_s1) % 65521;
        }
        return ((adler32_s2 << 16) | adler32_s1);
    }

    template <typename F>
    static constexpr void png_write_be(F &&char_out, uint32_t value) {
        std::forward<F>(char_out)(static_cast<char>((value >> 24) & 0xFF));
        std::forward<F>(char_out)(static_cast<char>((value >> 16) & 0xFF));
        std::forward<F>(char_out)(static_cast<char>((value >> 8) & 0xFF));
        std::forward<F>(char_out)(static_cast<char>((value >> 0) & 0xFF));
    }

    template <typename F, typename A>
    static constexpr void png_write_crc32(F &&char_out, const A &array, size_t bytes) {
        size_t idx = 0;
        uint32_t crc = 0xFFFFFFFF;
        for (size_t c = 0; c < bytes; c++) {
            auto d = static_cast<uint8_t>(array.at(idx++));
            crc ^= d;
            for (int i = 0; i < 8; ++i) {
                crc = (crc & 1) ? (crc >> 1) ^ 0xEDB88320u : crc >> 1;
            }
        }
        png_write_be(std::forward<F>(char_out), crc ^ 0xFFFFFFFF);
    }

    template <typename F, typename A>
    static constexpr void png_write_array(F &&char_out, const A &array, size_t bytes) {
        for (size_t c = 0; c < bytes; c++) {
            std::forward<F>(char_out)(static_cast<char>(array.at(c)));
        }
    }

    template <typename F>
    static constexpr void png_marker(F &&char_out) {
        std::forward<F>(char_out)(static_cast<char>(0x89));
        std::forward<F>(char_out)(0x50);
        std::forward<F>(char_out)(0x4E);
        std::forward<F>(char_out)(0x47);
        std::forward<F>(char_out)(0x0D);
        std::forward<F>(char_out)(0x0A);
        std::forward<F>(char_out)(0x1A);
        std::forward<F>(char_out)(0x0A);
    }

    template <typename F>
    static constexpr void png_header(F &&char_out, size_t w, size_t h, size_t depth) {
        const size_t chunkLength = 17;
        std::array<char, chunkLength> header{};
        size_t i = 0;
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
        if (depth <= 8) {
            header.at(i++) = static_cast<char>(depth);
            header.at(i++) = 3;
        } else if (depth == 24) {
            header.at(i++) = 8;
            header.at(i++) = 2;
        } else {
            header.at(i++) = 8;
            header.at(i++) = 6;
        }
        header.at(i++) = 0;
        header.at(i++) = 0;
        header.at(i++) = 0;
        png_write_be(std::forward<F>(char_out), static_cast<uint32_t>(i - 4));
        png_write_array(std::forward<F>(char_out), header, i);
        png_write_crc32(std::forward<F>(char_out), header, i);
    }

    template <typename F, typename P>
    static constexpr void png_palette(F &&char_out, const P &palette) {
        std::array<char, 256 * 3 + 4> header{};
        size_t i = 0;
        header.at(i++) = 'P';
        header.at(i++) = 'L';
        header.at(i++) = 'T';
        header.at(i++) = 'E';
        for (size_t c = 0; c < palette.size(); c++) {
            header.at(i++) = static_cast<char>((palette[c] >> 0) & 0xFF);
            header.at(i++) = static_cast<char>((palette[c] >> 8) & 0xFF);
            header.at(i++) = static_cast<char>((palette[c] >> 16) & 0xFF);
        }
        png_write_be(std::forward<F>(char_out), static_cast<uint32_t>(i - 4));
        png_write_array(std::forward<F>(char_out), header, i);
        png_write_crc32(std::forward<F>(char_out), header, i);
    }

    template <typename F>
    static constexpr void png_end(F &&char_out) {
        std::array<char, 4> header{};
        size_t i = 0;
        header.at(i++) = 'I';
        header.at(i++) = 'E';
        header.at(i++) = 'N';
        header.at(i++) = 'D';
        png_write_be(std::forward<F>(char_out), static_cast<uint32_t>(i - 4));
        png_write_array(std::forward<F>(char_out), header, i);
        png_write_crc32(std::forward<F>(char_out), header, i);
    }

    template <typename F>
    static constexpr void png_idat_zlib_header(F &&char_out) {
        std::array<char, 6> header{};
        size_t i = 0;
        header.at(i++) = 'I';
        header.at(i++) = 'D';
        header.at(i++) = 'A';
        header.at(i++) = 'T';
        header.at(i++) = 0x78;
        header.at(i++) = 0x01;
        png_write_be(std::forward<F>(char_out), static_cast<uint32_t>(i - 4));
        png_write_array(std::forward<F>(char_out), header, i);
        png_write_crc32(std::forward<F>(char_out), header, i);
    }

    template <typename F>
    [[nodiscard]] static constexpr uint32_t png_idat_zlib_stream(F &&char_out, const uint8_t *line, size_t bytes,
                                                                 uint32_t adler32_sum, bool last_line) {
        const auto max_data_use = size_t{1024};
        const auto extra_data = size_t{10};
        const size_t max_stack_use = max_data_use + extra_data;
        std::array<uint8_t, max_stack_use> header{};
        size_t filter_first = 1;
        while (bytes > 0) {
            size_t i = 0;
            header.at(i++) = 'I';
            header.at(i++) = 'D';
            header.at(i++) = 'A';
            header.at(i++) = 'T';

            const size_t bytes_to_copy = std::min(max_data_use, bytes);

            header.at(i++) = ((bytes_to_copy == bytes) && last_line) ? 0x01 : 0x00;

            header.at(i++) = (((bytes_to_copy + filter_first) >> 0) & 0xFF);
            header.at(i++) = (((bytes_to_copy + filter_first) >> 8) & 0xFF);
            header.at(i++) = ((((bytes_to_copy + filter_first) ^ 0xffff) >> 0) & 0xFF);
            header.at(i++) = ((((bytes_to_copy + filter_first) ^ 0xffff) >> 8) & 0xFF);

            const size_t adlersum32_start_pos = i;
            if (filter_first != 0) {
                filter_first = 0;
                header.at(i++) = 0;
            }
            for (size_t d = 0; d < bytes_to_copy; d++) {
                header.at(i++) = *line++;
            }
            adler32_sum = adler32(&header[adlersum32_start_pos], i - adlersum32_start_pos, adler32_sum);

            png_write_be(std::forward<F>(char_out), static_cast<uint32_t>(i - 4));
            png_write_array(std::forward<F>(char_out), header, i);
            png_write_crc32(std::forward<F>(char_out), header, i);

            bytes -= bytes_to_copy;
        }

        return adler32_sum;
    }

    template <typename F>
    static constexpr void png_idat_zlib_trailer(F &&char_out, uint32_t adler32_sum) {
        std::array<char, 8> header{};
        size_t i = 0;
        header.at(i++) = 'I';
        header.at(i++) = 'D';
        header.at(i++) = 'A';
        header.at(i++) = 'T';
        header.at(i++) = static_cast<char>((adler32_sum >> 24) & 0xFF);
        header.at(i++) = static_cast<char>((adler32_sum >> 16) & 0xFF);
        header.at(i++) = static_cast<char>((adler32_sum >> 8) & 0xFF);
        header.at(i++) = static_cast<char>((adler32_sum >> 0) & 0xFF);
        png_write_be(std::forward<F>(char_out), static_cast<uint32_t>(i - 4));
        png_write_array(std::forward<F>(char_out), header, i);
        png_write_crc32(std::forward<F>(char_out), header, i);
    }

    // Emit DCS sixel introducer `\eP;1q` — P2=1 marks unset pixels
    // as transparent, which avoids Microsoft Terminal #17887's
    // arbitrary-color cell-padding artifact. P2=1 is visually
    // equivalent to the default for fully-painted images, so it's
    // safe to use universally (compile-time + runtime palettes).
    template <typename F>
    static constexpr void sixel_header(F &&char_out) {
        std::forward<F>(char_out)(0x1b);
        std::forward<F>(char_out)('P');
        std::forward<F>(char_out)(';');
        std::forward<F>(char_out)('1');
        std::forward<F>(char_out)('q');
    }

    template <typename F>
    static constexpr void sixel_raster_attributes(F &&char_out, size_t w, size_t h, size_t s) {
        std::forward<F>(char_out)('\"');
        sixel_number(std::forward<F>(char_out), 1);
        std::forward<F>(char_out)(';');
        sixel_number(std::forward<F>(char_out), 1);
        std::forward<F>(char_out)(';');
        sixel_number(std::forward<F>(char_out), static_cast<uint16_t>(w * s));
        std::forward<F>(char_out)(';');
        sixel_number(std::forward<F>(char_out), static_cast<uint16_t>(h * s));
    }

    template <typename F>
    static constexpr void sixel_number(F &&char_out, uint16_t u) {
        if (u < 10) {
            std::forward<F>(char_out)(static_cast<char>('0' + u));
        } else if (u < 100) {
            std::forward<F>(char_out)(static_cast<char>('0' + (((u / 10) % 10))));
            std::forward<F>(char_out)(static_cast<char>('0' + (u % 10)));
        } else if (u < 1000) {
            std::forward<F>(char_out)(static_cast<char>('0' + ((u / 100) % 10)));
            std::forward<F>(char_out)(static_cast<char>('0' + ((u / 10) % 10)));
            std::forward<F>(char_out)(static_cast<char>('0' + (u % 10)));
        } else if (u < 10000) {
            std::forward<F>(char_out)(static_cast<char>('0' + ((u / 1000) % 10)));
            std::forward<F>(char_out)(static_cast<char>('0' + ((u / 100) % 10)));
            std::forward<F>(char_out)(static_cast<char>('0' + ((u / 10) % 10)));
            std::forward<F>(char_out)(static_cast<char>('0' + (u % 10)));
        } else {
            std::forward<F>(char_out)(static_cast<char>('0' + ((u / 10000) % 10)));
            std::forward<F>(char_out)(static_cast<char>('0' + ((u / 1000) % 10)));
            std::forward<F>(char_out)(static_cast<char>('0' + ((u / 100) % 10)));
            std::forward<F>(char_out)(static_cast<char>('0' + ((u / 10) % 10)));
            std::forward<F>(char_out)(static_cast<char>('0' + (u % 10)));
        }
    }

    template <typename F>
    static constexpr void sixel_color(F &&char_out, uint16_t i, uint32_t col) {
        std::forward<F>(char_out)('#');
        sixel_number(char_out, i);
        std::forward<F>(char_out)(';');
        std::forward<F>(char_out)('2');
        std::forward<F>(char_out)(';');
        for (size_t c = 0; c < 3; c++) {
            sixel_number(std::forward<F>(char_out), static_cast<uint16_t>((((col >> (8 * c)) & 0xFF) * 100) / 255));
            if (c < 2) {
                std::forward<F>(char_out)(';');
            }
        }
    }

    template <typename F>
    static constexpr void sixel_end(F &&char_out) {
        std::forward<F>(char_out)(0x1b);
        std::forward<F>(char_out)('\\');
    }

    // 24/32-bit sixel output caches the palette indices of the six lines of a band, column by column. A cache of
    // the image width quantizes every pixel once. A smaller one is a window which follows the column that is read
    // and quantizes once per color of a band.
    [[nodiscard]] static constexpr size_t sixel_cache_size(size_t columns) {
        return columns * 6;
    }

    // Columns of the cache a surface uses when its owner supplies none.
    static constexpr size_t sixel_cache_columns = 1024;

    // First cached column for a window of the given size which has to hold column x.
    [[nodiscard]] static constexpr size_t sixel_cache_start(size_t x, size_t width, size_t columns) {
        return std::min(x, width > columns ? width - columns : size_t{0});
    }

    template <typename PBT, size_t PBS>
    struct palette_bitset {
        constexpr void mark(PBT col) {
            const PBT idx = (col >> 5) & set_idx_mask;
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

    template <size_t PBS, typename P, typename F>
    static constexpr void png_image(const const_raster &pix, size_t line_bytes, const P &palette, F &&char_out) {
        png_marker(std::forward<F>(char_out));
        png_header(std::forward<F>(char_out), pix.width, pix.height, PBS);
        if constexpr (PBS <= 8) {
            png_palette(std::forward<F>(char_out), palette);
        }
        png_idat_zlib_header(std::forward<F>(char_out));
        uint32_t adler32_sum = 1;
        for (size_t y = 0; y < pix.height; y++) {
            adler32_sum = png_idat_zlib_stream(std::forward<F>(char_out), pix.data + y * pix.stride, line_bytes,
                                               adler32_sum, y == pix.height - 1);
        }
        png_idat_zlib_trailer(std::forward<F>(char_out), adler32_sum);
        png_end(std::forward<F>(char_out));
    }

    // Single sixel emission body — shared between compile-time call
    // sites and runtime-sized images. W/H/S are runtime parameters
    // so the same code drives both paths.
    template <typename PBT, size_t PBS, typename P, typename F, typename C, typename D>
    static constexpr void sixel_image(const uint8_t *data, const P &palette, F &&char_out, const rect<int32_t> &_r,
                                      size_t W, size_t H, size_t S, const C &collect6, const D &set6) {
        sixel_header(std::forward<F>(char_out));
        sixel_raster_attributes(std::forward<F>(char_out), W, H, S);
        for (size_t c = 0; c < palette.size(); c++) {
            sixel_color(std::forward<F>(char_out), static_cast<uint16_t>(c), palette[c]);
        }
        const auto r =
            rect<int32_t>{.x = _r.x * static_cast<int32_t>(S),
                          .y = _r.y * static_cast<int32_t>(S),
                          .w = _r.w * static_cast<int32_t>(S),
                          .h = _r.h * static_cast<int32_t>(S)} &
            rect<int32_t>{.x = 0, .y = 0, .w = static_cast<int32_t>(W * S), .h = static_cast<int32_t>(H * S)};
        std::array<PBT, PBS> stack{};
        palette_bitset<PBT, PBS> pset{};
        for (int32_t y = r.y; y < (r.y + r.h); y += 6) {
            pset.clear();
            set6(data, static_cast<size_t>(r.x), static_cast<size_t>(r.w), static_cast<size_t>(y), pset);
            const size_t stack_count = pset.genstack(stack);
            for (size_t s = 0; s < stack_count; s++) {
                const PBT col = stack[s];
                // Carriage-return-within-band ($) ONLY between layers,
                // never before the first. The previous `col != 0` test
                // assumed the built-in palette always has color 0 at
                // stack[0]; that's false for runtime palettes whose
                // index 0 is a reserved sentinel never marked into the
                // pset, so the wrong test emitted a stray $ at every
                // band start and Windows Terminal rendered black
                // artifacts for it.
                if (s != 0) {
                    std::forward<F>(char_out)('$');
                }
                std::forward<F>(char_out)('#');
                sixel_number(std::forward<F>(char_out), static_cast<uint16_t>(col));
                for (int32_t x = r.x; x < (r.x + r.w); x++) {
                    PBT bits6 = collect6(data, static_cast<size_t>(x), col, static_cast<size_t>(y));
                    uint16_t repeat_count = 0;
                    for (int32_t xr = x + 1; xr < (std::min(x + int32_t{255}, static_cast<int32_t>(W * S))); xr++) {
                        if (bits6 == collect6(data, static_cast<size_t>(xr), col, static_cast<size_t>(y))) {
                            repeat_count++;
                            continue;
                        }
                        break;
                    }
                    if (repeat_count > uint16_t{3}) {
                        std::forward<F>(char_out)('!');
                        sixel_number(std::forward<F>(char_out), static_cast<uint16_t>(repeat_count + 1));
                        x += repeat_count;
                    }
                    std::forward<F>(char_out)(static_cast<char>('?' + bits6));
                }
            }
            std::forward<F>(char_out)('-');
        }
        sixel_end(std::forward<F>(char_out));
    }
};
/// @endcond

/// @cond DOXYGEN_EXCLUDE
/// Runtime palette (<=255 caller entries, BGR-packed). Index 0 is reserved
/// as a sentinel because the DEC sixel spec treats color 0 as transparent
/// when the introducer's P2=1 — pixels mapped to #0 would render as
/// terminal background. Caller entries are stored at [1..count_].
struct runtime_palette {
    std::array<uint32_t, 256> pal{};
    size_t count_{0};
    constexpr runtime_palette() = default;
    explicit runtime_palette(std::span<const std::array<uint8_t, 3>> rgb) {
        const size_t n = std::min(rgb.size(), size_t{255});
        for (size_t i = 0; i < n; ++i) {
            pal[i + 1] = uint32_t(rgb[i][0]) | (uint32_t(rgb[i][1]) << 8) | (uint32_t(rgb[i][2]) << 16);
        }
        count_ = n + 1;
    }
    [[nodiscard]] constexpr size_t size() const noexcept {
        return count_;
    }
    [[nodiscard]] constexpr uint32_t at(size_t i) const noexcept {
        return pal[i];
    }
    [[nodiscard]] constexpr uint8_t nearest(int32_t r, int32_t g, int32_t b) const {
        int32_t best_d = std::numeric_limits<int32_t>::max();
        uint8_t best_i = 1;
        for (size_t i = 1; i < count_; ++i) {
            const int32_t dr = r - int32_t((pal[i] >> 0) & 0xFF);
            const int32_t dg = g - int32_t((pal[i] >> 8) & 0xFF);
            const int32_t db = b - int32_t((pal[i] >> 16) & 0xFF);
            const int32_t d = dr * dr + dg * dg + db * db;
            if (d < best_d) {
                best_d = d;
                best_i = uint8_t(i);
                if (d == 0)
                    break;
            }
        }
        return best_i;
    }
};
/// @endcond

/// @cond PRIVATE_CLASS
// Pixel operations of the palette formats with BITS (1, 2, 4 or 8) per pixel, independent of the image size.
template <size_t BITS, bool GRAYSCALE>
class pixel_packed : public format {
 public:
    /// @cond DOXYGEN_EXCLUDE
    static constexpr size_t sixel_bitset_size = BITS == 8 ? 256 : 32;
    static constexpr size_t bits_per_pixel = BITS;
    static constexpr size_t pixels_per_byte = 8 / BITS;
    static constexpr size_t color_mask = (1UL << BITS) - 1;
    [[nodiscard]] static constexpr size_t bytes_per_line(size_t w) {
        return (w * BITS + 7) / 8;
    }
    [[nodiscard]] static constexpr size_t image_size(size_t w, size_t h) {
        return h * bytes_per_line(w);
    }

    static consteval auto gen_palette_consteval() {
        if constexpr (BITS == 1) {
            return std::array<uint32_t, 2>({{0x00000000, 0x00ffffff}});
        } else if constexpr (BITS == 2) {
            if (GRAYSCALE) {
                return std::array<uint32_t, (1UL << bits_per_pixel)>({{0x000000, 0x444444, 0x888888, 0xffffff}});
            }
            return std::array<uint32_t, (1UL << bits_per_pixel)>({{0x000000, 0xffffff, 0xff0000, 0x0077ff}});
        } else if constexpr (BITS == 4) {
            if (GRAYSCALE) {
                return std::array<uint32_t, (1UL << bits_per_pixel)>(
                    {{0x000000, 0x111111, 0x222222, 0x333333, 0x444444, 0x555555, 0x666666, 0x777777, 0x888888,
                      0x999999, 0xaaaaaa, 0xbbbbbb, 0xcccccc, 0xdddddd, 0xeeeeee, 0xffffff}});
            }
            return std::array<uint32_t, (1UL << bits_per_pixel)>(
                {{0x000000, 0xffffff, 0xff0000, 0x00ff00, 0x0000ff, 0xffff00, 0x00ffff, 0xff00ff, 0x333333, 0x666666,
                  0x999999, 0xcccccc, 0x7f0000, 0x007f00, 0x00007f, 0x7f7f00}});
        } else {
            std::array<uint32_t, (1UL << bits_per_pixel)> pal{};
            if (GRAYSCALE) {
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
                    const uint32_t y = (0xff * static_cast<uint32_t>(c)) / 15;
                    pal[0x10 + c] = (y << 16) | (y << 8) | (y << 0);
                }
                for (size_t c = 0; c < 8; c++) {
                    const uint32_t y = (0xff * static_cast<uint32_t>(c)) / 7;
                    const uint32_t x = (0xff * (static_cast<uint32_t>(c) + 1)) / 8;
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
                    const hidden::oklab lft{.l = static_cast<double>(c) / 7 - 0.2, .a = 0.2, .b = 0.0};
                    const hidden::oklab rgh{.l = static_cast<double>(c) / 7 - 0.2, .a = 0.2, .b = 337.5};
                    for (size_t d = 0; d < 16; d++) {
                        auto res = hidden::oklab_to_srgb_consteval(hidden::oklch_to_oklab_consteval(hidden::oklch{
                            .l = std::lerp(lft.l, rgh.l, static_cast<double>(d) / 15.0),
                            .c = std::lerp(lft.a, rgh.a, static_cast<double>(d) / 15.0),
                            .h = std::lerp(lft.b, rgh.b, static_cast<double>(d) / 15.0),
                        }));
                        pal[0x80 + c * 16 + d] =
                            (static_cast<uint32_t>(std::max(0.0, std::min(1.0, res.r)) * 255.0) << 16) |
                            (static_cast<uint32_t>(std::max(0.0, std::min(1.0, res.g)) * 255.0) << 8) |
                            (static_cast<uint32_t>(std::max(0.0, std::min(1.0, res.b)) * 255.0) << 0);
                    }
                }
            }
            return pal;
        }
    }
    static constexpr const auto quant = hidden::quantize<1UL << BITS>(gen_palette_consteval());

    // Reverses the order of the pixels in a byte.
    static constexpr uint8_t reverse(uint8_t b) {
        if constexpr (BITS <= 4) {
            b = static_cast<uint8_t>((b & uint8_t{0xF0}) >> 4 | (b & uint8_t{0x0F}) << 4);
        }
        if constexpr (BITS <= 2) {
            b = static_cast<uint8_t>((b & uint8_t{0xCC}) >> 2 | (b & uint8_t{0x33}) << 2);
        }
        if constexpr (BITS == 1) {
            b = static_cast<uint8_t>((b & uint8_t{0xAA}) >> 1 | (b & uint8_t{0x55}) << 1);
        }
        return b;
    }

    static constexpr void transpose8x8(std::array<uint8_t, 8> &a) {
        // clang-format off
        uint32_t x = (static_cast<uint32_t>(a[0]) << 24) |
                     (static_cast<uint32_t>(a[1]) << 16) |
                     (static_cast<uint32_t>(a[2]) <<  8) |
                     (static_cast<uint32_t>(a[3]));
        uint32_t y = (static_cast<uint32_t>(a[4]) << 24) |
                     (static_cast<uint32_t>(a[5]) << 16) |
                     (static_cast<uint32_t>(a[6]) <<  8) |
                     (static_cast<uint32_t>(a[7]));

        uint32_t t = (x ^ (x >>  7)) & 0x00AA00AA; x = x ^ t ^ (t <<  7);
                 t = (y ^ (y >>  7)) & 0x00AA00AA; y = y ^ t ^ (t <<  7);
                 t = (x ^ (x >> 14)) & 0x0000CCCC; x = x ^ t ^ (t << 14);
                 t = (y ^ (y >> 14)) & 0x0000CCCC; y = y ^ t ^ (t << 14);

        t = ((y >> 4) & 0x0F0F0F0F) | (x & 0xF0F0F0F0);
        y = ((x << 4) & 0xF0F0F0F0) | (y & 0x0F0F0F0F);

        a[0] = static_cast<uint8_t>(t >> 24);
        a[1] = static_cast<uint8_t>(t >> 16);
        a[2] = static_cast<uint8_t>(t >>  8);
        a[3] = static_cast<uint8_t>(t);
        a[4] = static_cast<uint8_t>(y >> 24);
        a[5] = static_cast<uint8_t>(y >> 16);
        a[6] = static_cast<uint8_t>(y >>  8);
        a[7] = static_cast<uint8_t>(y);
        // clang-format on
    }

    template <bool FLIP_H = false, bool FLIP_V = false>
    static constexpr void transpose(const const_raster &pix, uint8_t *dst) {
        if constexpr (BITS == 1) {
            std::array<uint8_t, 8> tmp{};
            const size_t src_bytes = ((pix.width + 7) / 8);
            const size_t dst_stride = ((pix.height + 7) / 8);
            for (size_t y = 0; y < dst_stride; y++) {
                for (size_t x = 0; x < src_bytes; x++) {
                    size_t xl = std::min(size_t{8}, pix.height - (y * 8));
                    for (size_t c = 0; c < xl; c++) {
                        if constexpr (FLIP_V) {
                            tmp[c] = pix.data[(pix.height - 1 - (y * 8 + c)) * pix.stride + x];
                        } else {
                            tmp[c] = pix.data[(y * 8 + c) * pix.stride + x];
                        }
                    }
                    for (size_t c = xl; c < 8; c++) {
                        tmp[c] = 0;
                    }
                    transpose8x8(tmp);
                    size_t yl = std::min(size_t{8}, pix.width - (x * 8));
                    for (size_t c = 0; c < yl; c++) {
                        if constexpr (FLIP_H) {
                            dst[((src_bytes - x - 1) * 8 + 7 - c) * dst_stride + y] = tmp[c];
                        } else {
                            dst[(x * 8 + c) * dst_stride + y] = tmp[c];
                        }
                    }
                }
            }
        } else {
            // FLIP_H mirrors x on the 8-bit format and y on the 2 and 4-bit formats, FLIP_V the other axis.
            const bool flip_x = BITS == 8 ? FLIP_H : FLIP_V;
            const bool flip_y = BITS == 8 ? FLIP_V : FLIP_H;
            const raster out{dst, pix.height, pix.width, bytes_per_line(pix.height)};
            for (size_t y = 0; y < pix.height; y++) {
                for (size_t x = 0; x < pix.width; x++) {
                    plot(out, flip_y ? (pix.height - 1 - y) : y, flip_x ? (pix.width - 1 - x) : x, get_col(pix, x, y));
                }
            }
        }
    }

    template <bool SUPPORTED = (BITS >= 4)>
    static constexpr void compose(const raster &pix, size_t x, size_t y, float cola, float colr, float colg,
                                  float colb) {
#ifndef _MSC_VER
        static_assert(SUPPORTED, "composing not supported on 1-bit and 2-bit formats, use a mono font.");
#endif  // #ifndef _MSC_VER
        if constexpr (BITS == 8) {
            // clang-format off
#if defined(__ARM_NEON)
                if (!std::is_constant_evaluated()) {
                    auto bg = static_cast<size_t>(pix.data[y * pix.stride + x]);
                    const float32x4_t cola_v = vdupq_n_f32(cola);
                    const float32x4_t inv_cola_v = vdupq_n_f32(1.0f - cola);
                    const float32x4_t col_rgb = {colr, colg, colb, 0.0f};
                    // Note: this reads 1 float into linearpal_neon
                    const float32x4_t bg_rgb = vld1q_f32(&quant.linear_palette()[bg * 3]);
                    const float32x4_t result_rgb = vaddq_f32(vmulq_f32(col_rgb, cola_v), vmulq_f32(bg_rgb, inv_cola_v));
                    alignas(16) std::array<float, 4> result{};
                    vst1q_f32(result.data(), result_rgb);
                    plot(pix, x, y, quant.nearest_linear(result[0], result[1], result[2]));
                } else
#endif  // #if defined(__ARM_NEON)
#if defined(__AVX2__)
                if (!std::is_constant_evaluated()) {
                    __m128 cola_v = _mm_set1_ps(cola);
                    __m128 inv_cola_v = _mm_set1_ps(1.0f - cola);
                    __m128 col_rgb = _mm_set_ps(0.0f, colb, colg, colr);
                    auto bg = static_cast<size_t>(pix.data[y * pix.stride + x]);
                    // Note: this reads 1 float into linearpal_avx2
                    __m128 bg_rgb = _mm_loadu_ps(&quant.linear_palette()[bg * 3]);
                    __m128 result_rgb = _mm_add_ps(_mm_mul_ps(col_rgb, cola_v), _mm_mul_ps(bg_rgb, inv_cola_v));
                    alignas(16) float result[4];
                    _mm_store_ps(result, result_rgb);
                    plot(pix, x, y, quant.nearest_linear(result[0], result[1], result[2]));
                } else
#endif  // #if defined(__AVX2__)
                {
                    const auto bg = static_cast<size_t>(pix.data[y * pix.stride + x]);
                    const float Rl = colr * cola + quant.linear_palette().at(bg * 3 + 0) * (1.0f - cola);
                    const float Gl = colg * cola + quant.linear_palette().at(bg * 3 + 1) * (1.0f - cola);
                    const float Bl = colb * cola + quant.linear_palette().at(bg * 3 + 2) * (1.0f - cola);
                    plot(pix, x, y, quant.nearest_linear(Rl, Gl, Bl));
                }
            // clang-format on
        } else if constexpr (BITS == 4) {
            const auto bg = static_cast<size_t>(get_col(pix, x, y));
            const float Rl = colr * cola + quant.linear_palette().at(bg * 3 + 0) * (1.0f - cola);
            const float Gl = colg * cola + quant.linear_palette().at(bg * 3 + 1) * (1.0f - cola);
            const float Bl = colb * cola + quant.linear_palette().at(bg * 3 + 2) * (1.0f - cola);
            plot(pix, x, y, quant.nearest_linear(Rl, Gl, Bl));
        }
    }

    static constexpr void plot(const raster &pix, size_t x, size_t y, uint8_t col) {
        uint8_t *ptr = &pix.data[y * pix.stride + x / pixels_per_byte];
        if constexpr (BITS == 8) {
            *ptr = col;
        } else {
            const size_t shift = (pixels_per_byte - 1 - x % pixels_per_byte) * BITS;
            *ptr &= static_cast<uint8_t>(~(color_mask << shift));
            *ptr |= static_cast<uint8_t>((col & color_mask) << shift);
        }
    }

    static constexpr void extent(const raster &pix, size_t xl0, size_t xr0, size_t y, uint8_t col) {
        uint8_t *yptr = &pix.data[y * pix.stride];
        if constexpr (BITS == 8) {
            for (size_t x = xl0; x < xr0; x++) {
                yptr[x] = col;
            }
        } else {
            const size_t xl = xl0 / pixels_per_byte;
            const size_t xr = xr0 / pixels_per_byte;
            // col in every pixel of a byte, and the pixels of a byte right of / left of a position
            const auto fill = static_cast<uint8_t>((col & color_mask) * (0xFF / color_mask));
            const auto ml = static_cast<uint8_t>(0xFF >> ((xl0 % pixels_per_byte) * BITS));
            const auto mr = static_cast<uint8_t>(0xFF00 >> ((xr0 % pixels_per_byte) * BITS));
            if (xr > xl) {
                yptr[xl] = static_cast<uint8_t>((yptr[xl] & ~ml) | (ml & fill));
                for (size_t x = xl + 1; x < xr; x++) {
                    yptr[x] = fill;
                }
                if (mr != 0) {
                    yptr[xr] = static_cast<uint8_t>((yptr[xr] & ~mr) | (mr & fill));
                }
            } else {
                yptr[xl] = static_cast<uint8_t>((yptr[xl] & ~(ml & mr)) | (ml & mr & fill));
            }
        }
    }

    [[nodiscard]] static constexpr uint8_t get_col(const uint8_t *line, size_t x) {
        const size_t shift = (pixels_per_byte - 1 - x % pixels_per_byte) * BITS;
        return static_cast<uint8_t>((line[x / pixels_per_byte] >> shift) & color_mask);
    }

    [[nodiscard]] static constexpr uint8_t get_col(const const_raster &pix, size_t x, size_t y) {
        return get_col(pix.data + y * pix.stride, x);
    }

    static constexpr void RGBA_uint32(uint32_t *dst, const const_raster &pix) {
        for (size_t y = 0; y < pix.height; y++) {
            for (size_t x = 0; x < pix.width; x++) {
                const uint8_t col = get_col(pix, x, y);
                dst[y * pix.width + x] = quant.palette().at(col) | ((col != 0) ? 0xFF000000 : 0x00000000);
            }
        }
    }

    static constexpr void RGBA_uint8(uint8_t *dst, const const_raster &pix) {
        for (size_t y = 0; y < pix.height; y++) {
            for (size_t x = 0; x < pix.width; x++) {
                const uint8_t col = get_col(pix, x, y);
                const uint32_t rgb = quant.palette().at(col);
                dst[0] = static_cast<uint8_t>((rgb >> 0) & 0xFF);
                dst[1] = static_cast<uint8_t>((rgb >> 8) & 0xFF);
                dst[2] = static_cast<uint8_t>((rgb >> 16) & 0xFF);
                dst[3] = ((col != 0) ? 0xFF : 0x00);
                dst += 4;
            }
        }
    }

    // The 1-bit format picks black or white by luminance, the others the nearest palette color.
    static constexpr void blit_RGBA(const raster &pix, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        if constexpr (BITS == 1) {
            const auto r_w = static_cast<size_t>(r.w < 0 ? 0 : r.w);
            const auto r_h = static_cast<size_t>(r.h < 0 ? 0 : r.h);
            for (size_t y = 0; y < r_h; y++) {
                for (size_t x = 0; x < r_w; x++) {
                    const uint32_t R = ptr[y * static_cast<size_t>(stride) + x * 4 + 0];
                    const uint32_t G = ptr[y * static_cast<size_t>(stride) + x * 4 + 1];
                    const uint32_t B = ptr[y * static_cast<size_t>(stride) + x * 4 + 2];
                    plot(pix, x + static_cast<size_t>(r.x), y + static_cast<size_t>(r.y),
                         (R * 2 + G * 3 + B * 1) > 768 ? 1 : 0);
                }
            }
        } else {
            const auto r_w = static_cast<size_t>(r.w < 0 ? 0 : r.w);
            const auto r_h = static_cast<size_t>(r.h < 0 ? 0 : r.h);
            for (size_t y = 0; y < r_h; y++) {
                for (size_t x = 0; x < r_w; x++) {
                    plot(pix, (x + static_cast<uint32_t>(r.x)), (y + static_cast<uint32_t>(r.y)),
                         quant.nearest(ptr[y * static_cast<size_t>(stride) + x * 4 + 0],
                                       ptr[y * static_cast<size_t>(stride) + x * 4 + 1],
                                       ptr[y * static_cast<size_t>(stride) + x * 4 + 2]));
                }
            }
        }
    }

    static constexpr void blit_RGBA_diffused(const raster &pix, const rect<int32_t> &r, const uint8_t *ptr,
                                             int32_t stride) {
        if constexpr (BITS == 1) {
            const auto r_w = static_cast<size_t>(r.w < 0 ? 0 : r.w);
            const auto r_h = static_cast<size_t>(r.h < 0 ? 0 : r.h);
            for (size_t y = 0; y < r_h; y++) {
                int32_t err = 0;
                for (size_t x = 0; x < r_w; x++) {
                    const int32_t R = ptr[y * static_cast<size_t>(stride) + x * 4 + 0];
                    const int32_t G = ptr[y * static_cast<size_t>(stride) + x * 4 + 1];
                    const int32_t B = ptr[y * static_cast<size_t>(stride) + x * 4 + 2];
                    const int32_t V = (R * 2 + G * 3 + B * 1) + err;
                    const uint8_t n = V > 768 ? 1 : 0;
                    plot(pix, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                    err = std::clamp(V - ((n != 0) ? int32_t{0xFF * 6} : int32_t{0x00}), int32_t{-0xFF * 6},
                                     int32_t{0xFF * 6});
                }
            }
        } else {
            const auto r_w = static_cast<size_t>(r.w < 0 ? 0 : r.w);
            const auto r_h = static_cast<size_t>(r.h < 0 ? 0 : r.h);
            for (size_t y = 0; y < r_h; y++) {
                int32_t err_r = 0;
                int32_t err_g = 0;
                int32_t err_b = 0;
                for (size_t x = 0; x < r_w; x++) {
                    int32_t R = ptr[y * static_cast<size_t>(stride) + x * 4 + 0];
                    int32_t G = ptr[y * static_cast<size_t>(stride) + x * 4 + 1];
                    int32_t B = ptr[y * static_cast<size_t>(stride) + x * 4 + 2];
                    R = R + err_r;
                    G = G + err_g;
                    B = B + err_b;
                    const uint8_t n = quant.nearest(R, G, B);
                    plot(pix, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                    err_r = std::clamp(R - static_cast<int32_t>((quant.palette().at(n) >> 0) & 0xFF), int32_t{-255},
                                       int32_t{255});
                    err_g = std::clamp(G - static_cast<int32_t>((quant.palette().at(n) >> 8) & 0xFF), int32_t{-255},
                                       int32_t{255});
                    err_b = std::clamp(B - static_cast<int32_t>((quant.palette().at(n) >> 16) & 0xFF), int32_t{-255},
                                       int32_t{255});
                }
            }
        }
    }

    static constexpr void blit_RGBA_diffused_linear(const raster &pix, const rect<int32_t> &r, const uint8_t *ptr,
                                                    int32_t stride) {
        if constexpr (BITS == 1) {
            const auto r_w = static_cast<size_t>(r.w < 0 ? 0 : r.w);
            const auto r_h = static_cast<size_t>(r.h < 0 ? 0 : r.h);
            for (size_t y = 0; y < r_h; y++) {
                float err_r = 0;
                float err_g = 0;
                float err_b = 0;
                for (size_t x = 0; x < r_w; x++) {
                    const auto R =
                        static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 0]) * (1.0f / 255.0f);
                    const auto G =
                        static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 1]) * (1.0f / 255.0f);
                    const auto B =
                        static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 2]) * (1.0f / 255.0f);
                    float Rl = hidden::srgb_to_linear(R);
                    float Gl = hidden::srgb_to_linear(G);
                    float Bl = hidden::srgb_to_linear(B);
                    Rl = Rl + err_r;
                    Gl = Gl + err_g;
                    Bl = Bl + err_b;
                    const uint8_t n = (Rl * 2 + Gl * 3 + Bl * 1) > 3.0f ? 1 : 0;
                    plot(pix, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                    const float c = 0.75;
                    err_r = std::clamp(Rl - ((n != uint8_t{0}) ? 1.0f : 0.0f), -c, c);
                    err_g = std::clamp(Gl - ((n != uint8_t{0}) ? 1.0f : 0.0f), -c, c);
                    err_b = std::clamp(Bl - ((n != uint8_t{0}) ? 1.0f : 0.0f), -c, c);
                }
            }
        } else {
            const auto r_w = static_cast<size_t>(r.w < 0 ? 0 : r.w);
            const auto r_h = static_cast<size_t>(r.h < 0 ? 0 : r.h);
            for (size_t y = 0; y < r_h; y++) {
                float err_r = 0;
                float err_g = 0;
                float err_b = 0;
                for (size_t x = 0; x < r_w; x++) {
                    const auto R =
                        static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 0]) * (1.0f / 255.0f);
                    const auto G =
                        static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 1]) * (1.0f / 255.0f);
                    const auto B =
                        static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 2]) * (1.0f / 255.0f);
                    float Rl = hidden::srgb_to_linear(R);
                    float Gl = hidden::srgb_to_linear(G);
                    float Bl = hidden::srgb_to_linear(B);
                    Rl = Rl + err_r;
                    Gl = Gl + err_g;
                    Bl = Bl + err_b;
                    const uint8_t n = quant.nearest_linear(Rl, Gl, Bl);
                    plot(pix, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                    err_r = std::clamp(Rl - quant.linear_palette().at(n * size_t{3} + size_t{0}), -1.0f, 1.0f);
                    err_g = std::clamp(Gl - quant.linear_palette().at(n * size_t{3} + size_t{1}), -1.0f, 1.0f);
                    err_b = std::clamp(Bl - quant.linear_palette().at(n * size_t{3} + size_t{2}), -1.0f, 1.0f);
                }
            }
        }
    }

    static constexpr void blit_RGBA_with_palette(const raster &pix, const rect<int32_t> &r, const uint8_t *ptr,
                                                 int32_t stride, const runtime_palette &rp) {
        rect<int32_t> intersect_rect{
            .x = 0, .y = 0, .w = static_cast<int32_t>(pix.width), .h = static_cast<int32_t>(pix.height)};
        intersect_rect &= r;
        auto const r_x = static_cast<size_t>(intersect_rect.x);
        auto const r_y = static_cast<size_t>(intersect_rect.y);
        auto const r_w = static_cast<size_t>(intersect_rect.w);
        auto const r_h = static_cast<size_t>(intersect_rect.h);
        for (size_t y = 0; y < r_h; y++) {
            for (size_t x = 0; x < r_w; x++) {
                const uint8_t *src = &ptr[y * static_cast<size_t>(stride) + x * 4];
                pix.data[(r_y + y) * pix.stride + (r_x + x)] = rp.nearest(src[0], src[1], src[2]);
            }
        }
    }

    template <typename F>
    static constexpr void png(const const_raster &pix, F &&char_out) {
        png_image<bits_per_pixel>(pix, bytes_per_line(pix.width), quant.palette(), std::forward<F>(char_out));
    }

    // Gives a runtime palette the shape of quant.palette() for sixel_image. Only the 8-bit format has one.
    struct palette_view {
        const runtime_palette *rp;
        [[nodiscard]] constexpr size_t size() const {
            return (BITS == 8 && rp != nullptr) ? rp->size() : quant.palette().size();
        }
        [[nodiscard]] constexpr uint32_t operator[](size_t i) const {
            return (BITS == 8 && rp != nullptr) ? rp->at(i) : quant.palette()[i];
        }
    };

    template <size_t S, typename F>
    static constexpr void sixel(const const_raster &pix, F &&char_out, const rect<int32_t> &r,
                                const runtime_palette *rp = nullptr) {
        sixel_image<uint8_t, sixel_bitset_size>(
            pix.data, palette_view{rp}, std::forward<F>(char_out), r, pix.width, pix.height, S,
            [&pix](const uint8_t *data_raw, size_t x, size_t col, size_t y) {
                const uint8_t *ptr = &data_raw[(y / S) * pix.stride];
                uint8_t out = 0;
                size_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    out >>= 1;
                    if ((y + y6) < pix.height * S) {
                        out |= static_cast<uint8_t>((get_col(ptr, x / S) == col) ? (1UL << 5) : 0);
                        if ((y + y6) != ((pix.height * S) - 1)) {
                            if (++inc >= S) {
                                inc = 0;
                                ptr += pix.stride;
                            }
                        }
                    }
                }
                return out;
            },
            [&pix](const uint8_t *data_raw, size_t x, size_t w, size_t y,
                   palette_bitset<uint8_t, sixel_bitset_size> &set) {
                // The 8-bit format starts at column x / S + x, the others at x. Kept as it was.
                const size_t x0 = x + (BITS == 8 ? x / S : 0);
                for (size_t y6 = 0; y6 < 6 && (y + y6) < pix.height * S; y6++) {
                    const uint8_t *ptr = &data_raw[((y + y6) / S) * pix.stride];
                    for (size_t xx = 0; xx < (w + S - 1) / S; xx++) {
                        set.mark(get_col(ptr, xx + x0));
                    }
                }
            });
    }
    /// @endcond
};

using pixel_1bit = pixel_packed<1, false>;
template <bool GRAYSCALE>
using pixel_2bit = pixel_packed<2, GRAYSCALE>;
template <bool GRAYSCALE>
using pixel_4bit = pixel_packed<4, GRAYSCALE>;

// The 8-bit format is the one which can draw with a runtime palette.
template <bool GRAYSCALE>
class pixel_8bit : public pixel_packed<8, GRAYSCALE> {
 public:
    using runtime_palette = ::constixel::runtime_palette;
};

// Pixel operations of the RGB (BYTES == 3) and RGBA (BYTES == 4) formats, independent of the image size.
template <size_t BYTES, bool GRAYSCALE>
class pixel_rgb : public format {
 public:
    /// @cond DOXYGEN_EXCLUDE
    static constexpr size_t sixel_bitset_size = 256;
    static constexpr size_t bits_per_pixel = BYTES * 8;
    static constexpr size_t color_mask = 0xFF;
    [[nodiscard]] static constexpr size_t bytes_per_line(size_t w) {
        return w * BYTES;
    }
    [[nodiscard]] static constexpr size_t image_size(size_t w, size_t h) {
        return h * bytes_per_line(w);
    }

    // Colors are drawn from the palette of the 8-bit format.
    static consteval auto gen_palette_consteval() {
        return pixel_packed<8, GRAYSCALE>::gen_palette_consteval();
    }
    static constexpr const auto &quant = pixel_packed<8, GRAYSCALE>::quant;

    template <bool FLIP_H = false, bool FLIP_V = false>
    static constexpr void transpose(const const_raster &pix, uint8_t *dst) {
        for (size_t y = 0; y < pix.height; y++) {
            const uint8_t *src = pix.data + y * pix.stride;
            for (size_t x = 0; x < pix.width; x++) {
                uint8_t *out =
                    dst +
                    ((FLIP_H ? (pix.width - x - 1) : x) * pix.height + (FLIP_V ? (pix.height - y - 1) : y)) * BYTES;
                for (size_t c = 0; c < BYTES; c++) {
                    out[c] = *src++;
                }
            }
        }
    }

    // Stores one pixel. The RGBA format writes whole words at run time.
    static constexpr void store(uint8_t *ptr, uint32_t rgba) {
        if constexpr (BYTES == 4) {
            if (!std::is_constant_evaluated()) {
                *reinterpret_cast<uint32_t *>(ptr) = rgba;
                return;
            }
        }
        for (size_t c = 0; c < BYTES; c++) {
            ptr[c] = (rgba >> (8 * c)) & 0xFF;
        }
    }

    static constexpr void plot(const raster &pix, size_t x, size_t y, uint8_t col) {
        store(&pix.data[y * pix.stride + x * BYTES], quant.palette().at(col) | 0xFF000000);
    }

    static constexpr void extent(const raster &pix, size_t xl0, size_t xr0, size_t y, uint8_t col) {
        uint8_t *yptr = &pix.data[y * pix.stride + xl0 * BYTES];
        const uint32_t rgba = quant.palette().at(col) | 0xFF000000;
        for (size_t x = xl0; x < xr0; x++) {
            store(yptr, rgba);
            yptr += BYTES;
        }
    }

    static constexpr void compose(const raster &pix, size_t x, size_t y, float cola, float colr, float colg,
                                  float colb) {
        if constexpr (BYTES == 4) {
#if defined(__ARM_NEON)
            if (!std::is_constant_evaluated()) {
                const size_t off = y * pix.stride + x * 4;
                uint8x8_t px = vld1_u8(&pix.data[off]);
                float32x4_t src = {colr, colg, colb, cola};
                float32x4_t dst = {hidden::a2al_8bit[px[0]], hidden::a2al_8bit[px[1]], hidden::a2al_8bit[px[2]],
                                   static_cast<float>(px[3]) * (1.0f / 255.0f)};
                float32x4_t one = vdupq_n_f32(1.0f);
                float32x4_t inv_srca = vsubq_f32(one, vdupq_n_f32(cola));
                float32x4_t blended = vmlaq_f32(vmulq_n_f32(src, cola), dst, inv_srca);

                float outa = cola + dst[3] * (1.0f - cola);
                float32x4_t final = hidden::linear_to_srgb_approx_neon(blended);

                uint8x8_t out;
                out[0] = static_cast<uint8_t>(final[0] * 255.0f);
                out[1] = static_cast<uint8_t>(final[1] * 255.0f);
                out[2] = static_cast<uint8_t>(final[2] * 255.0f);
                out[3] = static_cast<uint8_t>(outa * 255.0f);

                vst1_lane_u32(&pix.data[off], vreinterpret_u32_u8(out), 0);
                return;
            }
#endif  // #if defined(__ARM_NEON)
#if defined(__AVX2__)
            if (!std::is_constant_evaluated()) {
                const size_t off = y * pix.stride + x * 4;

                const float lr = hidden::a2al_8bit[pix.data[off + 0]];
                const float lg = hidden::a2al_8bit[pix.data[off + 1]];
                const float lb = hidden::a2al_8bit[pix.data[off + 2]];
                const float la = pix.data[off + 3] * (1.0f / 255.0f);

                const float as = cola + la * (1.0f - cola);
                const float inv_cola = 1.0f - cola;

                __m128 dst = _mm_set_ps(0.0f, lb, lg, lr);
                __m128 src = _mm_set_ps(0.0f, colb, colg, colr);
                __m128 blended = _mm_add_ps(_mm_mul_ps(src, _mm_set1_ps(cola)), _mm_mul_ps(dst, _mm_set1_ps(inv_cola)));
                __m128 srgb = hidden::linear_to_srgb_approx_sse(blended);

                alignas(16) float out[4];
                _mm_store_ps(out, srgb);

                pix.data[off + 0] = static_cast<uint8_t>(out[0] * 255.0f);
                pix.data[off + 1] = static_cast<uint8_t>(out[1] * 255.0f);
                pix.data[off + 2] = static_cast<uint8_t>(out[2] * 255.0f);
                pix.data[off + 3] = static_cast<uint8_t>(as * 255.0f);
                return;
            }
#endif  // #if defined(__AVX2__)
            const size_t off = y * pix.stride + x * 4;

            const float lr = hidden::a2al_8bit[pix.data[off + 0]];
            const float lg = hidden::a2al_8bit[pix.data[off + 1]];
            const float lb = hidden::a2al_8bit[pix.data[off + 2]];
            const float la = pix.data[off + 3] * (1.0f / 255.0f);

            const float as = cola + la * (1.0f - cola);
            const float rs = hidden::linear_to_srgb(colr * cola + lr * (1.0f - cola));
            const float gs = hidden::linear_to_srgb(colg * cola + lg * (1.0f - cola));
            const float bs = hidden::linear_to_srgb(colb * cola + lb * (1.0f - cola));

            pix.data[off + 0] = static_cast<uint8_t>(rs * 255.0f);
            pix.data[off + 1] = static_cast<uint8_t>(gs * 255.0f);
            pix.data[off + 2] = static_cast<uint8_t>(bs * 255.0f);
            pix.data[off + 3] = static_cast<uint8_t>(as * 255.0f);
        } else {
            const size_t off = y * pix.stride + x * 3;

            const float lr = hidden::a2al_8bit[pix.data[off + 0]];
            const float lg = hidden::a2al_8bit[pix.data[off + 1]];
            const float lb = hidden::a2al_8bit[pix.data[off + 2]];

            const float rs = hidden::linear_to_srgb(colr * cola + lr * (1.0f - cola));
            const float gs = hidden::linear_to_srgb(colg * cola + lg * (1.0f - cola));
            const float bs = hidden::linear_to_srgb(colb * cola + lb * (1.0f - cola));

            pix.data[off + 0] = static_cast<uint8_t>(rs * 255.0f);
            pix.data[off + 1] = static_cast<uint8_t>(gs * 255.0f);
            pix.data[off + 2] = static_cast<uint8_t>(bs * 255.0f);
        }
    }

    static constexpr void RGBA_uint32(uint32_t *dst, const const_raster &pix) {
        if constexpr (BYTES == 4) {
            if (std::is_constant_evaluated()) {
                for (size_t y = 0; y < pix.height; y++) {
                    for (size_t x = 0; x < pix.width; x++) {
                        dst[y * pix.width + x] = static_cast<uint32_t>((pix.data[y * pix.stride + x * 4 + 0]) |
                                                                       (pix.data[y * pix.stride + x * 4 + 1] << 8) |
                                                                       (pix.data[y * pix.stride + x * 4 + 2] << 16) |
                                                                       (pix.data[y * pix.stride + x * 4 + 3] << 24));
                    }
                }
            } else {
                for (size_t y = 0; y < pix.height; y++) {
                    std::memcpy(dst + y * pix.width, pix.data + y * pix.stride, pix.width * 4);
                }
            }
        } else {
            for (size_t y = 0; y < pix.height; y++) {
                for (size_t x = 0; x < pix.width; x++) {
                    dst[y * pix.width + x] = static_cast<uint32_t>((pix.data[y * pix.stride + x * 3 + 0]) |
                                                                   (pix.data[y * pix.stride + x * 3 + 1] << 8) |
                                                                   (pix.data[y * pix.stride + x * 3 + 2] << 16)) |
                                             0xFF000000;
                }
            }
        }
    }

    static constexpr void RGBA_uint8(uint8_t *dst, const const_raster &pix) {
        if constexpr (BYTES == 4) {
            for (size_t y = 0; y < pix.height; y++) {
                const uint8_t *src = pix.data + y * pix.stride;
                if (std::is_constant_evaluated()) {
                    for (size_t c = 0; c < pix.width * 4; c++) {
                        dst[c] = src[c];
                    }
                } else {
                    std::memcpy(dst, src, pix.width * 4);
                }
                dst += pix.width * 4;
            }
        } else {
            const uint8_t *src_ptr = pix.data;
            uint8_t *dst_ptr = dst;
            for (size_t y = 0; y < pix.height; y++) {
                for (size_t x = 0; x < pix.width; x++) {
                    dst_ptr[x * 4 + 0] = src_ptr[x * 3 + 0];
                    dst_ptr[x * 4 + 1] = src_ptr[x * 3 + 1];
                    dst_ptr[x * 4 + 2] = src_ptr[x * 3 + 2];
                    dst_ptr[x * 4 + 3] = 0xFF;
                }
                src_ptr += pix.stride;
                dst_ptr += pix.width * 4;
            }
        }
    }

    static constexpr void blit_RGBA(const raster &pix, const rect<int32_t> &r, const uint8_t *ptr, int32_t stride) {
        if constexpr (BYTES == 4) {
            rect<int32_t> intersect_rect{
                .x = 0, .y = 0, .w = static_cast<int32_t>(pix.width), .h = static_cast<int32_t>(pix.height)};
            intersect_rect &= rect<int32_t>{.x = r.x, .y = r.y, .w = r.w, .h = r.h};
            auto const r_x = static_cast<size_t>(intersect_rect.x);
            auto const r_y = static_cast<size_t>(intersect_rect.y);
            auto const r_w = static_cast<size_t>(intersect_rect.w);
            auto const r_h = static_cast<size_t>(intersect_rect.h);
            const uint8_t *src = ptr;
            uint8_t *dst = pix.data + r_y * pix.stride + r_x * 4;
            for (size_t y = 0; y < r_h; y++) {
                std::memcpy(dst, src, r_w * 4);
                dst += pix.stride;
                src += stride;
            }
        } else {
            rect<int32_t> intersect_rect{
                .x = 0, .y = 0, .w = static_cast<int32_t>(pix.width), .h = static_cast<int32_t>(pix.height)};
            intersect_rect &= rect<int32_t>{.x = r.x, .y = r.y, .w = r.w, .h = r.h};
            auto const r_x = static_cast<size_t>(intersect_rect.x);
            auto const r_y = static_cast<size_t>(intersect_rect.y);
            auto const r_w = static_cast<size_t>(intersect_rect.w);
            auto const r_h = static_cast<size_t>(intersect_rect.h);
            const uint8_t *src = ptr;
            uint8_t *dst = pix.data + r_y * pix.stride + r_x * 3;
            for (size_t y = 0; y < r_h; y++) {
                for (size_t x = 0; x < r_w; x++) {
                    dst[x * 3 + 0] = src[x * 4 + 0];
                    dst[x * 3 + 1] = src[x * 4 + 1];
                    dst[x * 3 + 2] = src[x * 4 + 2];
                }
                dst += pix.stride;
                src += stride;
            }
        }
    }

    static constexpr void blit_RGBA_diffused(const raster &pix, const rect<int32_t> &r, const uint8_t *ptr,
                                             int32_t stride) {
        blit_RGBA(pix, r, ptr, stride);
    }

    static constexpr void blit_RGBA_diffused_linear(const raster &pix, const rect<int32_t> &r, const uint8_t *ptr,
                                                    int32_t stride) {
        blit_RGBA(pix, r, ptr, stride);
    }

    template <typename F>
    static constexpr void png(const const_raster &pix, F &&char_out) {
        png_image<bits_per_pixel>(pix, bytes_per_line(pix.width), quant.palette(), std::forward<F>(char_out));
    }

    // Marks the formats whose sixel output takes a line cache of sixel_cache_size(columns) bytes.
    static constexpr bool sixel_line_cache = true;

    template <size_t S, typename F>
    static constexpr void sixel(const const_raster &pix, F &&char_out, const rect<int32_t> &r,
                                std::span<uint8_t> line_cache) {
        const size_t cache_columns = line_cache.size() / sixel_cache_size(1);
        auto line_cache_line = std::numeric_limits<size_t>::max();
        size_t line_cache_x = 0;
        sixel_image<uint8_t, sixel_bitset_size>(
            pix.data, quant.palette(), std::forward<F>(char_out), r, pix.width, pix.height, S,
            [&pix, line_cache, cache_columns, &line_cache_line, &line_cache_x](const uint8_t *data_raw, size_t x,
                                                                               size_t col, size_t y) {
                const size_t xs = x / S;
                if (line_cache_line != (y / S) || (xs - line_cache_x) >= cache_columns) {
                    line_cache_line = (y / S);
                    line_cache_x = sixel_cache_start(xs, pix.width, cache_columns);
                    const size_t columns = std::min(cache_columns, pix.width - line_cache_x);
                    const uint8_t *ptr = &data_raw[(y / S) * pix.stride + line_cache_x * BYTES];
                    for (size_t y6 = 0; y6 < 6; y6++) {
                        if (((y / S) + y6) < pix.height) {
                            for (size_t xx = 0; xx < columns; xx++) {
                                line_cache[xx * 6 + y6] =
                                    quant.nearest(ptr[xx * BYTES + 0], ptr[xx * BYTES + 1], ptr[xx * BYTES + 2]);
                            }
                            if (((y / S) + y6) != (pix.height - 1)) {
                                ptr += pix.stride;
                            }
                        }
                    }
                }
                const uint8_t *column = &line_cache[(xs - line_cache_x) * 6];
                uint8_t out = 0;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    out >>= 1;
                    if ((y + y6) < (pix.height * S)) {
                        size_t y6_l = (((y + y6) % 6) / S);
                        const uint8_t *ptr = &column[y6_l];
                        out |= static_cast<uint8_t>((*ptr == col) ? (1UL << 5) : 0);
                    }
                }
                return out;
            },
            [&pix](const uint8_t *data_raw, size_t x, size_t w, size_t y,
                   palette_bitset<uint8_t, sixel_bitset_size> &set) {
                const uint8_t *ptr = &data_raw[(y / S) * pix.stride + (x / S) * BYTES];
                size_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    if ((y + y6) < pix.height * S) {
                        for (size_t xx = 0; xx < (w + S - 1) / S; xx++) {
                            uint8_t ncol = quant.nearest(ptr[(xx + x) * BYTES + 0], ptr[(xx + x) * BYTES + 1],
                                                         ptr[(xx + x) * BYTES + 2]);
                            set.mark(ncol);
                        }
                        if ((y + y6) != ((pix.height * S) - 1)) {
                            if (++inc >= S) {
                                inc = 0;
                                ptr += pix.stride;
                            }
                        }
                    }
                }
            });
    }
    /// @endcond
};

template <bool GRAYSCALE>
using pixel_24bit = pixel_rgb<3, GRAYSCALE>;
template <bool GRAYSCALE>
using pixel_32bit = pixel_rgb<4, GRAYSCALE>;
/// @endcond
/** @brief 1-bit format, just b/w. Use as template parameter for image. Example:
 *
 * \code{.cpp}
 * constixel::image<constixel::format_1bit, 640, 480> image;
 * \endcode
 *
 * @tparam W Width in pixels.
 * @tparam H Height in pixels.
 * @tparam GRAYSCALE Grayscale palette.
 */
template <size_t W, size_t H, bool GRAYSCALE, bool USE_SPAN>
class format_1bit : public pixel_1bit {
 public:
    /// @cond DOXYGEN_EXCLUDE
    using pixel_type = pixel_1bit;
    static constexpr size_t bytes_per_line = pixel_type::bytes_per_line(W);
    static constexpr size_t image_size = pixel_type::image_size(W, H);
    /// @endcond
};

/**
 * @brief 2-bit color format, 4 colors total. Use as template parameter for image.
 *
 * \code{.cpp}
 * constixel::image<constixel::format_2bit, 640, 480> image;
 * \endcode
 *
 * @tparam W Width in pixels.
 * @tparam H Height in pixels.
 * @tparam GRAYSCALE Grayscale palette.
 */
template <size_t W, size_t H, bool GRAYSCALE, bool USE_SPAN>
class format_2bit : public pixel_2bit<GRAYSCALE> {
 public:
    /// @cond DOXYGEN_EXCLUDE
    using pixel_type = pixel_2bit<GRAYSCALE>;
    static constexpr size_t bytes_per_line = pixel_type::bytes_per_line(W);
    static constexpr size_t image_size = pixel_type::image_size(W, H);
    /// @endcond
};

/**
 * @brief 4-bit color format, 16 colors total. Use as template parameter for image.
 *
 * \code{.cpp}
 * constixel::image<constixel::format_4bit, 640, 480> image;
 * \endcode
 *
 * @tparam W Width in pixels.
 * @tparam H Height in pixels.
 * @tparam GRAYSCALE Grayscale palette.
 */
template <size_t W, size_t H, bool GRAYSCALE, bool USE_SPAN>
class format_4bit : public pixel_4bit<GRAYSCALE> {
 public:
    /// @cond DOXYGEN_EXCLUDE
    using pixel_type = pixel_4bit<GRAYSCALE>;
    static constexpr size_t bytes_per_line = pixel_type::bytes_per_line(W);
    static constexpr size_t image_size = pixel_type::image_size(W, H);
    /// @endcond
};

/**
 * @brief 8-bit format, 256 colors total. Use as template parameter for image. Example:
 *
 * \code{.cpp}
 * constixel::image<constixel::format_8bit, 640, 480> image;
 * \endcode
 *
 * @tparam W Width in pixels.
 * @tparam H Height in pixels.
 * @tparam GRAYSCALE Grayscale palette.
 */
template <size_t W, size_t H, bool GRAYSCALE, bool USE_SPAN>
class format_8bit : public pixel_8bit<GRAYSCALE> {
 public:
    /// @cond DOXYGEN_EXCLUDE
    using pixel_type = pixel_8bit<GRAYSCALE>;
    static constexpr size_t bytes_per_line = pixel_type::bytes_per_line(W);
    static constexpr size_t image_size = pixel_type::image_size(W, H);
    /// @endcond
};

/**
 * @brief 24-bit format. Use as template parameter for image. Example:
 *
 * \code{.cpp}
 * constixel::image<constixel::format_24bit, 640, 480> image;
 * \endcode
 *
 * @tparam W Width in pixels.
 * @tparam H Height in pixels.
 * @tparam GRAYSCALE Grayscale palette.
 */
template <size_t W, size_t H, bool GRAYSCALE, bool USE_SPAN>
class format_24bit : public pixel_24bit<GRAYSCALE> {
 public:
    /// @cond DOXYGEN_EXCLUDE
    using pixel_type = pixel_24bit<GRAYSCALE>;
    static constexpr size_t bytes_per_line = pixel_type::bytes_per_line(W);
    static constexpr size_t image_size = pixel_type::image_size(W, H);
    /// @endcond
};

/**
 * @brief 32-bit format. Use as template parameter for image. Example:
 *
 * \code{.cpp}
 * constixel::image<constixel::format_32bit, 640, 480> image;
 * \endcode
 *
 * @tparam W Width in pixels.
 * @tparam H Height in pixels.
 * @tparam GRAYSCALE Grayscale palette.
 */
template <size_t W, size_t H, bool GRAYSCALE, bool USE_SPAN>
class format_32bit : public pixel_32bit<GRAYSCALE> {
 public:
    /// @cond DOXYGEN_EXCLUDE
    using pixel_type = pixel_32bit<GRAYSCALE>;
    static constexpr size_t bytes_per_line = pixel_type::bytes_per_line(W);
    static constexpr size_t image_size = pixel_type::image_size(W, H);
    /// @endcond
};

/// @cond DOXYGEN_EXCLUDE
// Former runtime-sized 8-bit format; the pixel classes are all runtime-sized.
using format_8bit_dyn = pixel_8bit<false>;

namespace hidden {
template <typename PF>
concept has_runtime_palette = requires { typename PF::runtime_palette; };
template <typename PF>
concept has_sixel_cache = requires { PF::sixel_line_cache; };
}  // namespace hidden
/// @endcond

/**
 * @enum color
 * @brief Extended color enum with all 256 colors from the 8-bit palette.
 * Provides named constants for common colors, gradients, and palette ranges.
 */
enum color : uint8_t {
    // Basic colors (0-15)
    BLACK = 0,       /*!< Black */
    TRANSPARENT = 0, /*!< Transparent (alias for BLACK) */
    WHITE = 1,       /*!< White */
    RED = 2,         /*!< Red */
    LIME = 3,        /*!< Lime/Bright Green */
    BLUE = 4,        /*!< Blue */
    YELLOW = 5,      /*!< Yellow */
    CYAN = 6,        /*!< Cyan */
    MAGENTA = 7,     /*!< Magenta */
    GRAY_20 = 8,     /*!< Gray 20% */
    GRAY_40 = 9,     /*!< Gray 40% */
    GRAY_60 = 10,    /*!< Gray 60% */
    GRAY_80 = 11,    /*!< Gray 80% */
    MAROON = 12,     /*!< Maroon (Dark Red) */
    GREEN = 13,      /*!< Green (Dark Green) */
    NAVY = 14,       /*!< Navy (Dark Blue) */
    OLIVE = 15,      /*!< Olive (Dark Yellow) */

    // Gray ramp (16-31)
    BLACK_OPAQUE = 16, /*!< Black (same as BLACK but in gray ramp) */
    GRAY_10 = 17,      /*!< Gray 10% */
    GRAY_11 = 18,      /*!< Gray 13% */
    GRAY_20_ALT = 19,  /*!< Gray 20% (alternative) */
    GRAY_30 = 20,      /*!< Gray 27% */
    GRAY_33 = 21,      /*!< Gray 33% */
    GRAY_40_ALT = 22,  /*!< Gray 40% (alternative) */
    GRAY_47 = 23,      /*!< Gray 47% */
    GRAY_53 = 24,      /*!< Gray 53% */
    GRAY_60_ALT = 25,  /*!< Gray 60% (alternative) */
    DARK_GRAY = 26,    /*!< Dark Gray (67%) */
    SILVER = 27,       /*!< Silver (73%) */
    GRAY_80_ALT = 28,  /*!< Gray 80% (alternative) */
    GAINSBORO = 29,    /*!< Gainsboro (87%) */
    LIGHT_GRAY = 30,   /*!< Light Gray (93%) */
    WHITE_ALT = 31,    /*!< White (alternative) */

    // Red luminance ramp (32-47)
    RED_DARK_1 = 32,   /*!< Very Dark Red */
    RED_DARK_2 = 33,   /*!< Dark Red */
    RED_DARK_3 = 34,   /*!< Medium Dark Red */
    RED_MEDIUM_1 = 35, /*!< Medium Red */
    RED_MEDIUM_2 = 36, /*!< Medium Bright Red */
    RED_BRIGHT_1 = 37, /*!< Bright Red */
    RED_BRIGHT_2 = 38, /*!< Very Bright Red */
    RED_FULL = 39,     /*!< Full Red */
    RED_PINK_1 = 40,   /*!< Light Red/Pink */
    RED_PINK_2 = 41,   /*!< Lighter Red/Pink */
    RED_PINK_3 = 42,   /*!< Light Pink */
    RED_PINK_4 = 43,   /*!< Very Light Pink */
    RED_PINK_5 = 44,   /*!< Pale Pink */
    RED_PINK_6 = 45,   /*!< Very Pale Pink */
    RED_PALE_1 = 46,   /*!< Extremely Pale Pink */
    RED_WHITE = 47,    /*!< Near White with Red tint */

    // Green luminance ramp (48-63)
    GREEN_DARK_1 = 48,   /*!< Very Dark Green */
    GREEN_DARK_2 = 49,   /*!< Dark Green */
    GREEN_DARK_3 = 50,   /*!< Medium Dark Green */
    GREEN_MEDIUM_1 = 51, /*!< Medium Green */
    GREEN_MEDIUM_2 = 52, /*!< Medium Bright Green */
    GREEN_BRIGHT_1 = 53, /*!< Bright Green */
    GREEN_BRIGHT_2 = 54, /*!< Very Bright Green */
    GREEN_FULL = 55,     /*!< Full Green */
    GREEN_LIGHT_1 = 56,  /*!< Light Green */
    GREEN_LIGHT_2 = 57,  /*!< Lighter Green */
    GREEN_LIGHT_3 = 58,  /*!< Very Light Green */
    GREEN_LIGHT_4 = 59,  /*!< Pale Green */
    GREEN_PALE_1 = 60,   /*!< Very Pale Green */
    GREEN_PALE_2 = 61,   /*!< Extremely Pale Green */
    GREEN_PALE_3 = 62,   /*!< Near White with Green tint */
    GREEN_WHITE = 63,    /*!< White with Green tint */

    // Blue luminance ramp (64-79)
    BLUE_DARK_1 = 64,   /*!< Very Dark Blue */
    BLUE_DARK_2 = 65,   /*!< Dark Blue */
    BLUE_DARK_3 = 66,   /*!< Medium Dark Blue */
    BLUE_MEDIUM_1 = 67, /*!< Medium Blue */
    BLUE_MEDIUM_2 = 68, /*!< Medium Bright Blue */
    BLUE_BRIGHT_1 = 69, /*!< Bright Blue */
    BLUE_BRIGHT_2 = 70, /*!< Very Bright Blue */
    BLUE_FULL = 71,     /*!< Full Blue */
    BLUE_LIGHT_1 = 72,  /*!< Light Blue */
    BLUE_LIGHT_2 = 73,  /*!< Lighter Blue */
    BLUE_LIGHT_3 = 74,  /*!< Very Light Blue */
    BLUE_LIGHT_4 = 75,  /*!< Pale Blue */
    BLUE_PALE_1 = 76,   /*!< Very Pale Blue */
    BLUE_PALE_2 = 77,   /*!< Extremely Pale Blue */
    BLUE_PALE_3 = 78,   /*!< Near White with Blue tint */
    BLUE_WHITE = 79,    /*!< White with Blue tint */

    // Yellow luminance ramp (80-95)
    YELLOW_DARK_1 = 80,   /*!< Very Dark Yellow */
    YELLOW_DARK_2 = 81,   /*!< Dark Yellow */
    YELLOW_DARK_3 = 82,   /*!< Medium Dark Yellow */
    YELLOW_MEDIUM_1 = 83, /*!< Medium Yellow */
    YELLOW_MEDIUM_2 = 84, /*!< Medium Bright Yellow */
    YELLOW_BRIGHT_1 = 85, /*!< Bright Yellow */
    YELLOW_BRIGHT_2 = 86, /*!< Very Bright Yellow */
    YELLOW_FULL = 87,     /*!< Full Yellow */
    YELLOW_LIGHT_1 = 88,  /*!< Light Yellow */
    YELLOW_LIGHT_2 = 89,  /*!< Lighter Yellow */
    YELLOW_LIGHT_3 = 90,  /*!< Very Light Yellow */
    YELLOW_LIGHT_4 = 91,  /*!< Pale Yellow */
    YELLOW_PALE_1 = 92,   /*!< Very Pale Yellow */
    YELLOW_PALE_2 = 93,   /*!< Extremely Pale Yellow */
    YELLOW_PALE_3 = 94,   /*!< Near White with Yellow tint */
    YELLOW_WHITE = 95,    /*!< White with Yellow tint */

    // Cyan luminance ramp (96-111)
    CYAN_DARK_1 = 96,    /*!< Very Dark Cyan */
    CYAN_DARK_2 = 97,    /*!< Dark Cyan */
    CYAN_DARK_3 = 98,    /*!< Medium Dark Cyan */
    CYAN_MEDIUM_1 = 99,  /*!< Medium Cyan */
    CYAN_MEDIUM_2 = 100, /*!< Medium Bright Cyan */
    CYAN_BRIGHT_1 = 101, /*!< Bright Cyan */
    CYAN_BRIGHT_2 = 102, /*!< Very Bright Cyan */
    CYAN_FULL = 103,     /*!< Full Cyan */
    CYAN_LIGHT_1 = 104,  /*!< Light Cyan */
    CYAN_LIGHT_2 = 105,  /*!< Lighter Cyan */
    CYAN_LIGHT_3 = 106,  /*!< Very Light Cyan */
    CYAN_LIGHT_4 = 107,  /*!< Pale Cyan */
    CYAN_PALE_1 = 108,   /*!< Very Pale Cyan */
    CYAN_PALE_2 = 109,   /*!< Extremely Pale Cyan */
    CYAN_PALE_3 = 110,   /*!< Near White with Cyan tint */
    CYAN_WHITE = 111,    /*!< White with Cyan tint */

    // Magenta luminance ramp (112-127)
    MAGENTA_DARK_1 = 112,   /*!< Very Dark Magenta */
    MAGENTA_DARK_2 = 113,   /*!< Dark Magenta */
    MAGENTA_DARK_3 = 114,   /*!< Medium Dark Magenta */
    MAGENTA_MEDIUM_1 = 115, /*!< Medium Magenta */
    MAGENTA_MEDIUM_2 = 116, /*!< Medium Bright Magenta */
    MAGENTA_BRIGHT_1 = 117, /*!< Bright Magenta */
    MAGENTA_BRIGHT_2 = 118, /*!< Very Bright Magenta */
    MAGENTA_FULL = 119,     /*!< Full Magenta */
    MAGENTA_LIGHT_1 = 120,  /*!< Light Magenta */
    MAGENTA_LIGHT_2 = 121,  /*!< Lighter Magenta */
    MAGENTA_LIGHT_3 = 122,  /*!< Very Light Magenta */
    MAGENTA_LIGHT_4 = 123,  /*!< Pale Magenta */
    MAGENTA_PALE_1 = 124,   /*!< Very Pale Magenta */
    MAGENTA_PALE_2 = 125,   /*!< Extremely Pale Magenta */
    MAGENTA_PALE_3 = 126,   /*!< Near White with Magenta tint */
    MAGENTA_WHITE = 127,    /*!< White with Magenta tint */

    DARK_RED = 12,   /*!< Dark Red (alias) */
    DARK_GREEN = 13, /*!< Dark Green (alias) */
    DARK_BLUE = 14,  /*!< Dark Blue (alias) */
    ORANGE = 38,     /*!< Orange (closest match) */
    PURPLE = 115,    /*!< Purple (closest match) */
    PINK = 125,      /*!< Pink (closest match) */
    BROWN = 83,      /*!< Brown (closest match) */
    TEAL = 99,       /*!< Teal (closest match) */

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
    MAGENTA_LUMA_RAMP_STOP = MAGENTA_LUMA_RAMP_START + 15,

    OKLCH_RAMP_START = 128,
    OKLCH_RAMP_COUNT = 128,
    OKLCH_RAMP_STOP = 255
};

/**
 * @enum text_rotation
 * @brief Text rotation. One of DEGREE_0, DEGREE_90, DEGREE_180 or DEGREE_270.
 */
enum text_rotation {
    DEGREE_0,
    DEGREE_90,
    DEGREE_180,
    DEGREE_270
};

/**
 * @enum text_align
 * @brief Where a string is placed relative to its x/y coordinate, along the direction of the text. One of ALIGN_LEFT
 * (the string starts there), ALIGN_CENTER or ALIGN_RIGHT (the string ends there).
 */
enum text_align {
    ALIGN_LEFT,
    ALIGN_CENTER,
    ALIGN_RIGHT
};

/**
 * @enum device_format
 * @brief Data formats for the convert function
 */
enum device_format {
    STRAIGHT_THROUGH,    //!< Just copy the data as is.
    RGB565_8BIT_SERIAL,  //!< RGB565 pixel data is stored from left to right, each two bytes containing 1 pixel value in
                         //!< the x direction.
                         //   Byte encoding: 0xRRRRRGGG 0xGGGBBBBB
    RGB666_8BIT_SERIAL_1,  //!< RGB565 pixel data is stored from left to right, each three bytes containing 1 pixel
                           //!< values in the x direction.
                           //   Byte encoding: 0x00RRRRRR 0x00GGGGGG 0x00BBBBBB
    RGB666_8BIT_SERIAL_2,  //!< RGB565 pixel data is stored from left to right, each three bytes containing 1 pixel
                           //!< values in the x direction.
                           //   Byte encoding: 0xRRRRRR00 0xGGGGGG00 0xBBBBBB00
    MONO_1BIT_VERTICAL_PAGES,  //!< 1-bit formats only. Pages of 8 rows from top to bottom, one byte per column from
                               //!< left to right, the least significant bit is the top row of the page. The page
                               //!< order of SSD1306 and SH1106 displays.
    MONO_1BIT_HORIZONTAL_LSB   //!< 1-bit formats only. Rows from top to bottom, 8 pixels per byte from left to
                               //!< right, the least significant bit is the leftmost pixel.
};

/**
 * @brief Position in the byte stream of convert_chunk(). Start with a default-constructed cursor and pass the same
 * one to every call for a frame or rectangle.
 */
struct convert_cursor {
    size_t line = 0;   /**< Output line (row, or page for MONO_1BIT_VERTICAL_PAGES) of the next byte. */
    size_t offset = 0; /**< Byte offset of the next byte in that line. */
};

/**
 * A struct which can be passed to: plot().
 */
struct plot {
    int32_t x = 0;              /**< X coordinate in pixels. */
    int32_t y = 0;              /**< Y coordinate in pixels. */
    uint8_t col = color::WHITE; /**< Color palette index to use. */
};

/**
 * A struct which can be passed to: draw_line(), draw_line_aa().
 */
struct draw_line {
    int32_t x0 = 0;             /**< First X coordinate in pixels. */
    int32_t y0 = 0;             /**< First Y coordinate in pixels. */
    int32_t x1 = 0;             /**< Second X coordinate in pixels. */
    int32_t y1 = 0;             /**< Second Y coordinate in pixels. */
    uint8_t col = color::WHITE; /**< Color palette index to use. */
    float sw = 1;               /**< Width of the stroke in pixels. */
};

/**
 * A struct which can be passed to: fill_rect(), stroke_rect(), .
 */
struct draw_rect {
    int32_t x = 0;              /**< X coordinate in pixels. */
    int32_t y = 0;              /**< Y coordinate in pixels. */
    int32_t w = 0;              /**< Width in pixels. */
    int32_t h = 0;              /**< Height in pixels. */
    uint8_t col = color::WHITE; /**< Color palette index to use. */
    int32_t sw = 1;             /**< Width of the stroke in pixels. */
};

/**
 * A struct which can be passed to: fill_circle(), stroke_circle(), fill_circle_aa(), stroke_circle_aa().
 */
struct draw_circle {
    int32_t cx = 0;             /**< Center X coordinate in pixels. */
    int32_t cy = 0;             /**< Center Y coordinate in pixels. */
    int32_t r = 0;              /**< Radius of the circle in pixels. */
    uint8_t col = color::WHITE; /**< Color palette index to use. */
    int32_t sw = 1;             /**< Width of the stroke in pixels. */
};

/**
 * A struct which can be passed to: fill_round_rect(), stroke_round_rect(), fill_round_rect_aa(),
 * stroke_round_rect_aa().
 */
struct draw_round_rect {
    int32_t x = 0;              /**< X coordinate in pixels. */
    int32_t y = 0;              /**< Y coordinate in pixels. */
    int32_t w = 0;              /**< Width in pixels. */
    int32_t h = 0;              /**< Height in pixels. */
    int32_t r = 0;              /**< Radius of the corners in pixels. */
    uint8_t col = color::WHITE; /**< Color palette index to use. */
    int32_t sw = 1;             /**< Width of the stroke in pixels. */
};

/**
 * A struct which can be passed to: draw_string(), draw_string_mono(), draw_string_aa().
 */
struct draw_string {
    int32_t x = 0;              /**< X coordinate in pixels. */
    int32_t y = 0;              /**< Y coordinate in pixels. */
    const char *str = nullptr;  /**< UTF8 string */
    uint8_t col = color::WHITE; /**< Color palette index to use. */
    text_align align = ALIGN_LEFT; /**< Placement relative to x/y. */
};

template <typename PF>
class surface;

/// @cond DOXYGEN_EXCLUDE
namespace hidden {
// View of a NUL-terminated string, empty for nullptr.
[[nodiscard]] constexpr std::string_view cstr_view(const char *str) {
    return str != nullptr ? std::string_view{str} : std::string_view{};
}
}  // namespace hidden
/// @endcond

/**
 * @brief Fluent shape API classes for method chaining
 */
namespace shapes {

/// @cond DOXYGEN_EXCLUDE
template <typename shader_func>
concept shader_function = std::is_invocable_r_v<std::array<float, 4>, shader_func, float, float, int32_t, int32_t>;
/// @endcond

/**
 * @brief Fluent API for drawing rectangles
 *
 * Provides a chainable interface for drawing filled and stroked rectangles.
 */
template <typename PF>
class rect {
    surface<PF> img;
    int32_t x, y, w, h;

 public:
    /**
     * @brief Construct a rectangle shape
     * @param image Target surface to draw on
     * @param x_ X coordinate of top-left corner
     * @param y_ Y coordinate of top-left corner
     * @param w_ Width of rectangle
     * @param h_ Height of rectangle
     */
    constexpr rect(const surface<PF> &image, int32_t x_, int32_t y_, int32_t w_, int32_t h_)
        : img(image), x(x_), y(y_), w(w_), h(h_) {
    }

    /**
     * @brief Fill the rectangle with a solid color
     * @param col Color value
     * @return Reference to this rectangle for chaining
     */
    constexpr rect &fill(uint8_t col) {
        img.fill_rect(x, y, w, h, col);
        return *this;
    }

    /**
     * @brief Fill the rectangle using a shader function
     * @param shader Function that returns RGBA values based on normalized coordinates
     * @return Reference to this rectangle for chaining
     */
    constexpr rect &fill_shader(const shader_function auto &shader) {
        img.fill_rect(x, y, w, h, shader);
        return *this;
    }

    /**
     * @brief Draw the rectangle outline
     * @param col Color value
     * @param stroke_width Width of the stroke (default: 1)
     * @return Reference to this rectangle for chaining
     */
    constexpr rect &stroke(uint8_t col, int32_t stroke_width = 1) {
        img.stroke_rect(x, y, w, h, col, stroke_width);
        return *this;
    }
};

/**
 * @brief Fluent API for drawing circles, plain (circle) or antialiased (circle_aa)
 *
 * Provides a chainable interface for drawing filled and stroked circles.
 * @tparam AA Draw with antialiasing.
 */
template <typename PF, bool AA>
class basic_circle {
    surface<PF> img;
    int32_t cx, cy, r;

 public:
    /**
     * @brief Construct a circle shape
     * @param image Target surface to draw on
     * @param cx_ X coordinate of center
     * @param cy_ Y coordinate of center
     * @param r_ Radius of circle
     */
    constexpr basic_circle(const surface<PF> &image, int32_t cx_, int32_t cy_, int32_t r_)
        : img(image), cx(cx_), cy(cy_), r(r_) {
    }

    /**
     * @brief Fill the circle with a solid color
     * @param col Color value
     * @return Reference to this circle for chaining
     */
    constexpr basic_circle &fill(uint8_t col) {
        if constexpr (AA) {
            img.fill_circle_aa(cx, cy, r, col);
        } else {
            img.fill_circle(cx, cy, r, col);
        }
        return *this;
    }

    /**
     * @brief Fill the circle using a shader function. Antialiased circles only.
     * @param shader Function that returns RGBA values based on normalized coordinates
     * @return Reference to this circle for chaining
     */
    constexpr basic_circle &fill_shader(const shader_function auto &shader)
        requires AA
    {
        img.fill_circle_aa(cx, cy, r, shader);
        return *this;
    }

    /**
     * @brief Draw the circle outline
     * @param col Color value
     * @param stroke_width Width of the stroke (default: 1)
     * @return Reference to this circle for chaining
     */
    constexpr basic_circle &stroke(uint8_t col, int32_t stroke_width = 1) {
        if constexpr (AA) {
            img.stroke_circle_aa(cx, cy, r, col, stroke_width);
        } else {
            img.stroke_circle(cx, cy, r, col, stroke_width);
        }
        return *this;
    }
};

/// @brief Fluent API for drawing circles
template <typename PF>
using circle = basic_circle<PF, false>;
/// @brief Fluent API for drawing anti-aliased circles
template <typename PF>
using circle_aa = basic_circle<PF, true>;

/**
 * @brief Fluent API for drawing rounded rectangles, plain (round_rect) or antialiased (round_rect_aa)
 *
 * Provides a chainable interface for drawing filled and stroked rounded rectangles.
 * @tparam AA Draw with antialiasing.
 */
template <typename PF, bool AA>
class basic_round_rect {
    surface<PF> img;
    int32_t x, y, w, h, radius;

 public:
    /**
     * @brief Construct a rounded rectangle shape
     * @param image Target surface to draw on
     * @param x_ X coordinate of top-left corner
     * @param y_ Y coordinate of top-left corner
     * @param w_ Width of rectangle
     * @param h_ Height of rectangle
     * @param radius_ Corner radius
     */
    constexpr basic_round_rect(const surface<PF> &image, int32_t x_, int32_t y_, int32_t w_, int32_t h_,
                               int32_t radius_)
        : img(image), x(x_), y(y_), w(w_), h(h_), radius(radius_) {
    }

    /**
     * @brief Fill the rounded rectangle with a solid color
     * @param col Color value
     * @return Reference to this rounded rectangle for chaining
     */
    constexpr basic_round_rect &fill(uint8_t col) {
        if constexpr (AA) {
            img.fill_round_rect_aa(x, y, w, h, radius, col);
        } else {
            img.fill_round_rect(x, y, w, h, radius, col);
        }
        return *this;
    }

    /**
     * @brief Fill the rounded rectangle using a shader function. Antialiased rounded rectangles only.
     * @param shader Function that returns RGBA values based on normalized coordinates
     * @return Reference to this rounded rectangle for chaining
     */
    constexpr basic_round_rect &fill_shader(const shader_function auto &shader)
        requires AA
    {
        img.fill_round_rect_aa(x, y, w, h, radius, shader);
        return *this;
    }

    /**
     * @brief Draw the rounded rectangle outline
     * @param col Color value
     * @param stroke_width Width of the stroke (default: 1)
     * @return Reference to this rounded rectangle for chaining
     */
    constexpr basic_round_rect &stroke(uint8_t col, int32_t stroke_width = 1) {
        if constexpr (AA) {
            img.stroke_round_rect_aa(x, y, w, h, radius, col, stroke_width);
        } else {
            img.stroke_round_rect(x, y, w, h, radius, col, stroke_width);
        }
        return *this;
    }
};

/// @brief Fluent API for drawing rounded rectangles
template <typename PF>
using round_rect = basic_round_rect<PF, false>;
/// @brief Fluent API for drawing anti-aliased rounded rectangles
template <typename PF>
using round_rect_aa = basic_round_rect<PF, true>;

/**
 * @brief Fluent API for drawing lines, plain (line) or antialiased (line_aa)
 *
 * Provides a chainable interface for drawing stroked lines.
 * @tparam AA Draw with antialiasing. The stroke width is a float then.
 */
template <typename PF, bool AA>
class basic_line {
    surface<PF> img;
    int32_t x0, y0, x1, y1;

 public:
    /**
     * @brief Construct a line shape
     * @param image Target surface to draw on
     * @param x0_ X coordinate of start point
     * @param y0_ Y coordinate of start point
     * @param x1_ X coordinate of end point
     * @param y1_ Y coordinate of end point
     */
    constexpr basic_line(const surface<PF> &image, int32_t x0_, int32_t y0_, int32_t x1_, int32_t y1_)
        : img(image), x0(x0_), y0(y0_), x1(x1_), y1(y1_) {
    }

    /**
     * @brief Draw the line with the specified color and width
     * @param col Color value
     * @param stroke_width Width of the line (default: 1)
     * @return Reference to this line for chaining
     */
    constexpr basic_line &stroke(uint8_t col, std::conditional_t<AA, float, int32_t> stroke_width = 1) {
        if constexpr (AA) {
            img.draw_line_aa(x0, y0, x1, y1, col, stroke_width);
        } else {
            img.draw_line(x0, y0, x1, y1, col, stroke_width);
        }
        return *this;
    }
};

/// @brief Fluent API for drawing lines
template <typename PF>
using line = basic_line<PF, false>;
/// @brief Fluent API for drawing anti-aliased lines
template <typename PF>
using line_aa = basic_line<PF, true>;

/**
 * @brief Fluent API for drawing points
 *
 * Provides a chainable interface for plotting individual pixels.
 */
template <typename PF>
class point {
    surface<PF> img;
    int32_t x, y;

 public:
    /**
     * @brief Construct a point shape
     * @param image Target surface to draw on
     * @param x_ X coordinate of the point
     * @param y_ Y coordinate of the point
     */
    constexpr point(const surface<PF> &image, int32_t x_, int32_t y_) : img(image), x(x_), y(y_) {
    }

    /**
     * @brief Plot the point with the specified color
     * @param col Color value
     * @return Reference to this point for chaining
     */
    constexpr point &plot(uint8_t col) {
        img.plot(x, y, col);
        return *this;
    }
};

/**
 * @brief Fluent API for drawing text: text_mono, text_aa, text_centered_mono and text_centered_aa
 *
 * Provides a chainable interface for drawing text using monospace or anti-aliased fonts.
 * @tparam AA The font is an anti-aliased font.
 * @tparam CENTERED The text is centered at its position.
 */
template <typename FONT, typename PF, bool KERNING, text_rotation ROTATION, bool AA, bool CENTERED>
class basic_text {
    static_assert(FONT::mono != AA, "text_mono needs a mono font, text_aa an antialiased font.");
    surface<PF> img;
    int32_t x, y;
    std::string_view str;

 public:
    /**
     * @brief Construct a text shape
     * @param image Target surface to draw on
     * @param x_ X coordinate of text position
     * @param y_ Y coordinate of text position
     * @param str_ Text string to draw
     */
    constexpr basic_text(const surface<PF> &image, int32_t x_, int32_t y_, std::string_view str_)
        : img(image), x(x_), y(y_), str(str_) {
    }

    /**
     * @brief Draw the text and return the width. Not for centered text.
     * @param col Color value
     * @param character_count Maximum number of characters to draw
     * @param character_actual Pointer to store actual characters drawn
     * @return Width of the drawn text in pixels
     */
    constexpr int32_t draw(uint8_t col, size_t character_count = std::numeric_limits<size_t>::max(),
                           size_t *character_actual = nullptr)
        requires(!CENTERED)
    {
        return img.template draw_string<FONT, KERNING, ROTATION>(x, y, str, col, ALIGN_LEFT, character_count,
                                                                 character_actual);
    }

    /**
     * @brief Draw the text with the specified color
     * @param col Color value
     * @return Reference to this text for chaining
     */
    constexpr basic_text &color(uint8_t col) {
        img.template draw_string<FONT, KERNING, ROTATION>(x, y, str, col, CENTERED ? ALIGN_CENTER : ALIGN_LEFT);
        return *this;
    }
};

/// @brief Fluent API for drawing monospace text
template <typename FONT, typename PF, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
using text_mono = basic_text<FONT, PF, KERNING, ROTATION, false, false>;
/// @brief Fluent API for drawing anti-aliased text
template <typename FONT, typename PF, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
using text_aa = basic_text<FONT, PF, KERNING, ROTATION, true, false>;
/// @brief Fluent API for drawing centered monospace text
template <typename FONT, typename PF, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
using text_centered_mono = basic_text<FONT, PF, KERNING, ROTATION, false, true>;
/// @brief Fluent API for drawing centered anti-aliased text
template <typename FONT, typename PF, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
using text_centered_aa = basic_text<FONT, PF, KERNING, ROTATION, true, true>;

}  // namespace shapes

/// @cond DOXYGEN_EXCLUDE
namespace hidden {

// Runtime palette pointer of a surface, only stored for formats which can use one.
template <typename PF>
struct palette_slot {
    constexpr explicit palette_slot(const runtime_palette * /*palette*/) {
    }
    [[nodiscard]] static constexpr const runtime_palette *get() {
        return nullptr;
    }
};

template <has_runtime_palette PF>
struct palette_slot<PF> {
    constexpr explicit palette_slot(const runtime_palette *palette) : slot(palette) {
    }
    [[nodiscard]] constexpr const runtime_palette *get() const {
        return slot;
    }
    const runtime_palette *slot;
};

// Runtime palette storage of image and dynamic_image, same condition. Empty for the other formats.
template <typename PF>
struct palette_store {
    [[nodiscard]] static constexpr const runtime_palette *get() {
        return nullptr;
    }
};

template <has_runtime_palette PF>
struct palette_store<PF> {
    [[nodiscard]] constexpr const runtime_palette *get() const {
        return custom ? &palette : nullptr;
    }
    runtime_palette palette{};
    bool custom{false};
};

}  // namespace hidden
/// @endcond

/**
 * @class basic_surface
 * @brief Non-owning view of a pixel buffer with the read-only operations. Use the aliases: const_surface for
 * read-only pixels, surface for drawing.
 *
 * @tparam PF Pixel class of the format: pixel_1bit, pixel_2bit<GRAYSCALE>, pixel_4bit<GRAYSCALE>,
 * pixel_8bit<GRAYSCALE>, pixel_24bit<GRAYSCALE> or pixel_32bit<GRAYSCALE>. image<...>::pixel_type names it.
 * @tparam BYTE uint8_t or const uint8_t.
 */
template <typename PF, typename BYTE>
class basic_surface : private hidden::palette_slot<PF> {
 public:
    /// @brief Pixel class of the format.
    using pixel_type = PF;

    /**
     * \brief Creates a view of existing pixel data.
     * \param data First byte of the first line.
     * \param width Width in pixels.
     * \param height Height in pixels.
     * \param stride Bytes from one line to the next, at least PF::bytes_per_line(width).
     * \param palette Optional runtime palette. Ignored by formats without runtime palette support.
     */
    constexpr basic_surface(BYTE *data, size_t width, size_t height, size_t stride,
                            const runtime_palette *palette = nullptr)
        : hidden::palette_slot<PF>(palette), pix{data, width, height, stride} {
    }

    /**
     * \brief Creates a read-only view from a drawing surface.
     * \param other Surface to view.
     */
    template <typename OTHER>
        requires(std::is_const_v<BYTE> && !std::is_const_v<OTHER>)
    constexpr basic_surface(const basic_surface<PF, OTHER> &other)
        : hidden::palette_slot<PF>(other.palette()),
          pix{other.data(), static_cast<size_t>(other.width()), static_cast<size_t>(other.height()), other.stride()} {
    }

    /**
     * \brief Bit depth of the pixel format.
     * \return Bit depth. 1 == 2 colors, 2 == 4 colors, 4 == 16 colors, 8 == 256 colors.
     */
    [[nodiscard]] static constexpr size_t bit_depth() {
        return PF::bits_per_pixel;
    }

    /**
     * \brief Width in pixels.
     * \return Width in pixels.
     */
    [[nodiscard]] constexpr int32_t width() const {
        return static_cast<int32_t>(pix.width);
    }

    /**
     * \brief Height in pixels.
     * \return Height in pixels.
     */
    [[nodiscard]] constexpr int32_t height() const {
        return static_cast<int32_t>(pix.height);
    }

    /**
     * \brief Bytes from one line to the next.
     * \return Stride in bytes.
     */
    [[nodiscard]] constexpr size_t stride() const {
        return pix.stride;
    }

    /**
     * \brief First byte of the first line.
     * \return Pointer to the pixel data.
     */
    [[nodiscard]] constexpr BYTE *data() const {
        return pix.data;
    }

    /**
     * \brief Runtime palette in use.
     * \return Runtime palette, or nullptr if the palette of the format is used.
     */
    [[nodiscard]] constexpr const runtime_palette *palette() const {
        return hidden::palette_slot<PF>::get();
    }

    /**
     * \brief Return closest match in the color palette based on the supplied red, green and blue values.
     * \param r Red value (0-255)
     * \param g Green value (0-255)
     * \param b Blue value (0-255)
     * \return The closest matching color palette index.
     */
    [[nodiscard]] constexpr uint8_t get_nearest_color(uint8_t r, uint8_t g, uint8_t b) const {
        return PF::quant.nearest(r, g, b);
    }

    /**
     * \brief Convert this surface to an RGBA8 buffer.
     * \param dst RGBA8 buffer of width() * height() entries.
     */
    constexpr void RGBA_uint32(uint32_t *dst) const {
        PF::RGBA_uint32(dst, pix);
    }

    /**
     * \brief Convert this surface to an RGBA8 buffer.
     * \param dst RGBA8 buffer of width() * height() * 4 bytes.
     */
    constexpr void RGBA_uint8(uint8_t *dst) const {
        PF::RGBA_uint8(dst, pix);
    }

    /**
     * \brief Convert this surface to an RGBA8 buffer, size-checked.
     * \param dst RGBA8 buffer of at least width() * height() entries.
     * \return false if dst is too small. Nothing is written then.
     */
    [[nodiscard]] constexpr bool RGBA_uint32(std::span<uint32_t> dst) const {
        if (dst.size() < pix.width * pix.height) {
            return false;
        }
        PF::RGBA_uint32(dst.data(), pix);
        return true;
    }

    /**
     * \brief Convert this surface to an RGBA8 buffer, size-checked.
     * \param dst RGBA8 buffer of at least width() * height() * 4 bytes.
     * \return false if dst is too small. Nothing is written then.
     */
    [[nodiscard]] constexpr bool RGBA_uint8(std::span<uint8_t> dst) const {
        if (dst.size() < pix.width * pix.height * 4) {
            return false;
        }
        PF::RGBA_uint8(dst.data(), pix);
        return true;
    }

    /**
     * \brief Write a transposed version of this surface. Data is copied.
     * \tparam FLIP_H Flip image horizontally.
     * \tparam FLIP_V Flip image vertically.
     * \param dst Pixel data of a height() x width() image of the same format, without line padding.
     */
    template <bool FLIP_H = false, bool FLIP_V = false>
    constexpr void transpose(uint8_t *dst) const {
        PF::template transpose<FLIP_H, FLIP_V>(pix, dst);
    }

    /**
     * \brief Convert the current instance into a png image. Typically an implementation looks like this:
     *
     * \code{.cpp}
     * auto some_container;
     * image.png([=](char ch) mutable {
     *    some_container.push_back(ch);
     * });
     * \endcode
     *
     * \param char_out A lambda function which consumes the png image data one byte at a time.
     */
    template <typename F>
    constexpr void png(F &&char_out) const {
        PF::png(pix, std::forward<F>(char_out));
    }

    /**
     * \brief Convert the current instance into a sixel stream. Typically an implementation looks like this:
     *
     * \code{.cpp}
     * auto some_container;
     * image.sixel([=](char ch) mutable {
     *    some_container.push_back(ch);
     * });
     * \endcode
     *
     * \tparam S scale of sixel output.
     * \param char_out A lambda function which consumes the sixel stream data one byte at a time.
     */
    template <size_t S = 1, typename F>
    constexpr void sixel(F &&char_out) const {
        sixel<S>(std::forward<F>(char_out), {0, 0, width(), height()});
    }

    /**
     * \brief Convert the current instance into a sixel stream. Typically an implementation looks like this:
     *
     * \code{.cpp}
     * auto some_container;
     * image.sixel({.x=0,.y=0,.w=100,.h=100},[=](char ch) mutable {
     *    some_container.push_back(ch);
     * });
     * \endcode
     *
     * \tparam S scale of sixel output.
     * \param char_out A lambda function which consumes the sixel stream data one byte at a time
     * \param rect Clipping rectangle; to only show a portion of the image.
     */
    template <size_t S = 1, typename F>
    constexpr void sixel(F &&char_out, const constixel::rect<int32_t> &rect) const {
        if constexpr (hidden::has_runtime_palette<PF>) {
            PF::template sixel<S>(pix, std::forward<F>(char_out), rect, palette());
        } else if constexpr (hidden::has_sixel_cache<PF>) {
            std::array<uint8_t, PF::sixel_cache_size(PF::sixel_cache_columns)> line_cache{};
            PF::template sixel<S>(pix, std::forward<F>(char_out), rect, std::span<uint8_t>{line_cache});
        } else {
            PF::template sixel<S>(pix, std::forward<F>(char_out), rect);
        }
    }

    /**
     * \brief Convert the current instance into a sixel stream, with a caller-supplied line cache. Only for the 24
     * and 32-bit formats. sixel() without a cache uses a window of 1024 columns on the stack, which is slow for
     * wider surfaces; a cache of the full width is not.
     *
     * \tparam S scale of sixel output.
     * \param char_out A lambda function which consumes the sixel stream data one byte at a time
     * \param rect Clipping rectangle; to only show a portion of the image.
     * \param line_cache Scratch memory, PF::sixel_cache_size(width()) bytes for full speed. It is overwritten.
     */
    template <size_t S = 1, typename F>
    constexpr void sixel(F &&char_out, const constixel::rect<int32_t> &rect, std::span<uint8_t> line_cache) const
        requires hidden::has_sixel_cache<PF>
    {
        if (line_cache.size() < PF::sixel_cache_size(1)) {
            sixel<S>(std::forward<F>(char_out), rect);
            return;
        }
        std::fill(line_cache.begin(), line_cache.end(), uint8_t{0});
        PF::template sixel<S>(pix, std::forward<F>(char_out), rect, line_cache);
    }

#ifdef CONSTIXEL_ENABLE_COUT
    /**
     * \brief Convert the current instance into a sixel stream and output it to std::cout.
     *
     * @tparam S scale of sixel output.
     */
    template <size_t S = 1>
    void sixel_to_cout() const {
        std::string out;
        sixel<S>([&out](char ch) mutable {
            out.push_back(ch);
        });
        std::cout << out << '\n';
    }

    /**
     * \brief Send a escape command to std::cout to clear the screen and scroll buffer of a vt100 compatible terminal.
     */
    void vt100_clear() const {
        std::cout << "\033[2J\033[3J\033[H";
    }

    /**
     * \brief Send a escape command to std::cout to clear the screen and scroll buffer of a vt100 compatible terminal.
     */
    void vt100_clear_scrollback() const {
        std::cout << "\033[3J";
    }

    /**
     * \brief Send a escape command to std::cout to home the cursor of a vt100 compatible terminal.
     */
    void vt100_home() const {
        std::cout << "\033[H";
    }

    /**
     * \brief Convert the current instance into a png and display it in iTerm.
     */
    void png_to_iterm() const {
        std::string output;
        output.append("\033]1337;File=inline=1:");
        append_png_as_base64(output);
        output.append("\07");
        std::cout << output << '\n';
    }

    /**
     * \brief Convert the current instance into a png and display it in a terminal with kitty graphics support.
     */
    void png_to_kitty() const {
        std::string base64{};
        std::string output{};
        append_png_as_base64(base64);
        bool first = true;
        for (; !base64.empty();) {
            if (first) {
                first = false;
                output.append("\033_Ga=T,f=100,");
            } else {
                output.append("\033_G");
            }
            output.append(base64.length() <= 4096 ? "m=0;" : "m=1;");
            const size_t bytes_to_append = std::min(base64.length(), static_cast<size_t>(4096));
            output.append(base64.substr(0, bytes_to_append));
            base64.erase(0, bytes_to_append);
            output.append("\033\\");
        }
        std::cout << output << '\n';
    }
#endif  // #ifdef CONSTIXEL_ENABLE_COUT
    /**
     * \brief Convert the current instance into a byte stream formatted for embedded displays.
     * \tparam dst_format The desired data format.
     * \param uint8_out A lambda function which consumes the data stream one byte at a time.
     */
    template <device_format dst_format, typename F>
    constexpr void convert(F &&uint8_out) const {
        convert<dst_format>(std::forward<F>(uint8_out), {0, 0, width(), height()});
    }

    /**
     * \brief Convert a rectangle of the current instance into a byte stream formatted for embedded displays, for
     * partial display updates. The stream holds the lines of the rectangle only.
     * \tparam dst_format The desired data format.
     * \param uint8_out A lambda function which consumes the data stream one byte at a time.
     * \param rect Rectangle to convert. It is clipped to the surface.
     */
    template <device_format dst_format, typename F>
    constexpr void convert(F &&uint8_out, const constixel::rect<int32_t> &rect) const {
        const convert_area area = convert_layout<dst_format>(rect);
        for (size_t line = 0; line < area.lines; line++) {
            for (size_t offset = 0; offset < area.line_bytes; offset++) {
                std::forward<F>(uint8_out)(convert_byte<dst_format>(area, line, offset));
            }
        }
    }

    /**
     * \brief Size of the byte stream convert() and convert_chunk() produce for a rectangle.
     * \tparam dst_format The desired data format.
     * \param rect Rectangle to convert. It is clipped to the surface.
     * \return Size in bytes.
     */
    template <device_format dst_format>
    [[nodiscard]] constexpr size_t convert_size(const constixel::rect<int32_t> &rect) const {
        const convert_area area = convert_layout<dst_format>(rect);
        return area.lines * area.line_bytes;
    }

    /**
     * \brief Size of the byte stream convert() and convert_chunk() produce for the whole surface.
     * \tparam dst_format The desired data format.
     * \return Size in bytes.
     */
    template <device_format dst_format>
    [[nodiscard]] constexpr size_t convert_size() const {
        return convert_size<dst_format>({0, 0, width(), height()});
    }

    /**
     * \brief Convert a rectangle of the current instance in chunks, the same byte stream as convert(). Call it
     * with the same cursor and rectangle until it returns 0. Example, a partial update over SPI:
     *
     * \code{.cpp}
     * std::array<uint8_t, 64> chunk;
     * constixel::convert_cursor cursor;
     * while (size_t n = image.convert_chunk<constixel::RGB565_8BIT_SERIAL>(chunk, cursor, {0, 0, 32, 16})) {
     *     spi_write(chunk.data(), n);
     * }
     * \endcode
     *
     * \tparam dst_format The desired data format.
     * \param dst Buffer to fill.
     * \param cursor Position in the stream. Start with convert_cursor{}; it is advanced by the bytes written.
     * \param rect Rectangle to convert. It is clipped to the surface.
     * \return Bytes written to dst. Less than dst.size() at the end of the stream, 0 when it is finished.
     */
    template <device_format dst_format>
    [[nodiscard]] constexpr size_t convert_chunk(std::span<uint8_t> dst, convert_cursor &cursor,
                                                 const constixel::rect<int32_t> &rect) const {
        const convert_area area = convert_layout<dst_format>(rect);
        size_t count = 0;
        while (count < dst.size() && cursor.line < area.lines && cursor.offset < area.line_bytes) {
            dst[count++] = convert_byte<dst_format>(area, cursor.line, cursor.offset);
            if (++cursor.offset >= area.line_bytes) {
                cursor.offset = 0;
                cursor.line++;
            }
        }
        return count;
    }

    /**
     * \brief Convert the current instance in chunks. See convert_chunk(std::span<uint8_t>, convert_cursor &, const
     * constixel::rect<int32_t> &).
     * \tparam dst_format The desired data format.
     * \param dst Buffer to fill.
     * \param cursor Position in the stream. Start with convert_cursor{}; it is advanced by the bytes written.
     * \return Bytes written to dst. Less than dst.size() at the end of the stream, 0 when it is finished.
     */
    template <device_format dst_format>
    [[nodiscard]] constexpr size_t convert_chunk(std::span<uint8_t> dst, convert_cursor &cursor) const {
        return convert_chunk<dst_format>(dst, cursor, {0, 0, width(), height()});
    }

 protected:
    /// @cond DOXYGEN_EXCLUDE
    hidden::raster<BYTE> pix;
    /// @endcond

 private:
    /// @cond DOXYGEN_EXCLUDE
    // Clipped source rectangle of a conversion and the shape of its byte stream.
    struct convert_area {
        size_t x = 0;
        size_t y = 0;
        size_t w = 0;
        size_t h = 0;
        size_t lines = 0;
        size_t line_bytes = 0;
    };

    template <device_format dst_format>
    [[nodiscard]] constexpr convert_area convert_layout(const constixel::rect<int32_t> &rect) const {
        const auto clipped = rect & constixel::rect<int32_t>{0, 0, width(), height()};
        convert_area area{static_cast<size_t>(clipped.x), static_cast<size_t>(clipped.y),
                          static_cast<size_t>(clipped.w), static_cast<size_t>(clipped.h)};
        area.lines = area.h;
        if constexpr (dst_format == STRAIGHT_THROUGH) {
            // the bytes which hold the pixels x .. x + w - 1
            area.line_bytes = area.w == 0 ? 0 : PF::bytes_per_line(area.x + area.w) - convert_first_byte(area.x);
        } else if constexpr (dst_format == RGB565_8BIT_SERIAL) {
            area.line_bytes = area.w * 2;
        } else if constexpr (dst_format == RGB666_8BIT_SERIAL_1 || dst_format == RGB666_8BIT_SERIAL_2) {
            area.line_bytes = area.w * 3;
        } else if constexpr (dst_format == MONO_1BIT_VERTICAL_PAGES) {
            area.lines = (area.h + 7) / 8;
            area.line_bytes = area.w;
        } else if constexpr (dst_format == MONO_1BIT_HORIZONTAL_LSB) {
            area.line_bytes = (area.w + 7) / 8;
        }
        return area;
    }

    [[nodiscard]] static constexpr size_t convert_first_byte(size_t x) {
        return (x * PF::bits_per_pixel) / 8;
    }

    // One byte of the stream. convert() and convert_chunk() both read the stream through this.
    template <device_format dst_format>
    [[nodiscard]] constexpr uint8_t convert_byte(const convert_area &area, size_t line, size_t offset) const {
        if constexpr (dst_format == STRAIGHT_THROUGH) {
            return pix.data[(area.y + line) * pix.stride + convert_first_byte(area.x) + offset];
        } else if constexpr (dst_format == RGB565_8BIT_SERIAL) {
            const uint32_t col = pixel_rgb(area.x + offset / 2, area.y + line);
            if ((offset & 1) == 0) {
                return static_cast<uint8_t>(((((col >> 0) & 0xff) >> 3) << 3) | ((((col >> 8) & 0xff) >> 2) >> 3));
            }
            return static_cast<uint8_t>((((((col >> 8) & 0xff) >> 2) & 0x7) << 5) | (((col >> 16) & 0xff) >> 3));
        } else if constexpr (dst_format == RGB666_8BIT_SERIAL_1) {
            const uint32_t col = pixel_rgb(area.x + offset / 3, area.y + line);
            return static_cast<uint8_t>(((col >> (8 * (offset % 3))) & 0xff) >> 2);
        } else if constexpr (dst_format == RGB666_8BIT_SERIAL_2) {
            const uint32_t col = pixel_rgb(area.x + offset / 3, area.y + line);
            return static_cast<uint8_t>((((col >> (8 * (offset % 3))) & 0xff) >> 2) << 2);
        } else if constexpr (dst_format == MONO_1BIT_VERTICAL_PAGES) {
            static_assert(PF::bits_per_pixel == 1, "MONO_1BIT_VERTICAL_PAGES needs a 1-bit format.");
            uint8_t out = 0;
            for (size_t bit = 0; bit < 8 && line * 8 + bit < area.h; bit++) {
                const uint8_t *row = pix.data + (area.y + line * 8 + bit) * pix.stride;
                out |= static_cast<uint8_t>((PF::get_col(row, area.x + offset) & 1) << bit);
            }
            return out;
        } else if constexpr (dst_format == MONO_1BIT_HORIZONTAL_LSB) {
            static_assert(PF::bits_per_pixel == 1, "MONO_1BIT_HORIZONTAL_LSB needs a 1-bit format.");
            const uint8_t *row = pix.data + (area.y + line) * pix.stride;
            uint8_t out = 0;
            for (size_t bit = 0; bit < 8 && offset * 8 + bit < area.w; bit++) {
                out |= static_cast<uint8_t>((PF::get_col(row, area.x + offset * 8 + bit) & 1) << bit);
            }
            return out;
        } else {
            return 0;
        }
    }

    // Color of a pixel as 0x00BBGGRR.
    [[nodiscard]] constexpr uint32_t pixel_rgb(size_t x, size_t y) const {
        const uint8_t *row = pix.data + y * pix.stride;
        if constexpr (PF::bits_per_pixel >= 24) {
            const uint8_t *ptr = row + x * (PF::bits_per_pixel / 8);
            return static_cast<uint32_t>(ptr[0]) | (static_cast<uint32_t>(ptr[1]) << 8) |
                   (static_cast<uint32_t>(ptr[2]) << 16);
        } else if constexpr (PF::bits_per_pixel == 8) {
            if (palette() != nullptr) {
                return palette()->at(row[x]);
            }
            return PF::quant.palette()[row[x]];
        } else {
            return PF::quant.palette()[PF::get_col(row, x)];
        }
    }
    /// @endcond

#ifdef CONSTIXEL_ENABLE_COUT
    /**
     * @private
     */
    void append_png_as_base64(std::string &output) const {
        size_t buffer = 0;
        size_t bits_collected = 0;
        PF::png(pix, [&buffer, &bits_collected, &output](char byte) mutable {
            static constexpr const char *base64_chars =
                "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
            buffer = (buffer << 8) | static_cast<uint8_t>(byte);
            bits_collected += 8;
            while (bits_collected >= 6) {
                bits_collected -= 6;
                output.push_back(base64_chars[(buffer >> bits_collected) & 0x3F]);
            }
            return [&output]() mutable {
                if (output.capacity() == 0) {
                    return;
                }
                if ((output.size() % 4) != 0) {
                    const size_t padding = 4 - output.size() % 4;
                    for (size_t c = 0; c < padding; c++) {
                        output.push_back('=');
                    }
                }
            };
        });
    }
#endif  // #ifdef CONSTIXEL_ENABLE_COUT
};

/**
 * @brief Read-only surface. See basic_surface.
 */
template <typename PF>
using const_surface = basic_surface<PF, const uint8_t>;

/**
 * @class surface
 * @brief Non-owning drawing surface: all drawing, text and output code of constixel works on this view of a pixel
 * buffer. image and dynamic_image forward to it; it can also wrap an external frame buffer, including one with
 * padded lines:
 *
 * \code{.cpp}
 * using pixel = constixel::pixel_8bit<false>;
 * constixel::surface<pixel> fb(memory, 320, 240, 512);
 * fb.fill_rect(0, 0, 10, 10, constixel::color::WHITE);
 * \endcode
 *
 * @tparam PF Pixel class of the format, see basic_surface.
 */
template <typename PF>
class surface : public basic_surface<PF, uint8_t> {
    using basic_surface<PF, uint8_t>::pix;

#if defined(__x86_64__) || defined(_M_X64) || defined(__aarch64__)
    static constexpr int32_t min_coord = -int32_t{1 << 28};
    static constexpr int32_t max_coord = +int32_t{1 << 28};
    using calc_square_type = int64_t;
#else   // #if defined(__x6_64__) || defined(_M_X64) || defined(__aarch64__)
    static constexpr int32_t min_coord = -int32_t{1 << 13} + int32_t{1};
    static constexpr int32_t max_coord = +int32_t{1 << 13} - int32_t{1};
    using calc_square_type = int32_t;
#endif  // #if defined(__x86_64__) || defined(_M_X64) || defined(__aarch64__)

 public:
    using basic_surface<PF, uint8_t>::basic_surface;

    surface() = delete;

    /**
     * \brief Plot a single pixel at the specified coordinates using the supplied color.
     * \param x X-coordinate in pixels.
     * \param y Y-coordinate in pixels.
     * \param col Color palette index to use.
     */
    constexpr void plot(int32_t x, int32_t y, uint8_t col) {
        plot_clipped(pix, x, y, col);
    }

    /**
     * \brief Plot a single pixel at the specified coordinates using the supplied color.
     *
     * \param p \ref plot initializer struct
     */
    constexpr void plot(const struct plot &p) {
        plot(p.x, p.y, p.col);
    }

    /**
     * \brief Draw a line with the specified color and thickness. Example:
     *
     * \code{.cpp}
     * image.draw_line(0, 0, 200, 100, constixel::color::WHITE, 2);
     * \endcode
     *
     * \param x0 Starting X-coordinate in pixels.
     * \param y0 Starting Y-coordinate in pixels.
     * \param x1 Ending X-coordinate in pixels.
     * \param y1 Ending Y-coordinate in pixels.
     * \param col Color palette index to use.
     * \param stroke_width Width of the stroke in pixels.
     */
    constexpr void draw_line(int32_t x0, int32_t y0, int32_t x1, int32_t y1, uint8_t col, int32_t stroke_width = 1) {
        const auto dst = pix;
        auto minmax_check = std::minmax({x0, y0, x1, y1, stroke_width});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }

        if (!clip_line(x0, y0, x1, y1, -stroke_width, -stroke_width, static_cast<int32_t>(pix.width) + stroke_width,
                       static_cast<int32_t>(pix.height) + stroke_width)) {
            return;
        }

        if (stroke_width <= 0) {
            return;
        }

        if (stroke_width == 1) {
            bool steep = abs(y1 - y0) > abs(x1 - x0);
            if (steep) {
                std::swap(x0, y0);
                std::swap(x1, y1);
            }
            if (x0 > x1) {
                std::swap(x0, x1);
                std::swap(y0, y1);
            }

            const int32_t dx = x1 - x0;
            const int32_t dy = abs(y1 - y0);
            int32_t err = dx / 2;
            int32_t ystep = 0;
            if (y0 < y1) {
                ystep = 1;
            } else {
                ystep = -1;
            }

            for (; x0 <= x1; x0++) {
                if (steep) {
                    plot_clipped(dst, y0, x0, col);
                } else {
                    plot_clipped(dst, x0, y0, col);
                }
                err -= dy;
                if (err < 0) {
                    y0 += ystep;
                    err += dx;
                }
            }
            return;
        }

        auto line_thick = [&]<typename I>() {
            const int32_t half_width = stroke_width / 2;

            const int32_t margin = half_width + 1;
            const int32_t min_x = std::max(int32_t{0}, std::min({x0, x1}) - margin);
            const int32_t max_x = std::min(static_cast<int32_t>(pix.width) - 1, std::max({x0, x1}) + margin);
            const int32_t min_y = std::max(int32_t{0}, std::min({y0, y1}) - margin);
            const int32_t max_y = std::min(static_cast<int32_t>(pix.height) - 1, std::max({y0, y1}) + margin);

            if (min_x > max_x || min_y > max_y) {
                return;
            }

            const I line_dx = x1 - x0;
            const I line_dy = y1 - y0;
            const I line_length_sq = line_dx * line_dx + line_dy * line_dy;

            if (line_length_sq <= 1) {
                if (x0 >= 0 && x0 < static_cast<int32_t>(pix.width) && y0 >= 0 &&
                    y0 < static_cast<int32_t>(pix.height)) {
                    fill_circle(x0, y0, half_width, col);
                }
                return;
            }

            for (int32_t py = min_y; py <= max_y; py++) {
                for (int32_t px = min_x; px <= max_x; px++) {
                    const I px_dx = px - x0;
                    const I px_dy = py - y0;

                    const I dot_product = px_dx * line_dx + px_dy * line_dy;
                    const I t_scaled = std::max(I{0}, std::min(line_length_sq, dot_product));

                    const I closest_x = x0 + (t_scaled * line_dx) / line_length_sq;
                    const I closest_y = y0 + (t_scaled * line_dy) / line_length_sq;

                    const I dist_x = px - closest_x;
                    const I dist_y = py - closest_y;
                    const I distance_sq = dist_x * dist_x + dist_y * dist_y;

                    if (distance_sq <= static_cast<I>(half_width) * static_cast<I>(half_width)) {
                        plot_clipped(pix, px, py, col);
                    }
                }
            }
        };
        line_thick.template operator()<calc_square_type>();
    }

    /**
     * \brief Draw a line with the specified color and thickness. Example:
     *
     * \code{.cpp}
     * image.draw_line({.x0=0, .y0=0, .x1=200, .y1=100, .col=color::WHITE, .sw=2});
     * \endcode
     *
     * \param d \ref draw_line initializer struct
     */
    constexpr void draw_line(const struct draw_line &d) {
        draw_line(d.x0, d.y0, d.x1, d.y1, d.col, static_cast<int32_t>(d.sw));
    }

    /**
     * \brief Draw an antialiased line with variable stroke width. Only format_8bit targets are supported.
     * Example:
     *
     * \code{.cpp}
     * image.draw_line_aa(0, 0, 200, 100, constixel::color::WHITE, 3.0f);
     * \endcode
     *
     * \param x0 Starting X-coordinate in pixels.
     * \param y0 Starting Y-coordinate in pixels.
     * \param x1 Ending X-coordinate in pixels.
     * \param y1 Ending Y-coordinate in pixels.
     * \param col Color palette index to use.
     * \param stroke_width Width of the line in pixels (can be fractional).
     */
    constexpr void draw_line_aa(int32_t x0, int32_t y0, int32_t x1, int32_t y1, uint8_t col,
                                float stroke_width = 1.0f) {
        const auto dst = pix;
        auto minmax_check = std::minmax({x0, y0, x1, y1});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        if (!std::isnormal(stroke_width) || stroke_width < static_cast<float>(min_coord) ||
            stroke_width > static_cast<float>(max_coord)) {
            return;
        }

        const auto [Rl, Gl, Bl] = linear_color(col);

        if (stroke_width <= 0.0f) {
            return;
        }

        if (stroke_width <= 1.0f) {
            if (!clip_line(x0, y0, x1, y1, 0, 0, static_cast<int32_t>(pix.width), static_cast<int32_t>(pix.height))) {
                return;
            }

            bool steep = abs(y1 - y0) > abs(x1 - x0);
            if (steep) {
                std::swap(x0, y0);
                std::swap(x1, y1);
            }
            if (x0 > x1) {
                std::swap(x0, x1);
                std::swap(y0, y1);
            }

            const auto dx = static_cast<float>(x1 - x0);
            const auto dy = static_cast<float>(y1 - y0);
            const float gradient = dx == 0.0f ? 1.0f : dy / dx;

            auto color_compose = [&](int32_t x, int32_t y, float a) {
                if (a < hidden::epsilon_low) {
                    return;
                }
                if (a >= hidden::epsilon_high) {
                    plot_clipped(pix, x, y, col);
                } else {
                    compose(pix, x, y, a, Rl, Gl, Bl);
                }
            };

            auto ipart = [](float x) {
                return std::floor(x);
            };
            auto fpart = [](float x) {
                return x - std::floor(x);
            };
            auto rfpart = [&](float x) {
                return 1.0f - fpart(x);
            };

            // first endpoint
            auto xend = static_cast<float>(x0);
            float yend = static_cast<float>(y0) + gradient * (xend - static_cast<float>(x0));
            float xgap = rfpart(static_cast<float>(x0) + 0.5f);
            const auto xpxl1 = static_cast<int32_t>(xend);
            const auto ypxl1 = static_cast<int32_t>(ipart(yend));
            if (steep) {
                color_compose(ypxl1, xpxl1, rfpart(yend) * xgap);
                color_compose(ypxl1 + 1, xpxl1, fpart(yend) * xgap);
            } else {
                color_compose(xpxl1, ypxl1, rfpart(yend) * xgap);
                color_compose(xpxl1, ypxl1 + 1, fpart(yend) * xgap);
            }
            float intery = yend + gradient;

            xend = static_cast<float>(x1);
            yend = static_cast<float>(y1) + gradient * (xend - static_cast<float>(x1));
            xgap = fpart(static_cast<float>(x1) + 0.5f);
            const auto xpxl2 = static_cast<int32_t>(xend);
            const auto ypxl2 = static_cast<int32_t>(ipart(yend));
            if (steep) {
                color_compose(ypxl2, xpxl2, rfpart(yend) * xgap);
                color_compose(ypxl2 + 1, xpxl2, fpart(yend) * xgap);
            } else {
                color_compose(xpxl2, ypxl2, rfpart(yend) * xgap);
                color_compose(xpxl2, ypxl2 + 1, fpart(yend) * xgap);
            }

            for (int32_t x = xpxl1 + 1; x < xpxl2; ++x) {
                auto y = static_cast<int32_t>(ipart(intery));
                if (steep) {
                    color_compose(y, x, rfpart(intery));
                    color_compose(y + 1, x, fpart(intery));
                } else {
                    color_compose(x, y, rfpart(intery));
                    color_compose(x, y + 1, fpart(intery));
                }
                intery += gradient;
            }
            return;
        }

        const float half_width = stroke_width * 0.5f;

        const auto margin = static_cast<int32_t>(std::ceil(half_width + 1.0f));
        const int32_t min_x = std::max(int32_t{0}, std::min({x0, x1}) - margin);
        const int32_t max_x = std::min(static_cast<int32_t>(pix.width) - 1, std::max({x0, x1}) + margin);
        const int32_t min_y = std::max(int32_t{0}, std::min({y0, y1}) - margin);
        const int32_t max_y = std::min(static_cast<int32_t>(pix.height) - 1, std::max({y0, y1}) + margin);

        if (min_x > max_x || min_y > max_y) {
            return;
        }

        const auto dx = static_cast<float>(x1 - x0);
        const auto dy = static_cast<float>(y1 - y0);
        const float line_length_sq = dx * dx + dy * dy;

        if (line_length_sq <= 1.0f) {
            const float radius = half_width;
            const auto r_ceil = static_cast<int32_t>(std::ceil(radius));
            for (int32_t py = y0 - r_ceil; py <= y0 + r_ceil; py++) {
                for (int32_t px = x0 - r_ceil; px <= x0 + r_ceil; px++) {
                    if (px >= 0 && px < static_cast<int32_t>(pix.width) && py >= 0 &&
                        py < static_cast<int32_t>(pix.height)) {
                        auto dist = static_cast<float>((px - x0) * (px - x0) + (py - y0) * (py - y0));
                        if (std::is_constant_evaluated()) {
                            dist = hidden::fast_sqrtf(dist);
                        } else {
                            dist = std::sqrt(dist);
                        }
                        const float coverage = std::max(0.0f, std::min(1.0f, radius + 0.5f - dist));
                        if (coverage > hidden::epsilon_low) {
                            if (coverage >= hidden::epsilon_high) {
                                plot_clipped(dst, px, py, col);
                            } else {
                                compose(dst, px, py, coverage, Rl, Gl, Bl);
                            }
                        }
                    }
                }
            }
            return;
        }

        for (int32_t py = min_y; py <= max_y; py++) {
            for (int32_t px = min_x; px <= max_x; px++) {
                const auto px_dx = static_cast<float>(px - x0);
                const auto px_dy = static_cast<float>(py - y0);
                const float t = std::max(0.0f, std::min(1.0f, (px_dx * dx + px_dy * dy) / line_length_sq));
                const float closest_x = static_cast<float>(x0) + t * dx;
                const float closest_y = static_cast<float>(y0) + t * dy;
                const float dist_x = static_cast<float>(px) - closest_x;
                const float dist_y = static_cast<float>(py) - closest_y;
                float dist = dist_x * dist_x + dist_y * dist_y;
                if (std::is_constant_evaluated()) {
                    dist = hidden::fast_sqrtf(dist);
                } else {
                    dist = std::sqrt(dist);
                }
                const float coverage = std::max(0.0f, std::min(1.0f, half_width + 0.5f - dist));
                if (coverage > hidden::epsilon_low) {
                    if (coverage >= hidden::epsilon_high) {
                        plot_clipped(dst, px, py, col);
                    } else {
                        compose(dst, px, py, coverage, Rl, Gl, Bl);
                    }
                }
            }
        }
    }

    /**
     * \brief Draw a 1-pixel wide antialiased line with the specified color. Only format_8bit targets are supported.
     * Example:
     *
     * \code{.cpp}
     * image.draw_line_aa({.x0=0, .y0=0, .x1=200, .y1=100, .col=color::WHITE, .sw=2.5f});
     * \endcode
     *
     * \param d \ref draw_line initializer struct
     */
    constexpr void draw_line_aa(const struct draw_line &d) {
        draw_line_aa(d.x0, d.y0, d.x1, d.y1, d.col, d.sw);
    }

    /**
     * \brief Return the width of a string using the specified font in the template parameter. Typical use:
     *
     * \code{.cpp}
     * #include "some_font_aa.h"
     * ...
     *     int32_t width = image.string_width<constixel::some_font_aa>();
     * ...
     * \endcode
     *
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available.
     * \param str UTF-8 string. Does not have to be NUL-terminated.
     * \param character_count How many utf32 characters in the string should be measured.
     * \param character_actual How many utf32 characters in the string were measured.
     * \return Width of the string in pixels.
     */
    template <typename FONT, bool KERNING = false>
    [[nodiscard]] constexpr int32_t string_width(std::string_view str,
                                                 size_t character_count = std::numeric_limits<size_t>::max(),
                                                 size_t *character_actual = nullptr) {
        const char *ptr = str.data();
        const char *end = ptr + str.size();
        int32_t x = 0;
        size_t count = 0;
        while (ptr != end && count++ < character_count) {
            uint32_t utf32 = 0;
            ptr = next_utf32(ptr, end, &utf32);
            size_t index = 0;
            if (lookup_glyph<FONT>(utf32, &index)) {
                const char_info<typename FONT::char_info_type> &ch_info = FONT::char_table.at(index);
                if (ptr == end) {
                    x += static_cast<int32_t>(ch_info.width);
                } else {
                    x += static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        x += get_kerning<FONT>(utf32, ptr, end);
                    }
                }
            }
        }
        if (character_actual) {
            *character_actual = count;
        }
        return x;
    }

    /**
     * \brief Return the width of a NUL-terminated string. See string_width(std::string_view, size_t, size_t *).
     * \param str UTF-8 string.
     * \param character_count How many utf32 characters in the string should be measured.
     * \param character_actual How many utf32 characters in the string were measured.
     * \return Width of the string in pixels.
     */
    template <typename FONT, bool KERNING = false>
    [[nodiscard]] constexpr int32_t string_width(const char *str,
                                                 size_t character_count = std::numeric_limits<size_t>::max(),
                                                 size_t *character_actual = nullptr) {
        return string_width<FONT, KERNING>(std::string_view{str}, character_count, character_actual);
    }

    /**
     * \brief Draw text at the specified coordinate. The font in the template parameter selects the kind of text: a
     * mono font draws pixelized text on every format, an antialiased font draws antialiased text on the formats
     * which can blend (4-bit and more). Typical use:
     *
     * \code{.cpp}
     * #include "some_font_aa.h"
     * ...
     *     image.draw_string<constixel::some_font_aa>(319, 0, "MyText", constixel::color::WHITE, ALIGN_RIGHT);
     * ...
     * \endcode
     *
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available. Default to false.
     * \tparam ROTATION Rotation around the x/y coordinate. Can be text_rotation::DEGREE_0, text_rotation::DEGREE_90,
     * text_rotation::DEGREE_180 or text_rotation::DEGREE_270
     * \param x X-coordinate in pixels.
     * \param y Y-coordinate in pixels.
     * \param str UTF-8 string. Does not have to be NUL-terminated.
     * \param col Color palette index to use.
     * \param align ALIGN_LEFT: the text starts at x/y. ALIGN_CENTER: it is centered there. ALIGN_RIGHT: it ends there.
     * \param character_count How many utf32 characters in the string should be drawn.
     * \param character_actual How many utf32 characters in the string were drawn.
     * \return Returns the new caret position.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string(int32_t x, int32_t y, std::string_view str, uint8_t col,
                                  text_align align = ALIGN_LEFT,
                                  size_t character_count = std::numeric_limits<size_t>::max(),
                                  size_t *character_actual = nullptr) {
        if (align != ALIGN_LEFT) {
            const int32_t width = string_width<FONT, KERNING>(str, character_count);
            const int32_t offset = align == ALIGN_CENTER ? width / 2 : width;
            if constexpr (ROTATION == DEGREE_0) {
                x -= offset;
            } else if constexpr (ROTATION == DEGREE_180) {
                x += offset;
            } else if constexpr (ROTATION == DEGREE_90) {
                y += offset;
            } else if constexpr (ROTATION == DEGREE_270) {
                y -= offset;
            }
        }
        return draw_text<FONT, KERNING, ROTATION>(x, y, str, col, character_count, character_actual);
    }

    /**
     * \brief Draw a NUL-terminated string. See draw_string(int32_t, int32_t, std::string_view, uint8_t, text_align,
     * size_t, size_t *).
     * \param x X-coordinate in pixels.
     * \param y Y-coordinate in pixels.
     * \param str UTF-8 string. nullptr draws nothing.
     * \param col Color palette index to use.
     * \param align Placement relative to x/y.
     * \param character_count How many utf32 characters in the string should be drawn.
     * \param character_actual How many utf32 characters in the string were drawn.
     * \return Returns the new caret position.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string(int32_t x, int32_t y, const char *str, uint8_t col, text_align align = ALIGN_LEFT,
                                  size_t character_count = std::numeric_limits<size_t>::max(),
                                  size_t *character_actual = nullptr) {
        if (str == nullptr) {
            return 0;
        }
        return draw_string<FONT, KERNING, ROTATION>(x, y, std::string_view{str}, col, align, character_count,
                                                    character_actual);
    }

    /**
     * \brief Draw text at the specified coordinate. Example:
     *
     * \code{.cpp}
     * image.draw_string<constixel::some_font_mono>({.x=191, .y=0, .str="MyText", .col=1, .align=ALIGN_RIGHT});
     * \endcode
     *
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available. Default to false.
     * \tparam ROTATION Rotation around the x/y coordinate.
     * \param d \ref draw_string initializer struct
     * \param character_count How many utf32 characters in the string should be drawn.
     * \param character_actual How many utf32 characters in the string were drawn.
     * \return Returns the new caret position.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string(const struct draw_string &d,
                                  size_t character_count = std::numeric_limits<size_t>::max(),
                                  size_t *character_actual = nullptr) {
        return draw_string<FONT, KERNING, ROTATION>(d.x, d.y, d.str, d.col, d.align, character_count,
                                                    character_actual);
    }

    /**
     * \brief Draw text at the specified coordinate. The template parameter selects which mono font to use. Only
     * format_8bit targets are supported.

     * \code{.cpp}
     * #include "some_font_mono.h"
     * ...
     *     image.draw_string_mono<constixel::some_font_mono>(0, 0, "MyText", constixel::color::WHITE);
     * ...
     * \endcode
     *
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available. Default to false.
     * \tparam ROTATION Rotation around the x/y coordinate. Defaults to text_rotation::DEGREE_0. Can be
     * text_rotation::DEGREE_0, text_rotation::DEGREE_90, text_rotation::DEGREE_180 or text_rotation::DEGREE_270
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param str UTF-8 string. Does not have to be NUL-terminated.
     * \param col Color palette index to use.
     * \param character_count How many utf32 characters in the string should be drawn.
     * \param character_actual How many utf32 characters in the string were drawn.
     * \return Returns the new caret X-coordinate position. Pass this value to the next draw_string call to get
     * continious text.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_mono(int32_t x, int32_t y, std::string_view str, uint8_t col,
                                       size_t character_count = std::numeric_limits<size_t>::max(),
                                       size_t *character_actual = nullptr) {
        static_assert(FONT::mono == true, "Can't use an antialiased font to draw mono/pixelized text.");
        return draw_text<FONT, KERNING, ROTATION>(x, y, str, col, character_count, character_actual);
    }

    /**
     * \brief Draw a NUL-terminated string. See draw_string_mono(int32_t, int32_t, std::string_view, uint8_t, size_t,
     * size_t *).
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param str UTF-8 string. nullptr draws nothing.
     * \param col Color palette index to use.
     * \param character_count How many utf32 characters in the string should be drawn.
     * \param character_actual How many utf32 characters in the string were drawn.
     * \return Returns the new caret position.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_mono(int32_t x, int32_t y, const char *str, uint8_t col,
                                       size_t character_count = std::numeric_limits<size_t>::max(),
                                       size_t *character_actual = nullptr) {
        if (str == nullptr) {
            return 0;
        }
        return draw_string_mono<FONT, KERNING, ROTATION>(x, y, std::string_view{str}, col, character_count,
                                                         character_actual);
    }

    /**
     * \brief Draw text at the specified coordinate. The template parameter selects which mono font to use. Only
     * format_8bit targets are supported.

     * \code{.cpp}
     * #include "some_font_mono.h"
     * ...
     *     image.draw_string_mono<constixel::some_font_mono>({.x=0, .y=0, .str="MyText", .col=constixel::color::WHITE});
     * ...
     * \endcode
     *
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available. Default to false.
     * \tparam ROTATION Rotation around the x/y coordinate. Defaults to text_rotation::DEGREE_0. Can be
     * text_rotation::DEGREE_0, text_rotation::DEGREE_90, text_rotation::DEGREE_180 or text_rotation::DEGREE_270
     * \param d \ref draw_string initializer struct.
     * \param character_count How many utf32 characters in the string should be drawn.
     * \param character_actual How many utf32 characters in the string were drawn.
     * \return Returns the new caret X-coordinate position. Pass this value to the next draw_string call to get
     * continious text.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_mono(const struct draw_string &d,
                                       size_t character_count = std::numeric_limits<size_t>::max(),
                                       size_t *character_actual = nullptr) {
        static_assert(FONT::mono == true, "Can't use an antialiased font to draw mono/pixelized text.");
        return draw_string<FONT, KERNING, ROTATION>(d, character_count, character_actual);
    }

    /**
     * \brief Draw text centered at the specified coordinate. The template parameter selects which mono font to use.
     * Only format_8bit targets are supported.
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available. Default to false.
     * \tparam ROTATION Rotation around the x/y coordinate. Can be text_rotation::DEGREE_0, text_rotation::DEGREE_90,
     * text_rotation::DEGREE_180 or text_rotation::DEGREE_270
     * \param x Center/Starting X-coordinate in pixels.
     * \param y Center/Starting Y-coordinate in pixels.
     * \param str UTF-8 string. Does not have to be NUL-terminated.
     * \param col Color palette index to use.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr void draw_string_centered_mono(int32_t x, int32_t y, std::string_view str, uint8_t col) {
        static_assert(FONT::mono == true, "Can't use an antialiased font to draw mono/pixelized text.");
        draw_string<FONT, KERNING, ROTATION>(x, y, str, col, ALIGN_CENTER);
    }

    /**
     * \brief Draw a centered NUL-terminated string. See draw_string_centered_mono(int32_t, int32_t,
     * std::string_view, uint8_t).
     * \param x Center/Starting X-coordinate in pixels.
     * \param y Center/Starting Y-coordinate in pixels.
     * \param str UTF-8 string. nullptr draws nothing.
     * \param col Color palette index to use.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr void draw_string_centered_mono(int32_t x, int32_t y, const char *str, uint8_t col) {
        draw_string_centered_mono<FONT, KERNING, ROTATION>(x, y, hidden::cstr_view(str), col);
    }

    /**
     * \brief Draw antialiased text at the specified coordinate. The template parameter selects which antialiased font
     * to use. Only format_8bit targets are supported. Typical use:
     *
     * \code{.cpp}
     * #include "some_font_aa.h"
     * ...
     *     image.draw_string_aa<constixel::some_font_aa>(0, 0, "MyText", constixel::color::WHITE);
     * ...
     * \endcode
     *
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available. Default to false.
     * \tparam ROTATION Rotation around the x/y coordinate. Can be text_rotation::DEGREE_0, text_rotation::DEGREE_90,
     * text_rotation::DEGREE_180 or text_rotation::DEGREE_270
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param str UTF-8 string. Does not have to be NUL-terminated.
     * \param col Color palette index to use.
     * \param character_count How many utf32 characters in the string should be drawn.
     * \param character_actual How many utf32 characters in the string were drawn.
     * \return Returns the new caret X-coordinate position. Pass this value to the next draw_string call to get
     * continious text.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_aa(int32_t x, int32_t y, std::string_view str, uint8_t col,
                                     size_t character_count = std::numeric_limits<size_t>::max(),
                                     size_t *character_actual = nullptr) {
        static_assert(FONT::mono == false, "Can't use a mono font to draw antialiased text.");
        return draw_text<FONT, KERNING, ROTATION>(x, y, str, col, character_count, character_actual);
    }

    /**
     * \brief Draw a NUL-terminated string. See draw_string_aa(int32_t, int32_t, std::string_view, uint8_t, size_t,
     * size_t *).
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param str UTF-8 string. nullptr draws nothing.
     * \param col Color palette index to use.
     * \param character_count How many utf32 characters in the string should be drawn.
     * \param character_actual How many utf32 characters in the string were drawn.
     * \return Returns the new caret position.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_aa(int32_t x, int32_t y, const char *str, uint8_t col,
                                       size_t character_count = std::numeric_limits<size_t>::max(),
                                     size_t *character_actual = nullptr) {
        if (str == nullptr) {
            return 0;
        }
        return draw_string_aa<FONT, KERNING, ROTATION>(x, y, std::string_view{str}, col, character_count,
                                                         character_actual);
    }

    /**
     * \brief Draw antialiased text at the specified coordinate. The template parameter selects which antialiased font
     * to use. Only format_8bit targets are supported. Typical use:
     *
     * \code{.cpp}
     * #include "some_font_aa.h"
     * ...
     *     image.draw_string_aa<constixel::some_font_aa>({.x=0, .y=0, .str="MyText", .col=constixel::color::WHITE});
     * ...
     * \endcode
     *
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available. Default to false.
     * \tparam ROTATION Rotation around the x/y coordinate. Can be text_rotation::DEGREE_0, text_rotation::DEGREE_90,
     * text_rotation::DEGREE_180 or text_rotation::DEGREE_270
     * \param d \ref draw_string initializer struct.
     * \param character_count How many utf32 characters in the string should be drawn.
     * \param character_actual How many utf32 characters in the string were drawn.
     * \return Returns the new caret X-coordinate position. Pass this value to the next draw_string call to get
     * continious text.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_aa(const struct draw_string &d,
                                     size_t character_count = std::numeric_limits<size_t>::max(),
                                     size_t *character_actual = nullptr) {
        static_assert(FONT::mono == false, "Can't use a mono font to draw antialiased text.");
        return draw_string<FONT, KERNING, ROTATION>(d, character_count, character_actual);
    }

    /**
     * \brief Draw antialiased text centered at the specified coordinate. The template parameter selects which
     * antialiased font to use. Only format_8bit targets are supported. Typical use
     *
     * \tparam FONT The font struct name.
     * \tparam KERNING Boolean, use kerning information if available. Default to false.
     * \tparam ROTATION Rotation around the x/y coordinate. Can be text_rotation::DEGREE_0, text_rotation::DEGREE_90,
     * text_rotation::DEGREE_180 or text_rotation::DEGREE_270
     * \param x Center/Starting X-coordinate in pixels.
     * \param y Center/Starting Y-coordinate in pixels.
     * \param str UTF-8 string. Does not have to be NUL-terminated.
     * \param col Color palette index to use.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr void draw_string_centered_aa(int32_t x, int32_t y, std::string_view str, uint8_t col) {
        static_assert(FONT::mono == false, "Can't use a mono font to draw antialiased text.");
        draw_string<FONT, KERNING, ROTATION>(x, y, str, col, ALIGN_CENTER);
    }

    /**
     * \brief Draw a centered NUL-terminated string. See draw_string_centered_aa(int32_t, int32_t,
     * std::string_view, uint8_t).
     * \param x Center/Starting X-coordinate in pixels.
     * \param y Center/Starting Y-coordinate in pixels.
     * \param str UTF-8 string. nullptr draws nothing.
     * \param col Color palette index to use.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr void draw_string_centered_aa(int32_t x, int32_t y, const char *str, uint8_t col) {
        draw_string_centered_aa<FONT, KERNING, ROTATION>(x, y, hidden::cstr_view(str), col);
    }

    /**
     * \brief Fill a rectangle with the specified color. Example:
     *
     * \code{.cpp}
     * image.fill_rect(0, 0, 320, 240, constixel::color::WHITE);
     * \endcode
     *
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param w Width of the rectangle in pixels.
     * \param h Height of the rectangle in pixels.
     * \param col Color palette index to use.
     */
    constexpr void fill_rect(int32_t x, int32_t y, int32_t w, int32_t h, uint8_t col) {
        const auto dst = pix;
        auto minmax_check = std::minmax({x, y, w, h});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        if (w <= 0 || h <= 0) {
            return;
        }
        if (check_not_in_bounds(x, y, w, h)) {
            return;
        }
        if (y < 0) {
            h += y;
            y = 0;
        }
        if (y + h >= static_cast<int32_t>(pix.height)) {
            h = static_cast<int32_t>(pix.height) - y;
        }
        h += y;
        for (; y < h; y++) {
            extent(dst, x, w, y, col);
        }
    }

    /**
     * \brief Fill a rectangle with the specified color. Example:
     *
     * \code{.cpp}
     * image.fill_rect({.x=0, .y=0, .w=320, .h=240, .col=constixel::color::WHITE});
     * \endcode
     *
     * \param d \ref draw_rect initializer struct
     */
    constexpr void fill_rect(const struct draw_rect &d) {
        fill_rect(d.x, d.y, d.w, d.h, d.col);
    }

    /**
     * \brief Fill a rectangle with a shader function that generates colors per pixel.
     *
     * \code{.cpp}
     * // Linear gradient from red to blue
     * image.fill_rect(0, 0, 320, 240, [](float u, float v, float au, float av) -> std::array<float, 4> {
     *     return {1.0f - u, 0.0f, u, 1.0f}; // R, G, B, A
     * });
     * \endcode
     *
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param w Width of the rectangle in pixels.
     * \param h Height of the rectangle in pixels.
     * \param shader Lambda function taking (u, v, x, y) normalized coordinates
     *                    and returning RGBA color as std::array<float, 4>
     */
    template <typename shader_func>
    constexpr auto fill_rect(int32_t x, int32_t y, int32_t w, int32_t h, const shader_func &shader) -> void
        requires std::is_invocable_r_v<std::array<float, 4>, shader_func, float, float, int32_t, int32_t>
    {
        fill_rect_int(x, y, w, h, shader, x, y, w, h);
    }

    /**
     * \brief Draw a stroked rectangle with the specified color and stroke width. Example:
     *
     * \code{.cpp}
     * image.stroke_rect(0, 0, 320, 240, constixel::color::WHITE);
     * \endcode
     *
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param w Width of the rectangle in pixels.
     * \param h Height of the rectangle in pixels.
     * \param col Color palette index to use.
     * \param stroke_width Width of the stroke in pixels.
     */
    constexpr void stroke_rect(int32_t x, int32_t y, int32_t w, int32_t h, uint8_t col, int32_t stroke_width = 1) {
        auto minmax_check = std::minmax({x, y, w, h, stroke_width});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        if (w <= 0 || h <= 0) {
            return;
        }
        stroke_width = std::min(abs(stroke_width), std::max(w, h) / 2);
        fill_rect(x, y, stroke_width, h, col);
        fill_rect(x + stroke_width, y, w - stroke_width * 2, stroke_width, col);
        fill_rect(x + w - stroke_width, y, stroke_width, h, col);
        fill_rect(x + stroke_width, y + h - stroke_width, w - stroke_width * 2, stroke_width, col);
    }

    /**
     * \brief Draw a stroked rectangle with the specified color and stroke width. Example:
     *
     * \code{.cpp}
     * image.stroke_rect({.x=0, .y=0, .w=320, .h=240, .col=constixel::color::WHITE, .sw=2});
     * \endcode
     *
     * \param d \ref draw_rect initializer struct
     */
    constexpr void stroke_rect(const struct draw_rect &d) {
        stroke_rect(d.x, d.y, d.w, d.h, d.col, d.sw);
    }

    /**
     * \brief Fill a circle with the specified radius and color. Example:
     *
     * \code{.cpp}
     * image.fill_circle(64, 64, 32, constixel::color::WHITE);
     * \endcode
     *
     * \param cx Center X-coordinate of the circle in pixels.
     * \param cy Center Y-coordinate of the circle in pixels.
     * \param radius radius of the circle in pixels.
     * \param col Color palette index to use.
     */
    constexpr void fill_circle(int32_t cx, int32_t cy, int32_t radius, uint8_t col) {
        auto minmax_check = std::minmax({cx, cy, radius});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        radius = abs(radius);
        if (radius == 1) {
            fill_rect(cx - 1, cy - 1, 2, 2, col);
            return;
        }
        circle_int<false, false>(cx, cy, radius, 0, 0, col, 0);
    }

    /**
     * \brief Fill a circle with the specified radius and color. Example:
     *
     * \code{.cpp}
     * image.fill_circle({.cx=64, .cy=64, .r=32, .col=constixel::color::WHITE});
     * \endcode
     *
     * \param d \ref draw_circle initializer struct
     */
    constexpr void fill_circle(const struct draw_circle &d) {
        fill_circle(d.cx, d.cy, d.r, d.col);
    }

    /**
     * \brief Stroke a circle with the specified radius and color. Example:
     *
     * \code{.cpp}
     * image.stroke_circle(64, 64, 32, 4, constixel::color::WHITE);
     * \endcode
     *
     * \param cx Center X-coordinate of the circle in pixels.
     * \param cy Center Y-coordinate of the circle in pixels.
     * \param radius radius of the circle in pixels.
     * \param col Color palette index to use.
     * \param stroke_width Width of the stroke in pixels.
     */
    constexpr void stroke_circle(int32_t cx, int32_t cy, int32_t radius, uint8_t col, int32_t stroke_width = 1) {
        auto minmax_check = std::minmax({cx, cy, radius, stroke_width});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        radius = abs(radius);
        stroke_width = abs(stroke_width);
        if (radius == 1) {
            fill_rect(cx - 1, cy - 1, 2, 2, col);
            return;
        }
        if (stroke_width >= radius) {
            circle_int<false, false>(cx, cy, radius, 0, 0, col, 0);
            return;
        }
        circle_int<false, true>(cx, cy, radius, 0, 0, col, stroke_width);
    }

    /**
     * \brief Stroke a circle with the specified radius and color. Example:
     *
     * \code{.cpp}
     * image.stroke_circle({.cx=64, .cy=64, .r=32, .col=constixel::color::WHITE, .sw=2});
     * \endcode
     *
     * \param d \ref draw_circle initializer struct
     */
    constexpr void stroke_circle(const struct draw_circle &d) {
        stroke_circle(d.cx, d.cy, d.r, d.col, d.sw);
    }

    /**
     * \brief Stroke a circle using antialiasing with the specified radius and color. Only format_8bit targets are
     * supported. Example:
     *
     * \code{.cpp}
     * image.stroke_circle_aa(64, 64, 32, constixel::color::WHITE, 2);
     * \endcode
     *
     * \param cx Center X-coordinate of the circle in pixels.
     * \param cy Center Y-coordinate of the circle in pixels.
     * \param radius radius of the circle in pixels.
     * \param col Color palette index to use.
     * \param stroke_width Width of the stroke in pixels.
     */
    constexpr void stroke_circle_aa(int32_t cx, int32_t cy, int32_t radius, uint8_t col, int32_t stroke_width = 1) {
        auto minmax_check = std::minmax({cx, cy, radius, stroke_width});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        radius = abs(radius);
        stroke_width = abs(stroke_width);
        if (radius == 1) {
            fill_rect(cx - 1, cy - 1, 2, 2, col);
            return;
        }
        if (stroke_width >= radius) {
            circle_int<true, false>(cx, cy, radius, 0, 0, col, 0);
            return;
        }
        circle_int<true, true>(cx, cy, radius, 0, 0, col, stroke_width);
    }

    /**
     * \brief Stroke a circle using antialiasing with the specified radius and color. Only format_8bit targets are
     * supported. Example:
     *
     * \code{.cpp}
     * image.stroke_circle_aa({.cx=64, .cy=64, .r=32, .col=constixel::color::WHITE, .sw=2});
     * \endcode
     *
     * \param d \ref draw_circle initializer struct
     */
    constexpr void stroke_circle_aa(const struct draw_circle &d) {
        stroke_circle_aa(d.cx, d.cy, d.r, d.col, d.sw);
    }

    /**
     * \brief Fill a circle using antialiasing with the specified radius and color. Only format_8bit targets are
     * supported. Example:
     *
     * \code{.cpp}
     * image.fill_circle_aa(64, 64, 32, constixel::color::WHITE);
     * \endcode
     *
     * \param cx Center X-coordinate of the circle in pixels.
     * \param cy Center Y-coordinate of the circle in pixels.
     * \param radius Radius of the circle in pixels.
     * \param col Color palette index to use.
     */
    constexpr void fill_circle_aa(int32_t cx, int32_t cy, int32_t radius, uint8_t col) {
        auto minmax_check = std::minmax({cx, cy, radius});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        radius = abs(radius);
        circle_int<true, false>(cx, cy, radius, 0, 0, col, 0);
    }

    /**
     * \brief Fill a circle using antialiasing with the specified radius and color. Only format_8bit targets are
     * supported. Example:
     *
     * \code{.cpp}
     * image.fill_circle_aa({.cx=64, .cy=64, .r=32, .col=constixel::color::WHITE});
     * \endcode
     *
     * \param d \ref draw_circle initializer struct
     */
    constexpr void fill_circle_aa(const struct draw_circle &d) {
        fill_circle_aa(d.cx, d.cy, d.r, d.col);
    }

    /**
     * \brief Fill a circle using antialiasing with a shader function that generates colors per pixel.
     *
     * \code{.cpp}
     * // Radial gradient circle
     * image.fill_circle_aa(100, 100, 50, [](float u, float v, float au, float av) -> std::array<float, 4> {
     *     float dist = std::sqrt((u - 0.5f) * (u - 0.5f) + (v - 0.5f) * (v - 0.5f)) * 2.0f;
     *     float intensity = 1.0f - std::clamp(dist, 0.0f, 1.0f);
     *     return {intensity, intensity * 0.8f, 1.0f, 1.0f}; // R, G, B, A
     * });
     * \endcode
     *
     * \param cx Center X-coordinate of the circle in pixels.
     * \param cy Center Y-coordinate of the circle in pixels.
     * \param radius Radius of the circle in pixels.
     * \param shader Lambda function taking (u, v, x, y) normalized coordinates
     *                    and returning RGBA color as std::array<float, 4>
     */
    template <typename shader_func>
    constexpr auto fill_circle_aa(int32_t cx, int32_t cy, int32_t radius, const shader_func &shader) -> void
        requires std::is_invocable_r_v<std::array<float, 4>, shader_func, float, float, int32_t, int32_t>
    {
        auto minmax_check = std::minmax({cx, cy, radius});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        radius = abs(radius);
        circle_int_shader(cx, cy, radius, 0, 0, shader);
    }

    /**
     * \brief Stroke a rounded rectangle with the specified color. Example:
     *
     * \code{.cpp}
     * image.stroke_round_rect(0, 0, 200, 100, 15, constixel::color::WHITE);
     * \endcode
     *
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param w Width of the rectangle in pixels.
     * \param h Height of the rectangle in pixels.
     * \param radius Radius of the rounded corners in pixels.
     * \param col Color palette index to use.
     * \param stroke_width Width of the stroke in pixels.
     */
    constexpr void stroke_round_rect(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius, uint8_t col,
                                     int32_t stroke_width = 1) {
        auto minmax_check = std::minmax({x, y, w, h, radius, stroke_width});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        if (w <= 0 || h <= 0) {
            return;
        }
        radius = abs(radius);
        stroke_width = abs(stroke_width);
        const int32_t cr = std::min({w / 2, h / 2, radius});
        const int32_t dx = w - cr * 2;
        const int32_t dy = h - cr * 2;
        if (radius == 0) {
            stroke_rect(x, y, w, h, col, stroke_width);
        } else {
            if (radius > stroke_width) {
                circle_int<false, true>(x + cr, y + cr, cr, dx, dy, col, stroke_width);
            } else {
                circle_int<false, false>(x + cr, y + cr, cr, dx, dy, col, 0);
            }
            fill_rect(x, y + cr, stroke_width, dy, col);
            fill_rect(x + w - stroke_width, y + cr, stroke_width, dy, col);
            fill_rect(x + cr, y, w - cr * 2, stroke_width, col);
            fill_rect(x + cr, y + h - stroke_width, w - cr * 2, stroke_width, col);
        }
    }

    /**
     * \brief Stroke a rounded rectangle with the specified color. Example:
     *
     * \code{.cpp}
     * image.stroke_round_rect({.x=0, .y=0, .w=200, .h=100, .r=15, .col=constixel::color::WHITE, .sw=8});
     * \endcode
     *
     * \param d \ref draw_round_rect initializer struct
     */
    constexpr void stroke_round_rect(const struct draw_round_rect &d) {
        stroke_round_rect(d.x, d.y, d.w, d.h, d.r, d.col, d.sw);
    }

    /**
     * \brief Fill a rounded rectangle with the specified color. Example:
     *
     * \code{.cpp}
     * image.fill_round_rect(0, 0, 200, 100, 15, constixel::color::WHITE);
     * \endcode
     *
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param w Width of the rectangle in pixels.
     * \param h Height of the rectangle in pixels.
     * \param radius Radius of the rounded corners in pixels.
     * \param col Color palette index to use.
     */
    constexpr void fill_round_rect(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius, uint8_t col) {
        auto minmax_check = std::minmax({x, y, w, h, radius});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        if (w <= 0 || h <= 0) {
            return;
        }
        const int32_t cr = std::min({w / 2, h / 2, radius});
        const int32_t dx = w - cr * 2;
        const int32_t dy = h - cr * 2;
        circle_int<false, false>(x + cr, y + cr, cr, dx, dy, col, 0);
        fill_rect(x, y + cr, cr, dy, col);
        fill_rect(x + w - cr, y + cr, cr, dy, col);
        fill_rect(x + cr, y, dx, h, col);
    }

    /**
     * \brief Fill a rounded rectangle with the specified color. Example:
     *
     * \code{.cpp}
     * image.fill_round_rect({.x=0, .y=0, .w=200, .h=100, .r=15, .col=constixel::color::WHITE});
     * \endcode
     *
     * \param d \ref draw_round_rect initializer struct
     */
    constexpr void fill_round_rect(const struct draw_round_rect &d) {
        fill_round_rect(d.x, d.y, d.w, d.h, d.r, d.col);
    }

    /**
     * \brief Fill a rounded rectangle using antialiasing with the specified color. Only format_8bit targets are
     * supported. Example:
     *
     * \code{.cpp}
     * image.fill_round_rect_aa(0, 0, 200, 100, 15, constixel::color::WHITE);
     * \endcode
     *
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param w Width of the rectangle.
     * \param h Height of the rectangle.
     * \param radius Radius of the rounded corners.
     * \param col Color palette index to use.
     */
    constexpr void fill_round_rect_aa(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius, uint8_t col) {
        auto minmax_check = std::minmax({x, y, w, h, radius});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        if (w <= 0 || h <= 0) {
            return;
        }
        radius = abs(radius);
        int32_t cr = std::min({w / 2, h / 2, radius});
        int32_t dx = w - cr * 2;
        int32_t dy = h - cr * 2;
        circle_int<true, false>(x + cr, y + cr, cr, dx, dy, col, 0);
        fill_rect(x, y + cr, cr, dy, col);
        fill_rect(x + w - cr, y + cr, cr, dy, col);
        fill_rect(x + cr, y, dx, h, col);
    }

    /**
     * \brief Fill a rounded rectangle using antialiasing with the specified color. Only format_8bit targets are
     * supported.
     *
     * \code{.cpp}
     * image.fill_round_rect_aa({.x=0, .y=0, .w=200, .h=100, .r=15, .col=constixel::color::WHITE, .sw=8});
     * \endcode
     *
     * \param d \ref draw_round_rect initializer struct
     */
    constexpr void fill_round_rect_aa(const struct draw_round_rect &d) {
        fill_round_rect_aa(d.x, d.y, d.w, d.h, d.r, d.col);
    }

    /**
     * \brief Fill a rounded rectangle using antialiasing with a shader function that generates colors per pixel.
     *
     * \code{.cpp}
     * // Gradient rounded rectangle
     * image.fill_round_rect_aa(10, 10, 200, 100, 15, [](float u, float v, float au, float av) -> std::array<float, 4> {
     *     return {u, v, 1.0f - u, 1.0f}; // Corner-to-corner gradient
     * });
     * \endcode
     *
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param w Width of the rectangle in pixels.
     * \param h Height of the rectangle in pixels.
     * \param radius Radius of the rounded corners in pixels.
     * \param shader Lambda function taking (u, v, x, y) normalized coordinates
     *                    and returning RGBA color as std::array<float, 4>
     */
    template <typename shader_func>
    constexpr auto fill_round_rect_aa(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius,
                                      const shader_func &shader) -> void
        requires std::is_invocable_r_v<std::array<float, 4>, shader_func, float, float, int32_t, int32_t>
    {
        auto minmax_check = std::minmax({x, y, w, h, radius});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        if (w <= 0 || h <= 0) {
            return;
        }
        radius = abs(radius);
        int32_t cr = std::min({w / 2, h / 2, radius});
        int32_t dx = w - cr * 2;
        int32_t dy = h - cr * 2;

        circle_int_shader(x + cr, y + cr, cr, dx, dy, shader, x, y, w, h);

        fill_rect_int(x, y + cr, cr, dy, shader, x, y, w, h);
        fill_rect_int(x + w - cr, y + cr, cr, dy, shader, x, y, w, h);
        fill_rect_int(x + cr, y, dx, h, shader, x, y, w, h);
    }

    /**
     * \brief Stroke a rounded rectangle using antialiasing with the specified color. Example:
     *
     * \code{.cpp}
     * image.stroke_round_rect_aa(0, 0, 200, 100, 15, constixel::color::WHITE);
     * \endcode
     *
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param w Width of the rectangle in pixels.
     * \param h Height of the rectangle in pixels.
     * \param radius Radius of the rounded corners in pixels.
     * \param col Color palette index to use.
     * \param stroke_width Width of the stroke in pixels.
     */
    constexpr void stroke_round_rect_aa(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius, uint8_t col,
                                        int32_t stroke_width = 1) {
        auto minmax_check = std::minmax({x, y, w, h, radius, stroke_width});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        if (w <= 0 || h <= 0) {
            return;
        }
        radius = abs(radius);
        stroke_width = abs(stroke_width);
        const int32_t cr = std::min({w / 2, h / 2, radius});
        const int32_t dx = w - cr * 2;
        const int32_t dy = h - cr * 2;
        if (radius == 0) {
            stroke_rect(x, y, w, h, col, stroke_width);
        } else {
            if (radius > stroke_width) {
                circle_int<true, true>(x + cr, y + cr, cr, dx, dy, col, stroke_width);
            } else {
                circle_int<true, false>(x + cr, y + cr, cr, dx, dy, col, 0);
            }
            fill_rect(x, y + cr, stroke_width, dy, col);
            fill_rect(x + w - stroke_width, y + cr, stroke_width, dy, col);
            fill_rect(x + cr, y, w - cr * 2, stroke_width, col);
            fill_rect(x + cr, y + h - stroke_width, w - cr * 2, stroke_width, col);
        }
    }

    /**
     * \brief Stroke a rounded rectangle using antialiasing with the specified color. Example:
     *
     * \code{.cpp}
     * image.stroke_round_rect_aa({.x=0, .y=0, .w=200, .h=100, .r=15, .col=constixel::color::WHITE, .sw=8});
     * \endcode
     *
     * \param d \ref draw_round_rect initializer struct
     */
    constexpr void stroke_round_rect_aa(const struct draw_round_rect &d) {
        stroke_round_rect_aa(d.x, d.y, d.w, d.h, d.r, d.col, d.sw);
    }

    /**
     * \brief Copy the pixels of a surface of the same size. No compositing occurs, replaces current content.
     * \param src Source surface.
     * \return false if the sizes differ. Nothing is copied then.
     */
    [[nodiscard]] constexpr bool copy(const const_surface<PF> &src) {
        if (src.width() != this->width() || src.height() != this->height()) {
            return false;
        }
        const size_t line_bytes = PF::bytes_per_line(pix.width);
        for (size_t y = 0; y < pix.height; y++) {
            std::copy_n(src.data() + y * src.stride(), line_bytes, pix.data + y * pix.stride);
        }
        return true;
    }

    /**
     * \brief Copy raw pixel data without line padding. No compositing occurs, replaces current content.
     * \param src Source data, height() * PF::bytes_per_line(width()) bytes.
     * \return false if the size of src differs. Nothing is copied then.
     */
    [[nodiscard]] constexpr bool copy(std::span<const uint8_t> src) {
        const size_t line_bytes = PF::bytes_per_line(pix.width);
        if (src.size() != line_bytes * pix.height) {
            return false;
        }
        return copy(const_surface<PF>{src.data(), pix.width, pix.height, line_bytes});
    }

    /**
     * \brief Flip the contents of this image horizontally.
     */
    constexpr void flip_h() {
        const size_t line_bytes = PF::bytes_per_line(pix.width);
        uint8_t *ptr = pix.data;
        for (size_t y = 0; y < pix.height; y++) {
            for (size_t x = 0; x < line_bytes / 2; x++) {
                const uint8_t a = PF::reverse(ptr[x]);
                const uint8_t b = PF::reverse(ptr[line_bytes - x - 1]);
                ptr[x] = b;
                ptr[line_bytes - x - 1] = a;
            }
            ptr += pix.stride;
        }
    }

    /**
     * \brief Flip the contents of this image vertically.
     */
    constexpr void flip_v() {
        const size_t line_bytes = PF::bytes_per_line(pix.width);
        uint8_t *ptr = pix.data;
        for (size_t x = 0; x < line_bytes; x++) {
            for (size_t y = 0; y < pix.height / 2; y++) {
                uint8_t *ptr_a = ptr + y * pix.stride;
                uint8_t *ptr_b = ptr + (pix.height - y - 1) * pix.stride;
                const uint8_t a = ptr_a[x];
                const uint8_t b = ptr_b[x];
                ptr_a[x] = b;
                ptr_b[x] = a;
            }
        }
    }

    /**
     * \brief Flip the contents of this image horizontally & vertically.
     */
    constexpr void flip_hv() {
        const size_t line_bytes = PF::bytes_per_line(pix.width);
        uint8_t *ptr = pix.data;
        for (size_t x = 0; x < line_bytes; x++) {
            for (size_t y = 0; y < pix.height / 2; y++) {
                uint8_t *ptr_a = ptr + y * pix.stride;
                uint8_t *ptr_b = ptr + (pix.height - y - 1) * pix.stride;
                uint8_t a = PF::reverse(ptr_a[x]);
                uint8_t b = PF::reverse(ptr_b[line_bytes - x - 1]);
                ptr_a[x] = b;
                ptr_b[line_bytes - x - 1] = a;
            }
        }
    }

    /**
     * \brief Blit an RGBA8 buffer into this instance using brute force color mapping.
     * \param x Starting X-coordinate position in the target instance in pixels.
     * \param y Starting Y-coordinate position in the target instance in pixels.
     * \param w Width of the rectangle. If the width is smaller than the RGBA8 buffer content will be clipped.
     * \param h Weight of the rectangle. If the height is smaller than the RGBA8 buffer content will be clipped.
     * \param ptr Pointer to the RGBA8 buffer.
     * \param iw Width in pixels of the RGBA8 buffer.
     * \param ih Height in pixels of the RGBA8 buffer.
     * \param stride Bytes per line of pixels of the RGBA8 buffer.
     */
    constexpr void blit_RGBA(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih,
                             int32_t stride) {
        constixel::rect<int32_t> blitrect{.x = x, .y = y, .w = w, .h = h};
        blitrect &= {.x = 0, .y = 0, .w = static_cast<int32_t>(pix.width), .h = static_cast<int32_t>(pix.height)};
        blitrect &= {.x = x, .y = y, .w = iw, .h = ih};
        if constexpr (hidden::has_runtime_palette<PF>) {
            if (this->palette() != nullptr) {
                PF::blit_RGBA_with_palette(pix, blitrect, ptr, stride, *this->palette());
                return;
            }
        }
        PF::blit_RGBA(pix, blitrect, ptr, stride);
    }

    /**
     * \brief Blit an RGBA8 buffer into this instance using brute force color mapping.
     * \param dstrect Rectangular area in the target buffer to blit into. If the rectangle is smaller than the RGBA8
     * buffer, clipping occurs.
     * \param ptr Pointer to the RGBA8 buffer.
     * \param iw Width in pixels of the RGBA8 buffer.
     * \param ih Weight in pixels of the RGBA8 buffer.
     * \param stride Bytes per line of pixels of the RGBA8 buffer.
     */
    constexpr void blit_RGBA(const constixel::rect<int32_t> &dstrect, const uint8_t *ptr, int32_t iw, int32_t ih,
                             int32_t stride) {
        constixel::rect<int32_t> blitrect{dstrect};
        blitrect &= {.x = 0, .y = 0, .w = static_cast<int32_t>(pix.width), .h = static_cast<int32_t>(pix.height)};
        blitrect &= {.x = dstrect.x, .y = dstrect.y, .w = iw, .h = ih};
        PF::blit_RGBA(pix, blitrect, ptr, stride);
    }

    /**
     * \brief Blit an RGBA8 buffer into this instance using brute force color mapping. Simple integer based diffusion is
     * applied.
     * \param x Starting X-coordinate position in the target instance in pixels.
     * \param y Starting Y-coordinate position in the target instance in pixels.
     * \param w Width of the rectangle. If the width is smaller than the RGBA8 buffer content will be clipped.
     * \param h Weight of the rectangle. If the height is smaller than the RGBA8 buffer content will be clipped.
     * \param ptr Pointer to the RGBA8 buffer.
     * \param iw Width in pixels of the RGBA8 buffer.
     * \param ih Height in pixels of the RGBA8 buffer.
     * \param stride Bytes per line of pixels of the RGBA8 buffer.
     */
    constexpr void blit_RGBA_diffused(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw,
                                      int32_t ih, int32_t stride) {
        constixel::rect<int32_t> blitrect{.x = x, .y = y, .w = w, .h = h};
        blitrect &= {.x = 0, .y = 0, .w = static_cast<int32_t>(pix.width), .h = static_cast<int32_t>(pix.height)};
        blitrect &= {.x = x, .y = y, .w = iw, .h = ih};
        PF::blit_RGBA_diffused(pix, blitrect, ptr, stride);
    }

    /**
     * \brief Blit an RGBA8 buffer into this instance using brute force color mapping. Simple integer based diffusion is
     * applied.
     * \param dstrect Rectangular area in the target buffer to blit into. If the rectangle is smaller than the RGBA8
     * buffer, clipping occurs.
     * \param ptr Pointer to the RGBA8 buffer.
     * \param iw Width in pixels of the RGBA8 buffer.
     * \param ih Weight in pixels of the RGBA8 buffer.
     * \param stride Bytes per line of pixels of the RGBA8 buffer.
     */
    constexpr void blit_RGBA_diffused(const constixel::rect<int32_t> &dstrect, const uint8_t *ptr, int32_t iw,
                                      int32_t ih, int32_t stride) {
        constixel::rect<int32_t> blitrect{dstrect};
        blitrect &= {.x = 0, .y = 0, .w = static_cast<int32_t>(pix.width), .h = static_cast<int32_t>(pix.height)};
        blitrect &= {.x = dstrect.x, .y = dstrect.y, .w = iw, .h = ih};
        PF::blit_RGBA_diffused(pix, blitrect, ptr, stride);
    }

    /**
     * \brief Blit an RGBA8 buffer into this instance using brute force color mapping. Diffusion in linear color space
     * is applied.
     * \param x Starting X-coordinate position in the target instance in pixels.
     * \param y Starting Y-coordinate position in the target instance in pixels.
     * \param w Width of the rectangle. If the width is smaller than the RGBA8 buffer content will be clipped.
     * \param h Weight of the rectangle. If the height is smaller than the RGBA8 buffer content will be clipped.
     * \param ptr Pointer to the RGBA8 buffer.
     * \param iw Width in pixels of the RGBA8 buffer.
     * \param ih Height in pixels of the RGBA8 buffer.
     * \param stride Bytes per line of pixels of the RGBA8 buffer.
     */
    constexpr void blit_RGBA_diffused_linear(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw,
                                             int32_t ih, int32_t stride) {
        constixel::rect<int32_t> blitrect{.x = x, .y = y, .w = w, .h = h};
        blitrect &= {.x = 0, .y = 0, .w = static_cast<int32_t>(pix.width), .h = static_cast<int32_t>(pix.height)};
        blitrect &= {.x = x, .y = y, .w = iw, .h = ih};
        PF::blit_RGBA_diffused_linear(pix, blitrect, ptr, stride);
    }

    /**
     * \brief Blit an RGBA8 buffer into this instance using brute force color mapping. Diffusion in linear color space
     * is applied.
     * \param dstrect Rectangular area in the target buffer to blit into. If the rectangle is smaller than the RGBA8
     * buffer, clipping occurs.
     * \param ptr Pointer to the RGBA8 buffer.
     * \param iw Width in pixels of the RGBA8 buffer.
     * \param ih Weight in pixels of the RGBA8 buffer.
     * \param stride Bytes per line of pixels of the RGBA8 buffer.
     */
    constexpr void blit_RGBA_diffused_linear(const constixel::rect<int32_t> &dstrect, const uint8_t *ptr, int32_t iw,
                                             int32_t ih, int32_t stride) {
        constixel::rect<int32_t> blitrect{dstrect};
        blitrect &= {.x = 0, .y = 0, .w = static_cast<int32_t>(pix.width), .h = static_cast<int32_t>(pix.height)};
        blitrect &= {.x = dstrect.x, .y = dstrect.y, .w = iw, .h = ih};
        PF::blit_RGBA_diffused_linear(pix, blitrect, ptr, stride);
    }

    /**
     * @brief Fluent shape API methods for method chaining
     * These provide a more expressive way to create shapes compared to struct-based calls.
     *
     * Example usage:
     * \code{.cpp}
     * image.rect(10, 10, 50, 30).fill(constixel::color::RED).stroke(constixel::color::BLACK, 2);
     * image.circle(100, 100, 20).fill(constixel::color::BLUE);
     * image.line(0, 0, 100, 100).stroke(constixel::color::WHITE, 3);
     * \endcode
     */

    /**
     * \brief Create a rectangle shape for fluent method chaining.
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param w Width of the rectangle in pixels.
     * \param h Height of the rectangle in pixels.
     * \return A rect shape object that supports .fill() and .stroke() methods.
     */
    constexpr auto rect(int32_t x, int32_t y, int32_t w, int32_t h) {
        return shapes::rect<PF>(*this, x, y, w, h);
    }

    /**
     * \brief Create a circle shape for fluent method chaining.
     * \param cx Center X-coordinate in pixels.
     * \param cy Center Y-coordinate in pixels.
     * \param r Radius of the circle in pixels.
     * \return A circle shape object that supports .fill() and .stroke() methods.
     */
    constexpr auto circle(int32_t cx, int32_t cy, int32_t r) {
        return shapes::circle<PF>(*this, cx, cy, r);
    }

    /**
     * \brief Create an antialiased circle shape for fluent method chaining.
     * Only format_8bit targets are supported.
     * \param cx Center X-coordinate in pixels.
     * \param cy Center Y-coordinate in pixels.
     * \param r Radius of the circle in pixels.
     * \return A circle_aa shape object that supports .fill() and .stroke() methods.
     */
    constexpr auto circle_aa(int32_t cx, int32_t cy, int32_t r) {
        return shapes::circle_aa<PF>(*this, cx, cy, r);
    }

    /**
     * \brief Create a rounded rectangle shape for fluent method chaining.
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param w Width of the rectangle in pixels.
     * \param h Height of the rectangle in pixels.
     * \param radius Radius of the rounded corners in pixels.
     * \return A round_rect shape object that supports .fill() and .stroke() methods.
     */
    constexpr auto round_rect(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius) {
        return shapes::round_rect<PF>(*this, x, y, w, h, radius);
    }

    /**
     * \brief Create an antialiased rounded rectangle shape for fluent method chaining.
     * Only format_8bit targets are supported.
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param w Width of the rectangle in pixels.
     * \param h Height of the rectangle in pixels.
     * \param radius Radius of the rounded corners in pixels.
     * \return A round_rect_aa shape object that supports .fill() and .stroke() methods.
     */
    constexpr auto round_rect_aa(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius) {
        return shapes::round_rect_aa<PF>(*this, x, y, w, h, radius);
    }

    /**
     * \brief Create a line shape for fluent method chaining.
     * \param x0 Starting X-coordinate in pixels.
     * \param y0 Starting Y-coordinate in pixels.
     * \param x1 Ending X-coordinate in pixels.
     * \param y1 Ending Y-coordinate in pixels.
     * \return A line shape object that supports .stroke() method.
     */
    constexpr auto line(int32_t x0, int32_t y0, int32_t x1, int32_t y1) {
        return shapes::line<PF>(*this, x0, y0, x1, y1);
    }

    /**
     * \brief Create an antialiased line shape for fluent method chaining.
     * Only format_8bit targets are supported.
     * \param x0 Starting X-coordinate in pixels.
     * \param y0 Starting Y-coordinate in pixels.
     * \param x1 Ending X-coordinate in pixels.
     * \param y1 Ending Y-coordinate in pixels.
     * \return A line_aa shape object that supports .stroke() method.
     */
    constexpr auto line_aa(int32_t x0, int32_t y0, int32_t x1, int32_t y1) {
        return shapes::line_aa<PF>(*this, x0, y0, x1, y1);
    }

    /**
     * \brief Create a point shape for fluent method chaining.
     * \param x X-coordinate in pixels.
     * \param y Y-coordinate in pixels.
     * \return A point shape object that supports .plot() method.
     */
    constexpr auto point(int32_t x, int32_t y) {
        return shapes::point<PF>(*this, x, y);
    }

    /**
     * \brief Create a text shape for fluent monospace string drawing with method chaining.
     * \tparam FONT The font type to use for rendering.
     * \tparam KERNING Enable kerning if true (default: false).
     * \tparam ROTATION Text rotation (default: DEGREE_0).
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param str The string to draw. Does not have to be NUL-terminated; it has to outlive the shape.
     * \return A text_mono shape object that supports .color() and .draw() methods.
     *
     * Example usage:
     * \code
     * img.text_mono<constixel::ibmplexmono_regular_12_mono>(10, 20, "Hello").color(constixel::color::WHITE);
     * \endcode
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_mono(int32_t x, int32_t y, std::string_view str) {
        return shapes::text_mono<FONT, PF, KERNING, ROTATION>(*this, x, y, str);
    }

    /**
     * \brief Create a text shape from a NUL-terminated string. See text_mono(int32_t, int32_t, std::string_view).
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param str The string to draw. nullptr draws nothing.
     * \return A text_mono shape object.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_mono(int32_t x, int32_t y, const char *str) {
        return text_mono<FONT, KERNING, ROTATION>(x, y, hidden::cstr_view(str));
    }

    /**
     * \brief Create a text shape for fluent antialiased string drawing with method chaining.
     * Only format_8bit targets are supported.
     * \tparam FONT The font type to use for rendering.
     * \tparam KERNING Enable kerning if true (default: false).
     * \tparam ROTATION Text rotation (default: DEGREE_0).
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param str The string to draw. Does not have to be NUL-terminated; it has to outlive the shape.
     * \return A text_aa shape object that supports .color() and .draw() methods.
     *
     * Example usage:
     * \code
     * img.text_aa<constixel::ibmplexmono_regular_12_aa>(10, 20, "Hello").color(constixel::color::WHITE);
     * \endcode
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_aa(int32_t x, int32_t y, std::string_view str) {
        return shapes::text_aa<FONT, PF, KERNING, ROTATION>(*this, x, y, str);
    }

    /**
     * \brief Create a text shape from a NUL-terminated string. See text_aa(int32_t, int32_t, std::string_view).
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param str The string to draw. nullptr draws nothing.
     * \return A text_aa shape object.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_aa(int32_t x, int32_t y, const char *str) {
        return text_aa<FONT, KERNING, ROTATION>(x, y, hidden::cstr_view(str));
    }

    /**
     * \brief Create a text shape for fluent centered monospace string drawing with method chaining.
     * \tparam FONT The font type to use for rendering.
     * \tparam KERNING Enable kerning if true (default: false).
     * \tparam ROTATION Text rotation (default: DEGREE_0).
     * \param x Center X-coordinate in pixels.
     * \param y Center Y-coordinate in pixels.
     * \param str The string to draw. Does not have to be NUL-terminated; it has to outlive the shape.
     * \return A text_centered_mono shape object that supports .color() method.
     *
     * Example usage:
     * \code
     * img.text_centered_mono<constixel::ibmplexmono_regular_12_mono>(100, 50,
     * "Centered").color(constixel::color::WHITE);
     * \endcode
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_centered_mono(int32_t x, int32_t y, std::string_view str) {
        return shapes::text_centered_mono<FONT, PF, KERNING, ROTATION>(*this, x, y, str);
    }

    /**
     * \brief Create a text shape from a NUL-terminated string. See text_centered_mono(int32_t, int32_t,
     * std::string_view).
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param str The string to draw. nullptr draws nothing.
     * \return A text_centered_mono shape object.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_centered_mono(int32_t x, int32_t y, const char *str) {
        return text_centered_mono<FONT, KERNING, ROTATION>(x, y, hidden::cstr_view(str));
    }

    /**
     * \brief Create a text shape for fluent centered antialiased string drawing with method chaining.
     * Only format_8bit targets are supported.
     * \tparam FONT The font type to use for rendering.
     * \tparam KERNING Enable kerning if true (default: false).
     * \tparam ROTATION Text rotation (default: DEGREE_0).
     * \param x Center X-coordinate in pixels.
     * \param y Center Y-coordinate in pixels.
     * \param str The string to draw. Does not have to be NUL-terminated; it has to outlive the shape.
     * \return A text_centered_aa shape object that supports .color() method.
     *
     * Example usage:
     * \code
     * img.text_centered_aa<constixel::ibmplexmono_regular_12_aa>(100, 50, "Centered").color(constixel::color::WHITE);
     * \endcode
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_centered_aa(int32_t x, int32_t y, std::string_view str) {
        return shapes::text_centered_aa<FONT, PF, KERNING, ROTATION>(*this, x, y, str);
    }

    /**
     * \brief Create a text shape from a NUL-terminated string. See text_centered_aa(int32_t, int32_t,
     * std::string_view).
     * \param x Starting X-coordinate in pixels.
     * \param y Starting Y-coordinate in pixels.
     * \param str The string to draw. nullptr draws nothing.
     * \return A text_centered_aa shape object.
     */
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_centered_aa(int32_t x, int32_t y, const char *str) {
        return text_centered_aa<FONT, KERNING, ROTATION>(x, y, hidden::cstr_view(str));
    }

 private:
#ifndef __INTELLISENSE__
    /// @cond DOXYGEN_EXCLUDE
    constexpr bool check_not_in_bounds(int32_t x, int32_t y, int32_t w, int32_t h) {
        constixel::rect<int32_t> intersect_rect{
            .x = 0, .y = 0, .w = static_cast<int32_t>(pix.width), .h = static_cast<int32_t>(pix.height)};
        intersect_rect &= constixel::rect<int32_t>{.x = x, .y = y, .w = w, .h = h};
        if (intersect_rect.w <= 0 || intersect_rect.h <= 0) {
            return true;
        }
        return false;
    }

    constexpr bool clip_line(int32_t &x0, int32_t &y0, int32_t &x1, int32_t &y1, int32_t xmin, int32_t ymin,
                             int32_t xmax, int32_t ymax) {
        enum : uint32_t {
            INSIDE = 0,
            XMIN = 1,
            XMAX = 2,
            YMIN = 4,
            YMAX = 8
        };

        auto calc_code = [&](int32_t x, int32_t y) {
            uint32_t code = INSIDE;
            if (x < xmin) {
                code |= XMIN;
            } else if (x > xmax) {
                code |= XMAX;
            }
            if (y < ymin) {
                code |= YMIN;
            } else if (y > ymax) {
                code |= YMAX;
            }
            return code;
        };

        uint32_t outcode0 = calc_code(x0, y0);
        uint32_t outcode1 = calc_code(x1, y1);
        auto clip_loop = [&]<typename I>() -> bool {
            for (size_t i = 0; i < 4; i++) {
                if ((outcode0 | outcode1) == INSIDE) {
                    return true;
                }
                if ((outcode0 & outcode1) != 0) {
                    return false;
                }
                uint32_t outcode_out = outcode1 > outcode0 ? outcode1 : outcode0;
                int32_t x = 0;
                int32_t y = 0;
                if ((outcode_out & YMAX) != 0) {
                    const auto x1x0 = static_cast<I>(x1 - x0);
                    const auto w1y0 = static_cast<I>(ymax - y0);
                    const auto y1y0 = static_cast<I>(y1 - y0);
                    x = x0 + static_cast<int32_t>((x1x0 * w1y0) / y1y0);
                    y = ymax;
                } else if ((outcode_out & YMIN) != 0) {
                    const auto x1x0 = static_cast<I>(x1 - x0);
                    const auto ymy0 = static_cast<I>(ymin - y0);
                    const auto y1y0 = static_cast<I>(y1 - y0);
                    x = x0 + static_cast<int32_t>((x1x0 * ymy0) / y1y0);
                    y = ymin;
                } else if ((outcode_out & XMAX) != 0) {
                    const auto y1y0 = static_cast<I>(y1 - y0);
                    const auto w1x0 = static_cast<I>(xmax - x0);
                    const auto x1x0 = static_cast<I>(x1 - x0);
                    y = y0 + static_cast<int32_t>((y1y0 * w1x0) / x1x0);
                    x = xmax;
                } else {
                    const auto y1y0 = static_cast<I>(y1 - y0);
                    const auto xmx0 = static_cast<I>(xmin - x0);
                    const auto x1x0 = static_cast<I>(x1 - x0);
                    y = y0 + static_cast<int32_t>((y1y0 * xmx0) / x1x0);
                    x = xmin;
                }
                if (outcode_out == outcode0) {
                    x0 = x;
                    y0 = y;
                    outcode0 = calc_code(x0, y0);
                } else {
                    x1 = x;
                    y1 = y;
                    outcode1 = calc_code(x1, y1);
                }
            }
            return false;
        };
        return clip_loop.template operator()<calc_square_type>();
    }

    /**
     * @private
     * Linear RGB of a palette color.
     */
    [[nodiscard]] static constexpr std::array<float, 3> linear_color(uint8_t col) {
        const size_t index = (col & PF::color_mask) * 3;
        const auto &linear = PF::quant.linear_palette();
        return {linear.at(index + 0), linear.at(index + 1), linear.at(index + 2)};
    }

    /**
     * @private
     */
    template <typename abs_T>
    [[nodiscard]] static constexpr abs_T abs(abs_T v) {
        return v < 0 ? -v : v;
    }

    /**
     * @private
     */
    template <typename FONT>
    constexpr bool lookup_glyph(uint32_t utf32, size_t *glyph_index) {
        if (utf32 >= static_cast<uint32_t>(std::numeric_limits<typename FONT::lookup_type>::max()) - 1) {
            return false;
        }
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

    /**
     * @private
     */
    template <typename FONT>
    constexpr int32_t get_kerning(uint32_t utf32, const char *str, const char *end) const {
        if constexpr (FONT::kerning_tree.byte_size() > 0) {
            uint32_t utf_l = utf32;
            uint32_t utf_r = 0;
            static_cast<void>(next_utf32(str, end, &utf_r));
            auto amount = FONT::kerning_tree.lookup(
                static_cast<FONT::kerning_lookup_type>(utf_l << FONT::kerning_code_shift | utf_r));
            if (amount != FONT::kerning_tree.invalid) {
                return static_cast<int32_t>(static_cast<FONT::kerning_amount_type>(amount) -
                                            static_cast<FONT::kerning_amount_type>(FONT::kerning_amount_offset));
            }
        }
        return 0;
    }

    /**
     * @private
     */
    // Decodes one code point of [str, end) and returns the position after it. Never reads at or past end: a
    // sequence cut off by the end decodes as 0.
    [[nodiscard]] static constexpr const char *next_utf32(const char *str, const char *end, uint32_t *utf32) {
        if (str == end) {
            *utf32 = 0;
            return str;
        }
        const uint32_t lead = static_cast<uint32_t>(str[0]) & 0xFF;
        const auto left = static_cast<size_t>(end - str);
        if (lead < 0x80) {
            *utf32 = lead;
            return str + 1;
        }
        if ((lead >> 5) == 0x06 && left >= 2) {
            *utf32 = ((lead & 0x1F) << 6) | (static_cast<uint32_t>(str[1]) & 0x3F);
            return str + 2;
        }
        if ((lead >> 4) == 0x0E && left >= 3) {
            *utf32 = ((lead & 0x0F) << 12) | ((static_cast<uint32_t>(str[1]) & 0x3F) << 6) |
                     (static_cast<uint32_t>(str[2]) & 0x3F);
            return str + 3;
        }
        if ((lead >> 3) == 0x1E && left >= 4) {
            *utf32 = ((lead & 0x07) << 18) | ((static_cast<uint32_t>(str[1]) & 0x3F) << 12) |
                     ((static_cast<uint32_t>(str[2]) & 0x3F) << 6) | (static_cast<uint32_t>(str[3]) & 0x3F);
            return str + 4;
        }
        *utf32 = 0;
        return str + 1;
    }

    /**
     * @brief Plot a pixel, clipped. The pixel helpers take the raster as a parameter: the drawing loops pass a local
     * copy of pix, which byte stores into the pixel data cannot alias.
     */
    static constexpr void plot_clipped(const hidden::raster<uint8_t> &dst, int32_t x, int32_t y, uint8_t col) {
        auto wS = static_cast<int32_t>(dst.width);
        auto hS = static_cast<int32_t>(dst.height);
        if (x >= wS || y >= hS || x < 0 || y < 0) {
            return;
        }
        PF::plot(dst, static_cast<uint32_t>(x), static_cast<uint32_t>(y), col);
    }

    /**
     * @private
     */
    static constexpr void compose(const hidden::raster<uint8_t> &dst, int32_t x, int32_t y, float cola, float colr,
                                  float colg, float colb) {
        auto wS = static_cast<int32_t>(dst.width);
        auto hS = static_cast<int32_t>(dst.height);
        if (x >= wS || y >= hS || x < 0 || y < 0) {
            return;
        }
        PF::compose(dst, static_cast<uint32_t>(x), static_cast<uint32_t>(y), cola, colr, colg, colb);
    }

    /**
     * @private
     */
    static constexpr void compose_unsafe(const hidden::raster<uint8_t> &dst, int32_t x, int32_t y, float cola,
                                         float colr, float colg, float colb) {
        PF::compose(dst, static_cast<uint32_t>(x), static_cast<uint32_t>(y), cola, colr, colg, colb);
    }

    /**
     * @private
     */
    static constexpr void plot_unsafe(const hidden::raster<uint8_t> &dst, int32_t x, int32_t y, uint8_t col) {
        PF::plot(dst, static_cast<uint32_t>(x), static_cast<uint32_t>(y), col);
    }

    /**
     * @private
     */
    constexpr auto calculate_circle_bounds(int32_t cx, int32_t cy, int32_t r) const {
        struct bounds {
            int32_t x0, y0, x0r, x0r2, y0r, y0r2;
        };

        const int32_t x0 = std::max(cx - r - int32_t{1}, int32_t{0});
        const int32_t y0 = std::max(cy - r - int32_t{1}, int32_t{0});
        const int32_t x0r = std::clamp(x0 + r, int32_t{0}, static_cast<int32_t>(pix.width) - int32_t{1});
        const int32_t x0r2 = std::clamp(x0 + r * int32_t{2}, int32_t{0}, static_cast<int32_t>(pix.width) - int32_t{1});
        const int32_t y0r = std::clamp(y0 + r, int32_t{0}, static_cast<int32_t>(pix.height) - int32_t{1});
        const int32_t y0r2 = std::clamp(y0 + r * int32_t{2}, int32_t{0}, static_cast<int32_t>(pix.height) - int32_t{1});

        return bounds{x0, y0, x0r, x0r2, y0r, y0r2};
    }

    /**
     * @private
     */
    template <typename shader_func>
    constexpr auto fill_rect_int(int32_t x, int32_t y, int32_t w, int32_t h, const shader_func &shader,
                                 int32_t parent_x, int32_t parent_y, int32_t parent_w, int32_t parent_h) -> void
        requires std::is_invocable_r_v<std::array<float, 4>, shader_func, float, float, int32_t, int32_t>
    {
        const auto dst = pix;
        auto minmax_check = std::minmax({x, y, w, h});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return;
        }
        if (w <= 0 || h <= 0) {
            return;
        }
        if (check_not_in_bounds(x, y, w, h)) {
            return;
        }

        int32_t x0 = std::max(x, int32_t{0});
        int32_t y0 = std::max(y, int32_t{0});
        int32_t x1 = std::min(x + w, static_cast<int32_t>(pix.width));
        int32_t y1 = std::min(y + h, static_cast<int32_t>(pix.height));

        const auto parent_wF = static_cast<float>(parent_w);
        const auto parent_hF = static_cast<float>(parent_h);

        for (int32_t py = y0; py < y1; py++) {
            for (int32_t px = x0; px < x1; px++) {
                float u = (static_cast<float>(px) - static_cast<float>(parent_x)) / parent_wF;
                float v = (static_cast<float>(py) - static_cast<float>(parent_y)) / parent_hF;

                auto rgba = shader(u, v, px, py);
                for (auto &p : rgba) {
                    p = std::clamp(p, 0.0f, 1.0f);
                }
                compose(dst, px, py, rgba[3], rgba[0], rgba[1], rgba[2]);
            }
        }
    }

    /**
     * @private
     */
    template <bool AA, bool STROKE>
    constexpr void circle_int(int32_t cx, int32_t cy, int32_t r, int32_t ox, int32_t oy, uint8_t col, int32_t s) {
        auto bounds = calculate_circle_bounds(cx, cy, r);

        if (check_not_in_bounds(bounds.x0, bounds.y0, r * int32_t{2} + ox + int32_t{1},
                                r * int32_t{2} + oy + int32_t{1})) {
            return;
        }

        auto for_each_quadrant = [&]<typename I>(auto &&plot_arc) {
            plot_arc.template operator()<I>(bounds.x0, bounds.y0, bounds.x0r, bounds.y0r, 0, 0);
            plot_arc.template operator()<I>(bounds.x0r, bounds.y0, bounds.x0r2, bounds.y0r, ox, 0);
            plot_arc.template operator()<I>(bounds.x0, bounds.y0r, bounds.x0r, bounds.y0r2, 0, oy);
            plot_arc.template operator()<I>(bounds.x0r, bounds.y0r, bounds.x0r2, bounds.y0r2, ox, oy);
        };

        auto limit_box = [this](int32_t &xmin, int32_t &ymin, int32_t &xmax, int32_t &ymax, int32_t x_off,
                                int32_t y_off) {
            xmin = std::max(xmin + x_off, int32_t{0}) - x_off;
            xmax = std::min(xmax + x_off, static_cast<int32_t>(pix.width) - int32_t{1}) - x_off;
            ymin = std::max(ymin + y_off, int32_t{0}) - y_off;
            ymax = std::min(ymax + y_off, static_cast<int32_t>(pix.height) - int32_t{1}) - y_off;
        };

        if constexpr (!AA) {
            if constexpr (!STROKE) {
                auto plot_arc = [&, this]<typename I>(int32_t xx0, int32_t yy0, int32_t xx1, int32_t yy1, int32_t x_off,
                                                      int32_t y_off) {
                    const auto dst = pix;
                    limit_box(xx0, yy0, xx1, yy1, x_off, y_off);
                    for (int32_t y = yy0; y <= yy1; y++) {
                        for (int32_t x = xx0; x <= xx1; x++) {
                            const I dx = (static_cast<I>(x) * I{2} + I{1}) - (static_cast<I>(cx) * I{2});
                            const I dy = (static_cast<I>(y) * I{2} + I{1}) - (static_cast<I>(cy) * I{2});
                            const I dist_sq = dx * dx + dy * dy;
                            if (dist_sq > (static_cast<I>(r) * static_cast<I>(r) * I{4} - I{3})) {
                                continue;
                            }
                            plot_unsafe(dst, x + x_off, y + y_off, col);
                        }
                    }
                };
                for_each_quadrant.template operator()<calc_square_type>(plot_arc);
            } else {
                auto plot_arc = [&, this]<typename I>(int32_t xx0, int32_t yy0, int32_t xx1, int32_t yy1, int32_t x_off,
                                                      int32_t y_off) {
                    const auto dst = pix;
                    limit_box(xx0, yy0, xx1, yy1, x_off, y_off);
                    for (int32_t y = yy0; y <= yy1; y++) {
                        for (int32_t x = xx0; x <= xx1; x++) {
                            const I dx = (static_cast<I>(x) * I{2} + I{1}) - (static_cast<I>(cx) * I{2});
                            const I dy = (static_cast<I>(y) * I{2} + I{1}) - (static_cast<I>(cy) * I{2});
                            const I dist_sq = dx * dx + dy * dy;
                            if (dist_sq > (static_cast<I>(r) * static_cast<I>(r) * I{4} - I{3})) {
                                continue;
                            }
                            if (dist_sq < (static_cast<I>(r - s) * static_cast<I>(r - s) * I{4} - I{3})) {
                                continue;
                            }
                            plot_unsafe(dst, x + x_off, y + y_off, col);
                        }
                    }
                };
                for_each_quadrant.template operator()<calc_square_type>(plot_arc);
            }
        } else {
            const auto rF = static_cast<float>(r);
            const auto [Rl, Gl, Bl] = linear_color(col);
            if constexpr (!STROKE) {
                auto plot_arc = [&, this]<typename I>(int32_t xx0, int32_t yy0, int32_t xx1, int32_t yy1, int32_t x_off,
                                                      int32_t y_off) {
                    const auto dst = pix;
                    limit_box(xx0, yy0, xx1, yy1, x_off, y_off);
                    for (int32_t y = yy0; y <= yy1; y++) {
                        for (int32_t x = xx0; x <= xx1; x++) {
                            const float dx = (static_cast<float>(x) + 0.5f) - static_cast<float>(cx);
                            const float dy = (static_cast<float>(y) + 0.5f) - static_cast<float>(cy);
                            const float dist_sq = dx * dx + dy * dy;
                            if (dist_sq > (rF + 0.5f) * (rF + 0.5f)) {
                                continue;
                            }
                            if (dist_sq < (rF - 0.5f) * (rF - 0.5f)) {
                                plot_unsafe(dst, x + x_off, y + y_off, col);
                                continue;
                            }
                            float a = rF;
                            if (std::is_constant_evaluated()) {
                                a -= hidden::fast_sqrtf(dist_sq);
                            } else {
                                a -= std::sqrt(dist_sq);
                            }
                            a = std::clamp(a + 0.5f, 0.0f, 1.0f);
                            if (a >= hidden::epsilon_low) {
                                if (a >= hidden::epsilon_high) {
                                    plot_unsafe(dst, x + x_off, y + y_off, col);
                                } else {
                                    compose_unsafe(dst, x + x_off, y + y_off, a, Rl, Gl, Bl);
                                }
                            }
                        }
                    }
                };
                for_each_quadrant.template operator()<float>(plot_arc);
            } else {
                const auto rsF = static_cast<float>(r - s);
                auto plot_arc = [&, this]<typename I>(int32_t xx0, int32_t yy0, int32_t xx1, int32_t yy1, int32_t x_off,
                                                      int32_t y_off) {
                    const auto dst = pix;
                    limit_box(xx0, yy0, xx1, yy1, x_off, y_off);
                    for (int32_t y = yy0; y <= yy1; y++) {
                        for (int32_t x = xx0; x <= xx1; x++) {
                            const float dx = (static_cast<float>(x) + 0.5f) - static_cast<float>(cx);
                            const float dy = (static_cast<float>(y) + 0.5f) - static_cast<float>(cy);
                            const float dist_sq = dx * dx + dy * dy;
                            if (dist_sq > ((rF + 0.5f) * (rF + 0.5f))) {
                                continue;
                            }
                            if (dist_sq < ((rsF - 0.5f) * (rsF - 0.5f))) {
                                continue;
                            }
                            if (dist_sq < ((rsF + 0.5f) * (rsF + 0.5f))) {
                                float a = rsF;
                                if (std::is_constant_evaluated()) {
                                    a -= hidden::fast_sqrtf(dist_sq);
                                } else {
                                    a -= std::sqrt(dist_sq);
                                }
                                a = 1.0f - std::clamp(a + 0.5f, 0.0f, 1.0f);
                                if (a >= hidden::epsilon_low) {
                                    if (a >= hidden::epsilon_high) {
                                        plot_unsafe(dst, x + x_off, y + y_off, col);
                                    } else {
                                        compose_unsafe(dst, x + x_off, y + y_off, a, Rl, Gl, Bl);
                                    }
                                } else {
                                    plot_unsafe(dst, x + x_off, y + y_off, col);
                                }
                            } else {
                                float a = rF;
                                if (std::is_constant_evaluated()) {
                                    a -= hidden::fast_sqrtf(dist_sq);
                                } else {
                                    a -= std::sqrt(dist_sq);
                                }
                                a = std::clamp(a + 0.5f, 0.0f, 1.0f);
                                if (a >= hidden::epsilon_low) {
                                    if (a >= hidden::epsilon_high) {
                                        plot_unsafe(dst, x + x_off, y + y_off, col);
                                    } else {
                                        compose_unsafe(dst, x + x_off, y + y_off, a, Rl, Gl, Bl);
                                    }
                                }
                            }
                        }
                    }
                };
                for_each_quadrant.template operator()<float>(plot_arc);
            }
        }
    }

    /**
     * @private
     */
    template <typename shader_func>
    constexpr void circle_int_shader(int32_t cx, int32_t cy, int32_t r, int32_t ox, int32_t oy,
                                     const shader_func &shader, int32_t parent_x = -1, int32_t parent_y = -1,
                                     int32_t parent_w = -1, int32_t parent_h = -1) {
        auto bounds = calculate_circle_bounds(cx, cy, r);

        if (check_not_in_bounds(bounds.x0, bounds.y0, r * int32_t{2} + ox + int32_t{1},
                                r * int32_t{2} + oy + int32_t{1})) {
            return;
        }

        auto for_each_quadrant = [&]<typename I>(auto &&plot_arc) {
            plot_arc.template operator()<I>(bounds.x0, bounds.y0, bounds.x0r, bounds.y0r, 0, 0);
            plot_arc.template operator()<I>(bounds.x0r, bounds.y0, bounds.x0r2, bounds.y0r, ox, 0);
            plot_arc.template operator()<I>(bounds.x0, bounds.y0r, bounds.x0r, bounds.y0r2, 0, oy);
            plot_arc.template operator()<I>(bounds.x0r, bounds.y0r, bounds.x0r2, bounds.y0r2, ox, oy);
        };

        auto limit_box = [this](int32_t &xmin, int32_t &ymin, int32_t &xmax, int32_t &ymax, int32_t x_off,
                                int32_t y_off) {
            xmin = std::max(xmin + x_off, int32_t{0}) - x_off;
            xmax = std::min(xmax + x_off, static_cast<int32_t>(pix.width) - int32_t{1}) - x_off;
            ymin = std::max(ymin + y_off, int32_t{0}) - y_off;
            ymax = std::min(ymax + y_off, static_cast<int32_t>(pix.height) - int32_t{1}) - y_off;
        };

        const auto rF = static_cast<float>(r);
        const int32_t actual_parent_x = (parent_x == -1) ? cx - r : parent_x;
        const int32_t actual_parent_y = (parent_y == -1) ? cy - r : parent_y;
        const int32_t actual_parent_w = (parent_w == -1) ? r * 2 : parent_w;
        const int32_t actual_parent_h = (parent_h == -1) ? r * 2 : parent_h;

        const auto parent_wF = static_cast<float>(actual_parent_w);
        const auto parent_hF = static_cast<float>(actual_parent_h);

        auto plot_arc = [&, this]<typename I>(int32_t xx0, int32_t yy0, int32_t xx1, int32_t yy1, int32_t x_off,
                                              int32_t y_off) {
            const auto dst = pix;
            limit_box(xx0, yy0, xx1, yy1, x_off, y_off);
            for (int32_t y = yy0; y <= yy1; y++) {
                for (int32_t x = xx0; x <= xx1; x++) {
                    const float dx = (static_cast<float>(x) + 0.5f) - static_cast<float>(cx);
                    const float dy = (static_cast<float>(y) + 0.5f) - static_cast<float>(cy);
                    const float dist_sq = dx * dx + dy * dy;
                    if (dist_sq > (rF + 0.5f) * (rF + 0.5f)) {
                        continue;
                    }

                    float u = (static_cast<float>(x + x_off) - static_cast<float>(actual_parent_x)) / parent_wF;
                    float v = (static_cast<float>(y + y_off) - static_cast<float>(actual_parent_y)) / parent_hF;

                    auto rgba = shader(u, v, x + x_off, y + y_off);
                    for (auto &p : rgba) {
                        p = std::clamp(p, 0.0f, 1.0f);
                    }

                    if (dist_sq < (rF - 0.5f) * (rF - 0.5f)) {
                        compose_unsafe(dst, x + x_off, y + y_off, rgba[3], rgba[0], rgba[1], rgba[2]);
                        continue;
                    }

                    float a = rF;
                    if (std::is_constant_evaluated()) {
                        a -= hidden::fast_sqrtf(dist_sq);
                    } else {
                        a -= std::sqrt(dist_sq);
                    }
                    a = std::clamp(a + 0.5f, 0.0f, 1.0f);
                    if (a >= hidden::epsilon_low) {
                        compose_unsafe(dst, x + x_off, y + y_off, a * rgba[3], rgba[0], rgba[1], rgba[2]);
                    }
                }
            }
        };
        for_each_quadrant.template operator()<float>(plot_arc);
    }

    /**
     * @private
     * Draws a string with a mono or an antialiased font, starting at x/y.
     */
    template <typename FONT, bool KERNING, text_rotation ROTATION>
    constexpr int32_t draw_text(int32_t x, int32_t y, std::string_view str, uint8_t col, size_t character_count,
                                size_t *character_actual) {
        static_assert(FONT::mono || PF::bits_per_pixel >= 4,
                      "An antialiased font needs a format which can blend, 4-bit or more. Use a mono font.");
        auto minmax_check = std::minmax({x, y});
        if (minmax_check.first < min_coord || minmax_check.second > max_coord) {
            return 0;
        }
        const char *ptr = str.data();
        const char *end = ptr + str.size();
        size_t count = 0;
        while (ptr != end && count++ < character_count) {
            uint32_t utf32 = 0;
            ptr = next_utf32(ptr, end, &utf32);
            size_t index = 0;
            if (lookup_glyph<FONT>(utf32, &index)) {
                const char_info<typename FONT::char_info_type> &ch_info = FONT::char_table.at(index);
                draw_char<FONT, ROTATION>(x, y, ch_info, col);
                if constexpr (ROTATION == DEGREE_0) {
                    x += static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        x += get_kerning<FONT>(utf32, ptr, end);
                    }
                } else if constexpr (ROTATION == DEGREE_90) {
                    y += static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        y += get_kerning<FONT>(utf32, ptr, end);
                    }
                } else if constexpr (ROTATION == DEGREE_180) {
                    x -= static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        x -= get_kerning<FONT>(utf32, ptr, end);
                    }
                } else if constexpr (ROTATION == DEGREE_270) {
                    y -= static_cast<int32_t>(ch_info.xadvance);
                    if constexpr (KERNING) {
                        y -= get_kerning<FONT>(utf32, ptr, end);
                    }
                }
            }
        }
        if (character_actual) {
            *character_actual = count;
        }
        if constexpr (ROTATION == DEGREE_90 || ROTATION == DEGREE_270) {
            return y;
        } else {
            return x;
        }
    }

    /**
     * @private
     * Draws one glyph of a mono (1 bit per pixel) or an antialiased (4 bits per pixel) font.
     */
    template <typename FONT, text_rotation ROTATION>
    constexpr void draw_char(int32_t x, int32_t y, const char_info<typename FONT::char_info_type> &ch, uint8_t col) {
        const auto dst = pix;
        constexpr int32_t ppb = FONT::mono ? 8 : 2;  // glyph pixels per byte
        int32_t ch_data_off = static_cast<int32_t>(ch.y) * static_cast<int32_t>(FONT::glyph_bitmap_stride) +
                              static_cast<int32_t>(ch.x) / ppb;
        std::array<float, 3> lin{};
        if constexpr (!FONT::mono) {
            lin = linear_color(col);
        }

        auto limit_box = [&](int32_t &xmin, int32_t &ymin, int32_t &xmax, int32_t &ymax) {
            ymin = std::max(ymin, int32_t{0});
            ymax = std::min(ymax, static_cast<int32_t>(pix.height) - int32_t{1});
            xmin = std::max(xmin, int32_t{0});
            xmax = std::min(xmax, static_cast<int32_t>(pix.width) - int32_t{1});
        };

        // Draws pixel x_off of the current glyph line at xx/yy.
        auto put = [&](int32_t xx, int32_t yy, int32_t x_off) {
            const int32_t bit_index = (ppb - 1 - (x_off % ppb)) * (8 / ppb);
            const auto byte_index = static_cast<size_t>(ch_data_off + x_off / ppb);
            if (byte_index < (FONT::glyph_bitmap_stride * FONT::glyph_bitmap_height)) {
                constexpr uint8_t full = FONT::mono ? 1 : 0xF;
                const auto a = static_cast<uint8_t>((FONT::glyph_bitmap[byte_index] >> bit_index) & full);
                if (a == full) {
                    plot_unsafe(dst, xx, yy, col);
                } else if constexpr (!FONT::mono) {
                    if (a != 0) {
                        compose_unsafe(dst, xx, yy, hidden::a2al_4bit[a], lin[0], lin[1], lin[2]);
                    }
                }
            }
        };

        if constexpr (ROTATION == DEGREE_0) {
            if (check_not_in_bounds(x + ch.xoffset, y + ch.yoffset, ch.width + 1, ch.height + 1)) {
                return;
            }

            x += ch.xoffset;
            y += ch.yoffset;
            int32_t x2 = x + static_cast<int32_t>(ch.width);
            int32_t y2 = y + static_cast<int32_t>(ch.height);
            limit_box(x, y, x2, y2);
            for (int32_t yy = y; yy < y2; yy++) {
                for (int32_t xx = x; xx < x2; xx++) {
                    put(xx, yy, (xx - x) + ch.x % ppb);
                }
                ch_data_off += FONT::glyph_bitmap_stride;
            }
        } else if constexpr (ROTATION == DEGREE_180) {
            if (check_not_in_bounds(x - ch.xoffset - ch.width, y - FONT::ascent, ch.width + 1, ch.height + 1)) {
                return;
            }

            x -= ch.xoffset;
            y -= FONT::ascent;
            int32_t x2 = x - static_cast<int32_t>(ch.width);
            int32_t y2 = y + static_cast<int32_t>(ch.height);
            limit_box(x2, y, x, y2);
            for (int32_t yy = y2 - 1; yy >= y; yy--) {
                for (int32_t xx = x2; xx < x; xx++) {
                    put(xx, yy, (ch.width - (xx - x2) - 1) + ch.x % ppb);
                }
                ch_data_off += FONT::glyph_bitmap_stride;
            }
        } else if constexpr (ROTATION == DEGREE_90) {
            if (check_not_in_bounds(x - FONT::ascent, y, ch.height + 1, ch.width + 1)) {
                return;
            }

            x -= FONT::ascent;
            y += ch.xoffset;
            int32_t x2 = x + static_cast<int32_t>(ch.height);
            int32_t y2 = y + static_cast<int32_t>(ch.width);
            limit_box(x, y, x2, y2);
            for (int32_t xx = x2 - 1; xx >= x; xx--) {
                for (int32_t yy = y; yy < y2; yy++) {
                    put(xx, yy, (yy - y) + ch.x % ppb);
                }
                ch_data_off += FONT::glyph_bitmap_stride;
            }
        } else if constexpr (ROTATION == DEGREE_270) {
            if (check_not_in_bounds(x + ch.yoffset, y - ch.xoffset - ch.width, ch.height + 1, ch.width + 1)) {
                return;
            }

            x += ch.yoffset;
            y -= ch.xoffset;
            int32_t x2 = x + static_cast<int32_t>(ch.height);
            int32_t y2 = y - static_cast<int32_t>(ch.width);
            limit_box(x, y2, x2, y);
            for (int32_t xx = x; xx < x2; xx++) {
                for (int32_t yy = y2; yy < y; yy++) {
                    put(xx, yy, (ch.width - (yy - y2) - 1) + ch.x % ppb);
                }
                ch_data_off += FONT::glyph_bitmap_stride;
            }
        }
    }


    /**
     * @private
     */
    static constexpr void extent(const hidden::raster<uint8_t> &dst, int32_t x, int32_t w, int32_t y, uint8_t col) {
        auto wS = static_cast<int32_t>(dst.width);
        auto hS = static_cast<int32_t>(dst.height);
        if (w <= 0 || y < 0 || y >= hS || x + w <= 0 || x >= wS) {
            return;
        }
        const int32_t x0 = std::max(x, int32_t{0});
        const int32_t x1 = std::min(x + w, static_cast<int32_t>(dst.width));
        PF::extent(dst, static_cast<size_t>(x0), static_cast<size_t>(x1), static_cast<size_t>(y), col);
    }
    /// @endcond  // DOXYGEN_EXCLUDE
#endif  // #ifndef __INTELLISENSE__
};

/**
 * @class surface_api
 * @brief Drawing API of image and dynamic_image. Every call builds a surface from view() on the fly and forwards to
 * it, so no pointer is stored in the image.
 *
 * @tparam D Storage class, provides view() and view() const.
 * @tparam PF Pixel class of the format, see basic_surface.
 */
template <typename D, typename PF>
class surface_api {
 public:
    /// @copydoc surface::plot(int32_t,int32_t,uint8_t)
    constexpr void plot(int32_t x, int32_t y, uint8_t col) {
        derived().view().plot(x, y, col);
    }

    /// @copydoc surface::plot(const struct plot &)
    constexpr void plot(const struct plot &p) {
        derived().view().plot(p);
    }

    /// @copydoc surface::draw_line(int32_t,int32_t,int32_t,int32_t,uint8_t,int32_t)
    constexpr void draw_line(int32_t x0, int32_t y0, int32_t x1, int32_t y1, uint8_t col, int32_t stroke_width = 1) {
        derived().view().draw_line(x0, y0, x1, y1, col, stroke_width);
    }

    /// @copydoc surface::draw_line(const struct draw_line &)
    constexpr void draw_line(const struct draw_line &d) {
        derived().view().draw_line(d);
    }

    /// @copydoc surface::draw_line_aa(int32_t,int32_t,int32_t,int32_t,uint8_t,float)
    constexpr void draw_line_aa(int32_t x0, int32_t y0, int32_t x1, int32_t y1, uint8_t col,
                                float stroke_width = 1.0f) {
        derived().view().draw_line_aa(x0, y0, x1, y1, col, stroke_width);
    }

    /// @copydoc surface::draw_line_aa(const struct draw_line &)
    constexpr void draw_line_aa(const struct draw_line &d) {
        derived().view().draw_line_aa(d);
    }

    /// @copydoc surface::string_width(const char *,size_t,size_t *)
    template <typename FONT, bool KERNING = false>
    [[nodiscard]] constexpr int32_t string_width(const char *str,
                                                 size_t character_count = std::numeric_limits<size_t>::max(),
                                                 size_t *character_actual = nullptr) {
        return derived().view().template string_width<FONT, KERNING>(str, character_count, character_actual);
    }

    /// @copydoc surface::string_width(std::string_view,size_t,size_t *)
    template <typename FONT, bool KERNING = false>
    [[nodiscard]] constexpr int32_t string_width(std::string_view str,
                                                 size_t character_count = std::numeric_limits<size_t>::max(),
                                                 size_t *character_actual = nullptr) {
        return derived().view().template string_width<FONT, KERNING>(str, character_count, character_actual);
    }

    /// @copydoc surface::draw_string(int32_t,int32_t,std::string_view,uint8_t,text_align,size_t,size_t *)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string(int32_t x, int32_t y, std::string_view str, uint8_t col,
                                  text_align align = ALIGN_LEFT,
                                  size_t character_count = std::numeric_limits<size_t>::max(),
                                  size_t *character_actual = nullptr) {
        return derived().view().template draw_string<FONT, KERNING, ROTATION>(x, y, str, col, align, character_count,
                                                                              character_actual);
    }

    /// @copydoc surface::draw_string(int32_t,int32_t,const char *,uint8_t,text_align,size_t,size_t *)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string(int32_t x, int32_t y, const char *str, uint8_t col, text_align align = ALIGN_LEFT,
                                  size_t character_count = std::numeric_limits<size_t>::max(),
                                  size_t *character_actual = nullptr) {
        return derived().view().template draw_string<FONT, KERNING, ROTATION>(x, y, str, col, align, character_count,
                                                                              character_actual);
    }

    /// @copydoc surface::draw_string(const struct draw_string &,size_t,size_t *)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string(const struct draw_string &d,
                                  size_t character_count = std::numeric_limits<size_t>::max(),
                                  size_t *character_actual = nullptr) {
        return derived().view().template draw_string<FONT, KERNING, ROTATION>(d, character_count, character_actual);
    }

    /// @copydoc surface::draw_string_mono(int32_t,int32_t,const char *,uint8_t,size_t,size_t *)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_mono(int32_t x, int32_t y, const char *str, uint8_t col,
                                       size_t character_count = std::numeric_limits<size_t>::max(),
                                       size_t *character_actual = nullptr) {
        return derived().view().template draw_string_mono<FONT, KERNING, ROTATION>(x, y, str, col, character_count,
                                                                                   character_actual);
    }

    /// @copydoc surface::draw_string_mono(int32_t,int32_t,std::string_view,uint8_t,size_t,size_t *)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_mono(int32_t x, int32_t y, std::string_view str, uint8_t col,
                                       size_t character_count = std::numeric_limits<size_t>::max(),
                                       size_t *character_actual = nullptr) {
        return derived().view().template draw_string_mono<FONT, KERNING, ROTATION>(x, y, str, col, character_count,
                                                                                   character_actual);
    }

    /// @copydoc surface::draw_string_mono(const struct draw_string &,size_t,size_t *)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_mono(const struct draw_string &d,
                                       size_t character_count = std::numeric_limits<size_t>::max(),
                                       size_t *character_actual = nullptr) {
        return derived().view().template draw_string_mono<FONT, KERNING, ROTATION>(d, character_count,
                                                                                   character_actual);
    }

    /// @copydoc surface::draw_string_centered_mono(int32_t,int32_t,const char *,uint8_t)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr void draw_string_centered_mono(int32_t x, int32_t y, const char *str, uint8_t col) {
        derived().view().template draw_string_centered_mono<FONT, KERNING, ROTATION>(x, y, str, col);
    }

    /// @copydoc surface::draw_string_centered_mono(int32_t,int32_t,std::string_view,uint8_t)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr void draw_string_centered_mono(int32_t x, int32_t y, std::string_view str, uint8_t col) {
        derived().view().template draw_string_centered_mono<FONT, KERNING, ROTATION>(x, y, str, col);
    }

    /// @copydoc surface::draw_string_aa(int32_t,int32_t,const char *,uint8_t,size_t,size_t *)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_aa(int32_t x, int32_t y, const char *str, uint8_t col,
                                     size_t character_count = std::numeric_limits<size_t>::max(),
                                     size_t *character_actual = nullptr) {
        return derived().view().template draw_string_aa<FONT, KERNING, ROTATION>(x, y, str, col, character_count,
                                                                                 character_actual);
    }

    /// @copydoc surface::draw_string_aa(int32_t,int32_t,std::string_view,uint8_t,size_t,size_t *)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_aa(int32_t x, int32_t y, std::string_view str, uint8_t col,
                                     size_t character_count = std::numeric_limits<size_t>::max(),
                                     size_t *character_actual = nullptr) {
        return derived().view().template draw_string_aa<FONT, KERNING, ROTATION>(x, y, str, col, character_count,
                                                                                 character_actual);
    }

    /// @copydoc surface::draw_string_aa(const struct draw_string &,size_t,size_t *)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr int32_t draw_string_aa(const struct draw_string &d,
                                     size_t character_count = std::numeric_limits<size_t>::max(),
                                     size_t *character_actual = nullptr) {
        return derived().view().template draw_string_aa<FONT, KERNING, ROTATION>(d, character_count, character_actual);
    }

    /// @copydoc surface::draw_string_centered_aa(int32_t,int32_t,const char *,uint8_t)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr void draw_string_centered_aa(int32_t x, int32_t y, const char *str, uint8_t col) {
        derived().view().template draw_string_centered_aa<FONT, KERNING, ROTATION>(x, y, str, col);
    }

    /// @copydoc surface::draw_string_centered_aa(int32_t,int32_t,std::string_view,uint8_t)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr void draw_string_centered_aa(int32_t x, int32_t y, std::string_view str, uint8_t col) {
        derived().view().template draw_string_centered_aa<FONT, KERNING, ROTATION>(x, y, str, col);
    }

    /// @copydoc surface::fill_rect(int32_t,int32_t,int32_t,int32_t,uint8_t)
    constexpr void fill_rect(int32_t x, int32_t y, int32_t w, int32_t h, uint8_t col) {
        derived().view().fill_rect(x, y, w, h, col);
    }

    /// @copydoc surface::fill_rect(const struct draw_rect &)
    constexpr void fill_rect(const struct draw_rect &d) {
        derived().view().fill_rect(d);
    }

    /// @copydoc surface::fill_rect(int32_t,int32_t,int32_t,int32_t,const shader_func &)
    template <typename shader_func>
    constexpr auto fill_rect(int32_t x, int32_t y, int32_t w, int32_t h, const shader_func &shader) -> void
        requires std::is_invocable_r_v<std::array<float, 4>, shader_func, float, float, int32_t, int32_t>
    {
        derived().view().fill_rect(x, y, w, h, shader);
    }

    /// @copydoc surface::stroke_rect(int32_t,int32_t,int32_t,int32_t,uint8_t,int32_t)
    constexpr void stroke_rect(int32_t x, int32_t y, int32_t w, int32_t h, uint8_t col, int32_t stroke_width = 1) {
        derived().view().stroke_rect(x, y, w, h, col, stroke_width);
    }

    /// @copydoc surface::stroke_rect(const struct draw_rect &)
    constexpr void stroke_rect(const struct draw_rect &d) {
        derived().view().stroke_rect(d);
    }

    /// @copydoc surface::fill_circle(int32_t,int32_t,int32_t,uint8_t)
    constexpr void fill_circle(int32_t cx, int32_t cy, int32_t radius, uint8_t col) {
        derived().view().fill_circle(cx, cy, radius, col);
    }

    /// @copydoc surface::fill_circle(const struct draw_circle &)
    constexpr void fill_circle(const struct draw_circle &d) {
        derived().view().fill_circle(d);
    }

    /// @copydoc surface::stroke_circle(int32_t,int32_t,int32_t,uint8_t,int32_t)
    constexpr void stroke_circle(int32_t cx, int32_t cy, int32_t radius, uint8_t col, int32_t stroke_width = 1) {
        derived().view().stroke_circle(cx, cy, radius, col, stroke_width);
    }

    /// @copydoc surface::stroke_circle(const struct draw_circle &)
    constexpr void stroke_circle(const struct draw_circle &d) {
        derived().view().stroke_circle(d);
    }

    /// @copydoc surface::stroke_circle_aa(int32_t,int32_t,int32_t,uint8_t,int32_t)
    constexpr void stroke_circle_aa(int32_t cx, int32_t cy, int32_t radius, uint8_t col, int32_t stroke_width = 1) {
        derived().view().stroke_circle_aa(cx, cy, radius, col, stroke_width);
    }

    /// @copydoc surface::stroke_circle_aa(const struct draw_circle &)
    constexpr void stroke_circle_aa(const struct draw_circle &d) {
        derived().view().stroke_circle_aa(d);
    }

    /// @copydoc surface::fill_circle_aa(int32_t,int32_t,int32_t,uint8_t)
    constexpr void fill_circle_aa(int32_t cx, int32_t cy, int32_t radius, uint8_t col) {
        derived().view().fill_circle_aa(cx, cy, radius, col);
    }

    /// @copydoc surface::fill_circle_aa(const struct draw_circle &)
    constexpr void fill_circle_aa(const struct draw_circle &d) {
        derived().view().fill_circle_aa(d);
    }

    /// @copydoc surface::fill_circle_aa(int32_t,int32_t,int32_t,const shader_func &)
    template <typename shader_func>
    constexpr auto fill_circle_aa(int32_t cx, int32_t cy, int32_t radius, const shader_func &shader) -> void
        requires std::is_invocable_r_v<std::array<float, 4>, shader_func, float, float, int32_t, int32_t>
    {
        derived().view().fill_circle_aa(cx, cy, radius, shader);
    }

    /// @copydoc surface::stroke_round_rect(int32_t,int32_t,int32_t,int32_t,int32_t,uint8_t,int32_t)
    constexpr void stroke_round_rect(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius, uint8_t col,
                                     int32_t stroke_width = 1) {
        derived().view().stroke_round_rect(x, y, w, h, radius, col, stroke_width);
    }

    /// @copydoc surface::stroke_round_rect(const struct draw_round_rect &)
    constexpr void stroke_round_rect(const struct draw_round_rect &d) {
        derived().view().stroke_round_rect(d);
    }

    /// @copydoc surface::fill_round_rect(int32_t,int32_t,int32_t,int32_t,int32_t,uint8_t)
    constexpr void fill_round_rect(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius, uint8_t col) {
        derived().view().fill_round_rect(x, y, w, h, radius, col);
    }

    /// @copydoc surface::fill_round_rect(const struct draw_round_rect &)
    constexpr void fill_round_rect(const struct draw_round_rect &d) {
        derived().view().fill_round_rect(d);
    }

    /// @copydoc surface::fill_round_rect_aa(int32_t,int32_t,int32_t,int32_t,int32_t,uint8_t)
    constexpr void fill_round_rect_aa(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius, uint8_t col) {
        derived().view().fill_round_rect_aa(x, y, w, h, radius, col);
    }

    /// @copydoc surface::fill_round_rect_aa(const struct draw_round_rect &)
    constexpr void fill_round_rect_aa(const struct draw_round_rect &d) {
        derived().view().fill_round_rect_aa(d);
    }

    /// @copydoc surface::fill_round_rect_aa(int32_t,int32_t,int32_t,int32_t,int32_t,const shader_func &)
    template <typename shader_func>
    constexpr auto fill_round_rect_aa(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius,
                                      const shader_func &shader) -> void
        requires std::is_invocable_r_v<std::array<float, 4>, shader_func, float, float, int32_t, int32_t>
    {
        derived().view().fill_round_rect_aa(x, y, w, h, radius, shader);
    }

    /// @copydoc surface::stroke_round_rect_aa(int32_t,int32_t,int32_t,int32_t,int32_t,uint8_t,int32_t)
    constexpr void stroke_round_rect_aa(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius, uint8_t col,
                                        int32_t stroke_width = 1) {
        derived().view().stroke_round_rect_aa(x, y, w, h, radius, col, stroke_width);
    }

    /// @copydoc surface::stroke_round_rect_aa(const struct draw_round_rect &)
    constexpr void stroke_round_rect_aa(const struct draw_round_rect &d) {
        derived().view().stroke_round_rect_aa(d);
    }

    /// @copydoc surface::flip_h()
    constexpr void flip_h() {
        derived().view().flip_h();
    }

    /// @copydoc surface::flip_v()
    constexpr void flip_v() {
        derived().view().flip_v();
    }

    /// @copydoc surface::flip_hv()
    constexpr void flip_hv() {
        derived().view().flip_hv();
    }

    /// @copydoc surface::blit_RGBA(int32_t,int32_t,int32_t,int32_t,const uint8_t *,int32_t,int32_t,int32_t)
    constexpr void blit_RGBA(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih,
                             int32_t stride) {
        derived().view().blit_RGBA(x, y, w, h, ptr, iw, ih, stride);
    }

    /// @copydoc surface::blit_RGBA(const constixel::rect<int32_t> &,const uint8_t *,int32_t,int32_t,int32_t)
    constexpr void blit_RGBA(const constixel::rect<int32_t> &dstrect, const uint8_t *ptr, int32_t iw, int32_t ih,
                             int32_t stride) {
        derived().view().blit_RGBA(dstrect, ptr, iw, ih, stride);
    }

    /// @copydoc surface::blit_RGBA_diffused(int32_t,int32_t,int32_t,int32_t,const uint8_t *,int32_t,int32_t,int32_t)
    constexpr void blit_RGBA_diffused(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw,
                                      int32_t ih, int32_t stride) {
        derived().view().blit_RGBA_diffused(x, y, w, h, ptr, iw, ih, stride);
    }

    /// @copydoc surface::blit_RGBA_diffused(const constixel::rect<int32_t> &,const uint8_t *,int32_t,int32_t,int32_t)
    constexpr void blit_RGBA_diffused(const constixel::rect<int32_t> &dstrect, const uint8_t *ptr, int32_t iw,
                                      int32_t ih, int32_t stride) {
        derived().view().blit_RGBA_diffused(dstrect, ptr, iw, ih, stride);
    }

    /// @copydoc surface::blit_RGBA_diffused_linear(int32_t,int32_t,int32_t,int32_t,const uint8_t *,int32_t,int32_t,int32_t)
    constexpr void blit_RGBA_diffused_linear(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw,
                                             int32_t ih, int32_t stride) {
        derived().view().blit_RGBA_diffused_linear(x, y, w, h, ptr, iw, ih, stride);
    }

    /// @copydoc surface::blit_RGBA_diffused_linear(const constixel::rect<int32_t> &,const uint8_t *,int32_t,int32_t,int32_t)
    constexpr void blit_RGBA_diffused_linear(const constixel::rect<int32_t> &dstrect, const uint8_t *ptr, int32_t iw,
                                             int32_t ih, int32_t stride) {
        derived().view().blit_RGBA_diffused_linear(dstrect, ptr, iw, ih, stride);
    }

    /// @copydoc surface::rect(int32_t,int32_t,int32_t,int32_t)
    constexpr auto rect(int32_t x, int32_t y, int32_t w, int32_t h) {
        return derived().view().rect(x, y, w, h);
    }

    /// @copydoc surface::circle(int32_t,int32_t,int32_t)
    constexpr auto circle(int32_t cx, int32_t cy, int32_t r) {
        return derived().view().circle(cx, cy, r);
    }

    /// @copydoc surface::circle_aa(int32_t,int32_t,int32_t)
    constexpr auto circle_aa(int32_t cx, int32_t cy, int32_t r) {
        return derived().view().circle_aa(cx, cy, r);
    }

    /// @copydoc surface::round_rect(int32_t,int32_t,int32_t,int32_t,int32_t)
    constexpr auto round_rect(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius) {
        return derived().view().round_rect(x, y, w, h, radius);
    }

    /// @copydoc surface::round_rect_aa(int32_t,int32_t,int32_t,int32_t,int32_t)
    constexpr auto round_rect_aa(int32_t x, int32_t y, int32_t w, int32_t h, int32_t radius) {
        return derived().view().round_rect_aa(x, y, w, h, radius);
    }

    /// @copydoc surface::line(int32_t,int32_t,int32_t,int32_t)
    constexpr auto line(int32_t x0, int32_t y0, int32_t x1, int32_t y1) {
        return derived().view().line(x0, y0, x1, y1);
    }

    /// @copydoc surface::line_aa(int32_t,int32_t,int32_t,int32_t)
    constexpr auto line_aa(int32_t x0, int32_t y0, int32_t x1, int32_t y1) {
        return derived().view().line_aa(x0, y0, x1, y1);
    }

    /// @copydoc surface::point(int32_t,int32_t)
    constexpr auto point(int32_t x, int32_t y) {
        return derived().view().point(x, y);
    }

    /// @copydoc surface::text_mono(int32_t,int32_t,const char *)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_mono(int32_t x, int32_t y, const char *str) {
        return derived().view().template text_mono<FONT, KERNING, ROTATION>(x, y, str);
    }

    /// @copydoc surface::text_mono(int32_t,int32_t,std::string_view)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_mono(int32_t x, int32_t y, std::string_view str) {
        return derived().view().template text_mono<FONT, KERNING, ROTATION>(x, y, str);
    }

    /// @copydoc surface::text_aa(int32_t,int32_t,const char *)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_aa(int32_t x, int32_t y, const char *str) {
        return derived().view().template text_aa<FONT, KERNING, ROTATION>(x, y, str);
    }

    /// @copydoc surface::text_aa(int32_t,int32_t,std::string_view)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_aa(int32_t x, int32_t y, std::string_view str) {
        return derived().view().template text_aa<FONT, KERNING, ROTATION>(x, y, str);
    }

    /// @copydoc surface::text_centered_mono(int32_t,int32_t,const char *)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_centered_mono(int32_t x, int32_t y, const char *str) {
        return derived().view().template text_centered_mono<FONT, KERNING, ROTATION>(x, y, str);
    }

    /// @copydoc surface::text_centered_mono(int32_t,int32_t,std::string_view)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_centered_mono(int32_t x, int32_t y, std::string_view str) {
        return derived().view().template text_centered_mono<FONT, KERNING, ROTATION>(x, y, str);
    }

    /// @copydoc surface::text_centered_aa(int32_t,int32_t,const char *)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_centered_aa(int32_t x, int32_t y, const char *str) {
        return derived().view().template text_centered_aa<FONT, KERNING, ROTATION>(x, y, str);
    }

    /// @copydoc surface::text_centered_aa(int32_t,int32_t,std::string_view)
    template <typename FONT, bool KERNING = false, text_rotation ROTATION = DEGREE_0>
    constexpr auto text_centered_aa(int32_t x, int32_t y, std::string_view str) {
        return derived().view().template text_centered_aa<FONT, KERNING, ROTATION>(x, y, str);
    }

    /// @copydoc basic_surface::get_nearest_color(uint8_t,uint8_t,uint8_t)
    [[nodiscard]] constexpr uint8_t get_nearest_color(uint8_t r, uint8_t g, uint8_t b) const {
        return derived().view().get_nearest_color(r, g, b);
    }

    /// @copydoc basic_surface::png(F &&)
    template <typename F>
    constexpr void png(F &&char_out) const {
        derived().view().png(std::forward<F>(char_out));
    }

    /// @copydoc basic_surface::sixel(F &&)
    template <size_t S = 1, typename F>
    constexpr void sixel(F &&char_out) const {
        const auto surface = derived().view();
        sixel<S>(std::forward<F>(char_out), {0, 0, surface.width(), surface.height()});
    }

    /// @copydoc basic_surface::sixel(F &&,const constixel::rect<int32_t> &)
    template <size_t S = 1, typename F>
    constexpr void sixel(F &&char_out, const constixel::rect<int32_t> &rect) const {
        if constexpr (hidden::has_sixel_cache<PF>) {
            // D::sixel_cache() returns the storage for a line cache of the full width.
            auto line_cache = derived().sixel_cache();
            derived().view().template sixel<S>(std::forward<F>(char_out), rect, std::span<uint8_t>{line_cache});
        } else {
            derived().view().template sixel<S>(std::forward<F>(char_out), rect);
        }
    }

#ifdef CONSTIXEL_ENABLE_COUT
    /// @copydoc basic_surface::sixel_to_cout()
    template <size_t S = 1>
    void sixel_to_cout() const {
        std::string out;
        sixel<S>([&out](char ch) mutable {
            out.push_back(ch);
        });
        std::cout << out << '\n';
    }

    /// @copydoc basic_surface::vt100_clear()
    void vt100_clear() const {
        derived().view().vt100_clear();
    }

    /// @copydoc basic_surface::vt100_clear_scrollback()
    void vt100_clear_scrollback() const {
        derived().view().vt100_clear_scrollback();
    }

    /// @copydoc basic_surface::vt100_home()
    void vt100_home() const {
        derived().view().vt100_home();
    }

    /// @copydoc basic_surface::png_to_iterm()
    void png_to_iterm() const {
        derived().view().png_to_iterm();
    }

    /// @copydoc basic_surface::png_to_kitty()
    void png_to_kitty() const {
        derived().view().png_to_kitty();
    }
#endif  // #ifdef CONSTIXEL_ENABLE_COUT
    /// @copydoc basic_surface::convert(F &&)
    template <device_format dst_format, typename F>
    constexpr void convert(F &&uint8_out) const {
        derived().view().template convert<dst_format>(std::forward<F>(uint8_out));
    }

    /// @copydoc basic_surface::convert(F &&,const constixel::rect<int32_t> &)
    template <device_format dst_format, typename F>
    constexpr void convert(F &&uint8_out, const constixel::rect<int32_t> &rect) const {
        derived().view().template convert<dst_format>(std::forward<F>(uint8_out), rect);
    }

    /// @copydoc basic_surface::convert_size(const constixel::rect<int32_t> &)
    template <device_format dst_format>
    [[nodiscard]] constexpr size_t convert_size(const constixel::rect<int32_t> &rect) const {
        return derived().view().template convert_size<dst_format>(rect);
    }

    /// @copydoc basic_surface::convert_size()
    template <device_format dst_format>
    [[nodiscard]] constexpr size_t convert_size() const {
        return derived().view().template convert_size<dst_format>();
    }

    /// @copydoc basic_surface::convert_chunk(std::span<uint8_t>,convert_cursor &,const constixel::rect<int32_t> &)
    template <device_format dst_format>
    [[nodiscard]] constexpr size_t convert_chunk(std::span<uint8_t> dst, convert_cursor &cursor,
                                                 const constixel::rect<int32_t> &rect) const {
        return derived().view().template convert_chunk<dst_format>(dst, cursor, rect);
    }

    /// @copydoc basic_surface::convert_chunk(std::span<uint8_t>,convert_cursor &)
    template <device_format dst_format>
    [[nodiscard]] constexpr size_t convert_chunk(std::span<uint8_t> dst, convert_cursor &cursor) const {
        return derived().view().template convert_chunk<dst_format>(dst, cursor);
    }

 private:
    [[nodiscard]] constexpr D &derived() {
        return static_cast<D &>(*this);
    }
    [[nodiscard]] constexpr const D &derived() const {
        return static_cast<const D &>(*this);
    }
};

/**
 * @class image
 * @brief Core class of constixel. Holds a buffer of an image width a certain size and format. Typical use:
 *
 * \code{.cpp}
 * constixel::image<constixel::format_8bit, 640, 480> image;
 * \endcode
 *
 * @tparam T Type of the image buffer. One of format_1bit, format_2bit, format_4bit, format_8bit, format_24bit or
 * format_32bit.
 * @tparam W Width in pixels.
 * @tparam H Height in pixels.
 * @tparam GRAYSCALE boolean to indicate if palette should be grayscale. Otherwise a colored palette will be used.
 * @tparam USE_SPAN Pass in your own std::span in the constructor to be used as a back buffer.
 * Defaults to false.
 */
template <template <size_t, size_t, bool, bool> class T, size_t W, size_t H, bool GRAYSCALE = false,
          bool USE_SPAN = false>
class image
    : public surface_api<image<T, W, H, GRAYSCALE, USE_SPAN>, typename T<W, H, GRAYSCALE, USE_SPAN>::pixel_type> {
    static_assert(sizeof(W) >= sizeof(uint32_t));
    static_assert(sizeof(H) >= sizeof(uint32_t));

    static_assert(W > 0 && H > 0);
    static_assert(W <= 65535 && H <= 65535);

 public:
    /// @brief Pixel class of the format, the PF of surface.
    using pixel_type = typename T<W, H, GRAYSCALE, USE_SPAN>::pixel_type;

    /**
     * \brief Creates a new image with internal storage.
     */
    image()
        requires(!USE_SPAN)
    = default;

    /**
     * \brief When USE_SPAN=true creates a new image with external storage based existing std::span.
     * \param other If USE_SPAN=true this constructor will accept a std::span.
     */
    explicit image(const std::span<uint8_t, T<W, H, GRAYSCALE, USE_SPAN>::image_size> &other)
        requires(USE_SPAN)
        : data(other) {
    }

    /**
     * \brief Creates an image that uses a caller-supplied runtime palette
     *        instead of the format's compile-time-baked palette.
     * \param palette Up to 256 sRGB byte triples.
     */
    explicit image(std::span<const std::array<uint8_t, 3>> palette)
        requires(!USE_SPAN && requires { typename T<W, H, GRAYSCALE, USE_SPAN>::runtime_palette; })
        : palette_{runtime_palette(palette), true} {
    }

    /**
     * \brief USE_SPAN + custom-palette variant.
     */
    image(const std::span<uint8_t, T<W, H, GRAYSCALE, USE_SPAN>::image_size> &other,
          std::span<const std::array<uint8_t, 3>> palette)
        requires(USE_SPAN && requires { typename T<W, H, GRAYSCALE, USE_SPAN>::runtime_palette; })
        : data(other), palette_{runtime_palette(palette), true} {
    }

    /**
     * \brief Boolean indicating that the palette is grayscale instead of color.
     * \return If true, the palette is grayscale. If false a colored palette is used.
     */
    [[nodiscard]] static constexpr bool grayscale() {
        return GRAYSCALE;
    }

    /**
     * \brief Bit depth of the image.
     * \return Bit depth of the image. 1 == 2 colors, 2 == 4 colors, 4 == 16 colors, 8 == 256 colors.
     */
    [[nodiscard]] static constexpr size_t bit_depth() {
        return T<W, H, GRAYSCALE, USE_SPAN>::bits_per_pixel;
    }

    /**
     * \brief Size in bytes of the image data. This does not include the image instance size.
     * \return Size in bytes of the image data.
     */
    [[nodiscard]] static constexpr size_t size() {
        return T<W, H, GRAYSCALE, USE_SPAN>::image_size;
    }

    /**
     * \brief Bytes per line in the image data. Also called stride.
     * \return Bytes per line in the image data.
     */
    [[nodiscard]] static constexpr size_t bytes_per_line() {
        return T<W, H, GRAYSCALE, USE_SPAN>::bytes_per_line;
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
    [[nodiscard]] constexpr std::array<uint8_t, T<W, H, GRAYSCALE, USE_SPAN>::image_size> &data_ref() {
        return data;
    }

    /**
     * \brief Returns a clone of this image. Data is copied.
     * \return Cloned image instance.
     */
    [[nodiscard]] constexpr image<T, W, H, GRAYSCALE> clone() const {
        return *this;
    }

    /**
     * \brief Copy source image into this instance. No compositing occurs, replaces current content.
     * \param src Source image.
     */
    constexpr void copy(const image<T, W, H, GRAYSCALE> &src) {
        static_cast<void>(view().copy(src.view()));
    }

    /**
     * \brief Copy raw source data into this instance. No compositing occurs.
     * \tparam BYTE_SIZE Amount data in the source data. Typically a sizeof() of an array. Must match image::size()
     * \param src Source data.
     */
    template <size_t BYTE_SIZE>
    constexpr void copy(const uint8_t *src) {
        static_assert(size() == BYTE_SIZE, "Copied length much match the image size");
        static_cast<void>(view().copy(std::span<const uint8_t>{src, BYTE_SIZE}));
    }

    /**
     * \brief Drawing surface over the pixels of this image. Valid while the image is alive and not moved.
     * \return surface of this image.
     */
    [[nodiscard]] constexpr surface<pixel_type> view() {
        return {data.data(), W, H, T<W, H, GRAYSCALE, USE_SPAN>::bytes_per_line, palette_.get()};
    }

    /**
     * \brief Read-only surface over the pixels of this image. Valid while the image is alive and not moved.
     * \return const_surface of this image.
     */
    [[nodiscard]] constexpr const_surface<pixel_type> view() const {
        return {data.data(), W, H, T<W, H, GRAYSCALE, USE_SPAN>::bytes_per_line, palette_.get()};
    }

    /**
     * \brief Convert this instance to an equivalent RGBA8 data array.
     * \param dst RGBA8 array made of image::width() * img::height() * uint32_t values.
     */
    constexpr void RGBA_uint32(std::array<uint32_t, W * H> &dst) const {
        view().RGBA_uint32(dst.data());
    }

    /**
     * \brief Convert this instance to an equivalent RGBA8 data array.
     * \param dst RGBA8 array made of image::width() * img::height() * 4 * uint8_t values.
     */
    constexpr void RGBA_uint8(std::array<uint8_t, W * H * 4> &dst) const {
        view().RGBA_uint8(dst.data());
    }

    /**
     * \brief Return a transposed version of this image.
     * \tparam FLIP_H Flip image horizontally
     * \tparam FLIP_V Flip image vertically
     */
    template <bool FLIP_H = false, bool FLIP_V = false>
    constexpr image<T, H, W, GRAYSCALE> transpose() const {
        image<T, H, W, GRAYSCALE> transposed;
        view().template transpose<FLIP_H, FLIP_V>(transposed.data_ref().data());
        return transposed;
    }

    /**
     * \brief Transpose this image into another.
     * \tparam FLIP_H Flip image horizontally
     * \tparam FLIP_V Flip image vertically
     */
    template <bool FLIP_H = false, bool FLIP_V = false>
    constexpr void transpose(image<T, H, W, GRAYSCALE> &dst) const {
        view().template transpose<FLIP_H, FLIP_V>(dst.data_ref().data());
    }

 private:
#ifndef __INTELLISENSE__
    /// @cond DOXYGEN_EXCLUDE
    /**
     * @private
     */
    using dataType = std::conditional_t<USE_SPAN, std::span<uint8_t, T<W, H, GRAYSCALE, USE_SPAN>::image_size>,
                                        std::array<uint8_t, T<W, H, GRAYSCALE, USE_SPAN>::image_size>>;

    friend surface_api<image<T, W, H, GRAYSCALE, USE_SPAN>, pixel_type>;

    /**
     * @private
     * Line cache for 24/32-bit sixel output, on the stack of the caller.
     */
    [[nodiscard]] static constexpr std::array<uint8_t, pixel_type::sixel_cache_size(W)> sixel_cache() {
        return {};
    }

    dataType data{};

    /**
     * @private
     */
    CONSTIXEL_NO_UNIQUE_ADDRESS hidden::palette_store<pixel_type> palette_{};
    /// @endcond  // DOXYGEN_EXCLUDE
#endif  // #ifndef __INTELLISENSE__
};

/**
 * @class dynamic_image
 * @brief Image with a size chosen at run time and heap storage. Same drawing API as image. Example:
 *
 * \code{.cpp}
 * constixel::dynamic_image<constixel::pixel_8bit<false>> image(640, 480);
 * \endcode
 *
 * @tparam PF Pixel class of the format, see basic_surface.
 */
template <typename PF>
class dynamic_image : public surface_api<dynamic_image<PF>, PF> {
 public:
    /// @brief Pixel class of the format, the PF of surface.
    using pixel_type = PF;

    /**
     * \brief Creates a new image, cleared to zero.
     * \param w Width in pixels.
     * \param h Height in pixels.
     */
    constexpr dynamic_image(size_t w, size_t h) : data_(PF::image_size(w, h), uint8_t{0}), w_(w), h_(h) {
    }

    /**
     * \brief Creates an image that uses a caller-supplied runtime palette.
     * \param w Width in pixels.
     * \param h Height in pixels.
     * \param palette Up to 255 sRGB byte triples.
     */
    dynamic_image(size_t w, size_t h, std::span<const std::array<uint8_t, 3>> palette)
        requires hidden::has_runtime_palette<PF>
        : data_(PF::image_size(w, h), uint8_t{0}), w_(w), h_(h), palette_{runtime_palette(palette), true} {
    }

    /**
     * \brief Bit depth of the image.
     * \return Bit depth of the image. 1 == 2 colors, 2 == 4 colors, 4 == 16 colors, 8 == 256 colors.
     */
    [[nodiscard]] static constexpr size_t bit_depth() {
        return PF::bits_per_pixel;
    }

    /**
     * \brief Size in bytes of the image data.
     * \return Size in bytes of the image data.
     */
    [[nodiscard]] constexpr size_t size() const noexcept {
        return data_.size();
    }

    /**
     * \brief Bytes per line in the image data. Also called stride.
     * \return Bytes per line in the image data.
     */
    [[nodiscard]] constexpr size_t bytes_per_line() const noexcept {
        return PF::bytes_per_line(w_);
    }

    /**
     * \brief Width in pixels of the image.
     * \return Width in pixels of the image.
     */
    [[nodiscard]] constexpr int32_t width() const noexcept {
        return static_cast<int32_t>(w_);
    }

    /**
     * \brief Height in pixels of the image.
     * \return Height in pixels of the image.
     */
    [[nodiscard]] constexpr int32_t height() const noexcept {
        return static_cast<int32_t>(h_);
    }

    /**
     * \brief Clear the image by setting everything to zero.
     */
    constexpr void clear() {
        std::fill(data_.begin(), data_.end(), uint8_t{0});
    }

    /**
     * \brief Get a reference to the underlying raw data of the image.
     * \return Reference to the vector which contains the raw image data.
     */
    [[nodiscard]] constexpr std::vector<uint8_t> &data_ref() {
        return data_;
    }

    /**
     * \brief Returns a clone of this image. Data is copied.
     * \return Cloned image instance.
     */
    [[nodiscard]] constexpr dynamic_image clone() const {
        return *this;
    }

    /**
     * \brief Copy source image into this instance. No compositing occurs, replaces current content.
     * \param src Source image of the same size.
     * \return false if the sizes differ. Nothing is copied then.
     */
    [[nodiscard]] constexpr bool copy(const dynamic_image &src) {
        return view().copy(src.view());
    }

    /**
     * \brief Copy raw source data into this instance. No compositing occurs.
     * \param src Source data. The size must match dynamic_image::size()
     * \return false if the size of src differs. Nothing is copied then.
     */
    [[nodiscard]] constexpr bool copy(std::span<const uint8_t> src) {
        return view().copy(src);
    }

    /**
     * \brief Convert this instance to equivalent RGBA8 data.
     * \param dst RGBA8 buffer of at least width() * height() uint32_t values.
     * \return false if dst is too small. Nothing is written then.
     */
    [[nodiscard]] constexpr bool RGBA_uint32(std::span<uint32_t> dst) const {
        return view().RGBA_uint32(dst);
    }

    /**
     * \brief Convert this instance to equivalent RGBA8 data.
     * \param dst RGBA8 buffer of at least width() * height() * 4 uint8_t values.
     * \return false if dst is too small. Nothing is written then.
     */
    [[nodiscard]] constexpr bool RGBA_uint8(std::span<uint8_t> dst) const {
        return view().RGBA_uint8(dst);
    }

    /**
     * \brief Return a transposed version of this image, width and height swapped.
     * \tparam FLIP_H Flip image horizontally
     * \tparam FLIP_V Flip image vertically
     */
    template <bool FLIP_H = false, bool FLIP_V = false>
    [[nodiscard]] constexpr dynamic_image transpose() const {
        dynamic_image transposed(h_, w_);
        view().template transpose<FLIP_H, FLIP_V>(transposed.data_.data());
        return transposed;
    }

    /**
     * \brief Transpose this image into another.
     * \tparam FLIP_H Flip image horizontally
     * \tparam FLIP_V Flip image vertically
     * \param dst Image of height() x width() pixels.
     * \return false if the size of dst does not match. Nothing is written then.
     */
    template <bool FLIP_H = false, bool FLIP_V = false>
    [[nodiscard]] constexpr bool transpose(dynamic_image &dst) const {
        if (dst.w_ != h_ || dst.h_ != w_) {
            return false;
        }
        view().template transpose<FLIP_H, FLIP_V>(dst.data_.data());
        return true;
    }

    /**
     * \brief Drawing surface over the pixels of this image. Valid while the image is alive and not moved.
     * \return surface of this image.
     */
    [[nodiscard]] constexpr surface<PF> view() {
        return {data_.data(), w_, h_, PF::bytes_per_line(w_), palette_.get()};
    }

    /**
     * \brief Read-only surface over the pixels of this image. Valid while the image is alive and not moved.
     * \return const_surface of this image.
     */
    [[nodiscard]] constexpr const_surface<PF> view() const {
        return {data_.data(), w_, h_, PF::bytes_per_line(w_), palette_.get()};
    }

 private:
    /// @cond DOXYGEN_EXCLUDE
    friend surface_api<dynamic_image<PF>, PF>;

    // Line cache for 24/32-bit sixel output.
    [[nodiscard]] constexpr std::vector<uint8_t> sixel_cache() const {
        return std::vector<uint8_t>(PF::sixel_cache_size(w_));
    }

    std::vector<uint8_t> data_;
    size_t w_;
    size_t h_;
    hidden::palette_store<PF> palette_{};
    /// @endcond
};

}  // namespace constixel

#undef CONSTIXEL_NO_UNIQUE_ADDRESS

#endif  // CONSTIXEL_HPP_
