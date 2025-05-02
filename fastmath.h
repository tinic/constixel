
/*
Copyright 2020 Tinic Uro

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#ifndef FASTMATH_H_
#define FASTMATH_H_

#include <math.h>
#include <memory.h>
#include <stdint.h>

#include <bit>

static constexpr float fast_rcp(const float x) {  // 23.8bits of accuracy
#ifdef GCC_BROKEN_BITCAST
    return 1.0f / x;
#else  // #ifdef GCC_BROKEN_BITCAST
#if (__FPU_PRESENT == 1) && (__FPU_USED == 1)
    return 1.0f / x;
#else   // #if (__FPU_PRESENT == 1) && (__FPU_USED == 1)
    // https://ieeexplore.ieee.org/document/8525803
    int32_t i = std::bit_cast<int>(x);
    i = 0x7eb53567 - i;
    float y = std::bit_cast<float>(i);
    y = 1.9395974f * y * fmaf(y, -x, 1.436142f);
    float r = fmaf(y, -x, 1.0f);
    y = fmaf(y, r, y);
    return y;
#endif  // #if (__FPU_PRESENT == 1) && (__FPU_USED == 1)
#endif  // #ifdef GCC_BROKEN_BITCAST
}

static constexpr float fast_cbrtf(const float x) {  // 22.84bits of accuracy
#ifdef GCC_BROKEN_BITCAST
    return cbrtf(x);
#else   // #ifdef GCC_BROKEN_BITCAST
    // https://www.mdpi.com/1996-1073/14/4/1058=
    float k1 = 1.7523196760f;
    float k2 = 1.2509524245f;
    float k3 = 0.5093818292f;
    uint32_t i = std::bit_cast<uint32_t>(x);
    i = 0x548c2b4b - i / 3;
    float y = std::bit_cast<float>(i);
    float c = x * y * y * y;
    y = y * (k1 - c * (k2 - k3 * c));
    float d = x * y * y;
    c = fmaf(-d, y, 1.0f);
    y = d * fmaf(1.0f / 3.0f, c, 1.0f);
    return y;
#endif  // #ifdef GCC_BROKEN_BITCAST
}

static constexpr float fast_rcbrtf(const float x) {  // 22.84bits of accuracy
#ifdef GCC_BROKEN_BITCAST
    return 1.0f / cbrtf(x);
#else  // #ifdef GCC_BROKEN_BITCAST
#if (__FPU_PRESENT == 1) && (__FPU_USED == 1)
    return 1.0f / cbrtf(x);
#else   // #if (__FPU_PRESENT == 1) && (__FPU_USED == 1)
    // https://www.mdpi.com/1996-1073/14/4/1058=
    float k1 = 1.7523196760f;
    float k2 = 1.2509524245f;
    float k3 = 0.5093818292f;
    uint32_t i = std::bit_cast<uint32_t>(x);
    i = 0x548c2b4b - i / 3;
    float y = std::bit_cast<float>(i);
    float c = x * y * y * y;
    y = y * (k1 - c * (k2 - k3 * c));
    c = 1.0f - x * y * y * y;
    y = y * fmaf(1.0f / 3.0f, c, 1.0f);
    return y;
#endif  // #if (__FPU_PRESENT == 1) && (__FPU_USED == 1)
#endif  // #ifdef GCC_BROKEN_BITCAST
}

static constexpr float fast_sqrtf(const float x) {  // 23.62bits of accuracy
#ifdef GCC_BROKEN_BITCAST
    return __builtin_sqrtf(x);
#else  // #ifdef GCC_BROKEN_BITCAST
#if (__FPU_PRESENT == 1) && (__FPU_USED == 1)
    return __builtin_sqrtf(x);
#else   // #if (__FPU_PRESENT == 1) && (__FPU_USED == 1)
    // https://www.mdpi.com/2079-3197/9/2/21
    int32_t i = std::bit_cast<int>(x);
    int k = i & 0x00800000;
    float y;
    if (k != 0) {
        i = 0x5ed9d098 - (i >> 1);
        y = std::bit_cast<float>(i);
        y = 2.33139729f * y * fmaf(-x, y * y, 1.07492042f);
    } else {
        i = 0x5f19d352 - (i >> 1);
        y = std::bit_cast<float>(i);
        y = 0.82420468f * y * fmaf(-x, y * y, 2.14996147f);
    }
    float c = x * y;
    float r = fmaf(y, -c, 1.0f);
    y = fmaf(0.5f * c, r, c);
    return y;
#endif  // #if (__FPU_PRESENT == 1) && (__FPU_USED == 1)
#endif  // #ifdef GCC_BROKEN_BITCAST
}

static constexpr float fast_rsqrtf(const float x) {  // 23.40bits of accuracy
#ifdef GCC_BROKEN_BITCAST
    return 1.0f / __builtin_sqrtf(x);
#else   // #ifdef GCC_BROKEN_BITCAST
    // https://www.mdpi.com/2079-3197/9/2/21
    int32_t i = std::bit_cast<int>(x);
    int k = i & 0x00800000;
    float y;
    if (k != 0) {
        i = 0x5ed9dbc6 - (i >> 1);
        y = std::bit_cast<float>(i);
        y = 2.33124018f * y * fmaf(-x, y * y, 1.07497406f);
    } else {
        i = 0x5f19d200 - (i >> 1);
        y = std::bit_cast<float>(i);
        y = 0.824212492f * y * fmaf(-x, y * y, 2.14996147f);
    }
    float c = x * y;
    float r = fmaf(y, -c, 1.0f);
    y = fmaf(0.5f * y, r, y);
    return y;
#endif  // #ifdef GCC_BROKEN_BITCAST
}

static constexpr float fast_exp2(const float p) {
#ifdef GCC_BROKEN_BITCAST
    return exp2f(p);
#else   // #ifdef GCC_BROKEN_BITCAST
    const float offset = (p < 0) ? 1.0f : 0.0f;
    const float clipp = (p < -126) ? -126.0f : p;
    const float z = clipp - static_cast<float>(static_cast<int32_t>(clipp)) + offset;
    return std::bit_cast<float>(static_cast<uint32_t>((1 << 23) * (clipp + 121.2740575f + 27.7280233f / (4.84252568f - z) - 1.49012907f * z)));
#endif  // #ifdef GCC_BROKEN_BITCAST
}

static constexpr float fast_log2(const float x) {
#ifdef GCC_BROKEN_BITCAST
    return log2f(x);
#else   // #ifdef GCC_BROKEN_BITCAST
    uint32_t xi = std::bit_cast<uint32_t>(x);
    float xf = std::bit_cast<float>((xi & 0x007FFFFF) | 0x3f000000);
    const float y = static_cast<float>(xi) * 1.1920928955078125e-7f;
    return y - 124.22551499f - 1.498030302f * xf - 1.72587999f / (0.3520887068f + xf);
#endif  // #ifdef GCC_BROKEN_BITCAST
}

static constexpr float fast_log(const float x) { return fast_log2(x) * 0.69314718f; }

static constexpr float fast_pow(const float x, const float p) {
#ifdef GCC_BROKEN_BITCAST
    return powf(x, p);
#else   // GCC_BROKEN_BITCAST
    return fast_exp2(p * fast_log2(x));
#endif  // GCC_BROKEN_BITCAST
}

static constexpr double frac(double v) {  // same as fmod(v, 1.0)
    return v - std::trunc(v);
}

static constexpr float fracf(float v) {  // same as fmodf(v, 1.0f)
    return v - std::truncf(v);
}

static constexpr float constexpr_pow(const float x, const float p) {
    return ::exp2f(p * ::log2f(x));
}

static constexpr float signf(float x) { return (x > 0.0f) ? 1.0f : ((x < 0.0f) ? -1.0f : 1.0f); }

#endif  // #ifndef FASTMATH_H_
