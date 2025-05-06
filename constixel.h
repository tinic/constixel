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
#include <vector>

namespace constixel {

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

static constexpr float fast_pow(const float x, const float p) {
#ifdef GCC_BROKEN_BITCAST
    return powf(x, p);
#else   // GCC_BROKEN_BITCAST
    return fast_exp2(p * fast_log2(x));
#endif  // GCC_BROKEN_BITCAST
}

static constexpr double m_pi_d = 3.14159265358979323846;

static consteval double cos(double x, int32_t terms = 10) {
    x = x - 6.283185307179586 * static_cast<int32_t>(x / 6.283185307179586);  // wrap x to [0, 2π)
    double res = 1.0, term = 1.0;
    double x2 = x * x;
    for (int32_t i = 1; i < terms; ++i) {
        term *= -x2 / ((2 * i - 1) * (2 * i));
        res += term;
    }
    return res;
}

static consteval double sin(double x, int32_t terms = 10) {
    x = x - 6.283185307179586 * static_cast<int32_t>(x / 6.283185307179586);  // wrap x to [0, 2π)
    double res = x, term = x;
    double x2 = x * x;
    for (int32_t i = 1; i < terms; ++i) {
        term *= -x2 / ((2 * i) * (2 * i + 1));
        res += term;
    }
    return res;
}

static consteval double pow(double base, double exp, int32_t terms = 10) {
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
    return {oklch.l, oklch.c * cos(oklch.h * m_pi_d / 180.0), oklch.c * sin(oklch.h * m_pi_d / 180.0)};
}

template <size_t S>
class quantize {
    static constexpr size_t palette_size = S;
    const std::array<uint32_t, palette_size> &pal;

   public:
    explicit constexpr quantize(const std::array<uint32_t, palette_size> &palette) : pal(palette) {
        for (size_t i = 0; i < pal.size(); ++i) {
            linearpal[i * 3 + 0] = srgb_to_linear(static_cast<float>((pal[i] >> 16) & 0xFF) * (1.0f / 255.0f));
            linearpal[i * 3 + 1] = srgb_to_linear(static_cast<float>((pal[i] >> 8) & 0xFF) * (1.0f / 255.0f));
            linearpal[i * 3 + 2] = srgb_to_linear(static_cast<float>((pal[i] >> 0) & 0xFF) * (1.0f / 255.0f));
        }
    }

    constexpr uint8_t nearest(int32_t r, int32_t g, int32_t b) const {
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

    constexpr uint8_t nearest_linear(float r, float g, float b) const {
        size_t best = 0;
        float bestd = 100.0f;
        for (size_t i = 0; i < pal.size(); ++i) {
            float dr = r - linearpal[i * 3 + 0];
            float dg = g - linearpal[i * 3 + 1];
            float db = b - linearpal[i * 3 + 2];
            float d = dr * dr + dg * dg + db * db;
            if (d < bestd) {
                bestd = d;
                best = i;
            }
        }
        return static_cast<uint8_t>(best);
    }
};

template <size_t N, typename T>
class hextree {
    static constexpr T bitslices = ((sizeof(T) * 8) / 4) - 1;
    static constexpr T invalid = std::numeric_limits<T>::max();

    struct node {
        T child[16]{};
        constexpr node() {
            for (auto &c : child) {
                c = invalid;
            }
        }
    };

   public:
    std::array<node, N> nodes;

    [[nodiscard]] size_t byte_size() const {
        return sizeof(node) * nodes.size();
    }

    hextree() = delete;
    hextree(const hextree &) = delete;
    hextree &operator=(const hextree &) = delete;

    template <std::size_t NS>
    explicit consteval hextree(const std::array<std::pair<T, T>, NS> &in) {
        nodes[0] = node{};
        T node_cnt = 1;
        for (auto [key, val] : in) {
            T idx = 0;
            for (T d = 0; d < bitslices; d++) {
                T nib = (key >> ((bitslices - d) * 4)) & 0xF;
                T next = nodes[idx].child[nib];
                if (next == invalid) {
                    next = node_cnt;
                    nodes[idx].child[nib] = next;
                    node_cnt++;
                }
                idx = next;
            }
            nodes[idx].child[key & 0xF] = val;
        }
    }

    template <std::size_t NS>
    [[nodiscard]] static consteval T size(const std::array<std::pair<T, T>, NS> &in) {
        std::vector<node> vnodes{1};
        vnodes.assign(1, node{});
        for (auto [key, val] : in) {
            T idx = 0;
            for (T d = 0; d < bitslices; d++) {
                T nib = (key >> ((bitslices - d) * 4)) & 0xF;
                T next = vnodes[idx].child[nib];
                if (next == invalid) {
                    next = static_cast<T>(vnodes.size());
                    vnodes[idx].child[nib] = next;
                    vnodes.emplace_back();
                }
                idx = next;
            }
            vnodes[idx].child[key & 0xF] = val;
        }
        return static_cast<T>(vnodes.size());
    }

    [[nodiscard]] constexpr T lookup(T key) const {
        T idx = 0;
        for (T d = 0; d < bitslices; ++d) {
            T nib = (key >> ((bitslices - d) * 4)) & 0xF;
            T next = nodes[idx].child[nib];
            if (next == invalid) {
                return 0;
            }
            idx = next;
        }
        if (nodes[idx].child[key & 0xF] == invalid) {
            return 0;
        }
        return nodes[idx].child[key & 0xF];
    }
};

struct char_info {
    int16_t x;
    int16_t y;
    int16_t width;
    int16_t height;
    int16_t xadvance;
    int16_t xoffset;
    int16_t yoffset;
};

template <typename T>
struct rect {
    T x = 0;
    T y = 0;
    T w = 0;
    T h = 0;

    constexpr rect operator&(const rect &other) const {
        T nax = std::max(x, other.x);
        T nay = std::max(y, other.y);
        T nix = std::min(x + w, other.x + other.w);
        T niy = std::min(x + h, other.y + other.h);
        return {nax, nay, nix - nax, niy - nay};
    }

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

class format {
   public:
    static constexpr uint32_t adler32(const uint8_t *data, std::size_t len, uint32_t adler32_sum) {
        uint32_t adler32_s1 = adler32_sum & 0xFFFF;
        uint32_t adler32_s2 = adler32_sum >> 16;
        for (size_t c = 0; c < len; c++) {
            adler32_s1 = (adler32_s1 + data[c]) % 65521;
            adler32_s2 = (adler32_s2 + adler32_s1) % 65521;
        }
        return ((adler32_s2 << 16) | adler32_s1);
    }

    template <typename F>
    static constexpr void png_write_be(F &&charOut, uint32_t value) {
        charOut(static_cast<char>((value >> 24) & 0xFF));
        charOut(static_cast<char>((value >> 16) & 0xFF));
        charOut(static_cast<char>((value >> 8) & 0xFF));
        charOut(static_cast<char>((value >> 0) & 0xFF));
    }

    template <typename F, typename A>
    static constexpr void png_write_crc32(F &&charOut, const A &array, size_t bytes) {
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
        png_write_be(charOut, crc ^ 0xFFFFFFFF);
    }

    template <typename F, typename A>
    static constexpr void png_write_array(F &&charOut, const A &array, size_t bytes) {
        for (size_t c = 0; c < bytes; c++) {
            charOut(static_cast<char>(array[c]));
        }
    }

    template <typename F>
    static constexpr void png_marker(F &&charOut) {
        charOut(0x89);
        charOut(0x50);
        charOut(0x4E);
        charOut(0x47);
        charOut(0x0D);
        charOut(0x0A);
        charOut(0x1A);
        charOut(0x0A);
    }

    template <typename F>
    static constexpr void png_header(F &&charOut, size_t w, size_t h, size_t depth) {
        const size_t chunkLength = 17;
        std::array<char, chunkLength> header;
        uint32_t i = 0;
        header[i++] = 'I';
        header[i++] = 'H';
        header[i++] = 'D';
        header[i++] = 'R';
        header[i++] = static_cast<char>((w >> 24) & 0xFF);
        header[i++] = static_cast<char>((w >> 16) & 0xFF);
        header[i++] = static_cast<char>((w >> 8) & 0xFF);
        header[i++] = static_cast<char>((w >> 0) & 0xFF);
        header[i++] = static_cast<char>((h >> 24) & 0xFF);
        header[i++] = static_cast<char>((h >> 16) & 0xFF);
        header[i++] = static_cast<char>((h >> 8) & 0xFF);
        header[i++] = static_cast<char>((h >> 0) & 0xFF);
        header[i++] = static_cast<char>(depth);
        header[i++] = 3;
        header[i++] = 0;
        header[i++] = 0;
        header[i++] = 0;
        png_write_be(charOut, i - 4);
        png_write_array(charOut, header, i);
        png_write_crc32(charOut, header, i);
    }

    template <typename F, typename P>
    static constexpr void png_palette(F &&charOut, const P &palette) {
        std::array<char, 256 * 3 + 4> header;
        uint32_t i = 0;
        header[i++] = 'P';
        header[i++] = 'L';
        header[i++] = 'T';
        header[i++] = 'E';
        for (size_t c = 0; c < palette.size(); c++) {
            header[i++] = static_cast<char>((palette[c] >> 16) & 0xFF);
            header[i++] = static_cast<char>((palette[c] >> 8) & 0xFF);
            header[i++] = static_cast<char>((palette[c] >> 0) & 0xFF);
        }
        png_write_be(charOut, i - 4);
        png_write_array(charOut, header, i);
        png_write_crc32(charOut, header, i);
    }

    template <typename F>
    static constexpr void png_end(F &&charOut) {
        std::array<char, 4> header;
        uint32_t i = 0;
        header[i++] = 'I';
        header[i++] = 'E';
        header[i++] = 'N';
        header[i++] = 'D';
        png_write_be(charOut, i - 4);
        png_write_array(charOut, header, i);
        png_write_crc32(charOut, header, i);
    }

    template <typename F>
    static constexpr void png_idat_zlib_header(F &&charOut) {
        std::array<char, 6> header;
        uint32_t i = 0;
        header[i++] = 'I';
        header[i++] = 'D';
        header[i++] = 'A';
        header[i++] = 'T';
        header[i++] = 0x78;
        header[i++] = 0x01;
        png_write_be(charOut, i - 4);
        png_write_array(charOut, header, i);
        png_write_crc32(charOut, header, i);
    }

    template <typename F>
    static constexpr void png_idat_zlib_stream(F &&charOut, const uint8_t *line, size_t bytes, uint32_t &adler32_sum) {
        while (bytes > 0) {
            static constexpr size_t max_data_use = 1024;
            static constexpr size_t extra_data = 24;
            static constexpr size_t max_stack_use = max_data_use + extra_data;
            std::array<uint8_t, max_stack_use> header;

            uint32_t i = 0;
            header[i++] = 'I';
            header[i++] = 'D';
            header[i++] = 'A';
            header[i++] = 'T';
            header[i++] = 0x00;

            size_t bytes_to_copy = std::min(static_cast<size_t>(max_data_use), bytes);
            header[i++] = (((bytes_to_copy + 1) >> 0) & 0xFF);
            header[i++] = (((bytes_to_copy + 1) >> 8) & 0xFF);
            header[i++] = ((((bytes_to_copy + 1) ^ 0xffff) >> 0) & 0xFF);
            header[i++] = ((((bytes_to_copy + 1) ^ 0xffff) >> 8) & 0xFF);

            uint32_t adlersum32_start_pos = i;
            header[i++] = 0;
            for (size_t c = 0; c < bytes_to_copy; c++) {
                header[i++] = line[c];
            }
            adler32_sum = adler32(&header[adlersum32_start_pos], i - adlersum32_start_pos, adler32_sum);

            png_write_be(charOut, i - 4);
            png_write_array(charOut, header, i);
            png_write_crc32(charOut, header, i);

            bytes -= bytes_to_copy;
        }
    }

    template <typename F>
    static constexpr void png_idat_zlib_trailer(F &&charOut, uint32_t adler32_sum) {
        std::array<char, 8> header;
        uint32_t i = 0;
        header[i++] = 'I';
        header[i++] = 'D';
        header[i++] = 'A';
        header[i++] = 'T';
        header[i++] = static_cast<char>((adler32_sum >> 24) & 0xFF);
        header[i++] = static_cast<char>((adler32_sum >> 16) & 0xFF);
        header[i++] = static_cast<char>((adler32_sum >> 8) & 0xFF);
        header[i++] = static_cast<char>((adler32_sum >> 0) & 0xFF);
        png_write_be(charOut, i - 4);
        png_write_array(charOut, header, i);
        png_write_crc32(charOut, header, i);
    }

    template <typename F>
    static constexpr void sixel_header(F &&charOut) {
        charOut(0x1b);
        charOut('P');
        charOut('0');
        charOut(';');
        charOut('1');
        charOut(';');
        charOut('0');
        charOut('q');
    }

    template <size_t W, size_t H, size_t S, typename F>
    static constexpr void sixel_raster_attributes(F &&charOut) {
        charOut('\"');
        sixel_number(charOut, 2);
        charOut(';');
        sixel_number(charOut, 2);
        charOut(';');
        sixel_number(charOut, W * S);
        charOut(';');
        sixel_number(charOut, H * S);
    }

    template <typename F>
    static constexpr void sixel_number(F &&charOut, uint16_t u) {
        if (u < 10) {
            charOut(static_cast<char>('0' + u));
        } else if (u < 100) {
            charOut(static_cast<char>('0' + (((u / 10) % 10))));
            charOut(static_cast<char>('0' + (u % 10)));
        } else if (u < 1000) {
            charOut(static_cast<char>('0' + ((u / 100) % 10)));
            charOut(static_cast<char>('0' + ((u / 10) % 10)));
            charOut(static_cast<char>('0' + (u % 10)));
        } else if (u < 10000) {
            charOut(static_cast<char>('0' + ((u / 1000) % 10)));
            charOut(static_cast<char>('0' + ((u / 100) % 10)));
            charOut(static_cast<char>('0' + ((u / 10) % 10)));
            charOut(static_cast<char>('0' + (u % 10)));
        } else {
            charOut(static_cast<char>('0' + ((u / 10000) % 10)));
            charOut(static_cast<char>('0' + ((u / 1000) % 10)));
            charOut(static_cast<char>('0' + ((u / 100) % 10)));
            charOut(static_cast<char>('0' + ((u / 10) % 10)));
            charOut(static_cast<char>('0' + (u % 10)));
        }
    }

    template <typename F>
    static constexpr void sixel_color(F &&charOut, uint16_t i, uint32_t col) {
        charOut('#');
        sixel_number(charOut, i);
        charOut(';');
        charOut('2');
        charOut(';');
        for (size_t c = 0; c < 3; c++) {
            sixel_number(charOut, static_cast<uint16_t>((((col >> (8 * (2 - c))) & 0xFF) * 100) / 255));
            charOut(';');
        }
    }

    template <typename F>
    static constexpr void sixel_end(F &&charOut) {
        charOut(0x1b);
        charOut('\\');
    }

    template <typename PBT, size_t PBS>
    struct palette_bitset {
        constexpr void mark(PBT col) {
            PBT idx = col >> 5;
            set[idx] |= 1UL << (col & 0x1F);
        }

        constexpr size_t genstack(std::array<PBT, PBS> &stack) const {
            size_t count = 0;
            for (size_t c = 0; c < PBS / 32; c++) {
                for (size_t d = 0; d < 32; d++) {
                    if ((set[c] >> d) & 1) {
                        stack[count++] = static_cast<PBT>(c * 32 + d);
                    }
                }
            }
            return count;
        }

        uint32_t set[PBS / 32]{};
    };

    template <size_t W, size_t H, int32_t S, typename PBT, size_t PBS, typename P, typename F, typename L>
    static constexpr void png_image(const uint8_t *data, const P &palette, F &&charOut, const L &linePtr) {
        png_marker(charOut);
        png_header(charOut, W, H, PBS);
        png_palette(charOut, palette);
        png_idat_zlib_header(charOut);
        uint32_t adler32_sum = 1;
        for (size_t y = 0; y < H; y++) {
            size_t bpl = 0;
            const uint8_t *ptr = linePtr(data, y, bpl);
            png_idat_zlib_stream(charOut, ptr, bpl, adler32_sum);
        }
        png_idat_zlib_trailer(charOut, adler32_sum);
        png_end(charOut);
    }

    template <size_t W, size_t H, int32_t S, typename PBT, size_t PBS, typename P, typename F, typename C, typename D>
    static constexpr void sixel_image(const uint8_t *data, const P &palette, F &&charOut, const rect<int32_t> &_r, const C &collect6, const D &set6) {
        sixel_header(charOut);
        sixel_raster_attributes<W, H, S>(charOut);
        for (size_t c = 0; c < palette.size(); c++) {
            sixel_color(charOut, static_cast<uint16_t>(c), palette.data()[c]);
        }
        const auto r = rect<int32_t>{_r.x * S, _r.y * S, _r.w * S, _r.h * S} & rect<int32_t>{0, 0, W * S, H * S};
        std::array<PBT, std::max(32UL, 1UL << PBS)> stack{};
        for (size_t y = static_cast<size_t>(r.y); y < static_cast<size_t>(r.y + r.h); y += 6) {
            palette_bitset<PBT, std::max(32UL, 1UL << PBS)> pset;
            set6(data, static_cast<size_t>(r.x), static_cast<size_t>(r.w), y, pset);
            size_t stackCount = pset.genstack(stack);
            for (size_t s = 0; s < stackCount; s++) {
                PBT col = stack[s];
                if (col != 0) {
                    charOut('$');
                }
                charOut('#');
                sixel_number(charOut, static_cast<uint16_t>(col));
                for (size_t x = static_cast<size_t>(r.x); x < static_cast<size_t>(r.x + r.w); x++) {
                    PBT bits6 = collect6(data, x, col, y);
                    uint16_t repeatCount = 0;
                    for (size_t xr = (x + 1); xr < (std::min(x + 255, W * S)); xr++) {
                        if (bits6 == collect6(data, xr, col, y)) {
                            repeatCount++;
                            continue;
                        }
                        break;
                    }
                    if (repeatCount > 3) {
                        charOut('!');
                        sixel_number(charOut, repeatCount + 1);
                        x += repeatCount;
                    }
                    charOut(static_cast<char>('?' + bits6));
                }
            }
            charOut('-');
        }
        sixel_end(charOut);
    }
};

template <size_t W, size_t H, int32_t S>
class format_1bit : public format {
   public:
    static constexpr size_t bits_per_pixel = 1;
    static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
    static constexpr size_t internal_height = ((H + 5) / 6) * 6;
    static constexpr size_t image_size = internal_height * bytes_per_line;
    static constexpr std::array<uint32_t, (1UL << bits_per_pixel)> palette = {0x00000000, 0x00ffffff};

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

    static constexpr uint8_t get(const uint8_t *line, size_t x) {
        size_t x8 = x / 8;
        size_t xb = x % 8;
        return (line[x8] >> (7 - xb)) & 1;
    }

    static constexpr auto RGBA_uint32(const std::array<uint8_t, image_size> &data) {
        std::array<uint32_t, W * H> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = get(ptr, x);
                rgba[y * W + x] = palette[col] | (col ? 0xFF000000 : 0x00000000);
            }
            ptr += bytes_per_line;
        }
        return rgba;
    }

    static constexpr auto RGBA_uint8(const std::array<uint8_t, image_size> &data) {
        std::array<uint8_t, W * H * 4> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = get(ptr, x);
                rgba[y * W * 4 + x * 4 + 0] = static_cast<uint8_t>((palette[col] >> 16) & 0xFF);
                rgba[y * W * 4 + x * 4 + 1] = static_cast<uint8_t>((palette[col] >> 8) & 0xFF);
                rgba[y * W * 4 + x * 4 + 2] = static_cast<uint8_t>((palette[col] >> 0) & 0xFF);
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
                err = std::clamp(V - (n ? 0xFF * 6 : 0x00), -0xFF * 6, 0xFF * 6);
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
    static constexpr void png(const std::array<uint8_t, image_size> &data, F &&charOut) {
        png_image<W, H, S, uint8_t, bits_per_pixel>(data.data(), palette, charOut, [](const uint8_t *dataRaw, size_t y, size_t &bpl) {
            bpl = bytes_per_line;
            return dataRaw + y * bytes_per_line;
        });
    }

    template <typename F>
    static constexpr void sixel(const std::array<uint8_t, image_size> &data, F &&charOut, const rect<int32_t> &r) {
        sixel_image<W, H, S, uint8_t, bits_per_pixel>(
            data.data(), palette, charOut, r,
            [](const uint8_t *dataRaw, size_t x, size_t col, size_t y) {
                const uint8_t *ptr = &dataRaw[(y / S) * bytes_per_line + (x / S) / 8];
                size_t x8 = (x / S) % 8;
                uint8_t out = 0;
                int32_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    out >>= 1;
                    if ((y + y6) < H * S) {
                        out |= ((((*ptr) >> (7 - x8)) & 1) == col) ? (1UL << 5) : 0;
                    }
                    if (y6 != 5) {
                        if (++inc >= S) {
                            inc = 0;
                            ptr += bytes_per_line;
                        }
                    }
                }
                return out;
            },
            [](const uint8_t *dataRaw, size_t x, size_t w, size_t y, palette_bitset<uint8_t, 32> &set) {
                int32_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    if ((y + y6) < H * S) {
                        for (size_t xx = 0; xx < (w + S - 1) / S; xx++) {
                            const uint8_t *ptr = &dataRaw[(y / S) * bytes_per_line + (xx + x) / 8];
                            size_t x8 = (xx + x) % 8;
                            set.mark((((*ptr) >> (7 - x8)) & 1));
                        }
                    }
                    if (y6 != 5) {
                        if (++inc >= S) {
                            inc = 0;
                        }
                    }
                }
            });
    }
};

template <size_t W, size_t H, int32_t S>
class format_2bit : public format {
   public:
    static constexpr size_t bits_per_pixel = 2;
    static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
    static constexpr size_t internal_height = ((H + 5) / 6) * 6;
    static constexpr size_t image_size = internal_height * bytes_per_line;
    static constexpr std::array<uint32_t, (1UL << bits_per_pixel)> palette = {0x000000, 0xffffff, 0xff0000, 0x0077ff};

    static constexpr const constixel::quantize<1UL << bits_per_pixel> gen_quant() {
        return constixel::quantize<1UL << bits_per_pixel>(palette);
    }

    static constexpr const constixel::quantize<1UL << bits_per_pixel> quant = gen_quant();

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

    static constexpr uint8_t get(const uint8_t *line, size_t x) {
        size_t x4 = x / 4;
        size_t xb = x % 4;
        return (line[x4] >> ((3 - xb) * 2)) & 0x3;
    }

    static constexpr auto RGBA_uint32(const std::array<uint8_t, image_size> &data) {
        std::array<uint32_t, W * H> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = get(ptr, x);
                rgba[y * W + x] = palette[col] | (col ? 0xFF000000 : 0x00000000);
            }
            ptr += bytes_per_line;
        }
        return rgba;
    }

    static constexpr auto RGBA_uint8(const std::array<uint8_t, image_size> &data) {
        std::array<uint8_t, W * H * 4> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = get(ptr, x);
                rgba[y * W * 4 + x * 4 + 0] = static_cast<uint8_t>((palette[col] >> 16) & 0xFF);
                rgba[y * W * 4 + x * 4 + 1] = static_cast<uint8_t>((palette[col] >> 8) & 0xFF);
                rgba[y * W * 4 + x * 4 + 2] = static_cast<uint8_t>((palette[col] >> 0) & 0xFF);
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
                err_r = std::clamp(R - static_cast<int32_t>((palette[n] >> 16) & 0xFF), -255, 255);
                err_g = std::clamp(G - static_cast<int32_t>((palette[n] >> 8) & 0xFF), -255, 255);
                err_b = std::clamp(B - static_cast<int32_t>((palette[n] >> 0) & 0xFF), -255, 255);
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
                err_r = std::clamp(Rl - quant.linearpal[n * 3 + 0], -1.0f, 1.0f);
                err_g = std::clamp(Gl - quant.linearpal[n * 3 + 1], -1.0f, 1.0f);
                err_b = std::clamp(Bl - quant.linearpal[n * 3 + 2], -1.0f, 1.0f);
            }
        }
    }

    template <typename F>
    static constexpr void png(const std::array<uint8_t, image_size> &data, F &&charOut) {
        png_image<W, H, S, uint8_t, bits_per_pixel>(data.data(), palette, charOut, [](const uint8_t *dataRaw, size_t y, size_t &bpl) {
            bpl = bytes_per_line;
            return dataRaw + y * bytes_per_line;
        });
    }

    template <typename F>
    static constexpr void sixel(const std::array<uint8_t, image_size> &data, F &&charOut, const rect<int32_t> &r) {
        sixel_image<W, H, S, uint8_t, bits_per_pixel>(
            data.data(), palette, charOut, r,
            [](const uint8_t *dataRaw, size_t x, size_t col, size_t y) {
                const uint8_t *ptr = &dataRaw[(y / S) * bytes_per_line + (x / S) / 4];
                size_t x4 = (x / S) % 4;
                uint8_t out = 0;
                int32_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    out >>= 1;
                    if ((y + y6) < H * S) {
                        out |= ((((*ptr) >> (6 - x4 * 2)) & 3) == col) ? (1UL << 5) : 0;
                    }
                    if (y6 != 5) {
                        if (++inc >= S) {
                            inc = 0;
                            ptr += bytes_per_line;
                        }
                    }
                }
                return out;
            },
            [](const uint8_t *dataRaw, size_t x, size_t w, size_t y, palette_bitset<uint8_t, 32> &set) {
                int32_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    if ((y + y6) < H * S) {
                        for (size_t xx = 0; xx < (w + S - 1) / S; xx++) {
                            const uint8_t *ptr = &dataRaw[(y / S) * bytes_per_line + (xx + x) / 4];
                            size_t x4 = (xx + x) % 4;
                            set.mark((((*ptr) >> (6 - x4 * 2)) & 3));
                        }
                    }
                    if (y6 != 5) {
                        if (++inc >= S) {
                            inc = 0;
                        }
                    }
                }
            });
    }
};

template <size_t W, size_t H, int32_t S>
class format_4bit : public format {
   public:
    static constexpr size_t bits_per_pixel = 4;
    static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
    static constexpr size_t internal_height = ((H + 5) / 6) * 6;
    static constexpr size_t image_size = internal_height * bytes_per_line;
    static constexpr std::array<uint32_t, (1UL << bits_per_pixel)> palette = {0x000000, 0xffffff, 0xff0000, 0x00ff00, 0x0000ff, 0xffff00, 0x00ffff, 0xff00ff,
                                                                              0x333333, 0x666666, 0x999999, 0xcccccc, 0x7f0000, 0x007f00, 0x00007f, 0x7f7f00};

    static constexpr const constixel::quantize<1UL << bits_per_pixel> gen_quant() {
        return constixel::quantize<1UL << bits_per_pixel>(palette);
    }

    static constexpr const constixel::quantize<1UL << bits_per_pixel> quant = gen_quant();

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

    static constexpr uint8_t get(const uint8_t *line, size_t x) {
        size_t x2 = x / 2;
        size_t xb = x % 2;
        return (line[x2] >> ((1 - xb) * 4)) & 0xF;
    }

    static constexpr auto RGBA_uint32(const std::array<uint8_t, image_size> &data) {
        std::array<uint32_t, W * H> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = get(ptr, x);
                rgba[y * W + x] = palette[col] | (col ? 0xFF000000 : 0x00000000);
            }
            ptr += bytes_per_line;
        }
        return rgba;
    }

    static constexpr auto RGBA_uint8(const std::array<uint8_t, image_size> &data) {
        std::array<uint8_t, W * H * 4> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = get(ptr, x);
                rgba[y * W * 4 + x * 4 + 0] = static_cast<uint8_t>((palette[col] >> 16) & 0xFF);
                rgba[y * W * 4 + x * 4 + 1] = static_cast<uint8_t>((palette[col] >> 8) & 0xFF);
                rgba[y * W * 4 + x * 4 + 2] = static_cast<uint8_t>((palette[col] >> 0) & 0xFF);
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
                err_r = std::clamp(R - static_cast<int32_t>((palette[n] >> 16) & 0xFF), -255, 255);
                err_g = std::clamp(G - static_cast<int32_t>((palette[n] >> 8) & 0xFF), -255, 255);
                err_b = std::clamp(B - static_cast<int32_t>((palette[n] >> 0) & 0xFF), -255, 255);
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
                err_r = std::clamp(Rl - quant.linearpal[n * 3 + 0], -1.0f, 1.0f);
                err_g = std::clamp(Gl - quant.linearpal[n * 3 + 1], -1.0f, 1.0f);
                err_b = std::clamp(Bl - quant.linearpal[n * 3 + 2], -1.0f, 1.0f);
            }
        }
    }

    template <typename F>
    static constexpr void png(const std::array<uint8_t, image_size> &data, F &&charOut) {
        png_image<W, H, S, uint8_t, bits_per_pixel>(data.data(), palette, charOut, [](const uint8_t *dataRaw, size_t y, size_t &bpl) {
            bpl = bytes_per_line;
            return dataRaw + y * bytes_per_line;
        });
    }

    template <typename F>
    static constexpr void sixel(const std::array<uint8_t, image_size> &data, F &&charOut, const rect<int32_t> &r) {
        sixel_image<W, H, S, uint8_t, bits_per_pixel>(
            data.data(), palette, charOut, r,
            [](const uint8_t *dataRaw, size_t x, size_t col, size_t y) {
                const uint8_t *ptr = &dataRaw[(y / S) * bytes_per_line + (x / S) / 2];
                size_t x2 = (x / S) % 2;
                uint8_t out = 0;
                int32_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    out >>= 1;
                    if ((y + y6) < H * S) {
                        out |= ((((*ptr) >> (4 - x2 * 4)) & 0xF) == col) ? (1UL << 5) : 0;
                    }
                    if (y6 != 5) {
                        if (++inc >= S) {
                            inc = 0;
                            ptr += bytes_per_line;
                        }
                    }
                }
                return out;
            },
            [](const uint8_t *dataRaw, size_t x, size_t w, size_t y, palette_bitset<uint8_t, 32> &set) {
                int32_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    if ((y + y6) < H * S) {
                        for (size_t xx = 0; xx < (w + S - 1) / S; xx++) {
                            const uint8_t *ptr = &dataRaw[(y / S) * bytes_per_line + (xx + x) / 2];
                            size_t x2 = (xx + x) % 2;
                            set.mark((((*ptr) >> (4 - x2 * 4)) & 0xF));
                        }
                    }
                    if (y6 != 5) {
                        if (++inc >= S) {
                            inc = 0;
                        }
                    }
                }
            });
    }
};

template <size_t W, size_t H, int32_t S>
class format_8bit : public format {
   public:
    static constexpr size_t bits_per_pixel = 8;
    static constexpr size_t bytes_per_line = W;
    static constexpr size_t internal_height = ((H + 5) / 6) * 6;
    static constexpr size_t image_size = internal_height * bytes_per_line;

    static consteval const std::array<uint32_t, (1UL << bits_per_pixel)> gen_palette() {
        std::array<uint32_t, (1UL << bits_per_pixel)> pal{};
        // clang-format off
        pal[ 0] = 0x000000; pal[ 1] = 0xffffff; pal[ 2] = 0xff0000; pal[ 3] = 0x00ff00;
        pal[ 4] = 0x0000ff; pal[ 5] = 0xffff00; pal[ 6] = 0x00ffff; pal[ 7] = 0xff00ff;
        pal[ 8] = 0x333333; pal[ 9] = 0x666666; pal[10] = 0x999999; pal[11] = 0xcccccc;
        pal[12] = 0x7f0000; pal[13] = 0x007f00; pal[14] = 0x00007f; pal[15] = 0x7f7f00;
        for (size_t c = 0; c < 16; c++) {
            uint32_t y = (0xff * static_cast<uint32_t>(c)) / 15;
            pal[0x10 + c] = (y << 16) | (y << 8) | (y << 0);
        }
        for (size_t c = 0; c < 8; c++) {
            uint32_t y = (0xff * static_cast<uint32_t>(c)) / 7;
            uint32_t x = (0xff * (static_cast<uint32_t>(c) + 1)) / 8;
            pal[0x20 + c + 0] = ( y  << 16) | ( 0  << 8) | ( 0  << 0); // Rxx
            pal[0x20 + c + 8] = (255 << 16) | ( x  << 8) | ( x  << 0); 
            pal[0x30 + c + 0] = ( 0  << 16) | ( y  << 8) | ( 0  << 0); // Gxx
            pal[0x30 + c + 8] = ( x  << 16) | (255 << 8) | ( x  << 0);
            pal[0x40 + c + 0] = ( 0  << 16) | ( 0  << 8) | ( y  << 0); // Bxx
            pal[0x40 + c + 8] = ( x  << 16) | ( x  << 8) | (255 << 0);
            pal[0x50 + c + 0] = ( y  << 16) | ( y  << 8) | ( 0  << 0); // RGx
            pal[0x50 + c + 8] = (255 << 16) | (255 << 8) | ( x  << 0);
            pal[0x60 + c + 0] = ( 0  << 16) | ( y  << 8) | ( y  << 0); // xGB
            pal[0x60 + c + 8] = ( x  << 16) | (255 << 8) | (255 << 0);
            pal[0x70 + c + 0] = ( y  << 16) | ( 0  << 8) | ( y  << 0); // RxB
            pal[0x70 + c + 8] = (255 << 16) | ( x  << 8) | (255 << 0);
        }
        // clang-format on
        for (size_t c = 0; c < 8; c++) {
            constixel::oklab lft{static_cast<double>(c) / 7 - 0.2, 0.2, 0.0};
            constixel::oklab rgh{static_cast<double>(c) / 7 - 0.2, 0.2, 360.0};
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
        return pal;
    }

    static constexpr std::array<uint32_t, 1UL << bits_per_pixel> palette = gen_palette();

    static constexpr const constixel::quantize<1UL << bits_per_pixel> gen_quant() {
        return constixel::quantize<1UL << bits_per_pixel>(palette);
    }

    static constexpr const constixel::quantize<1UL << bits_per_pixel> quant = gen_quant();

    static constexpr void plot(std::array<uint8_t, image_size> &data, size_t x0, size_t y, uint8_t col) {
        data.data()[y * bytes_per_line + x0] = static_cast<uint8_t>(col);
    }

    static constexpr void span(std::array<uint8_t, image_size> &data, size_t xl0, size_t xr0, size_t y, uint8_t col) {
        uint8_t *yptr = &data.data()[y * bytes_per_line];
        for (size_t x = xl0; x < xr0; x++) {
            yptr[x] = static_cast<uint8_t>(col);
        }
    }

    static constexpr auto RGBA_uint32(const std::array<uint8_t, image_size> &data) {
        std::array<uint32_t, W * H> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = ptr[x];
                rgba[y * W + x] = palette[col] | (col ? 0xFF000000 : 0x00000000);
            }
            ptr += bytes_per_line;
        }
        return rgba;
    }

    static constexpr auto RGBA_uint8(const std::array<uint8_t, image_size> &data) {
        std::array<uint8_t, W * H * 4> rgba{};
        const uint8_t *ptr = data.data();
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                uint8_t col = ptr[x];
                rgba[y * W * 4 + x * 4 + 0] = static_cast<uint8_t>((palette[col] >> 16) & 0xFF);
                rgba[y * W * 4 + x * 4 + 1] = static_cast<uint8_t>((palette[col] >> 8) & 0xFF);
                rgba[y * W * 4 + x * 4 + 2] = static_cast<uint8_t>((palette[col] >> 0) & 0xFF);
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
                err_r = std::clamp(R - static_cast<int32_t>((palette[n] >> 16) & 0xFF), -255, 255);
                err_g = std::clamp(G - static_cast<int32_t>((palette[n] >> 8) & 0xFF), -255, 255);
                err_b = std::clamp(B - static_cast<int32_t>((palette[n] >> 0) & 0xFF), -255, 255);
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
                err_r = std::clamp(Rl - quant.linearpal[n * 3 + 0], -1.0f, 1.0f);
                err_g = std::clamp(Gl - quant.linearpal[n * 3 + 1], -1.0f, 1.0f);
                err_b = std::clamp(Bl - quant.linearpal[n * 3 + 2], -1.0f, 1.0f);
            }
        }
    }

    template <typename F>
    static constexpr void png(const std::array<uint8_t, image_size> &data, F &&charOut) {
        png_image<W, H, S, uint8_t, bits_per_pixel>(data.data(), palette, charOut, [](const uint8_t *dataRaw, size_t y, size_t &bpl) {
            bpl = bytes_per_line;
            return dataRaw + y * bytes_per_line;
        });
    }

    template <typename F>
    static constexpr void sixel(const std::array<uint8_t, image_size> &data, F &&charOut, const rect<int32_t> &r) {
        sixel_image<W, H, S, uint8_t, bits_per_pixel>(
            data.data(), palette, charOut, r,
            [](const uint8_t *dataRaw, size_t x, size_t col, size_t y) {
                const uint8_t *ptr = &dataRaw[(y / S) * bytes_per_line + x / S];
                uint8_t out = 0;
                int32_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    out >>= 1;
                    if ((y + y6) < H * S) {
                        out |= (*ptr == col) ? (1UL << 5) : 0;
                    }
                    if (y6 != 5) {
                        if (++inc >= S) {
                            inc = 0;
                            ptr += bytes_per_line;
                        }
                    }
                }
                return out;
            },
            [](const uint8_t *dataRaw, size_t x, size_t w, size_t y, palette_bitset<uint8_t, 1UL << bits_per_pixel> &set) {
                const uint8_t *ptr = &dataRaw[(y / S) * bytes_per_line + x / S];
                int32_t inc = y % S;
                for (size_t y6 = 0; y6 < 6; y6++) {
                    if ((y + y6) < H * S) {
                        for (size_t xx = 0; xx < (w + S - 1) / S; xx++) {
                            set.mark(ptr[xx + x]);
                        }
                    }
                    if (y6 != 5) {
                        if (++inc >= S) {
                            inc = 0;
                            ptr += bytes_per_line;
                        }
                    }
                }
            });
    }
};

enum color : uint8_t {
    BLACK_TRANSPARENT = 0,
    BLACK_OPAQUE = 16,
    WHITE = 1,
    RED = 2,
    GREEN = 3,
    BLUE = 4,
    YELLOW = 5,
    CYAN = 6,
    MAGENTA = 7,
    GREY_80 = 8,
    GREY_60 = 9,
    GREY_40 = 10,
    GREY_20 = 11,
    DARK_RED = 12,
    DARK_GREEN = 13,
    DARK_BLUE = 14,
    DARK_YELLOW = 15,

    GREY_RAMP_START = 16,
    GREY_RAMP_COUNT = 16,
    GREY_RAMP_STOP = GREY_RAMP_START + 15,

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

template <template <size_t, size_t, int32_t> class T, size_t W, size_t H, int32_t S = 1>
class image {
    static_assert(sizeof(W) >= sizeof(uint32_t));
    static_assert(sizeof(H) >= sizeof(uint32_t));

    static_assert(W <= 16384 && H <= 16384);
    static_assert(S >= 1 && S <= 256);

   public:
    [[nodiscard]] static constexpr int32_t bit_depth() {
        return T<W, H, S>::bits_per_pixel;
    }

    [[nodiscard]] static constexpr int32_t size() {
        return T<W, H, S>::image_size;
    }

    [[nodiscard]] static constexpr size_t width() {
        return W;
    }

    [[nodiscard]] static constexpr size_t height() {
        return H;
    }

    constexpr void clear() {
        data.fill(0);
    }

    [[nodiscard]] static constexpr int32_t abs(int32_t v) {
        return v < 0 ? -v : v;
    }

    [[nodiscard]] constexpr std::array<uint8_t, T<W, H, S>::image_size> &data_ref() const {
        return data;
    }

    [[nodiscard]] constexpr std::array<uint8_t, T<W, H, S>::image_size> clone() const {
        return data;
    }

    constexpr void copy(const std::array<uint8_t, T<W, H, S>::image_size> &src) {
        data = src;
    }

    constexpr void line(int32_t x0, int32_t y0, int32_t x1, int32_t y1, uint8_t col, uint32_t width = 1) {
        int32_t steep = abs(y1 - y0) > abs(x1 - x0);

        if (steep) {
            x0 ^= y0;
            y0 ^= x0;
            x0 ^= y0;
            x1 ^= y1;
            y1 ^= x1;
            x1 ^= y1;
        }

        if (x0 > x1) {
            x0 ^= x1;
            x1 ^= x0;
            x0 ^= x1;
            y0 ^= y1;
            y1 ^= y0;
            y0 ^= y1;
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

        if (width == 1) {
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
        } else if (width > 1) {
            for (; x0 <= x1; x0++) {
                if (steep) {
                    fill_circle(y0, x0, (width + 1) / 2, col);
                } else {
                    fill_circle(x0, y0, (width + 1) / 2, col);
                }
                err -= dy;
                if (err < 0) {
                    y0 += ystep;
                    err += dx;
                }
            }
        }
    }

    constexpr void line(const rect<int32_t> &l, uint8_t col, uint32_t width = 1) {
        line(l.x, l.y, l.x + l.w, l.y + l.h, col, width);
    }

    constexpr void plot(int32_t x, int32_t y, uint8_t col) {
        if (x < 0 || x >= static_cast<int32_t>(W) || y < 0 || y >= static_cast<int32_t>(H)) {
            return;
        }
        T<W, H, S>::plot(data, static_cast<uint32_t>(x), static_cast<uint32_t>(y), col);
    }

    template <typename FONT>
    constexpr int32_t string_width(const char *str) {
        int32_t x = 0;
        while (*str != 0) {
            uint32_t utf32 = 0;
            uint32_t lead = static_cast<uint32_t>(*str);
            if (lead < 0x80) {
                utf32 = lead;
                str += 1;
            } else if ((lead >> 5) == 0x06) {
                utf32 = ((lead & 0x1F) << 6) | (static_cast<uint32_t>(str[1]) & 0x3F);
                str += 2;
            } else if ((lead >> 4) == 0x0E) {
                utf32 = ((lead & 0x0F) << 12) | ((static_cast<uint32_t>(str[1]) & 0x3F) << 6) | (static_cast<uint32_t>(str[3]) & 0x3F);
                str += 3;
            } else if ((lead >> 3) == 0x1E) {
                utf32 = ((lead & 0x07) << 18) | ((static_cast<uint32_t>(str[1]) & 0x3F) << 12) | ((static_cast<uint32_t>(str[2]) & 0x3F) << 6) |
                        (static_cast<uint32_t>(str[3]) & 0x3F);
                str += 4;
            } else {
                return x;
            }
#if __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wimplicit-int-conversion"
#endif  // #if __clang__
            const char_info &ch_info = FONT::char_table.at(FONT::glyph_tree.lookup(utf32));
#if __clang__
#pragma GCC diagnostic pop
#endif  // #if __clang__
            if (*str == 0) {
                x += ch_info.width;
            } else {
                x += ch_info.xadvance;
            }
        }
        return x;
    }

    template <typename FONT>
    constexpr int32_t draw_string_mono(int32_t x, int32_t y, const char *str, uint8_t col) {
        static_assert(FONT::mono == true);
        while (*str != 0) {
            uint32_t utf32 = 0;
            uint32_t lead = static_cast<uint32_t>(*str);
            if (lead < 0x80) {
                utf32 = lead;
                str += 1;
            } else if ((lead >> 5) == 0x06) {
                utf32 = ((lead & 0x1F) << 6) | (static_cast<uint32_t>(str[1]) & 0x3F);
                str += 2;
            } else if ((lead >> 4) == 0x0E) {
                utf32 = ((lead & 0x0F) << 12) | ((static_cast<uint32_t>(str[1]) & 0x3F) << 6) | (static_cast<uint32_t>(str[3]) & 0x3F);
                str += 3;
            } else if ((lead >> 3) == 0x1E) {
                utf32 = ((lead & 0x07) << 18) | ((static_cast<uint32_t>(str[1]) & 0x3F) << 12) | ((static_cast<uint32_t>(str[2]) & 0x3F) << 6) |
                        (static_cast<uint32_t>(str[3]) & 0x3F);
                str += 4;
            } else {
                return x;
            }
#if __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wimplicit-int-conversion"
#endif  // #if __clang__
            const char_info &ch_info = FONT::char_table.at(FONT::glyph_tree.lookup(utf32));
#if __clang__
#pragma GCC diagnostic pop
#endif  // #if __clang__
            draw_char_mono<FONT>(x, y, ch_info, col);
            x += ch_info.xadvance;
        }
        return x;
    }

    template <typename FONT>
    constexpr int32_t draw_string(int32_t x, int32_t y, const char *str, uint8_t col) {
        static_assert(FONT::mono == false);
        while (*str != 0) {
            uint32_t utf32 = 0;
            uint32_t lead = static_cast<uint32_t>(*str);
            if (lead < 0x80) {
                utf32 = lead;
                str += 1;
            } else if ((lead >> 5) == 0x06) {
                utf32 = ((lead & 0x1F) << 6) | (static_cast<uint32_t>(str[1]) & 0x3F);
                str += 2;
            } else if ((lead >> 4) == 0x0E) {
                utf32 = ((lead & 0x0F) << 12) | ((static_cast<uint32_t>(str[1]) & 0x3F) << 6) | (static_cast<uint32_t>(str[3]) & 0x3F);
                str += 3;
            } else if ((lead >> 3) == 0x1E) {
                utf32 = ((lead & 0x07) << 18) | ((static_cast<uint32_t>(str[1]) & 0x3F) << 12) | ((static_cast<uint32_t>(str[2]) & 0x3F) << 6) |
                        (static_cast<uint32_t>(str[3]) & 0x3F);
                str += 4;
            } else {
                return x;
            }
#if __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wimplicit-int-conversion"
#endif  // #if __clang__
            const char_info &ch_info = FONT::char_table.at(FONT::glyph_tree.lookup(utf32));
#if __clang__
#pragma GCC diagnostic pop
#endif  // #if __clang__
            draw_char<FONT>(x, y, ch_info, col);
            x += ch_info.xadvance;
        }
        return x;
    }

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

    constexpr void fill_rect(const rect<int32_t> &r, uint8_t col) {
        fill_rect(r.x, r.y, r.w, r.h, col);
    }

    constexpr void fill_circle(int32_t x, int32_t y, int32_t r, uint8_t col) {
        span(x - abs(r), 2 * abs(r) + 1, y, col);
        fill_arc(x, y, abs(r), 3, 0, col);
    }

    [[nodiscard]] constexpr std::array<uint32_t, W * H> RGBA_uint32() const {
        return T<W, H, S>::RGBA_uint32(data);
    }

    [[nodiscard]] constexpr std::array<uint8_t, W * H * 4> RGBA_uint8() const {
        return T<W, H, S>::RGBA_uint8(data);
    }

    constexpr void blit_RGBA(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{x, y, w, h};
        blitrect &= {0, 0, W, H};
        blitrect &= {x, y, iw, ih};
        T<W, H, S>::blit_RGBA(data, blitrect, ptr, stride);
    }

    constexpr void blit_RGBA(const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{r};
        blitrect &= {0, 0, W, H};
        blitrect &= {r.x, r.y, iw, ih};
        T<W, H, S>::blit_RGBA(data, blitrect, ptr, stride);
    }

    constexpr void blit_RGBA_diffused(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{x, y, w, h};
        blitrect &= {0, 0, W, H};
        blitrect &= {x, y, iw, ih};
        T<W, H, S>::blit_RGBA_diffused(data, blitrect, ptr, stride);
    }

    constexpr void blit_RGBA_diffused(const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{r};
        blitrect &= {0, 0, W, H};
        blitrect &= {r.x, r.y, iw, ih};
        T<W, H, S>::blit_RGBA_diffused(data, blitrect, ptr, stride);
    }

    constexpr void blit_RGBA_diffused_linear(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{x, y, w, h};
        blitrect &= {0, 0, W, H};
        blitrect &= {x, y, iw, ih};
        T<W, H, S>::blit_RGBA_diffused_linear(data, blitrect, ptr, stride);
    }

    constexpr void blit_RGBA_diffused_linear(const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{r};
        blitrect &= {0, 0, W, H};
        blitrect &= {r.x, r.y, iw, ih};
        T<W, H, S>::blit_RGBA_diffused_linear(data, blitrect, ptr, stride);
    }

    template <typename F>
    constexpr void png(F &&charOut) const {
        T<W, H, S>::png(data, charOut);
    }

    template <typename F>
    constexpr void sixel(F &&charOut) const {
        T<W, H, S>::sixel(data, charOut, {0, 0, W, H});
    }

    template <typename F>
    constexpr void sixel(F &&charOut, const rect<int32_t> &r) const {
        T<W, H, S>::sixel(data, charOut, r);
    }

    constexpr void sixel_to_cout() const {
        std::string out;
        T<W, H, S>::sixel(data,
                          [&out](char ch) mutable {
                              out.push_back(ch);
                          },
                          {0, 0, W, H});
        std::cout << out << std::endl;
    }

   private:
    template <typename FONT>
    constexpr void draw_char_mono(int32_t x, int32_t y, const char_info &ch, uint8_t col) {
        static_assert(FONT::mono == true);
        int32_t ch_data_off = static_cast<int32_t>(ch.y) * static_cast<int32_t>(FONT::glyph_bitmap_stride) + static_cast<int32_t>(ch.x) / 8;
        x += ch.xoffset;
        y += ch.yoffset;
        const int32_t x2 = x + static_cast<int32_t>(ch.width);
        const int32_t y2 = y + static_cast<int32_t>(ch.height);
        for (int32_t yy = y; yy < y2; yy++) {
            for (int32_t xx = x; xx < x2; xx++) {
                const int32_t x_off = (xx - x) + static_cast<int32_t>((ch.x % 8));
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

    template <typename FONT>
    constexpr void draw_char(int32_t x, int32_t y, const char_info &ch, uint8_t col) {
        static_assert(FONT::mono == false);
        int32_t ch_data_off = static_cast<int32_t>(ch.y) * static_cast<int32_t>(FONT::glyph_bitmap_stride) + static_cast<int32_t>(ch.x) / 2;
        x += ch.xoffset;
        y += ch.yoffset;
        const int32_t x2 = x + static_cast<int32_t>(ch.width);
        const int32_t y2 = y + static_cast<int32_t>(ch.height);
        for (int32_t yy = y; yy < y2; yy++) {
            for (int32_t xx = x; xx < x2; xx++) {
                const int32_t x_off = (xx - x) + static_cast<int32_t>((ch.x % 2));
                const int32_t bit_index = (1 - (x_off % 2)) * 4;
                const size_t byte_index = static_cast<size_t>(ch_data_off + x_off / 2);
                if (byte_index < (FONT::glyph_bitmap_stride * FONT::glyph_bitmap_height)) {
                    const uint8_t a = (FONT::glyph_bitmap[byte_index] >> bit_index) & 0x7;
                    plot(xx, yy, a ? col : 0);
                }
            }
            ch_data_off += FONT::glyph_bitmap_stride;
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
        T<W, H, S>::span(data, _xl, _xr, _y, col);
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
                    span(x0 - x, 2 * x + delta, y0 + y, col);
                if (corners & 2)
                    span(x0 - x, 2 * x + delta, y0 - y, col);
            }
            if (x != px) {
                if (corners & 1)
                    span(x0 - py, 2 * py + delta, y0 + px, col);
                if (corners & 2)
                    span(x0 - py, 2 * py + delta, y0 - px, col);
                px = x;
            }
            py = y;
        }
    }

    std::array<uint8_t, T<W, H, S>::image_size> data{};
    T<W, H, S> format{};
};

}  // namespace constixel

#endif  // CONSTIXEL_H_
