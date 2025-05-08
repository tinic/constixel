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
        return 1.055f * fast_pow(c, 1.0f / 2.4f) - 0.055f;
    }
}

[[nodiscard]] static constexpr float srgb_to_linear(float s) {
    if (s <= 0.040449936f) {
        return s / 12.92f;
    } else {
        return fast_pow((s + 0.055f) / 1.055f, 2.4f);
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
    std::array<node, N> nodes{};

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
            char_out(static_cast<char>(array[c]));
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
        png_write_be(char_out, i - 4);
        png_write_array(char_out, header, i);
        png_write_crc32(char_out, header, i);
    }

    template <typename F, typename P>
    static constexpr void png_palette(F &&char_out, const P &palette) {
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
        png_write_be(char_out, i - 4);
        png_write_array(char_out, header, i);
        png_write_crc32(char_out, header, i);
    }

    template <typename F>
    static constexpr void png_end(F &&char_out) {
        std::array<char, 4> header;
        uint32_t i = 0;
        header[i++] = 'I';
        header[i++] = 'E';
        header[i++] = 'N';
        header[i++] = 'D';
        png_write_be(char_out, i - 4);
        png_write_array(char_out, header, i);
        png_write_crc32(char_out, header, i);
    }

    template <typename F>
    static constexpr void png_idat_zlib_header(F &&char_out) {
        std::array<char, 6> header;
        uint32_t i = 0;
        header[i++] = 'I';
        header[i++] = 'D';
        header[i++] = 'A';
        header[i++] = 'T';
        header[i++] = 0x78;
        header[i++] = 0x01;
        png_write_be(char_out, i - 4);
        png_write_array(char_out, header, i);
        png_write_crc32(char_out, header, i);
    }

    template <typename F>
    [[nodiscard]] static constexpr uint32_t png_idat_zlib_stream(F &&char_out, const uint8_t *line, size_t bytes, uint32_t adler32_sum) {
        while (bytes > 0) {
            const size_t max_data_use = 1024;
            const size_t extra_data = 24;
            const size_t max_stack_use = max_data_use + extra_data;
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
        header[i++] = 'I';
        header[i++] = 'D';
        header[i++] = 'A';
        header[i++] = 'T';
        header[i++] = static_cast<char>((adler32_sum >> 24) & 0xFF);
        header[i++] = static_cast<char>((adler32_sum >> 16) & 0xFF);
        header[i++] = static_cast<char>((adler32_sum >> 8) & 0xFF);
        header[i++] = static_cast<char>((adler32_sum >> 0) & 0xFF);
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
            PBT idx = col >> 5;
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
                        stack[count++] = static_cast<PBT>(c * 32 + d);
                    }
                }
            }
            return count;
        }

        std::array<uint32_t, PBS / 32> set{};
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
            sixel_color(char_out, static_cast<uint16_t>(c), palette.data()[c]);
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
                    for (size_t xr = (x + 1); xr < (std::min(x + 255, W * S)); xr++) {
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

template <size_t W, size_t H, int32_t S, bool GR>
class format_1bit : public format {
 public:
    static constexpr bool grayscale = GR;
    static constexpr size_t bits_per_pixel = 1;
    static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
    static constexpr size_t internal_height = ((H + 5) / 6) * 6;
    static constexpr size_t image_size = internal_height * bytes_per_line;
    static constexpr std::array<uint32_t, (1UL << bits_per_pixel)> palette = {0x00000000, 0x00ffffff};

    static constexpr void compose(std::array<uint8_t, image_size> &, size_t, size_t, float, float, float, float) {
        static_assert(false, "composing not supported on 1-bit format, use a mono font.");
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
                rgba[y * W + x] = palette[col] | (col ? 0xFF000000 : 0x00000000);
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
                int32_t inc = y % S;
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
                int32_t inc = y % S;
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
};

template <size_t W, size_t H, int32_t S, bool GR>
class format_2bit : public format {
 public:
    static constexpr bool grayscale = GR;
    static constexpr size_t bits_per_pixel = 2;
    static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
    static constexpr size_t internal_height = ((H + 5) / 6) * 6;
    static constexpr size_t image_size = internal_height * bytes_per_line;
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

    static constexpr void compose(std::array<uint8_t, image_size> &, size_t, size_t, float, float, float, float) {
        static_assert(false, "composing not supported on 2-bit format, use a mono font.");
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
                rgba[y * W + x] = palette[col] | (col ? 0xFF000000 : 0x00000000);
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
                int32_t inc = y % S;
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
                int32_t inc = y % S;
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
};

template <size_t W, size_t H, int32_t S, bool GR>
class format_4bit : public format {
 public:
    static constexpr bool grayscale = GR;
    static constexpr size_t bits_per_pixel = 4;
    static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
    static constexpr size_t internal_height = ((H + 5) / 6) * 6;
    static constexpr size_t image_size = internal_height * bytes_per_line;

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

    static constexpr void compose(std::array<uint8_t, image_size> &data, size_t x, size_t y, float cola, float colr, float colg, float colb) {
        size_t bg = static_cast<size_t>(get_col(data, x,y));
        float Rl = colr + quant.linearpal[bg * 3 + 0] * (1.0f - cola);
        float Gl = colg + quant.linearpal[bg * 3 + 1] * (1.0f - cola);
        float Bl = colb + quant.linearpal[bg * 3 + 2] * (1.0f - cola);
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
                rgba[y * W + x] = palette[col] | (col ? 0xFF000000 : 0x00000000);
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
                int32_t inc = y % S;
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
                int32_t inc = y % S;
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
};

template <size_t W, size_t H, int32_t S, bool GR>
class format_8bit : public format {
 public:
    static constexpr bool grayscale = GR;
    static constexpr size_t bits_per_pixel = 8;
    static constexpr size_t bytes_per_line = W;
    static constexpr size_t internal_height = ((H + 5) / 6) * 6;
    static constexpr size_t image_size = internal_height * bytes_per_line;

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
        float Rl = colr + quant.linearpal[bg * 3 + 0] * (1.0f - cola);
        float Gl = colg + quant.linearpal[bg * 3 + 1] * (1.0f - cola);
        float Bl = colb + quant.linearpal[bg * 3 + 2] * (1.0f - cola);
        plot(data, x, y, quant.nearest_linear(Rl, Gl, Bl));
    }

    [[nodiscard]] static constexpr auto RGBA_uint32(const std::array<uint8_t, image_size> &data) {
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

    [[nodiscard]] static constexpr auto RGBA_uint8(const std::array<uint8_t, image_size> &data) {
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
                int32_t inc = y % S;
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
                int32_t inc = y % S;
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

template <template <size_t, size_t, int32_t, bool> class T, size_t W, size_t H, int32_t S = 1, bool GR = false>
class image {
    static_assert(sizeof(W) >= sizeof(uint32_t));
    static_assert(sizeof(H) >= sizeof(uint32_t));

    static_assert(W > 0 && H > 0);
    static_assert((W * S) <= 16384 && (H * S) <= 16384);
    static_assert(S >= 1 && S <= 256);

 public:
    [[nodiscard]] static constexpr size_t grayscale() {
        return T<W, H, S, GR>::grayscale;
    }

    [[nodiscard]] static constexpr size_t bit_depth() {
        return T<W, H, S, GR>::bits_per_pixel;
    }

    [[nodiscard]] static constexpr size_t size() {
        return T<W, H, S, GR>::image_size;
    }

    [[nodiscard]] static constexpr int32_t width() {
        return static_cast<int32_t>(W);
    }

    [[nodiscard]] static constexpr int32_t height() {
        return static_cast<int32_t>(H);
    }

    constexpr void clear() {
        data.fill(0);
    }

    template <typename abs_T>
    [[nodiscard]] static constexpr abs_T abs(abs_T v) {
        return v < 0 ? -v : v;
    }

    [[nodiscard]] constexpr std::array<uint8_t, T<W, H, S, GR>::image_size> &data_ref() const {
        return data;
    }

    [[nodiscard]] constexpr image<T, W, H, S, GR> clone() const {
        return *this;
    }

    constexpr void copy(const image<T, W, H, S, GR> &src) {
        data = src.data;
    }

    constexpr void line(int32_t x0, int32_t y0, int32_t x1, int32_t y1, uint8_t col, uint32_t stroke_width = 1) {
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

    constexpr void line(const rect<int32_t> &l, uint8_t col, uint32_t stroke_width = 1) {
        line(l.x, l.y, l.x + l.w, l.y + l.h, col, stroke_width);
    }

    constexpr void compose(int32_t x, int32_t y, float cola, float colr, float colg, float colb) {
        if (x < 0 || x >= static_cast<int32_t>(W) || y < 0 || y >= static_cast<int32_t>(H)) {
            return;
        }
        T<W, H, S, GR>::compose(data, static_cast<uint32_t>(x), static_cast<uint32_t>(y), cola, colr, colg, colb);
    }

    constexpr void plot(int32_t x, int32_t y, uint8_t col) {
        if (x < 0 || x >= static_cast<int32_t>(W) || y < 0 || y >= static_cast<int32_t>(H)) {
            return;
        }
        T<W, H, S, GR>::plot(data, static_cast<uint32_t>(x), static_cast<uint32_t>(y), col);
    }

    [[nodiscard]] constexpr uint8_t get_nearest_color(uint8_t r, uint8_t g, uint8_t b) const {
        return format.quant.nearest(r, g, b);
    }

    template <typename FONT>
    [[nodiscard]] constexpr int32_t string_width(const char *str) {
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
            const char_info &ch_info = FONT::char_table.at(FONT::glyph_tree.lookup(static_cast<FONT::lookup_type>(utf32)));
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
        static_assert(FONT::mono == true, "Can't use a antialiased font to draw mono/pixelized text.");
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
            const char_info &ch_info = FONT::char_table.at(FONT::glyph_tree.lookup(static_cast<FONT::lookup_type>(utf32)));
            draw_char_mono<FONT>(x, y, ch_info, col);
            x += ch_info.xadvance;
        }
        return x;
    }

    template <typename FONT>
    constexpr int32_t draw_string_aa(int32_t x, int32_t y, const char *str, uint8_t col) {
        static_assert(FONT::mono == false, "Can't use a mono font to draw antialiased text.");
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
            const char_info &ch_info = FONT::char_table.at(FONT::glyph_tree.lookup(static_cast<FONT::lookup_type>(utf32)));
            draw_char_aa<FONT>(x, y, ch_info, col);
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

    constexpr void stroke_rect(int32_t x, int32_t y, int32_t w, int32_t h, uint8_t col, uint32_t stroke_width = 1) {
        line(x, y, x + w, y, col, stroke_width);
        line(x + w, y, x + w, y + h, col, stroke_width);
        line(x + w, y + h, x, y + h, col, stroke_width);
        line(x, y + h, x, y, col, stroke_width);
    }

    constexpr void stroke_rect(const rect<int32_t> &r, uint8_t col, uint32_t stroke_width = 1) {
        stroke_rect(r.x, r.y, r.w, r.h, col, stroke_width);
    }

    constexpr void fill_circle(int32_t cx, int32_t cy, int32_t r, uint8_t col) {
        if (r == 1) {
            fill_rect(cx - 1, cy - 1, 2, 2, col);
            return;
        }
        fill_arc(cx, cy - 1, r, 9, 0, col);
        fill_arc(cx, cy, r, 10, 0, col);
        fill_arc(cx, cy - 1, r, 1, -1, col);
        fill_arc(cx, cy, r, 2, -1, col);
    }

    constexpr void fill_round_rect(int32_t x, int32_t y, int32_t w, int32_t h, int32_t r, uint8_t col) {
        int32_t cr = std::min((w) / 2, std::min((w) / 2, r));
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

    constexpr void fill_circle_aa(int32_t cx, int32_t cy, int32_t r, uint8_t col) {
        fill_circle_aa(cx, cy, r, 0, 0, col);
    }

    constexpr void fill_round_rect_aa(int32_t x, int32_t y, int32_t w, int32_t h, int32_t r, uint8_t col) {
        int32_t cr = std::min((w) / 2, std::min((w) / 2, r));
        int32_t dx = w - cr * 2;
        int32_t dy = h - cr * 2;
        fill_circle_aa(x + cr, y + cr, cr, dx, dy, col);
        fill_rect(x, y + cr, cr, dy, col);
        fill_rect(x + w - cr, y + cr, cr, dy, col);
        fill_rect(x + cr, y, dx, h, col);
    }

    [[nodiscard]] constexpr std::array<uint32_t, W * H> RGBA_uint32() const {
        return T<W, H, S, GR>::RGBA_uint32(data);
    }

    [[nodiscard]] constexpr std::array<uint8_t, W * H * 4> RGBA_uint8() const {
        return T<W, H, S, GR>::RGBA_uint8(data);
    }

    constexpr void blit_RGBA(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{x, y, w, h};
        blitrect &= {0, 0, W, H};
        blitrect &= {x, y, iw, ih};
        T<W, H, S, GR>::blit_RGBA(data, blitrect, ptr, stride);
    }

    constexpr void blit_RGBA(const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{r};
        blitrect &= {0, 0, W, H};
        blitrect &= {r.x, r.y, iw, ih};
        T<W, H, S, GR>::blit_RGBA(data, blitrect, ptr, stride);
    }

    constexpr void blit_RGBA_diffused(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{x, y, w, h};
        blitrect &= {0, 0, W, H};
        blitrect &= {x, y, iw, ih};
        T<W, H, S, GR>::blit_RGBA_diffused(data, blitrect, ptr, stride);
    }

    constexpr void blit_RGBA_diffused(const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{r};
        blitrect &= {0, 0, W, H};
        blitrect &= {r.x, r.y, iw, ih};
        T<W, H, S, GR>::blit_RGBA_diffused(data, blitrect, ptr, stride);
    }

    constexpr void blit_RGBA_diffused_linear(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{x, y, w, h};
        blitrect &= {0, 0, W, H};
        blitrect &= {x, y, iw, ih};
        T<W, H, S, GR>::blit_RGBA_diffused_linear(data, blitrect, ptr, stride);
    }

    constexpr void blit_RGBA_diffused_linear(const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{r};
        blitrect &= {0, 0, W, H};
        blitrect &= {r.x, r.y, iw, ih};
        T<W, H, S, GR>::blit_RGBA_diffused_linear(data, blitrect, ptr, stride);
    }

    template <typename F>
    constexpr void png(F &&char_out) const {
        T<W, H, S, GR>::png(data, char_out);
    }

    template <typename F>
    constexpr void sixel(F &&char_out) const {
        T<W, H, S, GR>::sixel(data, char_out, {0, 0, W, H});
    }

    template <typename F>
    constexpr void sixel(F &&char_out, const rect<int32_t> &r) const {
        T<W, H, S, GR>::sixel(data, char_out, r);
    }

    constexpr void sixel_to_cout() const {
        std::string out;
        T<W, H, S, GR>::sixel(data,
                              [&out](char ch) mutable {
                                  out.push_back(ch);
                              },
                              {0, 0, W, H});
        std::cout << out << std::endl;
    }

    void png_to_iterm() const {
        size_t buffer = 0;
        size_t bits_collected = 0;
        std::string output;
        output.append("\033]1337;File=inline=1:");
        T<W, H, S, GR>::png(data, [&buffer, &bits_collected, &output](char byte) mutable {
            static constexpr char base64_chars[] =
                "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        
            buffer = (buffer << 8) | static_cast<uint8_t>(byte);
            bits_collected += 8;
            while (bits_collected >= 6) {
                bits_collected -= 6;
                output.push_back(base64_chars[(buffer >> bits_collected) & 0x3F]);
            }
            return [&output]() mutable {
                if (output.capacity() == 0) return;
                if (output.size() % 4) {
                    size_t padding = 4 - output.size() % 4;
                    while (padding--) output.push_back('=');
                }
            };
        });
        output.append("\07");
        std::cout << output << std::endl;
    }

 private:
    constexpr void fill_circle_aa(int32_t cx, int32_t cy, int32_t r, int32_t ox, int32_t oy, uint8_t col) {
        int32_t x0 = cx - r - 1;
        int32_t y0 = cy - r - 1;
        float Rl = format.quant.linearpal[col * 3 + 0];
        float Gl = format.quant.linearpal[col * 3 + 1];
        float Bl = format.quant.linearpal[col * 3 + 2];
        for (int y = y0; y <= cy; ++y) {
            for (int x = x0; x <= cx; ++x) {
                float dx = (x + 0.5f) - cx;
                float dy = (y + 0.5f) - cy;
                float dist_sq = dx * dx + dy * dy;
                float a = r - fast_sqrtf(dist_sq);
                a = std::clamp(a + 0.5f, 0.0f, 1.0f);
                if (a != 0.0f) {
                    int32_t lx = x;
                    int32_t ly = y;
                    int32_t rx = cx + (x0 - x) + r + ox;
                    int32_t ry = cy + (y0 - y) + r + oy;
                    if (a >= 1.0f) {
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

    template <typename FONT>
    constexpr void draw_char_mono(int32_t x, int32_t y, const char_info &ch, uint8_t col) {
        static_assert(FONT::mono == true, "Can't use a antialiased font to draw mono/pixelized text.");
        int32_t ch_data_off = static_cast<int32_t>(ch.y) * static_cast<int32_t>(FONT::glyph_bitmap_stride) + static_cast<int32_t>(ch.x) / 8;
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
    }

    static consteval auto gen_a2al() {
        std::array<float, 16> a2alg{};
        for (size_t c = 0; c < 16; c++) {
            a2alg[c] = constixel::srgb_to_linear(static_cast<float>(c) * (1.0f / 15.0f));
        }
        return a2alg;
    }

    static constexpr std::array<float, 16> a2al = gen_a2al();

    template <typename FONT>
    constexpr void draw_char_aa(int32_t x, int32_t y, const char_info &ch, uint8_t col) {
        static_assert(FONT::mono == false, "Can't use a mono font to draw antialiased text.");
        int32_t ch_data_off = static_cast<int32_t>(ch.y) * static_cast<int32_t>(FONT::glyph_bitmap_stride) + static_cast<int32_t>(ch.x) / 2;
        x += ch.xoffset;
        y += ch.yoffset;
        const int32_t x2 = x + static_cast<int32_t>(ch.width);
        const int32_t y2 = y + static_cast<int32_t>(ch.height);
        float Rl = format.quant.linearpal[col * 3 + 0];
        float Gl = format.quant.linearpal[col * 3 + 1];
        float Bl = format.quant.linearpal[col * 3 + 2];
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
