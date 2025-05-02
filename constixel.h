#ifndef _SIXEL_H_
#define _SIXEL_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#include <bit>

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

static constexpr float fast_log(const float x) {
    return fast_log2(x) * 0.69314718f;
}

static constexpr float fast_pow(const float x, const float p) {
#ifdef GCC_BROKEN_BITCAST
    return powf(x, p);
#else   // GCC_BROKEN_BITCAST
    return fast_exp2(p * fast_log2(x));
#endif  // GCC_BROKEN_BITCAST
}

struct node {
    int16_t child[8]{-1, -1, -1, -1, -1, -1, -1, -1};
    int16_t palette = -1;
};

template <size_t S, size_t N>
class octree_impl {
    size_t node_idx = 0;
    static constexpr size_t palette_size = S;
    static constexpr size_t node_size = N;
    std::array<node, node_size> nodes {};
    const std::array<uint32_t, palette_size> &pal;

    static constexpr size_t child_index(uint8_t r, uint8_t g, uint8_t b, int32_t level) {
        return static_cast<size_t>(((r >> level) & 1) << 2 | ((g >> level) & 1) << 1 | ((b >> level) & 1));
    }

    constexpr void insert(uint32_t c, size_t idx) {
        size_t n = 0;
        for (int32_t lvl = 7; lvl >= 0; --lvl) {
            size_t ci = child_index((c >> 16) & 0xFF, (c >> 8) & 0xFF, (c >> 0) & 0xFF, lvl);
            if (nodes[n].child[ci] == -1) {
                nodes[n].child[ci] = static_cast<int16_t>(node_idx);
                node_idx++;
            }
            n = static_cast<size_t>(nodes[n].child[ci]);
        }
        nodes[n].palette = static_cast<int16_t>(idx);
    }

   public:
    constexpr octree_impl(const std::array<uint32_t, palette_size> &palette) : pal(palette) {
        for (size_t i = 0; i < palette_size; ++i) {
            insert(palette[i], i);
        }
    }

    constexpr uint8_t nearest(uint8_t r, uint8_t g, uint8_t b) const {
        size_t n = 0;
        for (int32_t lvl = 7; lvl >= 0; --lvl) {
            size_t ci = child_index(r, g, b, lvl);
            size_t next = static_cast<size_t>(nodes[n].child[ci]);
            if (next == -1)
                break;
            n = next;
            if (nodes[n].palette != -1)
                return static_cast<uint8_t>(nodes[n].palette);
        }
        // fallback: brute-force
        int32_t best = 0, bestd = 1UL<<30;
        for (size_t i = 0; i < pal.size(); ++i) {
            int32_t dr = static_cast<int32_t>(r) - static_cast<int32_t>((pal[i] >> 16) & 0xFF);
            int32_t dg = static_cast<int32_t>(g) - static_cast<int32_t>((pal[i] >> 8) & 0xFF);
            int32_t db = static_cast<int32_t>(b) - static_cast<int32_t>((pal[i] >> 0) & 0xFF);
            int32_t d = dr * dr + dg * dg + db * db;
            if (d < bestd) {
                bestd = d;
                best = static_cast<int32_t>(i);
            }
        }
        return static_cast<uint8_t>(best);
    }

    constexpr size_t nodes_used_length() const {
        return node_idx;
    }
    constexpr size_t nodes_memory_length() const {
        return node_size;
    }
};

static constexpr double m_pi_d = 3.14159265358979323846;

static consteval double cos(double x, int32_t terms = 10) {
    x = x - 6.283185307179586 * int(x / 6.283185307179586);  // wrap x to [0, 2π)
    double res = 1.0, term = 1.0;
    double x2 = x * x;
    for (int32_t i = 1; i < terms; ++i) {
        term *= -x2 / ((2 * i - 1) * (2 * i));
        res += term;
    }
    return res;
}

static consteval double sin(double x, int32_t terms = 10) {
    x = x - 6.283185307179586 * int(x / 6.283185307179586);  // wrap x to [0, 2π)
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

template <std::size_t N, typename K, typename V>
class map {
   public:
    using key_type = K;
    using value_type = V;
    using size_type = decltype(N);
    using data_type = std::array<std::pair<K, V>, N>;
    using const_iterator = typename data_type::const_iterator;

    template <typename... E>
    explicit constexpr map(E &&...elements) noexcept : data{std::forward<E>(elements)...} {
        static_assert(N > 0, "N should be positive");
        static_assert(N == sizeof...(elements), "Elements size doesn't match expected size of a hash-map");
    }

    [[nodiscard]] constexpr const_iterator find(const K &key) const noexcept {
        return search<0, N>(key);
    }

    [[nodiscard]] constexpr std::pair<bool, const V &> at(const K &key) const noexcept {
        const auto it = find(key);

        if (it != cend()) {
            return {true, it->second};
        }

        return {false, {}};
    }

    [[nodiscard]] constexpr const V &operator[](const K &key) const noexcept {
        return find(key)->second;
    }

    [[nodiscard]] constexpr bool contains(const K &key) const noexcept {
        return search<0, N>(key) != cend();
    }

    [[nodiscard]] constexpr size_type size() const noexcept {
        return data.size();
    }

    [[nodiscard]] constexpr const_iterator begin() const noexcept {
        return cbegin();
    }

    [[nodiscard]] constexpr const_iterator cbegin() const noexcept {
        return std::cbegin(data);
    }

    [[nodiscard]] constexpr const_iterator end() const noexcept {
        return cend();
    }

    [[nodiscard]] constexpr const_iterator cend() const noexcept {
        return std::cend(data);
    }

    [[nodiscard]] constexpr bool empty() const noexcept {
        return false;
    }

   protected:
    using index_type = size_type;

    template <index_type L, index_type R>
    [[nodiscard]] constexpr const_iterator search(const K &key) const noexcept {
        if constexpr (L < R) {
            if (equal(data[L].first, key)) {
                return std::next(cbegin(), L);
            }

            return search<L + 1, R>(key);
        }

        return cend();
    }

    template <typename T = K>
    [[nodiscard]] constexpr bool equal(const T &lhs, const T &rhs) const noexcept {
        return lhs == rhs;
    }

    [[nodiscard]] constexpr bool equal(const char *lhs, const char *rhs) const noexcept {
        return *lhs == *rhs && (*lhs == '\0' || equal(lhs + 1, rhs + 1));
    }

   private:
    data_type data;
};

struct character {
    uint32_t id;
    int16_t x;
    int16_t y;
    int16_t width;
    int16_t height;
    int16_t xadvance;
    int16_t xoffset;
    int16_t yoffset;
};

struct font {
    std::string face;
    int8_t size;
    int8_t bold : 1;
    int8_t italic : 1;
    int8_t outline : 1;
};

struct kerning {
    int16_t amount;
    int16_t first;
    int16_t second;
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
    static constexpr void sixel_color(F &&charOut, uint16_t index, uint32_t col) {
        charOut('#');
        sixel_number(charOut, index);
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

    template <size_t W, size_t H, int32_t S, typename P, typename F, typename C>
    static constexpr void sixel_image(const uint8_t *data, const P &palette, F &&charOut, const rect<int32_t> &_r, bool preserveBackground, const C &collect6) {
        sixel_header(charOut);
        if (!preserveBackground) {
            sixel_raster_attributes<W, H, S>(charOut);
        }
        for (size_t c = 0; c < palette.size(); c++) {
            sixel_color(charOut, static_cast<uint16_t>(c), palette.data()[c]);
        }
        const auto r = rect<int32_t>{_r.x * S, _r.y * S, _r.w * S, _r.h * S} & rect<int32_t>{0, 0, W * S, H * S};
        for (size_t y = static_cast<size_t>(r.y); y < static_cast<size_t>(r.y + r.h); y += 6) {
            for (size_t c = 0; c < palette.size(); c++) {
                uint8_t test6 = 0;
                for (size_t x = static_cast<size_t>(r.x); x < static_cast<size_t>(r.x + r.w); x++) {
                    test6 |= collect6(data, x, c, y);
                }
                if (!test6) {
                    continue;
                }
                if (c != 0) {
                    charOut('$');
                }
                charOut('#');
                sixel_number(charOut, static_cast<uint16_t>(c));
                for (size_t x = static_cast<size_t>(r.x); x < static_cast<size_t>(r.x + r.w); x++) {
                    uint8_t bits6 = collect6(data, x, c, y);
                    uint16_t repeatCount = 0;
                    for (size_t xr = (x + 1); xr < (std::min(x + 255, W * S)); xr++) {
                        if (bits6 == collect6(data, xr, c, y)) {
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

    static constexpr void blitRGBA(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                uint32_t R = ptr[y * static_cast<size_t>(stride) + x * 4 + 0];
                uint32_t G = ptr[y * static_cast<size_t>(stride) + x * 4 + 1];
                uint32_t B = ptr[y * static_cast<size_t>(stride) + x * 4 + 2];
                plot(data, x + static_cast<size_t>(r.x), y + static_cast<size_t>(r.y), (R * 2 + G * 3 + B * 1) > 768 ? 1 : 0);
            }
        }
    }

    static constexpr void blitRGBADiffused(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih,
                                           int32_t stride) {
        int32_t err_r = 0;
        int32_t err_g = 0;
        int32_t err_b = 0;
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                int32_t R = ptr[y * static_cast<size_t>(stride) + x * 4 + 0];
                int32_t G = ptr[y * static_cast<size_t>(stride) + x * 4 + 1];
                int32_t B = ptr[y * static_cast<size_t>(stride) + x * 4 + 2];
                R = std::clamp(R + err_r, 0, 255);
                G = std::clamp(G + err_g, 0, 255);
                B = std::clamp(B + err_b, 0, 255);
                uint8_t n = (R * 2 + G * 3 + B * 1) > 768 ? 1 : 0;
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                err_r = R - (n ? 0xFF : 0x00);
                err_g = G - (n ? 0xFF : 0x00);
                err_b = B - (n ? 0xFF : 0x00);
            }
        }
    }

    static constexpr void blitRGBADiffusedLinear(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih,
                                                 int32_t stride) {
        float err_r = 0;
        float err_g = 0;
        float err_b = 0;
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                float R = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 0]) * (1.0f / 255.0f);
                float G = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 1]) * (1.0f / 255.0f);
                float B = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 2]) * (1.0f / 255.0f);
                float Rl = constixel::srgb_to_linear(R);
                float Gl = constixel::srgb_to_linear(G);
                float Bl = constixel::srgb_to_linear(B);
                Rl = std::clamp(Rl + err_r, 0.0f, 1.0f);
                Gl = std::clamp(Gl + err_g, 0.0f, 1.0f);
                Bl = std::clamp(Bl + err_b, 0.0f, 1.0f);
                R = constixel::linear_to_srgb(Rl);
                G = constixel::linear_to_srgb(Gl);
                B = constixel::linear_to_srgb(Bl);
                uint8_t n = (R * 2 * +G * 3 + B * 1) > 3 ? 1 : 0;
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                err_r = Rl - constixel::srgb_to_linear(n ? 1.0f : 0.0f);
                err_g = Gl - constixel::srgb_to_linear(n ? 1.0f : 0.0f);
                err_b = Bl - constixel::srgb_to_linear(n ? 1.0f : 0.0f);
            }
        }
    }

    template <typename F>
    static constexpr void sixel(std::array<uint8_t, image_size> &data, F &&charOut, const rect<int32_t> &r, bool preserveBackground) {
        sixel_image<W, H, S>(data.data(), palette, charOut, r, preserveBackground, [](const uint8_t *dataRaw, size_t x, size_t col, size_t y) {
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
        });
    }
};

template <size_t W, size_t H, int32_t S>
class format_2bit : public format {
   public:
    static constexpr size_t octree_size = 25;  // strongly depends on fixed palette below.
    static constexpr size_t bits_per_pixel = 2;
    static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
    static constexpr size_t internal_height = ((H + 5) / 6) * 6;
    static constexpr size_t image_size = internal_height * bytes_per_line;
    static constexpr std::array<uint32_t, (1UL << bits_per_pixel)> palette = {0x000000, 0xffffff, 0xff0000, 0x0077ff};

    static consteval const constixel::octree_impl<1UL << bits_per_pixel, octree_size> gen_octree() {
        return constixel::octree_impl<1UL << bits_per_pixel, octree_size>(palette);
    }

    static constexpr const constixel::octree_impl<1UL << bits_per_pixel, octree_size> octree = gen_octree();

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

    static constexpr void blitRGBA(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                plot(data, (x + static_cast<uint32_t>(r.x)), (y + static_cast<uint32_t>(r.y)), octree.nearest(ptr[y * static_cast<size_t>(stride) + x * 4 + 0], ptr[y * static_cast<size_t>(stride) + x * 4 + 1], ptr[y * static_cast<size_t>(stride) + x * 4 + 2]));
            }
        }
    }

    static constexpr void blitRGBADiffused(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih,
                                           int32_t stride) {
        int32_t err_r = 0;
        int32_t err_g = 0;
        int32_t err_b = 0;
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                int32_t R = ptr[y * static_cast<size_t>(stride) + x * 4 + 0];
                int32_t G = ptr[y * static_cast<size_t>(stride) + x * 4 + 1];
                int32_t B = ptr[y * static_cast<size_t>(stride) + x * 4 + 2];
                R = std::clamp(R + err_r, 0, 255);
                G = std::clamp(G + err_g, 0, 255);
                B = std::clamp(B + err_b, 0, 255);
                uint8_t n = octree.nearest(static_cast<uint8_t>(R), static_cast<uint8_t>(G), static_cast<uint8_t>(B));
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                err_r = R - static_cast<int32_t>((palette[n] >> 16) & 0xFF);
                err_g = G - static_cast<int32_t>((palette[n] >> 8) & 0xFF);
                err_b = B - static_cast<int32_t>((palette[n] >> 0) & 0xFF);
            }
        }
    }

    static constexpr void blitRGBADiffusedLinear(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih,
                                                 int32_t stride) {
        float err_r = 0;
        float err_g = 0;
        float err_b = 0;
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                float R = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 0]) * (1.0f / 255.0f);
                float G = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 1]) * (1.0f / 255.0f);
                float B = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 2]) * (1.0f / 255.0f);
                float Rl = constixel::srgb_to_linear(R);
                float Gl = constixel::srgb_to_linear(G);
                float Bl = constixel::srgb_to_linear(B);
                Rl = std::clamp(Rl + err_r, 0.0f, 1.0f);
                Gl = std::clamp(Gl + err_g, 0.0f, 1.0f);
                Bl = std::clamp(Bl + err_b, 0.0f, 1.0f);
                R = constixel::linear_to_srgb(Rl);
                G = constixel::linear_to_srgb(Gl);
                B = constixel::linear_to_srgb(Bl);
                uint8_t n = octree.nearest(static_cast<uint8_t>(R * 255.0f), static_cast<uint8_t>(G * 255.0f), static_cast<uint8_t>(B * 255.0f));
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                err_r = Rl - constixel::srgb_to_linear(static_cast<float>((palette[n] >> 16) & 0xFF) * (1.0f / 255.0f));
                err_g = Gl - constixel::srgb_to_linear(static_cast<float>((palette[n] >> 8) & 0xFF) * (1.0f / 255.0f));
                err_b = Bl - constixel::srgb_to_linear(static_cast<float>((palette[n] >> 0) & 0xFF) * (1.0f / 255.0f));
            }
        }
    }

    template <typename F>
    static constexpr void sixel(std::array<uint8_t, image_size> &data, F &&charOut, const rect<int32_t> &r, bool preserveBackground) {
        sixel_image<W, H, S>(data.data(), palette, charOut, r, preserveBackground, [](const uint8_t *dataRaw, size_t x, size_t col, size_t y) {
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
        });
    }
};

template <size_t W, size_t H, int32_t S>
class format_4bit : public format {
   public:
    static constexpr size_t octree_size = 70;  // strongly depends on fixed palette below.
    static constexpr size_t bits_per_pixel = 4;
    static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
    static constexpr size_t internal_height = ((H + 5) / 6) * 6;
    static constexpr size_t image_size = internal_height * bytes_per_line;
    static constexpr std::array<uint32_t, (1UL << bits_per_pixel)> palette = {0x000000, 0xffffff, 0xff0000, 0x00ff00, 0x0000ff, 0xffff00, 0x00ffff, 0xff00ff,
                                                                              0x333333, 0x666666, 0x999999, 0xcccccc, 0x7f0000, 0x007f00, 0x00007f, 0x7f7f00};

    static consteval const constixel::octree_impl<1UL << bits_per_pixel, octree_size> gen_octree() {
        return constixel::octree_impl<1UL << bits_per_pixel, octree_size>(palette);
    }

    static constexpr const constixel::octree_impl<1UL << bits_per_pixel, octree_size> octree = gen_octree();

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

    static constexpr void blitRGBA(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), octree.nearest(ptr[y * static_cast<size_t>(stride) + x * 4 + 0], ptr[y * static_cast<size_t>(stride) + x * 4 + 1], ptr[y * static_cast<size_t>(stride) + x * 4 + 2]));
            }
        }
    }

    static constexpr void blitRGBADiffused(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih,
                                           int32_t stride) {
        int32_t err_r = 0;
        int32_t err_g = 0;
        int32_t err_b = 0;
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                int32_t R = ptr[y * static_cast<size_t>(stride) + x * 4 + 0];
                int32_t G = ptr[y * static_cast<size_t>(stride) + x * 4 + 1];
                int32_t B = ptr[y * static_cast<size_t>(stride) + x * 4 + 2];
                R = std::clamp(R + err_r, 0, 255);
                G = std::clamp(G + err_g, 0, 255);
                B = std::clamp(B + err_b, 0, 255);
                uint8_t n = octree.nearest(static_cast<uint8_t>(R), static_cast<uint8_t>(G), static_cast<uint8_t>(B));
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                err_r = R - static_cast<int32_t>((palette[n] >> 16) & 0xFF);
                err_g = G - static_cast<int32_t>((palette[n] >>  8) & 0xFF);
                err_b = B - static_cast<int32_t>((palette[n] >>  0) & 0xFF);
            }
        }
    }

    static constexpr void blitRGBADiffusedLinear(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih,
                                                 int32_t stride) {
        float err_r = 0;
        float err_g = 0;
        float err_b = 0;
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                float R = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 0]) * (1.0f / 255.0f);
                float G = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 1]) * (1.0f / 255.0f);
                float B = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 2]) * (1.0f / 255.0f);
                float Rl = constixel::srgb_to_linear(R);
                float Gl = constixel::srgb_to_linear(G);
                float Bl = constixel::srgb_to_linear(B);
                Rl = std::clamp(Rl + err_r, 0.0f, 1.0f);
                Gl = std::clamp(Gl + err_g, 0.0f, 1.0f);
                Bl = std::clamp(Bl + err_b, 0.0f, 1.0f);
                R = constixel::linear_to_srgb(Rl);
                G = constixel::linear_to_srgb(Gl);
                B = constixel::linear_to_srgb(Bl);
                uint8_t n = octree.nearest(static_cast<uint8_t>(R * 255.0f), static_cast<uint8_t>(G * 255.0f), static_cast<uint8_t>(B * 255.0f));
                plot(data, (x + static_cast<size_t>(r.x)), (y + static_cast<size_t>(r.y)), n);
                err_r = Rl - constixel::srgb_to_linear(static_cast<float>((palette[n] >> 16) & 0xFF) * (1.0f / 255.0f));
                err_g = Gl - constixel::srgb_to_linear(static_cast<float>((palette[n] >> 8) & 0xFF) * (1.0f / 255.0f));
                err_b = Bl - constixel::srgb_to_linear(static_cast<float>((palette[n] >> 0) & 0xFF) * (1.0f / 255.0f));
            }
        }
    }

    template <typename F>
    static constexpr void sixel(std::array<uint8_t, image_size> &data, F &&charOut, const rect<int32_t> &r, bool preserveBackground) {
        sixel_image<W, H, S>(data.data(), palette, charOut, r, preserveBackground, [](const uint8_t *dataRaw, size_t x, size_t col, size_t y) {
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
        });
    }
};

template <size_t W, size_t H, int32_t S>
class format_8bit : public format {
   public:
    static constexpr size_t octree_size = 1061;  // strongly depends on calculated palette below.
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

    static consteval const constixel::octree_impl<1UL << bits_per_pixel, octree_size> gen_octree() {
        return constixel::octree_impl<1UL << bits_per_pixel, octree_size>(palette);
    }

    static constexpr const constixel::octree_impl<1UL << bits_per_pixel, octree_size> octree = gen_octree();

    static constexpr void plot(std::array<uint8_t, image_size> &data, size_t x0, size_t y, uint8_t col) {
        data.data()[y * bytes_per_line + x0] = col;
    }

    static constexpr void span(std::array<uint8_t, image_size> &data, size_t xl0, size_t xr0, size_t y, uint8_t col) {
        uint8_t *yptr = &data.data()[y * bytes_per_line];
        for (size_t x = xl0; x < xr0; x++) {
            yptr[x] = col;
        }
    }

    static constexpr void blitRGBA(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                data.data()[(y + static_cast<size_t>(r.y)) * bytes_per_line + (x + static_cast<size_t>(r.x))] =
                    octree.nearest(ptr[y * static_cast<size_t>(stride) + x * 4 + 0], ptr[y * static_cast<size_t>(stride) + x * 4 + 1], ptr[y * static_cast<size_t>(stride) + x * 4 + 2]);
            }
        }
    }

    static constexpr void blitRGBADiffused(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih,
                                           int32_t stride) {
        int32_t err_r = 0;
        int32_t err_g = 0;
        int32_t err_b = 0;
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                int32_t R = ptr[y * static_cast<size_t>(stride) + x * 4 + 0];
                int32_t G = ptr[y * static_cast<size_t>(stride) + x * 4 + 1];
                int32_t B = ptr[y * static_cast<size_t>(stride) + x * 4 + 2];
                R = std::clamp(R + err_r, 0, 255);
                G = std::clamp(G + err_g, 0, 255);
                B = std::clamp(B + err_b, 0, 255);
                uint8_t n = octree.nearest(static_cast<uint8_t>(R), static_cast<uint8_t>(G), static_cast<uint8_t>(B));
                data.data()[(y + static_cast<size_t>(r.y)) * bytes_per_line + (x + static_cast<size_t>(r.x))] = n;
                err_r = R - static_cast<int32_t>((palette[n] >> 16) & 0xFF);
                err_g = G - static_cast<int32_t>((palette[n] >> 8) & 0xFF);
                err_b = B - static_cast<int32_t>((palette[n] >> 0) & 0xFF);
            }
        }
    }

    static constexpr void blitRGBADiffusedLinear(std::array<uint8_t, image_size> &data, const rect<int32_t> &r, const uint8_t *ptr, int32_t iw, int32_t ih,
                                                 int32_t stride) {
        float err_r = 0;
        float err_g = 0;
        float err_b = 0;
        for (size_t y = 0; y < static_cast<size_t>(r.h); y++) {
            for (size_t x = 0; x < static_cast<size_t>(r.w); x++) {
                float R = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 0]) * (1.0f / 255.0f);
                float G = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 1]) * (1.0f / 255.0f);
                float B = static_cast<float>(ptr[y * static_cast<size_t>(stride) + x * 4 + 2]) * (1.0f / 255.0f);
                float Rl = constixel::srgb_to_linear(R);
                float Gl = constixel::srgb_to_linear(G);
                float Bl = constixel::srgb_to_linear(B);
                Rl = std::clamp(Rl + err_r, 0.0f, 1.0f);
                Gl = std::clamp(Gl + err_g, 0.0f, 1.0f);
                Bl = std::clamp(Bl + err_b, 0.0f, 1.0f);
                R = constixel::linear_to_srgb(Rl);
                G = constixel::linear_to_srgb(Gl);
                B = constixel::linear_to_srgb(Bl);
                uint8_t n = octree.nearest(static_cast<uint8_t>(R * 255.0f), static_cast<uint8_t>(G * 255.0f), static_cast<uint8_t>(B * 255.0f));
                data.data()[(y + static_cast<size_t>(r.y)) * bytes_per_line + (x + static_cast<size_t>(r.x))] = n;
                err_r = Rl - constixel::srgb_to_linear(static_cast<float>((palette[n] >> 16) & 0xFF) * (1.0f / 255.0f));
                err_g = Gl - constixel::srgb_to_linear(static_cast<float>((palette[n] >> 8) & 0xFF) * (1.0f / 255.0f));
                err_b = Bl - constixel::srgb_to_linear(static_cast<float>((palette[n] >> 0) & 0xFF) * (1.0f / 255.0f));
            }
        }
    }

    template <typename F>
    static constexpr void sixel(std::array<uint8_t, image_size> &data, F &&charOut, const rect<int32_t> &r, bool preserveBackground) {
        sixel_image<W, H, S>(data.data(), palette, charOut, r, preserveBackground, [](const uint8_t *dataRaw, size_t x, size_t col, size_t y) {
            const uint8_t *ptr = &dataRaw[(y / S) * bytes_per_line + (x / S)];
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
        });
    }
};

template <template <size_t, size_t, int32_t> class T, size_t W, size_t H, int32_t S = 1>
class image {
    static_assert(sizeof(W) >= sizeof(uint32_t));
    static_assert(sizeof(H) >= sizeof(uint32_t));

    static_assert(W <= 16384 && H <= 16384);
    static_assert(S >= 1 && S <= 256);

   public:
    [[nodiscard]] constexpr int32_t Size() const {
        return T<W, H, S>::image_size;
    }

    [[nodiscard]] constexpr size_t width() const {
        return W;
    }

    [[nodiscard]] constexpr size_t height() const {
        return H;
    }

    constexpr void clear() {
        data.fill(0);
    }

    [[nodiscard]] constexpr int32_t abs(int32_t v) const {
        return v < 0 ? -v : v;
    }

    [[nodiscard]] constexpr std::array<uint8_t, T<W, H, S>::image_size> &dataRef() const {
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
                    fillcircle(y0, x0, (width + 1) / 2, col);
                } else {
                    fillcircle(x0, y0, (width + 1) / 2, col);
                }
                err -= dy;
                if (err < 0) {
                    y0 += ystep;
                    err += dx;
                }
            }
        }
    }

    constexpr void fillrect(int32_t x, int32_t y, int32_t w, int32_t h, uint8_t col) {
        h += y;
        for (; y < h; y++) {
            span(x, w, y, col);
        }
    }

    constexpr void fillcircle(int32_t x, int32_t y, int32_t r, uint8_t col) {
        span(x - abs(r), 2 * abs(r) + 1, y, col);
        fillarc(x, y, abs(r), 3, 0, col);
    }

    constexpr void blitRGBA(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{x, y, w, h};
        blitrect &= {0, 0, W, H};
        blitrect &= {x, y, iw, ih};
        T<W, H, S>::blitRGBA(data, blitrect, ptr, iw, ih, stride);
    }

    constexpr void blitRGBADiffused(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{x, y, w, h};
        blitrect &= {0, 0, W, H};
        blitrect &= {x, y, iw, ih};
        T<W, H, S>::blitRGBADiffused(data, blitrect, ptr, iw, ih, stride);
    }

    constexpr void blitRGBADiffusedLinear(int32_t x, int32_t y, int32_t w, int32_t h, const uint8_t *ptr, int32_t iw, int32_t ih, int32_t stride) {
        rect<int32_t> blitrect{x, y, w, h};
        blitrect &= {0, 0, W, H};
        blitrect &= {x, y, iw, ih};
        T<W, H, S>::blitRGBADiffusedLinear(data, blitrect, ptr, iw, ih, stride);
    }

    template <typename F>
    constexpr void sixel(F &&charOut, bool preserveBackground = false) {
        T<W, H, S>::sixel(data, charOut, {0, 0, W, H}, preserveBackground);
    }

    template <typename F>
    constexpr void sixel(F &&charOut, const rect<int32_t> &r, bool preserveBackground = true) {
        T<W, H, S>::sixel(data, charOut, r, preserveBackground);
    }

    constexpr size_t octree_used_length() const {
        return format.octree.nodes_used_length();
    }
    constexpr size_t octree_memory_length() const {
        return format.octree.nodes_memory_length();
    }

   private:
    constexpr void plot(int32_t x, int32_t y, uint8_t col) {
        size_t _x = static_cast<size_t>(x);
        _x %= W;
        size_t _y = static_cast<size_t>(y);
        _y %= H;
        T<W, H, S>::plot(data, _x, _y, col);
    }

    constexpr void span(int32_t x, int32_t w, int32_t y, uint8_t col) {
        while (y < 0) {
            y += H;
        }
        while (x < 0) {
            x += W;
        }
        size_t _xl = static_cast<size_t>(x);
        _xl %= W;
        size_t _xr = static_cast<size_t>(x + w);
        _xr %= W;
        size_t _y = static_cast<size_t>(y);
        _y %= H;
        if (_xl + static_cast<size_t>(w) < W) {
            T<W, H, S>::span(data, _xl, _xr, _y, col);
        } else {
            T<W, H, S>::span(data, _xl, W - 1, _y, col);
            T<W, H, S>::span(data, 0, _xr, _y, col);
        }
    }

    constexpr void fillarc(int32_t x0, int32_t y0, int32_t r, uint8_t corners, int32_t delta, uint8_t col) {
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

static_assert(constixel::image<constixel::format_2bit, 1, 1>().octree_memory_length() == constixel::image<constixel::format_2bit, 1, 1>().octree_used_length(), "Octree memory size does not match 2bit palette contents.");
static_assert(constixel::image<constixel::format_4bit, 1, 1>().octree_memory_length() == constixel::image<constixel::format_4bit, 1, 1>().octree_used_length(), "Octree memory size does not match 4bit palette contents.");
static_assert(constixel::image<constixel::format_8bit, 1, 1>().octree_memory_length() == constixel::image<constixel::format_8bit, 1, 1>().octree_used_length(), "Octree memory size does not match 8bit palette contents.");

}  // namespace constixel
#endif  // #ifndef _SIXEL_H_
